#!/bin/bash
#SBATCH --job-name="gatk"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mail-user=Tyrone.Chen@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --partition=rhel_short
#SBATCH --output="gatk.out"
#SBATCH --mem-per-cpu=4000

# Script: 4.run_gatk.sh
# Purpose: Perform germline variant calling on a cohort of BAM files using GATK HaplotypeCaller and GenotypeGVCFs

# Load gatk 4.5.0
module load gatk/4.5.0.0-gcc-13.2.0
module load samtools/1.19.2-gcc-13.2.0
THREAD=16
MEMORY=64

# Setup paths
SAMPLE_METADATA="../data/samplesheet_new.tsv"
OUT_DIR="../results/"
REF_GENOME=$(realpath "../data/Reference/GCA_000002035.2_Zv9_genomic.fna")

# List of BAM files
BAM_FILES=$(cut -f2 ${SAMPLE_METADATA} | tail -n +2)

# Create output directories, GVCF is temporary
VCF_GVCF_DIR="${OUT_DIR}/VCF_GVCF"
VCF_Original="${OUT_DIR}/VCF_Original"
mkdir -p ${VCF_Original}
mkdir -p ${VCF_GVCF_DIR}

# GenomicsDBImport sample map
SAMPLE_MAP="${OUT_DIR}/sample_map.txt"
GDB_DIR="${OUT_DIR}/GDB"
mkdir -p ${GDB_DIR}

create_index() {
    samtools faidx ${REF_GENOME}
    gatk CreateSequenceDictionary -R "${REF_GENOME}" -O "${REF_GENOME/fna/dict}"
}

make_gvcfs() {
    # Step 1: Generate GVCFs for each sample
    # TODO: 
    #   - compare against old pipeline
    #   - formulate as pooling ref genotypes
    for BAM in $(echo "${BAM_FILES}"); do
        OUT_VCF=$(echo $BAM | cut -d '/' -f4)
        gatk --java-options "-Xmx${MEMORY}g" HaplotypeCaller \
            -R "${REF_GENOME}" \
            -I "${BAM}" \
            -O "${VCF_GVCF_DIR}/${OUT_VCF/bam}g.vcf.gz" \
            -ERC GVCF
    done
}

make_samplemap() {
    # Step 2: Create a sample map for GenomicsDBImport
    # TODO:
    #   - as above
    rm -f "${SAMPLE_MAP}"
    for GVCF in "${VCF_GVCF_DIR}/"*.g.vcf.gz; do
        SAMPLE=$(basename "${GVCF}" .g.vcf.gz)
        echo -e "${SAMPLE}\t${GVCF}" >> "${SAMPLE_MAP}"
    done
}

make_genomicsdb() {
    # Step 3: Import GVCFs into GenomicsDB
    # TODO:
    #   - as above
    gatk --java-options "-Xmx${MEMORY}g" GenomicsDBImport \
        --sample-name-map "${SAMPLE_MAP}" \
        --genomicsdb-workspace-path "${GDB_DIR}" \
        --reader-threads "${THREAD}"
}

run_genotyping() {
    # Step 4: Joint genotyping
    # TODO:
    #   - as above
    gatk --java-options "-Xmx${MEMORY}g" GenotypeGVCFs \
        -R "${REF_GENOME}" \
        -V gendb://${GDB_DIR} \
        -O "${VCF_GVCF_DIR}/cohort.vcf.gz"
}

main() {
    create_index
    make_gvcfs
    make_samplemap
    make_genomicsdb
    run_genotyping
    echo "Variant calling pipeline complete. Output VCF: ${VCF_GVCF_DIR}/cohort.vcf.gz"
}

main