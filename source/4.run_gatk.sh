#!/bin/bash
# Script: 4.run_gatk.sh
# Purpose: Perform germline variant calling on a cohort of BAM files using GATK HaplotypeCaller and GenotypeGVCFs

# Load gatk 4.5.0
setup_params() {
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
    GDB_DIR="../results/GDB"
}

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
        bash 4.1.make_gvcfs.slurm \
            $(basename ${BAM/.bam}) \
            ${MEMORY} \
            ${REF_GENOME} \
            ${BAM} \
            ${VCF_GVCF_DIR}/${OUT_VCF/bam}g.vcf.gz \
            ${THREAD}
    # this is a batch submission
    # slurm scripts are excluded from the repository
    # we show the command that would be run
    # gatk --java-options "-Xmx${MEMORY}g" HaplotypeCaller \
    #     -R "${REF_GENOME}" \
    #     -I "${BAM}" \
    #     -O "${VCF_GVCF_DIR}/${OUT_VCF/bam}g.vcf.gz" \
    #     -ERC GVCF
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
        --genomicsdb-workspace-path "${GDB_DIR}" \
        $(cut -f2 $SAMPLE_MAP | sed -e "s|^|\-V |") \
        $(cut -f1 ${REF_GENOME}.fai | head -n 25 | sed "s|^|--intervals |") \
        --reader-threads "${THREAD}"
        # --sample-name-map "${SAMPLE_MAP}" \    
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
    echo "note:"
    echo "  - do not run this directly"
    echo "  - each of the steps below submits one or more slurm jobs"
    echo "  - wait until each completes successfully (can take days)"
    echo "  - then run the next step"

    echo "Step 4.1: bash 4.1.make_gvcfs.sh"
    # bash 4.1.make_gvcfs.sh
    echo "Step 4.2: bash 4.2.make_genomicsdb.slurm"
    echo "Outside a cluster env, this is equivalent to running the following 3 functions:"
    echo "  setup_params"
    echo "  make_samplemap [OPTIONAL]"
    echo "  make_genomicsdb"
    # setup_params
    # make_samplemap [OPTIONAL]
    # make_genomicsdb
    echo "Step 4.3: bash 4.3.joint_genotyping.slurm"
    echo "  coming soon"
    # echo "Variant calling pipeline complete. Output VCF: ${VCF_GVCF_DIR}/cohort.vcf.gz"
}

main