#!/bin/bash
# Run gatk HaplotypeCaller per individual file concurrently

# Load gatk 4.5.0
setup_params() {
    module load gatk/4.5.0.0-gcc-13.2.0
    module load samtools/1.19.2-gcc-13.2.0
    THREAD=4
    MEMORY=16

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

main() {
    setup_params
    create_index
    make_gvcfs
}

main