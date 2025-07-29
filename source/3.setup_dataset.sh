#!/bin/bash
# setup new IO directory structure based on original input samplesheet
SAMPLESHEET_ORIGINAL="../data/samplesheet_original.tsv"
SAMPLESHEET_NEW="../data/samplesheet_new.tsv"
SAMPLESHEET_VCF="../data/samplesheet_postprocess.tsv"
ALN_DIR="../data/Alignment_File/"
VCF_DIR="../data/VCF_Original/"
SRC_DIR="../source/"

copy_bams () {
    # input goes here
    mkdir -p ${ALN_DIR}

    # copy files over to specified data dir
    for i in $(cut -f2 ${SAMPLESHEET_ORIGINAL} | grep -v 'Alignment_File'); do
        cp ${i} ${ALN_DIR};
        cp ${i/bam/bai} ${ALN_DIR};
    done
}

copy_vcfs () {
    # input goes here
    mkdir -p ${VCF_DIR}

    # copy files over to specified data dir
    for i in $(cut -f3 ${SAMPLESHEET_ORIGINAL} | grep -v 'VCF_Original'); do
        cp ${i} ${VCF_DIR};
    done
}

index_vcfs () {
    for i in ${VCF_DIR}/*vcf; do
        bgzip ${i};
    done
    for i in ${VCF_DIR}/*vcf.gz; do
        bcftools index -t ${i}
    done
}

rename_vcfs() {
    # long names get truncated by some pipeline modules, shorten first
    cd ${VCF_DIR}
    for vcf in *vcf.gz; do
        ln -s $(basename ${vcf}) $(basename $vcf | cut -d '-' -f1-3).vcf.gz
    done
    for tbi in *vcf.gz.tbi; do
        ln -s $(basename ${tbi}) $(basename $tbi | cut -d '-' -f1-3).vcf.gz.tbi
    done
    cd ${SRC_DIR}
}

rename_bams() {
    # long names get truncated by some pipeline modules, shorten first
    cd ${ALN_DIR}
    for bam in *bam; do
        ln -s $(basename ${bam}) $(basename $bam | cut -d '-' -f1-3).bam
    done
    for bai in *bai; do
        ln -s $(basename ${bai}) $(basename $bai | cut -d '-' -f1-3).bai
    done
    cd ${SRC_DIR}
}

make_outdirs() {
    # results go here
    mkdir -p $(echo \
        $(head -n 1 ${SAMPLESHEET_ORIGINAL} | \
            tr '\t' '\n' | \
            sed "s|^|..\/results\/|" | \
            tail -n +3))
}

make_samplesheet() {
    if [ ! -f "${SAMPLESHEET_NEW}" ]; then
        printf "Sample_Identity\tAlignment_File\tVCF_Original\tVCF_Merged\tVCF_Norm\tVCF_Snps\tVCF_ChrFixed\tVCF_Annotated\tVCF_Candidates_Strict\tVCF_Candidates_NonStrict\tSnzl_NoGaps_NBases\tSnzl_NoGaps_NSnps\tSnzl_WithGaps_NBases\tSnzl_WithGaps_NSnps\tJson\tSample_Short\n" > ${SAMPLESHEET_NEW}
        for i in $(cut -f13 ${SAMPLESHEET_ORIGINAL} | grep -v 'Sample_Short'); do
            printf "${i}\t../data/Alignment_File/${i}.bam\t\t\t\t\t\t\t\t\t\t\t\t\t\t${i}\n" >> ${SAMPLESHEET_NEW}
        done    
    else
        for i in $(cut -f13 ${SAMPLESHEET_ORIGINAL} | grep -v 'Sample_Short'); do
            printf "${i}\t../data/Alignment_File/${i}.bam\t\t\t\t\t\t\t\t\t\t\t\t\t\t${i}\n" >> ${SAMPLESHEET_NEW}
        done
    fi
}

make_samplesheet_postprocess() {
    if [ ! -f "${SAMPLESHEET_VCF}" ]; then
        printf "Sample_Identity\tAlignment_File\tVCF_Original\tVCF_Merged\tVCF_Norm\tVCF_Snps\tVCF_ChrFixed\tVCF_Annotated\tVCF_Candidates_Strict\tVCF_Candidates_NonStrict\tSnzl_NoGaps_NBases\tSnzl_NoGaps_NSnps\tSnzl_WithGaps_NBases\tSnzl_WithGaps_NSnps\tJson\tSample_Short\n" > ${SAMPLESHEET_VCF}
        for i in $(cut -f13 ${SAMPLESHEET_ORIGINAL} | grep -v 'Sample_Short'); do
            printf "${i}\t../data/Alignment_File/${i}.bam\t../data/VCF_Original/${i}.vcf.gz\t\t\t\t\t\t\t\t\t\t\t\t\t${i}\n" >> ${SAMPLESHEET_VCF}
        done    
    else
        for i in $(cut -f13 ${SAMPLESHEET_ORIGINAL} | grep -v 'Sample_Short'); do
            printf "${i}\t../data/Alignment_File/${i}.bam\t../data/VCF_Original/${i}.vcf.gz\t\t\t\t\t\t\t\t\t\t\t\t\t${i}\n" >> ${SAMPLESHEET_VCF}
        done
    fi    
}

main() {
    copy_bams
    rename_bams
    copy_vcfs
    index_vcfs
    rename_vcfs
    make_outdirs
    make_samplesheet
    make_samplesheet_postprocess 
}

main