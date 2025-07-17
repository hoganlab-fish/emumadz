#!/bin/bash
# setup new IO directory structure based on original input samplesheet
SAMPLESHEET_ORIGINAL="../data/samplesheet_original.tsv"
SAMPLESHEET_NEW="../data/samplesheet_new.tsv"
ALN_DIR="../data/Alignment_File/"
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
        printf "Sample_Identity\tAlignment_File\tVCF_Original\tVCF_Merged\tVCF_ChrFixed\tVCF_Annotated\tVCF_Candidates\tSnzl_NoGaps_NBases\tSnzl_NoGaps_NSnps\tSnzl_WithGaps_NBases\tSnzl_WithGaps_NSnps\tJson\tSample_Short\n" > ${SAMPLESHEET_NEW}
        for i in $(cut -f13 ${SAMPLESHEET_ORIGINAL} | grep -v 'Sample_Short'); do
            printf "${i}\t../data/Alignment_File/${i}.bam\t\t\t\t\t\t\t\t\t\t\t${i}\n" >> ${SAMPLESHEET_NEW}
        done    
    else
        for i in $(cut -f13 ${SAMPLESHEET_ORIGINAL} | grep -v 'Sample_Short'); do
            printf "${i}\t../data/Alignment_File/${i}.bam\t\t\t\t\t\t\t\t\t\t\t${i}\n" >> ${SAMPLESHEET_NEW}
        done
    fi
}

main() {
    copy_bams
    rename_bams
    make_outdirs
    make_samplesheet    
}

main