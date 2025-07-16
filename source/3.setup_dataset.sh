#!/bin/bash
# setup new IO directory structure based on original input samplesheet
SAMPLESHEET_ORIGINAL="../data/samplesheet_original.tsv"
ALN_DIR="../data/Alignment_File/"

# input goes here
mkdir -p ${ALN_DIR}

# copy files over to specified data dir
for i in $(cut -f2 ../data/samplesheet_original.tsv | grep -v 'Alignment_File'); do
    echo ${i}
    cp ${i} ${ALN_DIR};
    cp ${i/bam/bai} ${ALN_DIR};
done

# results go here
mkdir -p $(echo \
    $(head -n 1 ../data/samplesheet_original.tsv | \
        tr '\t' '\n' | \
        sed "s|^|..\/results\/|" | \
        tail -n +3))

# TODO: create new samplesheet pointing to the relevant directories for output