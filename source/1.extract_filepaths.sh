#!/bin/bash
# setup file paths per sample
BASE_DIR="/hogan_lab/Hogan_Lab_Shared/Zebrafish_Homozygosity/"
DATA_DIR="../data/"
EXCLUDE_STR="candidates_old|chr[0-9]+.vcf|debug|experimental|IMBScorer|log|Paired|PostProcess|noindels|QC|remRef|scripts|viewer_dev"

generate_sample_sheet() {
    local Sample_Identity="$1"
    local Samplesheet_Out="$2"
    echo "Appending ${Sample_Identity} to ${Samplesheet_Out}"

    if [ -n "$3" ]; then 
        local exclude_str_original="${EXCLUDE_STR}"
        EXCLUDE_STR="${3}|${exclude_str_original}"
    fi

    local Sample_Identity_Pattern="${Sample_Identity:0:20}"
    
    local filepaths=$(find ${BASE_DIR} -type f -name "*${Sample_Identity_Pattern}*" | sort | grep -v -P "${EXCLUDE_STR}")

    local Fastq_File="NA"
    local Alignment_File=$(echo ${filepaths} | tr ' ' '\n' | grep ".bam$")
    local VCF_Original=$(echo ${filepaths}   | tr ' ' '\n' | grep "UG.vcf.gz$")
    local VCF_Merged=$(echo ${filepaths}     | tr ' ' '\n' | grep "merged.vcf.gz$")
    local VCF_ChrFixed=$(echo ${filepaths}   | tr ' ' '\n' | grep "ChromFixed.vcf.gz$")
    local VCF_Annotated=$(echo ${filepaths}  | tr ' ' '\n' | grep "annotated.vcf.gz$")
    local VCF_Candidates=$(dirname $(find ${BASE_DIR} -type f -name "*${Sample_Identity_Pattern}*" | sort | grep -v -P "${EXCLUDE_STR/chr\[0-9\]\+\.vcf|}" | grep "chr.*vcf" | sed "s|chr.*vcf||" | uniq)).vcf

    local Snzl_NoGaps_NBases=$(echo ${filepaths}   | tr ' ' '\n' | grep "nogap_nbases.bedgraph" | sed "s|\.bedgraph||")
    local Snzl_NoGaps_NSnps=$(echo ${filepaths}    | tr ' ' '\n' | grep "nogap_nsnps.bedgraph" | sed "s|\.bedgraph||")
    local Snzl_WithGaps_NBases=$(echo ${filepaths} | tr ' ' '\n' | grep "withgap_nbases.bedgraph" | sed "s|\.bedgraph||")
    local Snzl_WithGaps_NSnps=$(echo ${filepaths}  | tr ' ' '\n' | grep "withgap_nsnps.bedgraph" | sed "s|\.bedgraph||")

    local Json=$(echo ${filepaths} | tr ' ' '\n' | grep "_candidates.json$")

    echo "Column:Data assignment for ${Sample_Identity}"
    paste <(printf '${Sample_Identity}\t${Fastq_File}\t${Alignment_File}\t${VCF_Original}\t${VCF_Merged}\t${VCF_ChrFixed}\t${VCF_Annotated}\t${VCF_Candidates}\t${Snzl_NoGaps_NBases}\t${Snzl_NoGaps_NSnps}\t${Snzl_WithGaps_NBases}\t${Snzl_WithGaps_NSnps}\t${Json}\n' | tr '\t' '\n' | tr -d '${}' ) <(echo "${Sample_Identity} ${Fastq_File} ${Alignment_File} ${VCF_Original} ${VCF_Merged} ${VCF_ChrFixed} ${VCF_Annotated} ${VCF_Candidates} ${Snzl_NoGaps_NBases} ${Snzl_NoGaps_NSnps} ${Snzl_WithGaps_NBases} ${Snzl_WithGaps_NSnps} ${Json}" | tr ' ' '\n')

    if [ ! -f "${Samplesheet_Out}" ]; then
        printf "Sample_Identity\tFastq_File\tAlignment_File\tVCF_Original\tVCF_Merged\tVCF_ChrFixed\tVCF_Annotated\tVCF_Candidates\tSnzl_NoGaps_NBases\tSnzl_NoGaps_NSnps\tSnzl_WithGaps_NBases\tSnzl_WithGaps_NSnps\tJson\n" > ${Samplesheet_Out}
    else
        printf "${Sample_Identity}\t${Fastq_File}\t${Alignment_File}\t${VCF_Original}\t${VCF_Merged}\t${VCF_ChrFixed}\t${VCF_Annotated}\t${VCF_Candidates}\t${Snzl_NoGaps_NBases}\t${Snzl_NoGaps_NSnps}\t${Snzl_WithGaps_NBases}\t${Snzl_WithGaps_NSnps}\t${Json}\n" >> ${Samplesheet_Out}
    fi

    if [ -n "$3" ]; then EXCLUDE_STR=${exclude_str_original}; fi
}

parse_sample_names() {
    local dir_path=${BASE_DIR}/Alignment_VariantCall/"$1"
    for i in $(find "$dir_path" -maxdepth 1 -type d | sort | tail -n +2 | grep -Pv "logs|scripts"); do
        echo $(basename ${i})
    done
}

# generate_sample_sheet TL2312073-163-4L-MAN-20231116 samplesheet.tsv

main() {
    mkdir -p ${DATA_DIR}

    # f0 reference generation
    SEQRUN_220701="220701_A01221_0125_BHCCHNDMXY"
    SAMPLES=$(for i in ${SEQRUN_220701}; do parse_sample_names ${i}; done | tr ' ' '\n' | sort | uniq)
    
    Samplesheet_Ref=${DATA_DIR}/samplesheet.${SEQRUN_220701}.tsv
    if [ ! -f "${Samplesheet_Ref}" ]; then
        printf "Sample_Identity\tFastq_File\tAlignment_File\tVCF_Original\tVCF_Merged\tVCF_ChrFixed\tVCF_Annotated\tVCF_Candidates\tSnzl_NoGaps_NBases\tSnzl_NoGaps_NSnps\tSnzl_WithGaps_NBases\tSnzl_WithGaps_NSnps\tJson\n" > ${Samplesheet_Ref}
    fi
    for i in ${SAMPLES}; do 
        bam=$(find ${BASE_DIR}/Alignment_VariantCall/${SEQRUN_220701}/${i} -type f | grep "bam$")
        vcf=$(find ${BASE_DIR}/Alignment_VariantCall/${SEQRUN_220701}/${i} -type f | grep "vcf$")
        printf "${i}\tNA\t${bam}\t${vcf}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n" >> \
            ${Samplesheet_Ref}
    done

    # all normal exclusion patterns should apply here
    SEQRUN_220930="220930_A00692_316_AHVMJFDSX3"
    SAMPLES=$(for i in ${SEQRUN_220930}; do parse_sample_names ${i}; done | tr ' ' '\n' | sort | uniq)
    for i in ${SAMPLES}; do generate_sample_sheet $i ${DATA_DIR}/samplesheet.${SEQRUN_220930}.tsv; done
    
    # two versions of bedgraphs exist, pick the latest ones
    SEQRUN_221014="221014_A00692_0319_BHGJ7CDMXY"
    SAMPLES=$(for i in ${SEQRUN_221014}; do parse_sample_names ${i}; done | tr ' ' '\n' | sort | uniq)
    for i in ${SAMPLES}; do generate_sample_sheet $i ${DATA_DIR}/samplesheet.${SEQRUN_221014}.tsv bedgraph/; done
    
    # exclude two specific incomplete/irrelevant samples
    SEQRUN_231201="231201_A00692_0394_231208_A00692_0396"    
    SAMPLES=$(for i in ${SEQRUN_231201}; do parse_sample_names ${i}; done | tr ' ' '\n' | sort | uniq | \
    grep -Pv '175-1-Sib|TL233872-22-1-I-MAN-20231124')
    for i in ${SAMPLES}; do generate_sample_sheet $i ${DATA_DIR}/samplesheet.${SEQRUN_231201}.tsv; done
    
    (cat ${Samplesheet_Ref}; 
        tail -n +2 ${DATA_DIR}/samplesheet.${SEQRUN_220930}.tsv; \
        tail -n +2 ${DATA_DIR}/samplesheet.${SEQRUN_221014}.tsv; \
        tail -n +2 ${DATA_DIR}/samplesheet.${SEQRUN_231201}.tsv) \
        > ${DATA_DIR}/samplesheet.tmp

    python 2.validate_samplesheet.py \
        ${DATA_DIR}/samplesheet.tmp \
        -o ${DATA_DIR}/samplesheet_original.tsv
}

main
