#!/bin/bash
# setup file paths per sample
BASE_DIR="/hogan_lab/Hogan_Lab_Shared/Zebrafish_Homozygosity/"
EXCLUDE_STR="chr[0-9]+.vcf|experimental|IMBScorer|Paired|PostProcess|noindels|QC|remRef|viewer_dev"

generate_sample_sheet() {
    local Sample_Identity="$1"
    local Samplesheet_Out="$2"
    echo "Appending ${Sample_Identity} to ${Samplesheet_Out}"

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

    if [ ! -f "${Samplesheet_Out}" ]; then
        printf "Sample_Identity\tFastq_File\tAlignment_File\tVCF_Original\tVCF_Merged\tVCF_ChrFixed\tVCF_Annotated\tVCF_Candidates\tSnzl_NoGaps_NBases\tSnzl_NoGaps_NSnps\tSnzl_WithGaps_NBases\tSnzl_WithGaps_NSnps\tJson\n" > ${Samplesheet_Out}
    else
        printf "${Sample_Identity}\t${Fastq_File}\t${Alignment_File}\t${VCF_Original}\t${VCF_Merged}\t${VCF_ChrFixed}\t${VCF_Annotated}\t${VCF_Candidates}\t${Snzl_NoGaps_NBases}\t${Snzl_NoGaps_NSnps}\t${Snzl_WithGaps_NBases}\t${Snzl_WithGaps_NSnps}\t${Json}\n" >> ${Samplesheet_Out}
    fi
}

parse_sample_names() {
    local dir_path=${BASE_DIR}/Alignment_VariantCall/"$1"
    for i in $(find "$dir_path" -maxdepth 1 -type d | sort | tail -n +2 | grep -Pv "logs|scripts"); do
        echo $(basename ${i})
    done
}

# generate_sample_sheet TL2312073-163-4L-MAN-20231116 samplesheet.tsv

main() {
    SEQRUNS="220930_A00692_316_AHVMJFDSX3 221014_A00692_0319_BHGJ7CDMXY 231201_A00692_0394_231208_A00692_0396"
    
    SEQRUNS="220930_A00692_316_AHVMJFDSX3"
    SAMPLES=$(for i in ${SEQRUNS}; do parse_sample_names ${i}; done | tr ' ' '\n' | sort | uniq)
    for i in ${SAMPLES}; do generate_sample_sheet $i samplesheet.${SEQRUNS}.tmp; done
    SEQRUNS="221014_A00692_0319_BHGJ7CDMXY"
    SAMPLES=$(for i in ${SEQRUNS}; do parse_sample_names ${i}; done | tr ' ' '\n' | sort | uniq)
    for i in ${SAMPLES}; do generate_sample_sheet $i samplesheet.${SEQRUNS}.tmp; done
    SEQRUNS="231201_A00692_0394_231208_A00692_0396"    
    SAMPLES=$(for i in ${SEQRUNS}; do parse_sample_names ${i}; done | tr ' ' '\n' | sort | uniq)
    for i in ${SAMPLES}; do generate_sample_sheet $i samplesheet.${SEQRUNS}.tmp; done
}

main

# # Encapsulate in loop later
# Sample_Identity="TL2312073-163-4L-MAN-20231116"

# # Create a regex pattern to match truncated identifiers (e.g., start of string)
# Sample_Identity_Pattern="${Sample_Identity:0:20}"

# # Use grep -E to match either full or truncated pattern
# exclude_str="chr[0-9]+.vcf|experimental|IMBScorer|Paired|PostProcess|noindels|QC|remRef|viewer_dev"

# # search both full sample identity and truncated sample identity
# # filepaths=$(find . -type f \( -name "*${Sample_Identity}*" -o -name "*${Sample_Identity_Pattern}*" \) | sort | grep -v -P "${exclude_str}")

# # search truncated sample identity only
# filepaths=$(find ${BASE_DIR} -type f -name "*${Sample_Identity_Pattern}*" | sort | grep -v -P "${exclude_str}")

# Fastq_File="NA"

# Alignment_File=$(echo ${filepaths} | tr ' ' '\n' | grep ".bam")

# VCF_Original=$(echo ${filepaths}   | tr ' ' '\n' | grep "UG.vcf.gz$")
# VCF_Merged=$(echo ${filepaths}     | tr ' ' '\n' | grep "merged.vcf.gz$")
# VCF_ChrFixed=$(echo ${filepaths}   | tr ' ' '\n' | grep "ChromFixed.vcf.gz$")
# VCF_Annotated=$(echo ${filepaths}  | tr ' ' '\n' | grep "annotated.vcf.gz$")
# VCF_Candidates=$(dirname $(find ${BASE_DIR} -type f -name "*${Sample_Identity_Pattern}*" | sort | grep -v -P "${exclude_str/chr\[0-9\]\+\.vcf|}" | grep "chr.*vcf" | sed "s|chr.*vcf||" | uniq)).vcf

# Snzl_NoGaps_NBases=$(echo ${filepaths}   | tr ' ' '\n' | grep "nogap_nbases.bedgraph" | sed "s|\.bedgraph||")
# Snzl_NoGaps_NSnps=$(echo ${filepaths}    | tr ' ' '\n' | grep "nogap_nsnps.bedgraph" | sed "s|\.bedgraph||")
# Snzl_WithGaps_NBases=$(echo ${filepaths} | tr ' ' '\n' | grep "withgap_nbases.bedgraph" | sed "s|\.bedgraph||")
# Snzl_WithGaps_NSnps=$(echo ${filepaths}  | tr ' ' '\n' | grep "withgap_nsnps.bedgraph" | sed "s|\.bedgraph||")

# Json=$(echo ${filepaths} | tr ' ' '\n' | grep "_candidates.json$")

# printf "Sample_Identity\tFastq_File\tAlignment_File\tVCF_Original\tVCF_Merged\tVCF_ChrFixed\tVCF_Annotated\tVCF_Candidates\tSnzl_NoGaps_NBases\tSnzl_NoGaps_NSnps\tSnzl_WithGaps_NBases\tSnzl_WithGaps_NSnps\tJson\n" > ${SAMPLESHEET_OUT}
# printf "${Sample_Identity}\t${Fastq_File}\t${Alignment_File}\t${VCF_Original}\t${VCF_Merged}\t${VCF_ChrFixed}\t${VCF_Annotated}\t${VCF_Candidates}\t${Snzl_NoGaps_NBases}\t${Snzl_NoGaps_NSnps}\t${Snzl_WithGaps_NBases}\t${Snzl_WithGaps_NSnps}\t${Json}\n" >> ${SAMPLESHEET_OUT}