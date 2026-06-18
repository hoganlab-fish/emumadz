#!/bin/bash

# Usage: ./rename_bam_chroms.sh mapping.txt input.bam output.bam
# mapping.txt: tab-separated file with old_name (col1) and new_name (col2)

MAPPING_FILE=$1
INPUT_BAM=$2
OUTPUT_BAM=$3

if [ $# -ne 3 ]; then
    echo "Usage: $0 <mapping.txt> <input.bam> <output.bam>"
    exit 1
fi

# Extract header and rename chromosomes using awk
samtools view -H "$INPUT_BAM" | awk -v mapfile="$MAPPING_FILE" '
BEGIN {
    FS = OFS = "\t"
    # Load mapping into array
    while ((getline line < mapfile) > 0) {
        split(line, arr, "\t")
        map[arr[1]] = arr[2]
    }
    close(mapfile)
}
{
    # Process @SQ lines
    if ($1 == "@SQ") {
        for (i = 2; i <= NF; i++) {
            if ($i ~ /^SN:/) {
                old_name = substr($i, 4)
                if (old_name in map) {
                    $i = "SN:" map[old_name]
                }
            }
        }
    }
    print
}
' > new_header.sam

# Reheader the BAM file
samtools reheader new_header.sam "$INPUT_BAM" > "$OUTPUT_BAM"

# Index the output BAM
samtools index "$OUTPUT_BAM"

# Clean up
rm -f new_header.sam

echo "Done! Output written to $OUTPUT_BAM"
