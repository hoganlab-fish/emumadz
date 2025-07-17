#!/usr/bin/python
# check all filepaths for samplesheet
import argparse
import os
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
        description="Load samplesheet file into pandas DataFrame"
        )
    parser.add_argument("samplesheet",
                        help="Path to the TSV samplesheet file")
    parser.add_argument("--outfile_path", "-o", default="samplesheet_original.tsv",
                        help="Write out the metadata")
    parser.add_argument("--dropna", action="store_true",
                        help="Drop rows with null values in the samplesheet")
    args = parser.parse_args()

    data = pd.read_csv(args.samplesheet, sep='\t', index_col=0)
    data.drop(["Fastq_File"], axis=1, inplace=True)
    print(f"Drop rows without file paths: {data.isna().any(axis=1).sum()}")
    if args.dropna:
        data.dropna(inplace=True)
    missing_files = list()
    columns = [
        'Alignment_File', 'VCF_Original', 'VCF_Merged', 'VCF_ChrFixed',
        'VCF_Annotated', 'VCF_Candidates', 'Json'
        ]
    for col in columns:
        for filepath in data[col]:
            if type(filepath) == float:
                pass
            elif not os.path.isfile(filepath):
                missing_files.append(filepath)
    colsnzl = [
        'Snzl_NoGaps_NBases', 'Snzl_NoGaps_NSnps', 
        'Snzl_WithGaps_NBases', 'Snzl_WithGaps_NSnps'
        ]                
    for col in colsnzl:
        for filepath in data[col]:
            if type(filepath) == float:
                pass            
            elif not os.path.isfile(".".join([filepath, "bedgraph"])):
                missing_files.append(filepath)
            elif not os.path.isfile(".".join([filepath, "tdf"])):
                missing_files.append(filepath)

    if missing_files:
        print("Missing files or invalid file paths:")
        for f in missing_files:
            print(f)
    else:
        print("All remaining file paths valid.")

    # add short sample id
    data["Sample_Short"] = data.index.str.replace(
        r'^((?:[^-]*-){2}[^-]*)-.*', r'\1', regex=True
        )
    data.to_csv(args.outfile_path, sep="\t")

if __name__ == "__main__":
    main()