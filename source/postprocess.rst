Quickstart
==========

.. whole_genome_sequencing documentation master file, created by
   sphinx-quickstart on Mon Jun 30 12:23:50 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. raw:: html

  Copyright (c) 2025 <a href="https://orcid.org/0000-0002-9207-0385">Tyrone Chen, <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a><a href="https://orcid.org/0000-0003-3746-0695">Oliver Yu, <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a><a href="https://orcid.org/0000-0002-0651-7065">Benjamin Hogan <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>.

Code in this repository is provided under a `MIT license`_. This documentation is provided under a `CC-BY-3.0 AU license`_.

.. _MIT license: https://opensource.org/licenses/MIT

.. _CC-BY-3.0 AU license: https://creativecommons.org/licenses/by/3.0/au/

`Visit our lab website here.`_ Contact Benjamin Hogan at `ben.hogan@petermac.org`_.

.. _Visit our lab website here.: https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/hogan-laboratory-vascular-developmental-genetics-and-cell-biology

.. _ben.hogan@petermac.org: mailto:ben.hogan@petermac.org


Usage
-----

.. note::
    All steps are run in directory ``source``.

Data setup
++++++++++

Original data
#############

.. caution::
   The original data setup was sample-centric as opposed to more conventional process-centric directory structure. In addition, code had hardcoded paths which limited file movement. This impacted all steps of analysis and reduced reproducibility. To lower the impact of the original file layout, the data directory now contains samplesheets with updated file paths to symlinks and their attributes. However, it is possible that some errors remain.

Run ``1.extract_filepaths.sh`` with the required command line arguments to generate the samplesheets with paths to old files. The samplesheets are tab separated files with the following fields:

.. csv-table:: Samplesheet column information
   :file: tables/samplesheet_fields.csv
   :header-rows: 1

.. warning::
   A downstream step silently truncates file names if it exceeds a certain limit. No explanation for this behaviour was available. The ``1.extract_filepaths.sh`` script will try to "guess" file paths by performing an arbitrary truncation and subsequent substring matching.

.. note::
   ``Vis_`` fields should contain the file name but not the file extensions, both ``bedgraph`` and ``tdf`` files will be generated. Note that both files contain the same information but ``tdf`` is optimised for viewing with the ``IGV`` Genome Browser.

.. csv-table:: Example samplesheet showing one sample
   :file: tables/samplesheet_example.csv
   :header-rows: 1

The ``2.validate_samples.py`` script checks the validity of each file in the samplesheet per sample. This also drops the ``fastq`` column. Automatically runs in ``1.extract_filepaths.sh``.

.. code:: shell

   python 2.validate_samples.py ../data/samplesheet.tmp -o ../data/samplesheet.tsv

.. caution::
   Custom directories and/or files exist since the data passed through multiple iterations. To some extent this is accommodated in the setup and validation scripts, but make sure to double-check everything.

For the purposes of rerunning the postprocessing, we only want the ``bam`` alignment files and original ``vcf`` files, since we will be recreating everything after this step. Running ``3.setup_dataset.sh`` will create the corresponding input and output directories::

   ../data/Alignment_File
   ../data/VCF_Original
   ../results/Sample_Identity
   ../results/Fastq_File
   ../results/Alignment_File
   ../results/VCF_Merged
   ../results/VCF_ChrFixed
   ../results/VCF_Annotated
   ../results/VCF_Candidates
   ../results/Snzl_NoGaps_NBases
   ../results/Snzl_NoGaps_NSnps
   ../results/Snzl_WithGaps_NBases
   ../results/Snzl_WithGaps_NSnps
   ../results/Json
   
While the ``bam`` and ``vcf`` files are not directly included in this repository due to size constraints, ``md5`` sums are preserved in ``data/Alignment_File/bam.md5`` and ``data/VCF_Original/vcf.md5``.

Four samplesheet files are generated:
- ``samplesheet.220701_A01221_0125_BHCCHNDMXY.tsv``
- ``samplesheet.220930_A00692_316_AHVMJFDSX3.tsv``
- ``samplesheet.221014_A00692_0319_BHGJ7CDMXY.tsv``
- ``samplesheet.231201_A00692_0394_231208_A00692_0396.tsv``

The ``samplesheet.220701_A01221_0125_BHCCHNDMXY.tsv`` contains ``F0`` grandparent reference generation metadata. All other samplesheets contain ``F2`` generation metadata. ``samplesheet_original.tsv`` contains these legacy file paths.

For postprocessing (steps not carried out)
##########################################

``3.setup_dataset.sh`` copies ``bam`` alignment files and ``vcf`` variant calling files over from their specified locations in ``samplesheet_original.tsv`` to ``data/Alignment_File`` and ``data/VCF_Original`` respectively. We generate ``samplesheet_postprocess.tsv`` for use in this postprocessing analysis.

Whole zebrafish genome assembly
+++++++++++++++++++++++++++++++

*The ``seqliner`` assembly step is not covered in this documentation. The pipeline starts from the variant calling files*

Variant calling
+++++++++++++++

*The ``gatk`` variant calling step is not covered in this documentation. The pipeline starts from the variant calling files*

For postprocessing (steps carried out)
######################################

Rename chromosomes in fasta file with map file
++++++++++++++++++++++++++++++++++++++++++++++

``data/chrom_table.tsv`` contains the chromosome name mappings. Use this to remap the chromosomes.

.. code-block:: shell

    module load samtools

    REF_GEN="../data/Reference/GCA_000002035.2_Zv9_genomic.fna"
    CHR_MAP="../data/chrom_table.tsv"
    OUT_GEN="../data/Reference/GCA_000002035.2_Zv9_genomic.ChrFixed.fna"

    {
        while read seqname length rest; do
            if grep -q "^${seqname}" ${CHR_MAP}; then
                newname=$(awk -v seq="$seqname" '$1==seq {print $2}' ${CHR_MAP})
                samtools faidx ${REF_GEN} "$seqname" | sed "s/^>$seqname/>$newname/"
            else
                samtools faidx ${REF_GEN} "$seqname"
            fi
        done < ${REF_GEN}.fai
    } > ${OUT_GEN}

Rename chromosomes in vcf with map file, and discard non-chromosomes
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

``data/chrom_table.tsv`` contains the chromosome name mappings. Use this to remap the chromosomes.

At the same time, we are only interested in the chromosomes, which start with ``chr`` and go from ``chr1`` to ``chr25``. So we discard all other contigs (including from the header).

For one example:

.. code-block:: shell

    THREADS=16
    # 33550336
    for i in {1..25}; do 
        echo -e "chr${i}\t1\t1111111111111000000000000"
    done > ../data/chromosomes_only.txt

    bcftools annotate --threads ${THREADS} \
        --rename-chr "../data/chrom_table.tsv" \
        ../data/VCF_Original/TL2312073-163-4L.vcf.gz | \
    grep -v '##contig=<ID=FR' | \
    bgzip -c ../results/VCF_ChrFixed/TL2312073-163-4L.vcf.gz
    bcftools index -t ../results/VCF_ChrFixed/TL2312073-163-4L.vcf.gz

For an example on all samples:

.. code-block:: shell

    THREADS=16
    metadata="../data/samplesheet_postprocess.tsv"
    CHR_MAP="../data/chrom_table.tsv"
    CHR_FIL="../data/chromosomes_only.txt"
    # 33550336
    for i in {1..25}; do 
        echo -e "chr${i}\t1\t1111111111111000000000000"
    done > ../data/chromosomes_only.txt

    tail -n +2 ${metadata} | \
        while IFS="\t" read -r line; do
            in=$(echo $line | cut -d ' ' -f3)
            out="../results/VCF_ChrFixed/"$(basename ${in})
            bcftools annotate \
                --threads ${THREADS} \
                --rename-chrs ${CHR_MAP} \
                ${in} | \
            grep -v '##contig=<ID=FR' | \
            bgzip -c > ${out}
            bcftools index -t ${out}
            chgrp -R hogan_lab_bioinf ${out}*
        done     
        
Normalise VCF files
+++++++++++++++++++

For one example:

.. code-block:: shell

    # $REF_GEN is the previous $OUT_GEN, with fixed chromosome names
    REF_GEN="../data/Reference/GCA_000002035.2_Zv9_genomic.ChrFixed.fna"
    THREADS=16
    bcftools norm --threads ${THREADS} -m-both -f ${REF_GEN} \
        ../results/VCF_ChrFixed/TL2312073-163-4L.vcf.gz \
        -Ob -o ../results/VCF_Norm/TL2312073-163-4L.vcf.gz
    bcftools index --threads ${THREADS} -t \
        ../results/VCF_Norm/TL2312073-163-4L.vcf.gz

For an example on all samples:

.. code-block:: shell

    # $REF_GEN is the previous $OUT_GEN, with fixed chromosome names
    REF_GEN="../data/Reference/GCA_000002035.2_Zv9_genomic.ChrFixed.fna"
    THREADS=16
    INFILE_DIR="../results/VCF_ChrFixed/*gz"
    OUTFILE_DIR="../results/VCF_Norm/"

    find ${INFILE_DIR} | sort | \
        while IFS="\t" read -r line; do
            in=$(echo $line)
            out=${OUTFILE_DIR}$(basename ${in})
            bcftools norm --threads ${THREADS} -m-both \
                -f ${REF_GEN} ${in} -Ob -o ${out}
            bcftools index --threads ${THREADS} -t ${out}
        done


Merge files with refs
+++++++++++++++++++++

For each sample, merge the four references:
- TL2209397-ENUref-Female.vcf.gz
- TL2209398-ENUref-Male.vcf.gz
- TL2209399-TUref-Female.vcf.gz
- TL2209400-TUref-Male.vcf.gz

For one example:

.. code-block:: shell

    THREADS=16
    REF_PATHS="../results/VCF_ChrFixed/TL2209397-ENUref-Female.vcf.gz
    ../results/VCF_ChrFixed/TL2209398-ENUref-Male.vcf.gz
    ../results/VCF_ChrFixed/TL2209399-TUref-Female.vcf.gz
    ../results/VCF_ChrFixed/TL2209400-TUref-Male.vcf.gz"

    bcftools merge --threads ${THREADS} \
        ../results/VCF_ChrFixed/TL2312073-163-4L.vcf.gz ${REF_PATHS} \
        -Ob -o ../results/VCF_Merged/TL2312073-163-4L.vcf.gz
    bcftools index --threads ${THREADS} -t \
        ../results/VCF_Merged/TL2312073-163-4L.merged.vcf.gz

For an example on all samples:

.. code-block:: shell

    THREADS=16
    INFILE_DIR="../results/VCF_ChrFixed/*gz"
    OUTFILE_DIR="../results/VCF_Merged/"
    REF_PATHS="../results/VCF_ChrFixed/TL2209397-ENUref-Female.vcf.gz
    ../results/VCF_ChrFixed/TL2209398-ENUref-Male.vcf.gz
    ../results/VCF_ChrFixed/TL2209399-TUref-Female.vcf.gz
    ../results/VCF_ChrFixed/TL2209400-TUref-Male.vcf.gz"

    find ${INFILE_DIR} | grep -v 'ref' | \
        while IFS="\t" read -r line; do
            in=$(echo $line)
            out=${OUTFILE_DIR}$(basename ${in})
            bcftools merge --threads ${THREADS} \
                ${in} ${REF_PATHS} -Ob -o ${out}
            bcftools index --threads ${THREADS} -t ${out}
        done

Filter out non-SNPs
+++++++++++++++++++

For one example:

.. code-block:: shell

    THREADS=16
    bcftools view --threads ${THREADS} -i 'TYPE="snp"' \
        ../results/VCF_Merged/TL2312073-163-4L.vcf.gz \
        -Ob -o ../results/VCF_Snps/TL2312073-163-4L.vcf.gz
    bcftools index --threads ${THREADS} -t \
        ../results/VCF_SNPs/TL2312073-163-4L.vcf.gz

For an example on all samples:

.. code-block:: shell

    THREADS=16
    INFILE_DIR="../results/VCF_Merged/*gz"
    OUTFILE_DIR="../results/VCF_Snps/"

    find ${INFILE_DIR} | sort | \
        while IFS="\t" read -r line; do
            in=$(echo $line)
            out=${OUTFILE_DIR}$(basename ${in})
            bcftools view --threads ${THREADS} \
                -i 'TYPE="snp"' ${in} -Ob -o ${out}
            bcftools index --threads ${THREADS} -t ${out}
        done

Obtain SNPs which occur in the sample but not the references
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Unrefined implementation (works only if grandparents have the reference allele)
*******************************************************************************

Ideally we only should include samples homozygous for the SNPs. However, false heterozygotes appear in the data due to reads with low mapping quality introducing false heterozygoisty. Therefore, this filter is intentionally relaxed to accommodate heterozygotes. The information can be verified by inspecting the ``bam`` file directly.

For one example:

.. code-block:: shell

    THREADS=16
    bcftools view --threads ${THREADS} \
        -i 'GT[0]="hom" && F_MISSING=0.8 || GT[0]="het" && F_MISSING=0.8' \
        ../results/VCF_Snps/TL2312073-163-4L.vcf.gz \
        -Ob -o ../results/VCF_Candidates/TL2312073-163-4L.vcf.gz

To verify that the conditional statements are formulated correctly and that this returns expected results:

.. code-block:: shell

    # both return the same line count of 1481032
    THREADS=16
    bcftools view --threads ${THREADS} \
        -i 'GT[0]="hom" && F_MISSING=0.8 || GT[0]="het" && F_MISSING=0.8' \
        ../results/VCF_Snps/TL2312073-163-4L.vcf.gz | \
        grep -v '#' | wc -l
    bcftools view --threads ${THREADS}\
        -i 'GT[0]="hom" && F_MISSING=0.8 || GT[0]="het" && F_MISSING=0.8' \
        ../results/VCF_Snps/TL2312073-163-4L.vcf.gz | \
        grep -Pc './.:.:.:.:.\t./.:.:.:.:.\t./.:.:.:.:.\t./.:.:.:.:.$'

For an example on all samples:

.. code-block:: shell

    THREADS=16
    INFILE_DIR="../results/VCF_Snps/*gz"
    OUTFILE_DIR="../results/VCF_Candidates/"

    find ${INFILE_DIR} | sort | \
        while IFS="\t" read -r line; do
            in=$(echo $line)
            out=${OUTFILE_DIR}$(basename ${in})
            bcftools view --threads ${THREADS} \
                -i 'GT[0]="hom" && F_MISSING=0.8 || GT[0]="het" && F_MISSING=0.8' \
                ${in} -Ob -o ${out}
            bcftools index --threads ${THREADS} -t ${out}
        done

.. note::
    To be refined. Cannot just obtain all missing snps in refs. e.g:
        - Reference allele may be A
        - Grandparents may be C
        - Mutant may be G

Refined implementation (works if grandparents have the reference allele)
************************************************************************

``vcf`` file specifications:
- We have a multisample vcf file
- The first sample is the mutant
- The remaining 4 samples are the grandparent `F0` samples used as a reference
- Each sample is the result of variant calling against the primary reference genome

Original claimed criteria
^^^^^^^^^^^^^^^^^^^^^^^^^

Candidate criteria:
    1. The mutant must be covered to a depth of >= 1 read
    2. At least one of the unaffected sibling or other samples must be covered to a depth of >= 1 read
    3. If the mutant and sibling have a single allele, it must not be the same in both
    4. The mutant must have a majority allele that is:
        a. not the majority allele in any of the 'other' samples'''

Strict candidate criteria:
    1. Mutant is homozygous
    2. All of the 'other' samples are called
    3. Mutant allele is not present in any of the 'other' samples

Revised criteria
^^^^^^^^^^^^^^^^

Candidate criteria:
    1. The mutant must be covered to a depth of >= 1 read
    2. At least one of the other samples must be covered to a depth of >= 1 read
    3. If the mutant and grandparent reference have a single allele, it must not be the same in both
    4. The mutant must have a majority allele that is not the majority allele in any of the F0 references

Strict candidate criteria:
    1. Mutant is homozygous
    2. All of the other F0 references are called
    3. Mutant allele is not present in any of the other F0 references

Modular bcftools-native implementation (pseudocode):

.. code-block:: shell

    THREADS=16
    INFILE_DIR="../results/VCF_Snps/*gz"
    OUTFILE_DIR="../results/VCF_Candidates/"

    bcftools view  --threads ${THREADS} -i \
        'FORMAT/DP[0]>=1 && (FORMAT/DP[1]>=1 || FORMAT/DP[2]>=1 || FORMAT/DP[3]>=1 || FORMAT/DP[4]>=1)' \
        ../results/VCF_Snps/TL2312073-163-4L.vcf.gz | \        
    bcftools view --threads ${THREADS} -i \
        
    # 3. If the mutant and grandparent reference have a single allele, it must not be the same in both
    # 4. The mutant must have a majority allele that is not the majority allele in any of the references

.. code-block:: python

    # 1. Check that mutant is called (covered by at least one read)
    if not GT[mutant_id].called:
        return False

    # 2. Check at least one of the unaffected sibling or other samples is called
    if sibling_id is not None:
        if (not GT[sibling_id].called) and (no other GTs are called):
        return False
    else:
        if (no other GTs are called):
        return False

    # 3. If both mutant and sibling are called and homozygous, their alleles must differ
    if sibling_id is not None and GT[sibling_id].called:
        if is_homozygous(GT[mutant_id]) and is_homozygous(GT[sibling_id]):
        if most_abundant_allele(GT[mutant_id]) == most_abundant_allele(GT[sibling_id]):
            return False

    # 4. The mutant's majority allele must not be the majority allele in any 'other' sample
    mutant_major_allele = most_abundant_allele(GT[mutant_id])
    for sample in samples:
        if sample not in [mutant_id, sibling_id] + exclude_ids:
        if GT[sample].called:
            if most_abundant_allele(GT[sample]) == mutant_major_allele:
            return False

    return True

    # Helper functions (bcftools plugin C API):
    # - is_homozygous: check if GT is hom (both alleles same and called)
    # - most_abundant_allele: get allele with max AD for sample
    # Use bcf_get_genotypes, bcf_get_format_int32, etc.

Strict candidate modular bcftools-native implementation (pseudocode):

.. code-block:: python

    # 1. Check if mutant is homozygous
    if (GT[mutant_id] is not called) or (not is_homozygous(GT[mutant_id])):
        return False

    # 2. Check all 'other' samples are called
    for sample in samples:
        if sample not in [mutant_id, unaffected_sibling_id] + exclude_ids:
        if GT[sample] is not called:
            return False

    # 3. Check mutant allele is not present in any 'other' sample
    mutant_allele = most_abundant_allele(GT[mutant_id])
    for sample in samples:
        if sample not in [mutant_id, unaffected_sibling_id] + exclude_ids:
        if mutant_allele in alleles_with_depth(GT[sample]):
            return False

    return True

    # Helper functions (bcftools plugin C API):
    # - is_homozygous: check if GT is hom (both alleles same and called)
    # - most_abundant_allele: get allele with max AD for sample
    # - alleles_with_depth: alleles with AD > 0
    # Use bcf_get_genotypes, bcf_get_format_int32, etc.    

Filter by read depth
++++++++++++++++++++

*Not intended to be carried out for this experiment*

Filter by snp impact
++++++++++++++++++++

    can just grep, dont need SnpSift

TEMP
++++

Quick, unrefined way to pull out strict candidates for reference only. Not to be considered as final results!

.. code-block:: shell

    outdir="../results/tmp/"
    pattern='#|./.:.:.:.:.\t./.:.:.:.:.\t./.:.:.:.:.\t./.:.:.:.:.$'
    metadata="../data/samplesheet_original.tsv"
    mkdir -p ${outdir}
    cut -f6 ${metadata} | tail -n +6 | \
        while IFS="\t" read -r line; do
            in=$(echo $line | cut -d ' ' -f3)
            out="${outdir}$(basename ${in} | cut -d '-' -f1-3).vcf"
            bgzip -cd ${in} | grep -P ${pattern} > ${out}
        done