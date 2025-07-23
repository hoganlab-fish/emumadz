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

For postprocessing
##################

``3.setup_dataset.sh`` copies ``bam`` alignment files and ``vcf`` variant calling files over from their specified locations in ``samplesheet_original.tsv`` to ``data/Alignment_File`` and ``data/VCF_Original`` respectively. We generate ``samplesheet_postprocess.tsv`` for use in this postprocessing analysis.

Whole zebrafish genome assembly
+++++++++++++++++++++++++++++++

*The ``seqliner`` assembly step is not covered in this documentation. The pipeline starts from the variant calling files*

Variant calling
+++++++++++++++

*The ``gatk`` variant calling step is not covered in this documentation. The pipeline starts from the variant calling files*

Postprocessing VCF
++++++++++++++++++

Rename chromosomes in vcf with map file
+++++++++++++++++++++++++++++++++++++++

``data/chrom_table.tsv`` contains the chromosome name mappings. Use this to remap the chromosomes.

For one example:

.. code-block:: shell

    bcftools annotate --rename-chr "../data/chrom_table.tsv" \
        ../data/VCF_Original/TL2312073-163-4L.vcf.gz \
        -Ob -o ../results/VCF_ChrFixed/TL2312073-163-4L.vcf.gz
    bcftools index -t ../results/VCF_ChrFixed/TL2312073-163-4L.vcf.gz

For an example on all samples:

.. code-block:: shell

    metadata="../data/samplesheet_postprocess.tsv"
    map="../data/chrom_table.tsv"
    tail -n +2 ${metadata} | \
        while IFS="\t" read -r line; do
            in=$(echo $line | cut -d ' ' -f3)
            out="../results/VCF_ChrFixed/"$(basename ${in})
            bcftools annotate --rename-chrs ${map} \
                ${in} -Ob -o ${out}
            bcftools index -t ${out}
        done


Merge files with refs
+++++++++++++++++++++

For each sample, merge the four references:
- TL2209397-ENUref-Female.vcf.gz
- TL2209398-ENUref-Male.vcf.gz
- TL2209399-TUref-Female.vcf.gz
- TL2209400-TUref-Male.vcf.gz

.. code-block:: shell

    bcftools merge


Filter out non-snps
+++++++++++++++++++

bcftools

Get strict candidates
+++++++++++++++++++++

.. code-block:: shell

    grep -P './.:.:.:.:.\t./.:.:.:.:.\t./.:.:.:.:.\t./.:.:.:.:.$'


Filter by read depth
++++++++++++++++++++

    awk DP with threshold filter

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