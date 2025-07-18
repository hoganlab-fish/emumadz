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

High-level automated approach
+++++++++++++++++++++++++++++

`A nextflow pipeline developed by the Peter MacCallum Cancer Centre Bioinformatics Core is also available.`_ This ingests raw `fastq` files as input, performs alignments to generate `bam` files, and performs variant calling to generate `vcf` files. Note that the paths are hardcoded at time of writing. Using this pipeline incorporates all steps below.

.. _A nextflow pipeline developed by the Peter MacCallum Cancer Centre Bioinformatics Core is also available.: https://github.com/PMCC-BioinformaticsCore/zebrafish_postprocess/tree/main

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

For the purposes of rerunning this experiment, we only want the ``bam`` alignment files since we will be recreating everything starting from the variant calling step. Running ``3.setup_dataset.sh`` will create the corresponding input and output directories::

   ../data/Alignment_File
   ../results/Sample_Identity
   ../results/Fastq_File
   ../results/Alignment_File
   ../results/VCF_Original
   ../results/VCF_Merged
   ../results/VCF_ChrFixed
   ../results/VCF_Annotated
   ../results/VCF_Candidates
   ../results/Snzl_NoGaps_NBases
   ../results/Snzl_NoGaps_NSnps
   ../results/Snzl_WithGaps_NBases
   ../results/Snzl_WithGaps_NSnps
   ../results/Json
   
While the ``bam`` files are not directly included in this repository due to size constraints, ``md5`` sums are preserved in ``data/Alignment_File/bam.md5``.

Four samplesheet files are generated:
- ``samplesheet.220701_A01221_0125_BHCCHNDMXY.tsv``
- ``samplesheet.220930_A00692_316_AHVMJFDSX3.tsv``
- ``samplesheet.221014_A00692_0319_BHGJ7CDMXY.tsv``
- ``samplesheet.231201_A00692_0394_231208_A00692_0396.tsv``

The ``samplesheet.220701_A01221_0125_BHCCHNDMXY.tsv`` contains ``F0`` grandparent reference generation metadata. All other samplesheets contain ``F2`` generation metadata. ``samplesheet_original.tsv`` contains these legacy file paths.

This dataset
############

``3.setup_dataset.sh`` copies ``bam`` alignment files over from their specified locations in ``samplesheet_original.tsv`` to ``data/Alignment_File``. A new ``samplesheet_new.tsv`` is created for the latest analysis.

Whole zebrafish genome assembly
+++++++++++++++++++++++++++++++

*The ``seqliner`` assembly step is not covered in this documentation. The pipeline starts from the aligned bam files*

Variant calling
+++++++++++++++

NCBI reference genome
#####################

`Download the reference genome here.`_ Note that the chromosome names may not be the same as the newly assembled ones.

.. _Download the reference genome here.: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000002035.2/

F0 genomes
##########

There are 4 F0 references.
- ``TL2209397-ENUref-Female-6-MAN-20220655``
- ``TL2209398-ENUref-Male-6-MAN-20220656``
- ``TL2209399-TUref-Female-4-MAN-20220657``
- ``TL2209400-TUref-Male-4-MAN-20220658``

F2 genomes
##########

There are 44 F2 samples.

Call references against the NCBI reference
##########################################

The ``vcf`` files are generated by calling variants from the:
- ``NCBI reference genome`` against the ``F0 genomes``
- ``NCBI reference genome`` against the ``F2 genomes``

We want to find SNPs in the F2 generation which are not in the F1 generation.

F0 pooling
##########

Perform the variant calling for the F0, pooling samples.

.. attention::
   In the first iteration, samples were aggregated only but not pooled during variant calling.

Run the ``4.run_gatk.sh`` script and follow the instructions.

.. caution::
   The script is not designed to be run in one go. Each step submits a series of slurm jobs, which generates the files used in the next stage of the pipeline.

.. attention::
   The original iteration of the pipeline used an older variant caller ``UnifiedGenotyper``. This method is now obsolete and documentation is sparse. We use the updated ``HaplotypeCaller``, which is functionally similar and has a higher performance. `We follow the pipeline described here.`_

.. _We follow the pipeline described here.: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels


3. Generate homozygosity peaks for visualisation
++++++++++++++++++++++++++++++++++++++++++++++++


4. Quantify SNP impact
++++++++++++++++++++++


5. Identify strict candidate SNPs
+++++++++++++++++++++++++++++++++


6. Filter out predicted high impact SNPs
++++++++++++++++++++++++++++++++++++++++


7. Filter out only SNPs (excluding CNV, INDELS, etc)
++++++++++++++++++++++++++++++++++++++++++++++++++++

.. raw:: html

   <details>
   <summary><a>Example hyperparameter.json file</a></summary>

.. code-block:: shell

    echo TEST

.. raw:: html

   </details>