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

.. _A nextflow pipeline developed by the Peter MacCallum Cancer Centre Bioinformatics Core is also available.: https://github.com/PMCC-BioinformaticsCore/zebrafish_postprocess/tree/main: https://github.com/PMCC-BioinformaticsCore/zebrafish_postprocess/tree/main

0. Data setup
+++++++++++++

.. caution::
   The original data setup was sample-centric as opposed to more conventional process-centric directory structure. In addition, code had hardcoded paths which limited file movement. This impacted all steps of analysis and reduced reproducibility. To lower the impact of the original file layout, the data directory now contains samplesheets with updated file paths to symlinks and their attributes. However, it is possible that some errors remain.

The samplesheets are tab separated files with the following fields:

.. Sample_Identity   Sample_Identity_Intermediate   Fastq_File  Alignment_File VCF_Original  VCF_Merged  VCF_ChrFixed  VCF_Annotated   VCF_Candidates Snzl_Gaps_NBases   Snzl_NoGaps_NBases  Snzl_Gaps_NSnps   Snzl_NoGaps_NSnps Json

.. csv-table:: Samplesheet column information
   :file: tables/samplesheet_fields.csv
   :header-rows: 1

.. warning::
   ``Sample_Identity_Intermediate`` is required. A downstream step silently truncates file names if it exceeds a certain limit. No explanation for this behaviour was available.

.. note::
   ``Vis_`` fields should contain the file name but not the file extensions, both ``bedgraph`` and ``tdf`` files will be generated. Note that both files contain the same information but ``tdf`` is optimised for viewing with the ``IGV`` Genome Browser.

Run the ``validate_samples.py`` script if unsure. This checks the validity of each file in the samplesheet per sample. *To be written*

1. Whole zebrafish genome assembly
++++++++++++++++++++++++++++++++++


2. Variant calling
++++++++++++++++++


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