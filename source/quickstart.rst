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

0. High-level automated approach
++++++++++++++++++++++++++++++++

`A nextflow pipeline developed by the Peter MacCallum Cancer Centre Bioinformatics Core is also available.`_ This ingests raw `fastq` files as input, performs alignments to generate `bam` files, and performs variant calling to generate `vcf` files. Note that the paths are hardcoded at time of writing. Using this pipeline incorporates all steps below.

.. _https://github.com/PMCC-BioinformaticsCore/zebrafish_postprocess/tree/main: https://github.com/PMCC-BioinformaticsCore/zebrafish_postprocess/tree/main

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