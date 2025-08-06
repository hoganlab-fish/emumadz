Introduction
============

.. raw:: html

  Copyright (c) 2025 <a href="https://orcid.org/0000-0002-9207-0385">Tyrone Chen, <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a><a href="https://orcid.org/0000-0003-3746-0695">Oliver Yu, <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a><a href="https://orcid.org/0000-0002-0651-7065">Benjamin Hogan <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>.

Code in this repository is provided under a `MIT license`_. This documentation is provided under a `CC-BY-3.0 AU license`_.

.. _MIT license: https://opensource.org/licenses/MIT

.. _CC-BY-3.0 AU license: https://creativecommons.org/licenses/by/3.0/au/

`Visit our lab website here.`_ Contact Benjamin Hogan at `ben.hogan@petermac.org`_.

.. _Visit our lab website here.: https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/hogan-laboratory-vascular-developmental-genetics-and-cell-biology

.. _ben.hogan@petermac.org: mailto:ben.hogan@petermac.org


Highlights
----------

We provide a pipeline which:
1. Performs whole zebrafish genome assembly
2. Variant calling
3. Generates homozygosity peaks for visualisation
4. Predicts the impact of identified SNPs
5. Identify strict candidate SNPs
6. Filter out predicted high impact SNPs
7. Filter out only SNPs (excluding CNV, INDELS, etc)

*put a workflow figure here*

Changelog
---------
- Uses updated versions of ``gatk``, ``bcftools``, ``snpEff`` to replace obsolete functionality.
- Uses ``csi`` instead of ``tbi`` indexing for stability with long chromosomes.
- Reworked logic of the original pipeline to handle missing values.
- Increased threshold for read mapping quality which was causing false negatives.
- Lowered threshold for homozygous calls to handle leakage which was causing false negatives.
- All-in-one pipeline which can handle any combination of mutant and control data.
- Now handles mulitallelic SNPs correctly (i.e. if a reference has a SNP but this SNP is different to the mutant, it is correctly identified as a SNP event in the mutant.)
- Concept of candidate strictness is removed and replaced with allele frequencies with user-defined filters.
- Modern visualisation methods independent of ``IGV`` genome browser.

Cite us with
------------

*eventual manuscript and/or zenodo citation*
