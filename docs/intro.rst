Introduction
============

.. raw:: html

   <a href="https://opensource.org/licenses/MIT">
       <img src="https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge" alt="MIT License">
   </a>

   <a href="https://creativecommons.org/licenses/by/4.0/">
       <img src="https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg?style=for-the-badge" alt="CC BY 4.0">
   </a>

   <a href="https://github.com/hoganlab-fish/emumadz">
       <img src="https://img.shields.io/badge/GitHub-hoganlab--fish%2Femumadz-181717?style=for-the-badge&logo=github&logoColor=white" alt="GitHub">
   </a>

   <a href="https://github.com/hoganlab-fish/emumadz/actions/workflows/docker.yml">
       <img src="https://img.shields.io/github/actions/workflow/status/hoganlab-fish/emumadz/docker.yml?style=for-the-badge&logo=github-actions&logoColor=white" alt="Build Status">
   </a>

   <a href="https://hub.docker.com/r/tyronechen/pysam_pandas">
       <img src="https://img.shields.io/docker/pulls/tyronechen/pysam_pandas?style=for-the-badge&logo=docker&logoColor=white" alt="Docker Pulls">
   </a>

   <a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/hogan-laboratory-vascular-developmental-genetics-and-cell-biology">
       <img src="https://img.shields.io/badge/Website-hoganlab-4285F4?style=for-the-badge&logo=google-chrome&logoColor=white" alt="Website">
   </a>

.. raw:: html

    Copyright © 2025 <a href="https://orcid.org/0000-0002-9207-0385">Tyrone Chen <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-6435-7100">Richard Lupat <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0009-0005-5595-3882">Michelle Meier <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-6528-891X">Maia Zethoven <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-6130-250X">Gregory Baillie <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0000-0000-0000">Scott Paterson <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0009-0001-5651-1331">Oguzhan Baltaci <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0003-3147-8042">Cas Simons <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-1150-3549">Jason Li <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0003-4189-9422">Andrew Cox <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-8283-9760">Kelly Smith <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-0651-7065">Benjamin Hogan<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>

Code in this repository is provided under a `MIT license`_. This documentation is provided under a `CC-BY-4.0 license`_.

.. _MIT license: https://opensource.org/licenses/MIT

.. _CC-BY-4.0 license: https://creativecommons.org/licenses/by/4.0/

`Visit our lab website here.`_ Contact Benjamin Hogan at `ben.hogan@petermac.org`_.

.. _Visit our lab website here.: https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/hogan-laboratory-vascular-developmental-genetics-and-cell-biology

.. _ben.hogan@petermac.org: mailto:ben.hogan@petermac.org


Highlights
----------

We provide a pipeline with the following features:

- Takes fastq files, samplesheet as input
- Returns candidate SNPs for further investigation
- Interactive visualisation is a component
- Detailed documentation on multiple case studies


Graphical abstract
------------------

For more information, `visit our poster here`_:

.. _visit our poster here: https://f1000research.com/posters/14-1261

Abstract
--------

Forward genetic screening in model organisms remains a powerful tool for identifying unexpected new genes involved in any biological process of interest. While biological and computational protocols for identifying genetic variation are mature, many workflows are designed for human and mouse data. Therefore, we developed a species-agnostic pipeline which compares any combination of mutant and reference samples. We demonstrate the capabilities of our pipeline with a case study using the common laboratory model organism Danio rerio: the zebrafish. A large scale forward genetic screen was performed using the chemical mutagen N-ethyl-N-nitrosourea (ENU) to discover new autosomal recessive mutations in biological processes of interest. To identify candidate causative mutations for screened mutants, a low coverage (~10x) whole genome sequencing (WGS) was performed. To map homozygous mutants to the genome for subsequent identification of candidate mutations, WGS data was generated from mutant animals and control samples. We optimised a pipeline that (a) mapped mutations to a region of genomic linkage and (b) identified candidate SNPs predicted to be damaging to gene function, matching (c) essential criteria of homozygous mutations which are ENU-induced and absent in control reference samples. Other features of our pipeline include modern, interactive visualisations. Future plans include AlphaFold integration as an additional predictor for SNP impact as well as improved scalability with nextflow. Our pipeline is documented with a case study for reference and is hosted in an open-source software repository. The tool has minimal dependencies, with apptainer and conda instances available for usability. 

Cite us with
------------

*Manuscript and zenodo citation to follow*
