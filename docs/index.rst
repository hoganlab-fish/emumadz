EMUMADZ
=======

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

    <a href="https://doi.org/10.5281/zenodo.17810567">
        <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.17810567.svg?style=for-the-badge" alt="Zenodo">
    </a>

Enhanced MUtation MApping and Detection in Zebrafish
----------------------------------------------------
.. whole_genome_sequencing documentation master file, created by
   sphinx-quickstart on Mon Jun 30 12:23:50 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. raw:: html

   Copyright © 2025 <a href="https://orcid.org/0000-0002-9207-0385">Tyrone Chen <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-6435-7100">Richard Lupat <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0009-0005-5595-3882">Michelle Meier <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-6528-891X">Maia Zethoven <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-6130-250X">Gregory Baillie <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0000-0000-0000">Scott Paterson <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0009-0001-5651-1331">Oguzhan Baltaci <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0003-3147-8042">Cas Simons <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-1150-3549">Jason Li <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0003-4189-9422">Andrew Cox <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-8283-9760">Kelly Smith <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-0651-7065">Benjamin Hogan <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>

Code in this repository is provided under a `MIT license`_. This documentation is provided under a `CC-BY-4.0 license`_.

.. _MIT license: https://opensource.org/licenses/MIT

.. _CC-BY-4.0 license: https://creativecommons.org/licenses/by/4.0/

`Visit our lab website here.`_ Contact Benjamin Hogan at `ben.hogan@petermac.org`_.

.. _Visit our lab website here.: https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/hogan-laboratory-vascular-developmental-genetics-and-cell-biology

.. _ben.hogan@petermac.org: mailto:ben.hogan@petermac.org

.. raw:: html

    We thank <a href="https://orcid.org/0000-0002-9207-0385">Tyrone Chen <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-6435-7100">Richard Lupat <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0009-0005-5595-3882">Michelle Meier <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-6528-891X">Maia Zethoven <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0000-0000-0000">Song Li <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-3089-0871">Niko Thio <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, and <a href="https://orcid.org/0000-0002-1150-3549">Jason Li <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a> from the Bioinformatics Core (RRID: `SCR_025901`) at the Peter MacCallum Cancer Centre for their helpful comments and contributions to the pipeline.

Contains the documentation for our custom zebrafish whole genome sequencing pipeline.

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   intro
   install
   case_study_1
   case_study_2
   case_study_3
   api

Appendix
========

.. toctree::
   :maxdepth: 1
   :caption: Appendix:
   
   appendix_internal_ref_only

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
