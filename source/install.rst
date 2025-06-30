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


Install
-------

Apptainer
+++++++++

First, install `apptainer`.

Run the provided `apptainer` image file.

*eventual apptainer image setup*

Conda
+++++

First, install `conda`, `mamba` or `micromamba`. You can find install instructions at `https://mamba.readthedocs.io/en/latest/`_.

.. _https://mamba.readthedocs.io/en/latest/: https://mamba.readthedocs.io/en/latest/

Install the provided `conda` environment.

*eventual conda package (optional, priority is above)*

Manual
++++++

The critical packages and their version numbers are listed below for reference::

    bcftools==1.19-gcc-13.2.0
    samtools==1.19.2-gcc-13.2.0
    nextflow==25.04.2.5947
