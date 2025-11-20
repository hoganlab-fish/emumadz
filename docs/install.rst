Installation
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


Install
-------

Apptainer
+++++++++

First, `install apptainer`_.

.. _install apptainer: https://apptainer.org/get-started/

Run the provided `apptainer` image file.

.. code-block:: shell

    apptainer pull docker://tyronechen/emumadz:latest
    singularity pull docker://tyronechen/emumadz:latest

Conda
+++++

First, install `conda`, `mamba` or `micromamba`. You can find install instructions at `https://mamba.readthedocs.io/en/latest/`_.

.. _https://mamba.readthedocs.io/en/latest/: https://mamba.readthedocs.io/en/latest/

Install the provided `conda` environment.

*to be written.*

Manual
++++++

Base software
*************

The critical packages and their version numbers are listed below for reference::

    bcftools==1.19-gcc-13.2.0
    gatk/4.5.0.0-gcc-13.2.0
    pandas==2.3.0
    python==3.13.3
    pysam==0.16.0
    samtools==1.19.2-gcc-13.2.0
    snpeff==5.2
    ensembl-vep==109.3

Variant annotators
******************

.. caution::
    Make sure to install the correct genome assemblies for your use case so you have the right coordinates. Having the exact version number is a lower priority, since it affects genome annotations only.

VEP
###

For a manual install, `you can follow the install instructions on their github`_:, or use the ``conda`` / ``docker`` environments provided.

.. note::

    The ``conda`` environment is maintained by a third party unaffiliated with the authors. If that does not work, try building from source or using the official container instead.

.. code-block:: shell

    conda install bioconda::ensembl-vep

    apptainer pull --name vep.sif docker://ensemblorg/ensembl-vep
    singularity pull --name vep.sif docker://ensemblorg/ensembl-vep

.. _you can follow the install instructions on their github: https://github.com/Ensembl/ensembl-vep

Manually install the following libraries if you are running their ``INSTALL.pl`` script:

.. code-block:: shell

    perl==5.32.1
    perl-dbi==1.643
    perl-archive-zip==1.6.8
    perl-dbd-mysql==4.050
    perl-set-intervaltree==0.12
    perl-json==4.10
    perl-perlio-gzip==0.20
    perl-bio-bigfile==1.07
    perl-list-moreutils==0.430

Regardless of container or install method used, additional setup is required for the offline cache to work.
1. First, download the cache corresponding to the latest zebrafish ``Zv9`` assembly ``ENSEMBL 79``.
2. Unpack, rename and move the cache to its required location.
During use, you will need to provide specific options.

.. note::
  
    We want the regulatory regions also, so we get the merged tarball specifically.

.. caution::

    The cache directory defaults to ``~/.vep``, but suggest saving your downloaded cache files elsewhere as home directory usually has limited storage. You can then symlink the cache files to the directory. For example:

        .. code-block:: shell

            cd /place/with/storage/
            mkdir .vep
            mv my_cache_dir /place/with/storage/.vep/
            cd ~ && ln -s /place/with/storage/.vep
            

.. code-block:: shell

    # get cache
    wget 'https://ftp.ensembl.org/pub/release-79/variation/VEP/danio_rerio_merged_vep_79_Zv9.tar.gz'
    tar -xzvf danio_rerio_merged_vep_79_Zv9.tar.gz

    # if a cache dir isnt generated, this will default to ~/.vep
    # the directory must be renamed to drop the trailing _merged
    mv danio_rerio_merged ~/.vep/danio_rerio

    # then this should work, note that these specific options are required
    # [--cache --dir_cache --species --assembly --cache_version --offline]
    # this will be covered in detail in the relevant documentation section
    vep \
      --cache \
      --dir_cache ~/.vep/ \
      --species danio_rerio \
      --assembly Zv9 \
      --cache_version 79 \
      --offline \
      --regulatory \
      --vcf \
      ... \
      --input_file /path/to/input.vcf \
      --output_file /path/to/output.vcf
    

.. caution::
    There may be problems if there is conflict between VEP version, genome version and chromosome nomenclature. If you get any errors please check your files and `follow the guidelines on the ENSEMBL-VEP website`_.

.. _follow the guidelines on the ENSEMBL-VEP website: https://asia.ensembl.org/info/docs/tools/vep/script/vep_download.html

snpEff
######

Recommend installing through ``conda``.

.. code-block:: shell

    conda install bioconda::snpeff

.. warning::

    Manually installing this software is not recommended.

`If necessary, follow the instructions on the website`_.

.. _If necessary, follow the instructions on the website: https://pcingola.github.io/SnpEff/snpeff/introduction/

You will need to download the corresponding zebrafish genome database. If install is successful, run the code below.

.. code-block:: shell
  
    snpEff download Zv9.75

Visualisation module (OPTIONAL)
*******************************

.. caution::
    This part is more involved and is intended for developers who want to host a web server to view the data. This will also work on a local machine though.

Install ``npm`` `following the instructions for your own operating system`_. A few ``linux`` examples are provided.

.. _following the instructions for your own operating system: https://docs.npmjs.com/downloading-and-installing-node-js-and-npm

.. code-block:: shell

    # Ubuntu/Debian
    sudo apt install nodejs npm

    # CentOS/RHEL
    sudo yum install nodejs npm

.. hint::
    If you get an error saying the host name cannot be resolved, try the following (at your own risk).

    .. code-block:: shell

        sudo sed -i 's/mirrorlist/#mirrorlist/g' /etc/yum.repos.d/CentOS-*
        sudo sed -i 's|#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|g' /etc/yum.repos.d/CentOS-*


Then install ``igv-dist``.

.. code-block:: shell

    npm install express

    mkdir igv-dist
    curl -o igv-dist/igv.min.js https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js
    curl -o igv-dist/igv.css https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.css