Installation
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

Base software
*************

The critical packages and their version numbers are listed below for reference::

    bcftools==1.19-gcc-13.2.0
    gatk/4.5.0.0-gcc-13.2.0
    pandas==2.3.0
    python==3.13.3
    pyvcf==0.6.8
    samtools==1.19.2-gcc-13.2.0

    snpeff==5.2
    ensembl-vep==114.2




Variant annotators
******************

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
