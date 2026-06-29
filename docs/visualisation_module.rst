Visualisation module
====================

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

   <a href="https://hub.docker.com/r/tyronechen/pysam_pandas">
       <img src="https://img.shields.io/badge/Docker-tyronechen%pysam_pandas-2496ED?style=for-the-badge&logo=docker&logoColor=white" alt="Docker">
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

.. raw:: html

   Copyright © 2025 <a href="https://orcid.org/0000-0002-9207-0385">Tyrone Chen <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-6435-7100">Richard Lupat <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0009-0005-5595-3882">Michelle Meier <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-6528-891X">Maia Zethoven <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-6130-250X">Gregory Baillie <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0000-0000-0000">Scott Paterson <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0009-0001-5651-1331">Oguzhan Baltaci <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0003-3147-8042">Cas Simons <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-1150-3549">Jason Li <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0003-4189-9422">Andrew Cox <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-8283-9760">Kelly Smith <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0002-0651-7065">Benjamin Hogan <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>

Code in this repository is provided under a `MIT license`_. This documentation is provided under a `CC-BY-4.0 license`_.

.. _MIT license: https://opensource.org/licenses/MIT

.. _CC-BY-4.0 license: https://creativecommons.org/licenses/by/4.0/

`Visit our lab website here.`_ Contact Benjamin Hogan at `ben.hogan@petermac.org`_.

.. _Visit our lab website here.: https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/hogan-laboratory-vascular-developmental-genetics-and-cell-biology

.. _ben.hogan@petermac.org: mailto:ben.hogan@petermac.org

Overview
--------

This is an independent visualisation module distinct from the modules shown in the existing case studies. The key differences between the old version and this version are:

- Independent of pipeline, will run with any combination of files and genome annotations
- Advanced navigation tools including high-level variant maps
- Lazy-loading to handle massive files
- Sort and filter functions for regions of interest

.. warning::
    While theoretically possible to load entire genome ``bam`` and ``vcf`` files, this is not recommended due to practical data volume and memory constraints. At most, a chromosome of interest should be surveyed.

.. caution::
    The old visualisation module is functional but is no longer maintained. It requires a very specific setup. The new version is redesigned to be independent of the broader pipeline and more modular.

.. tip::
    Can also be set up on localhost to view.

Dependencies
------------

No major dependencies are required other than ``python`` and its base libraries. This should be pre-installed in most operating systems and available by default. ``igv.js`` will be retrieved automatically.

Usage
-----

.. hint:: 
    Most of the input data is optional and the interactive report will compile in their absence. Its possible to include only one piece of data of interest, while skipping others depending on circumstances.

Input directory structure
*************************

Here is a sample input directory structure::

    data
    ├── SOME_SAMPLE
    │   ├── chr12.vcf.gz
    │   ├── chr12.vcf.gz.csi
    │   ├── genome_overview
    │   │   ├── SOME_SAMPLE_MUT.nogap_nbases.bedgraph
    │   │   ├── SOME_SAMPLE_MUT.nogap_nsnps.bedgraph
    │   │   ├── SOME_SAMPLE_MUT.withgap_nbases.bedgraph
    │   │   ├── SOME_SAMPLE_MUT.withgap_nsnps.bedgraph
    │   │   ├── genome_0.png
    │   │   ├── genome_chr12.png
    │   │   └── SOME_SAMPLE_MUT.tdf
    │   ├── nucleotide_view
    │   │   ├── chr12_9404103-9406014.png
    │   │   ├── chr18_4727041-4727633.png
    │   │   ├── chr18_4730933_Dut.png
    │   │   ├── mut
    │   │   │   ├── SOME_SAMPLE_MUT.chr12.0-100Mb.bam
    │   │   │   └── SOME_SAMPLE_MUT.chr12.0-100Mb.bam.bai
    │   │   └── wt
    │   │       ├── SOME_SAMPLE_REF.chr12.0-100Mb.bam
    │   │       └── SOME_SAMPLE_REF.chr12.0-100Mb.bam.bai
    │   └── sample.md
    └── index.md

Overall metadata
++++++++++++++++

The root dataset directory can contain any number of subdirectories by sample name. Here one example ``SOME_SAMPLE`` is shown. An ``index.md`` file can be optionally provided. The index file has a flexible format allowing developers to add notes for users, while its main contents are derived automatically from parsing samples in each subdirectory.

Individual sample metadata
++++++++++++++++++++++++++

The special file ``sample.md`` can be included in the root of ``SOME_SAMPLE``::

    ---
    sample_id: some_id
    date: YYYY-MM-DD
    description: Lab X
    colour: green
    ---
    Some extra info

This will be rendered into a table of all samples on the index landing page.

Variant files
+++++++++++++

Any number of indexed, ``bgzipped``, ``vcf`` files can be added in the root of ``SOME_SAMPLE``. In this example we extract chromosome 12 only. 

.. caution::
    Ensure that chromosome names in all cases match annotation ``gff`` file.

Homozygosity maps
+++++++++++++++++

Place any number of ``bedgraph``, ``tdf`` and ``png`` files in ``genome_overview``. Image files will be rendered in collapsible blocks, and track files will be rendered in ``IGV``. 

These will appear above the single nucleotide resolution plots.

Single nucleotide resolution
++++++++++++++++++++++++++++

Place any number of ``png`` files in ``nucleotide_view`` which be rendered in image files. In ``mut`` and ``wt`` subdirectories, place indexed ``bam`` files which will be rendered in ``IGV``. Files in ``mut`` subdirectory will be rendered first in the hierarchy, followed by ``wt`` reference. In this example, we extract the whole chromosome 12.

These will appear below the homozygosity maps.

Genome file and annotations
+++++++++++++++++++++++++++

We use ``danRer7 Zv9`` in our example, but the user can insert any genome file and annotation corresponding to their study of interest.

.. warning::
    Chromosome names must match across data and annotations. A common error source is mismatches across different databases, eg ``1`` vs ``chr1``.

In this case we use files from NCBI and fix chromosomes manually (not shown) for our use case::

    wget 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.4_Zv9/GCF_000002035.4_Zv9_genomic.fna.gz' -o danRer7.fna.gz
    wget 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.4_Zv9/GCF_000002035.4_Zv9_genomic.gff.gz' -o danRer7.gff.gz

Files need to be indexed::

    samtools faidx danRer7.fna.gz
    tabix -p gff danRer7.gff.gz

If needed, sort the ``gff`` file before ``bgzip`` and index.
    
Usage
*****

Generating the files
++++++++++++++++++++

Once all the above is in place::

    python generate_report_json.py data -o report -g danRer7.gff.gz -f danRer7.fna.gz

You can then view the output in ``report``.

.. caution::
    The ``gff`` and ``fna`` files are symlinked into the output directory, in this case ``report``, to avoid large file duplication. If files are moved, the file paths will have to be updated accordingly for the web server to visualise genome annotations and tracks correctly.

Visualising the data
++++++++++++++++++++

Run the following, specifying the path to your output directory, in this case ``report``::

    python -m http.server 8000 --directory report

Open ``localhost:8000`` in a web browser of choice to view the interactive report.

.. caution::
    This is insecure and should be used with caution if hosted remotely on a web server. Standard cybersecurity measures should be present, but details are beyond the scope of this document.

Detailed information
++++++++++++++++++++

::
    usage: generate_report_json.py [-h] [-o OUTPUT_DIR] [-r MAX_VCF_ROWS] [-c NCPU] [-g GFF] [-f FASTA] [--skip-json] data

    Genomics report generator (JSON sidecar mode).

    positional arguments:
    data                  Root data directory

    options:
    -h, --help            show this help message and exit
    -o, --output-dir OUTPUT_DIR
                            Output directory (default: report/)
    -r, --max-vcf-rows MAX_VCF_ROWS
                            Max VCF rows per sample, 0=unlimited (default: 0)
    -c, --ncpu NCPU       Workers for parallel HTML rendering (default: 1)
    -g, --gff GFF         Sorted bgzipped+tabix-indexed GFF/GFF3 file (optional)
    -f, --fasta FASTA     Optional reference genome assembly FASTA (.fna/.fa)
    --skip-json           Skip writing JSON sidecars (re-use existing files for fast UI testing)

[OPTIONAL] Secure hosting
+++++++++++++++++++++++++

To set up a ``https`` secure host:

1. Obtain a free domain name at ``duckdns``
2. Apply ``certbot`` to generate a SSL certificate
3. Setup a ``nginx`` configuration.

.. note::
    This is highly user-specific, and will change depending on your level of access to the host server, file location and open port(s). For this reason, no specific instructions are provided. However, following the instructions for each of the above three services should work in most scenarios.
