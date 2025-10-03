---
title: 'EMUMADZ: Enhanced MUtation MApping and Detection in Zebrafish'
tags:
  - nextflow
  - forward genetics screen
  - variant calling
  - single nucleotide polymorphism
  - variant annotation
  - genome browser visualisation
authors:
  - name: Tyrone Chen
    orcid: 0000-0002-9207-0385
    affiliation: "1, 2"
  - name: Richard Lupat
    orcid: 0000-0002-6435-7100
    affiliation: 1
  - name: Michelle Meier
    orcid: 0009-0005-5595-3882
    affiliation: "1, 2"
  - name: Maia Zethoven
    orcid: 0000-0002-6528-891X
    affiliation: "1, 2"
  - name: Gregory Baillie
    orcid: 0000-0002-6130-250X
    affiliation: 3
  - name: Scott Paterson
    affiliation: 2
  - name: Oguzhan Baltaci
    orcid: 0009-0001-5651-1331
    affiliation: 2
  - name: Cas Simons
    orcid: 0000-0003-3147-8042
    affiliation: "1, 2"
  - name: Jason Li
    orcid: 0000-0002-1150-3549
    affiliation: 1
  - name: Benjamin Hogan
    orcid: 0000-0002-0651-7065
    affiliation: 2
    corresponding: true

affiliations:
 - name: Bioinformatics Core, Peter MacCallum Cancer Centre, Australia
   index: 1
 - name: Organogenesis and Cancer Program, Peter MacCallum Cancer Centre, Australia
   index: 2
 - name: Institute for Molecular Bioscience, University of Queensland, Australia
   index: 3
date: 1 November 2025
bibliography: paper.bib
---

# Summary

`EMUMADZ` is a species-agnostic tool incorporating multiple command line tools for forward genetics screening. We optimised a pipeline that (a) mapped mutations to a region of genomic linkage and (b) identified candidate Single Nucleotide Polymorphisms (SNPs) predicted to be damaging to gene function, matching (c) essential criteria of homozygous mutations which are _N-ethyl-N-nitrosourea_(ENU)-induced and absent in control reference samples. `Nextflow` was used to improve the scalability and portability of the pipeline, along with other quality of life features such as process checkpointing [@di2017nextflow]. At the same time, the high-level API is designed to provide easy usage for individuals running the pipeline. Code is designed in a modular way, allowing users to reuse or modify specific processes to suit their own needs. Extensive documentation exists, with two case studies on _Danio rerio_: the zebrafish. 

The primary motivation for devleoping this pipeline stems from the fact that forward genetic screening in model organisms remains a powerful tool for identifying unexpected new genes involved in any biological process of interest. Characterising mutants and their candidate causative mutations, often in the form of SNPs, frequently lead to new insights on gene function in previously unknown genes or other genomic regions of interest. 

# Statement of need

Individual software tools have been designed for forward genetic screens, but most tools and higher-level pipelines have (a) an organism bias towards humans and mice and (b) are designed for addressing specific aspects of SNPs instead of the SNP characterisation process as a whole.

Therefore, `EMUMADZ` incorporates a multi-step process encompassing a complete forward genetic screening pipeline for SNPs:

1. Variant calling
2. SNP identification
3. Filtering strong candidate SNPs
4. Annotating SNP impacts
5. Visualising data as tracks on genome browsers

Additional important features are also performed in the background as needed:

1. Chromosome name alignment
2. Decomposition of multiallelic variants
3. Merging of files to accommodate multi-sample comparisons

`EMUMADZ` is designed to be used by bioinformaticians interested in analysing forward genetics screening data, for example with the chemical mutagen ENU to discover new autosomal recessive mutations in biological processes of interest. It incorporates state of the art tools such as `gatk` [@mckenna2010genome], `bcftools` [@danecek2021twelve], `samtools` [@danecek2021twelve], `snpEff` [@cingolani2012variant] and `ensembl-vep` [@mclaren2016ensembl], all which will be familiar to those in the field. However, we note that no prior knowledge of these tools are required. Case studies demonstrating a run of the pipeline exist. These were intentionally written to be informative to both biologists and bioinformaticians, where individual steps are explained in their biological context with their corresponding code block.

A large scale forward genetic screen was performed using the chemical mutagen _N-ethyl-N-nitrosourea_ (ENU) to discover new autosomal recessive mutations in biological processes of interest. To identify candidate causative mutations for screened mutants, a low coverage (~10x) whole genome sequencing (WGS) was performed. To map homozygous mutants to the genome for subsequent identification of candidate mutations, WGS data was generated from mutant animals and control samples. We optimised a pipeline that (a) mapped mutations to a region of genomic linkage and (b) identified candidate SNPs predicted to be damaging to gene function, matching (c) essential criteria of homozygous mutations which are ENU-induced and absent in control reference samples.

<!-- # Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% } -->

# Acknowledgements

We thank Tyrone Chen, Richard Lupat, Song Li and Jason Li from the Bioinformatics Core (`RRID: SCR_025901`) at the Peter MacCallum Cancer Centre for their helpful comments and contributions to the pipeline.

# References
