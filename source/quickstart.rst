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

.. _A nextflow pipeline developed by the Peter MacCallum Cancer Centre Bioinformatics Core is also available.: https://github.com/PMCC-BioinformaticsCore/zebrafish_postprocess/tree/main

Data setup
++++++++++

Original data
#############

.. caution::
   The original data setup was sample-centric as opposed to more conventional process-centric directory structure. In addition, code had hardcoded paths which limited file movement. This impacted all steps of analysis and reduced reproducibility. To lower the impact of the original file layout, the data directory now contains samplesheets with updated file paths to symlinks and their attributes. However, it is possible that some errors remain.

Run ``1.extract_filepaths.sh`` with the required command line arguments to generate the samplesheets with paths to old files. The samplesheets are tab separated files with the following fields:

.. csv-table:: Samplesheet column information
   :file: tables/samplesheet_fields.csv
   :header-rows: 1

.. warning::
   A downstream step silently truncates file names if it exceeds a certain limit. No explanation for this behaviour was available. The ``1.extract_filepaths.sh`` script will try to "guess" file paths by performing an arbitrary truncation and subsequent substring matching.

.. note::
   ``Vis_`` fields should contain the file name but not the file extensions, both ``bedgraph`` and ``tdf`` files will be generated. Note that both files contain the same information but ``tdf`` is optimised for viewing with the ``IGV`` Genome Browser.

.. csv-table:: Example samplesheet showing one sample
   :file: tables/samplesheet_example.csv
   :header-rows: 1

The ``2.validate_samples.py`` script checks the validity of each file in the samplesheet per sample. This also drops the ``fastq`` column. Automatically runs in ``1.extract_filepaths.sh``.

.. code:: shell

   python 2.validate_samples.py ../data/samplesheet.tmp -o ../data/samplesheet.tsv

.. caution::
   Custom directories and/or files exist since the data passed through multiple iterations. To some extent this is accommodated in the setup and validation scripts, but make sure to double-check everything.

For the purposes of rerunning this experiment, we only want the ``bam`` alignment files since we will be recreating everything starting from the variant calling step. Running ``3.setup_dataset.sh`` will create the corresponding input and output directories::

   ../data/Alignment_File
   ../results/Sample_Identity
   ../results/Fastq_File
   ../results/Alignment_File
   ../results/VCF_Original
   ../results/VCF_Merged
   ../results/VCF_ChrFixed
   ../results/VCF_Annotated
   ../results/VCF_Candidates
   ../results/Snzl_NoGaps_NBases
   ../results/Snzl_NoGaps_NSnps
   ../results/Snzl_WithGaps_NBases
   ../results/Snzl_WithGaps_NSnps
   ../results/Json
   
While the ``bam`` files are not directly included in this repository due to size constraints, ``md5`` sums are preserved in ``data/Alignment_File/bam.md5``.

Four samplesheet files are generated:
- ``samplesheet.220701_A01221_0125_BHCCHNDMXY.tsv``
- ``samplesheet.220930_A00692_316_AHVMJFDSX3.tsv``
- ``samplesheet.221014_A00692_0319_BHGJ7CDMXY.tsv``
- ``samplesheet.231201_A00692_0394_231208_A00692_0396.tsv``

The ``samplesheet.220701_A01221_0125_BHCCHNDMXY.tsv`` contains ``F0`` grandparent reference generation metadata. All other samplesheets contain ``F2`` generation metadata. ``samplesheet_original.tsv`` contains these legacy file paths.

This dataset
############

``3.setup_dataset.sh`` copies ``bam`` alignment files over from their specified locations in ``samplesheet_original.tsv`` to ``data/Alignment_File``. A new ``samplesheet_new.tsv`` is created for the latest analysis.

Whole zebrafish genome assembly
+++++++++++++++++++++++++++++++

*The ``seqliner`` assembly step is not covered in this documentation. The pipeline starts from the aligned bam files*

Variant calling
+++++++++++++++

NCBI reference genome
#####################

`Download the reference genome here.`_ Note that the chromosome names may not be the same as the newly assembled ones.

.. _Download the reference genome here.: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000002035.2/

F0 genomes
##########

There are 4 F0 references.
- ``TL2209397-ENUref-Female-6-MAN-20220655``
- ``TL2209398-ENUref-Male-6-MAN-20220656``
- ``TL2209399-TUref-Female-4-MAN-20220657``
- ``TL2209400-TUref-Male-4-MAN-20220658``

F2 genomes
##########

There are 44 F2 samples.

Call references against the NCBI reference
##########################################

The ``vcf`` files are generated by calling variants from the:
- ``NCBI reference genome`` against the ``F0 genomes``
- ``NCBI reference genome`` against the ``F2 genomes``

We want to find SNPs in the F2 generation which are not in the F1 generation.

F0 pooling
##########

Perform the variant calling for the F0, pooling samples.

.. attention::
   In the first iteration, samples were aggregated only but not pooled during variant calling.

Run the ``4.run_gatk.sh`` script and follow the instructions.

.. caution::
   The script is not designed to be run in one go. Each step submits a series of slurm jobs, which generates the files used in the next stage of the pipeline.

.. attention::
   The original iteration of the pipeline used an older variant caller ``UnifiedGenotyper``. This method is now obsolete and documentation is sparse. We use the updated ``HaplotypeCaller``, which is functionally similar and has a higher performance. `We follow the pipeline described here.`_

.. _We follow the pipeline described here.: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels


3. Generate homozygosity peaks for visualisation
++++++++++++++++++++++++++++++++++++++++++++++++


4. Quantify SNP impact
++++++++++++++++++++++


5. Identify strict candidate SNPs
+++++++++++++++++++++++++++++++++

Description of claims
#####################

The ``snzl`` module claims to identify strict candidate SNPs using the following criteria:
1. The mutant must be covered to a depth of \>= 1 read
2. At least one of the unaffected sibling or other samples must be covered to a depth of \>= 1 read
3. If the mutant and sibling have a single allele, it must not be the same in both
4. The mutant must have a majority allele that is:
   a. not the majority allele in any of the 'other' samples


.. raw:: html

   <details>
   <summary><a>Unit tests of ``snzl``</a></summary>


Setup test environment
######################

.. caution::
   These steps are specific to the Peter Mac compute cluster

Load singularity container. You can either load it manually::
   
   apptainer shell --bind \
      /team_folders/hogan_lab/Hogan_Lab_Shared/Zebrafish_Homozygosity/,/etc/profile.d/ --home /hogan_lab/Hogan_Lab_People/Tyrone_Chen/repos/hogan-lab-shared-repository/wgs_vbm/ /config/spack/containers/centos7/container.sif
   source /etc/profile.d/modules.sh
   source /team_folders/hogan_lab/Hogan_Lab_Shared/Zebrafish_Homozygosity/gjbzebrafishtools/venv/bin/activate
   module load tabix
   module load bcftools

Or specify this as an option in the ``sbatch`` script::

   #SBATCH --container /config/spack/containers/centos7/container.sif
   source /etc/profile.d/modules.sh
   export LD_LIBRARY_PATH=/config/binaries/gsl/2.7.1/lib:/config/binaries/gcc/12.2.0/lib64:/config/binaries/R/4.2.0.Core/lib64/R/lib:/config/binaries/python/3.8.1/lib:/config/binaries/hdf5/1.10.5/lib
   source /team_folders/hogan_lab/Hogan_Lab_Shared/Zebrafish_Homozygosity/gjbzebrafishtools/venv/bin/activate

``sinteractive`` option::

   sinteractive -p rhel_short --time 0-08:00:00 --mem 64G --container /config/spack/containers/centos7/container.sif
   source /etc/profile.d/modules.sh
   source /team_folders/hogan_lab/Hogan_Lab_Shared/Zebrafish_Homozygosity/gjbzebrafishtools/venv/bin/activate

.. warning::
   Library ``snzl`` only works on CentOS7.

.. warning::
   Library ``snzl`` no longer exists in public.

Test design
###########

Check the original assumptions by manually editing the VCF, varying SNP and DP, and see if ``snzl`` picks it up.

To test these claims, in each case a ``vcf`` file of a target sample was subsampled to generate a small test sample. Most headers were removed.

All of the steps below are carried out in the ``test`` directory.

Controls (SNP quantity)
***********************

.. note::
    ``./.:.:.:.:.`` indicates a non-SNP event in the corresponding sample.

Five fields with varying numbers of SNP across the sample and four references::

   cd test
   # a file from the original sample set
   head -n1172 TL2312073-163-4L-MAN-20231_Ref_merged_ChromFixed.vcf.annotated.vcf > snps.1.4.vcf
   
   # take 5 regions where the sample and all four references have SNPs
   grep -v '##' TL2312073-163-4L-MAN-20231_Ref_merged_ChromFixed.vcf.annotated.vcf | head -n20001 | grep -v './.:.:.:.' >> snps.1.4.vcf

   # bgzip    
   bgzip snps.1.4.vcf
   bcftools index -t snps.1.4.vcf.gz

There are 8 SNPs in the vcf file. Each SNP has SNPs removed from one sample to end up with the following design, where the first column is the sample and all trailing columns are the references::

   1 1 1 1 1
   1 1 1 1 1
   0 1 1 1 1
   1 1 1 1 0
   1 1 1 0 0
   1 1 0 0 0
   1 0 0 0 0
   0 0 0 0 0

DP values are also varied.

Controls (DP quantity)
**********************

.. note::
    ``DP`` indicates read depth. Discussion on an ideal read depth is outside the scope of this document. 10-30 is generally considered OK.

Five entries with varying levels of *read depth* {0,10,20,30,40}. In this sample, SNPs are not detected in any references, making this a strict candidate. SNP effect was predicted to be ``MODERATE``, matching filter threshold. For the purposes of this test, SNP positions are artificial.

.. code-block:: shell

   head -n1172 TL2312073-163-4L-MAN-20231_Ref_merged_ChromFixed.vcf.annotated.vcf > header.vcf
   (cat header.vcf; for i in {1..5}; do echo $(grep chr24 TL2312073-163-4L-MAN-20231_Ref_merged_ChromFixed.vcf.annotated.vcf | grep -m 1 423800); done) > test_dp.tmp
   tac test_dp.tmp | sed -e "1,5s| |\t|g" -e "1s/DP=10/DP=40/" -e "1s/423800/5/" -e "2s/DP=10/DP=30/" -e "2s/423800/4/" -e "3s/DP=10/DP=20/" -e "3s/423800/3/" -e "4s/423800/2/" -e "5s/423800/1/" -e "5s/DP=10/DP=0/" | tac > test_dp.vcf
   rm test_dp.tmp
   bgzip test_dp.vcf
   bcftools index -t test_dp.vcf.gz

Here is an example of one entry in the file. Only the ``DP`` value is modified in each iteration::

   chr24   5       .       G       T       269.98  .       BaseQRankSum=0;Dels=0;ExcessHet=3.0103;FS=0;HaplotypeScore=0.9997;MQ=60;MQ0=0;MQRankSum=0;QD=27;ReadPosRankSum=1.036;SOR=0.307;DP=40;AF=0.5;MLEAC=1;MLEAF=0.5;AN=2;AC=1;EFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Gtt/Ttt|V7F|128|RECK|protein_coding|CODING|ENSDART00000129135|1|T|WARNING_TRANSCRIPT_NO_START_CODON),DOWNSTREAM(MODIFIER||2776||157|INSC|protein_coding|CODING|ENSDART00000131091||T|WARNING_TRANSCRIPT_NO_STOP_CODON)      GT:AD:DP:GQ:PL  0/1:1,9:10:13:298,0,13  ./.:.:.:.:.    ./.:.:.:.:.      ./.:.:.:.:.     ./.:.:.:.:.

The ``snzl`` strict candidate pipeline is run.

.. code:: shell

   infile_path="test_dp.vcf.gz"
   outfile_path="test_dp_out.vcf"
   sample_name="TL2312073-163-4L-MAN-20231116"
   chromosome="chr24"

   # see only 2 alleles (not sure what cases would result in >2 in zebrafish)
   bcftools view --max-alleles 2 ${infile_path} ${chromosome} | \
      snzl --no-filtered --output ${outfile_path} \
         - candidate_snp -csm ${sample_name} \
         snpeff -sem ${sample_name} -see EFF -semi MODERATE

Although this was a strict candidate, no candidates were detected::

   grep -v '#' test_dp_out.vcf | wc -l
   # get zero lines detected

Controls (SNP quantity)
***********************

.. note::
    ``./.:.:.:.:.`` indicates a non-SNP event in the corresponding sample.

Sample with varying numbers of SNPs across the four references {0,1,2,3,4}. In this sample, *read depth* is set to 100. SNP effect was predicted to be ``MODERATE``, matching filter threshold. For the purposes of this test, SNP positions are artificial.

This sample was identified as a region of high interest by the original pipeline.

.. code-block:: shell

   OUT="test_snp.vcf"
   SNP='1/1:0,16:16:15:158,15,0'
   NAN='./.:.:.:.:.'

   head -n1172 TL2312073-163-4L-MAN-20231_Ref_merged_ChromFixed.vcf.annotated.vcf > header.vcf
   mv header.vcf ${OUT}
   paste -d '\t' \
      <(grep -m 1 5996897 chr24.vcf | cut -f1-10 | sed "s|5996897|1|") \
      <(printf "${NAN}\t${NAN}\t${NAN}\t${NAN}\n") >> ${OUT}
   paste -d '\t' \
      <(grep -m 1 5996897 chr24.vcf | cut -f1-10 | sed "s|5996897|2|") \
      <(printf "${SNP}\t${NAN}\t${NAN}\t${NAN}\n") >> ${OUT}
   paste -d '\t' \
      <(grep -m 1 5996897 chr24.vcf | cut -f1-10 | sed "s|5996897|3|") \
      <(printf "${SNP}\t${SNP}\t${NAN}\t${NAN}\n") >> ${OUT}
   paste -d '\t' \
      <(grep -m 1 5996897 chr24.vcf | cut -f1-10 | sed "s|5996897|4|") \
      <(printf "${SNP}\t${SNP}\t${SNP}\t${NAN}\n") >> ${OUT}
   paste -d '\t' \
      <(grep -m 1 5996897 chr24.vcf | cut -f1-10 | sed "s|5996897|5|") \
      <(printf "${SNP}\t${SNP}\t${SNP}\t${SNP}\n") >> ${OUT}
   sed -i "s|DP=25|DP=100|" $OUT
   bgzip ${OUT}
   bcftools index -t ${OUT}.gz
   
Here is an example of one entry in the file. Only the ``SNP`` field is modified in each iteration, with an increasing number of blanks: ``./.:.:.:.:.`` ::

   chr24   5       .       A       T       129.9   .       BaseQRankSum=0;Dels=0;ExcessHet=3.0103;FS=0;HaplotypeScore=0;MQ=12.77;MQ0=6;MQRankSum=0.967;QD=3.97;ReadPosRankSum=0.967;SOR=1.179;DP=100;AF=1;MLEAC=2;MLEAF=1;AN=4;AC=4;EFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gaA/gaT|E32D|310|si:ch211-193e5.4|protein_coding|CODING|ENSDART00000132686|2|T|WARNING_TRANSCRIPT_INCOMPLETE),UPSTREAM(MODIFIER||1545||310|si:ch211-193e5.3|protein_coding|CODING|ENSDART00000077933||T|WARNING_TRANSCRIPT_INCOMPLETE),UPSTREAM(MODIFIER||2423||279|si:ch211-193e5.3|protein_coding|CODING|ENSDART00000153736||T|WARNING_TRANSCRIPT_NO_START_CODON),UPSTREAM(MODIFIER||237||278|si:ch211-193e5.4|protein_coding|CODING|ENSDART00000077922||T|WARNING_TRANSCRIPT_NO_START_CODON),DOWNSTREAM(MODIFIER||2894|||CT573382.1|miRNA|NON_CODING|ENSDART00000119433||T),DOWNSTREAM(MODIFIER||705|||CT573382.2|miRNA|NON_CODING|ENSDART00000119567||T),DOWNSTREAM(MODIFIER||2377||503|acbd5a|protein_coding|CODING|ENSDART00000135124||T),DOWNSTREAM(MODIFIER||2377||502|acbd5a|protein_coding|CODING|ENSDART00000122018||T),DOWNSTREAM(MODIFIER||2377||501|acbd5a|protein_coding|CODING|ENSDART00000007373||T)      GT:AD:DP:GQ:PL  1/1:7,2:9:6:63,6,0      1/1:0,16:16:15:158,15,01/1:0,16:16:15:158,15,0  1/1:0,16:16:15:158,15,0 1/1:0,16:16:15:158,15,0

The ``snzl`` strict candidate pipeline is run.

.. code:: shell

   infile_path="test_snp.vcf.gz"
   outfile_path="test_snp_out.vcf"
   sample_name="TL2312073-163-4L-MAN-20231116"
   chromosome="chr24"

   # see only 2 alleles (not sure what cases would result in >2 in zebrafish)
   bcftools view --max-alleles 2 ${infile_path} ${chromosome} | \
      snzl --no-filtered --output ${outfile_path} \
         - candidate_snp -csm ${sample_name} \
         snpeff -sem ${sample_name} -see EFF -semi MODERATE

All non-strict candidates were retained and all strict candidates were discarded::

   grep -v '#' test_snp_out.vcf | grep -P './.:.:.:.:.\t./.:.:.:.:.\t./.:.:.:.:.\t./.:.:.:.:.' | wc -l
   # 0

Controls (SNP presence and DP quantity)
***************************************

Unit test of old code
#####################

.. code:: shell

   # we know chr24:423800 has info we want
   infile_path="TL2312073-163-4L-MAN-20231_Ref_merged_ChromFixed.vcf.annotated.vcf.gz"
   outfile_path="strict_candidates.vcf"
   sample_name="TL2312073-163-4L-MAN-20231116"
   chromosome="chr24"

   # see only 2 alleles (not sure what cases would result in >2 in zebrafish)
   bcftools view --max-alleles 2 ${infile_path} ${chromosome} | \
      snzl --no-filtered --output ${outfile_path} \
         - candidate_snp -csm ${sample_name} \
         snpeff -sem ${sample_name} -see EFF -semi MODERATE

.. raw:: html

   </details>


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