EMUMADZ pipeline case study 1
=============================
BAM to SNP for F0-F2 comparisons
--------------------------------

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

.. note::
    All steps are run in directory ``source``.

Data setup
++++++++++

Software requirements:

- bcftools
- gatk
- samtools
- snpEff

.. note::
    Here we assume that ``bcftools,gatk,samtools,snpEff`` are all installed and the user has completed the steps in the install instructions.

Data requirements:
- Reference genome fasta file
- Samplesheet with sample information
- Chromosome name mappings

Here is the directory structure::

    project_root/
    ├── source/         # Scripts are run from here
    ├── data/           # Input data (BAM files, reference genome, samplesheet)
    └── results/        # All output files organized by processing step
        ├── 01_variant_calling/
        ├── 02_chromosome_filtering/
        ├── 03_chromosome_filtering/
        ├── 04_mutant_analysis/
        ├── 05_annotation/
        └── 06_final_results/

Create a samplesheet with sample identifiers, paths to bam file and the sample type::

    sample_identity alignment_file  sample_type
    some_sample_1   /path/to/1.bam  reference
    some_sample_2   /path/to/2.bam  reference
    some_sample_3   /path/to/3.bam  reference
    some_sample_4   /path/to/4.bam  reference
    some_sample_5   /path/to/5.bam  mutant
    ...             ...             ...     

In our case, this is how our table looks like:

.. csv-table:: Example samplesheet for F0-F2 comparison
   :file: tables/samplesheet.F0-F2.csv
   :header-rows: 1

Here we setup the file structure and some other metadata.

.. code-block:: shell

    # convenience function since it will be used in most downstream operations
    setup() {
        # navigate to source directory
        if [[ ! $(basename $PWD) -eq source ]]; then
            echo "All scripts are set to run from directory source"
            echo "See directory structure above for more information"
            exit 1
        fi

        # various configuration variables
        THREADS=16

        # GATK
        GATK_MEM="64g"
        GATK_JAVA_OPTS="-Xmx${GATK_MEM} -XX:+UseParallelGC"

        # ENSEMBL VEP
        VEP_ASSEMBLY="Zv9"
        VEP_VERSION="79"
        VEP_CACHE="${HOME}/.vep"
        VEP_BUFFER="8192"

        # snpEff
        SNPEFF_CONFIG="/home/tchen/.local/share/mamba/envs/snps/share/snpeff-5.2-1/snpEff.config"
        SNPEFF_GENOME="Zv9.75"

        # paths relative to source directory
        DATA_DIR="../data/"
        RESULTS_DIR="../results/"
        SAMPLESHEET="${DATA_DIR}/samplesheet.F0-F2.tsv"
        REFERENCE_FA="${DATA_DIR}/Reference/GCA_000002035.2_Zv9_genomic.fna"
        REF_FIXED_FA="${DATA_DIR}/Reference/GCA_000002035.2_Zv9_genomic.ChrFixed.fna"
        CHROM_MAP="${DATA_DIR}/chrom_table.tsv"

        # parse samplesheet and load samples into array
        MUT_SAMPLES=($(awk '$3=="mutant" {print $1}' ${SAMPLESHEET}))
        REF_SAMPLES=($(awk '$3=="reference" {print $1}' ${SAMPLESHEET}))
        echo "${#MUT_SAMPLES[@]} mutants: ${MUT_SAMPLES[@]}"
        echo "${#REF_SAMPLES[@]} references: ${REF_SAMPLES[@]}"

        # create results directory structure
        mkdir -p ${RESULTS_DIR}/{01_variant_called,02_chromosome_filtered,03_stringent_filtered,04_normalised,05_samples_merged,06_snps_filtered,07_mutant_candidates,08_annot_eff,08_annot_vep,08_annot_all,09_finalised,10_visualisations}

        # set downstream directories and suffixes
        # if step 3 is actioned, edit these as needed
        PROCESSED_VCF_DIR="${RESULTS_DIR}/02_chromosome_filtered"
    }

Create some metadata files. The chromosome file maps the chromosome IDs onto the ``chr{1..25}`` more human-readable format. Have the ``fasta`` genome reference file on hand.

.. code-block:: shell

    setup_metadata() {
        # verify reference genome preparation
        if [[ ! -f "${REFERENCE_FA}.fai" ]]; then
            echo "Creating FASTA index..."
            samtools faidx ${REFERENCE_FA}
        fi
        if [[ ! -f "${REFERENCE_FA%.*}.dict" ]]; then
            echo "Creating sequence dictionary..."
            gatk CreateSequenceDictionary -R ${REFERENCE_FA}
        fi

        # same but for the fixed chromosomes
        if [[ ! -f "${REF_FIXED_FA}.fai" ]]; then
            echo "Creating FASTA index..."
            samtools faidx ${REFERENCE_FA}
        fi

        # create chromosome list (adjust range as needed)
        grep -P "^chr[0-9]+\t" ${REF_FIXED_FA}.fai > \
            ${DATA_DIR}/chromosomes.txt
    }

    setup
    setup_metadata

Whole zebrafish genome assembly
+++++++++++++++++++++++++++++++

*The ``seqliner`` assembly step is not covered in this documentation.*

Variant calling and read filtering
++++++++++++++++++++++++++++++++++

We filter the ``bam`` files for low quality reads and perform variant calling simultaneously.

.. note::
    We apply intentionally minimal filtering at the early stage and apply custom filters later. Here, we are mainly interested in filtering poorly mapped reads.

.. code-block:: shell

    # this helper function retrieves bam file path given sample id
    get_bam_path() {
        local sample_id=$1
        awk -v id="$sample_id" '$1==id {print $2}' ${SAMPLESHEET}
    }

    call_variants() {
        local sample=$1
        local sample_type=$2
        local bam_path=$(get_bam_path $sample)
        
        echo "Processing ${sample_type}: $sample ($bam_path)"
        
        # Verify BAM file and index exist
        if [[ ! -f "${bam_path}" ]]; then
            echo "No BAM: ${bam_path}"
            return 1
        fi
        if [[ ! -f "${bam_path/bam/bai}" ]]; then
            echo "No index: ${bam_path/bam/bai}"
            return 1
        fi
        
        # gatk HaplotypeCaller with quality filters and threading
        gatk --java-options "${GATK_JAVA_OPTS}" HaplotypeCaller \
            -R ${REFERENCE_FA} \
            -I ${bam_path} \
            -O "${RESULTS_DIR}/01_variant_calling/${sample}.vcf.gz" \
            --minimum-mapping-quality 20 \
            --min-base-quality-score 20 \
            --stand-call-conf 30 \
            --max-reads-per-alignment-start 50 \
            --max-mnp-distance 1 \
            --native-pair-hmm-threads ${THREADS} \
            --verbosity INFO
                
        echo "Completed variant calling for ${sample_type}: $sample"        
    }

    setup
    export REFERENCE_FA RESULTS_DIR GATK_JAVA_OPTS THREADS SAMPLESHEET

    for ref in "${REF_SAMPLES[@]}"; do call_variants ${ref} reference; done
    for mut in "${MUT_SAMPLES[@]}"; do call_variants ${mut} mutant; done

    echo "Variant calling completed for all samples"

.. caution:: 
    The above was tested given an allocation of 40 cpus and 80 GB of memory. You may want to adjust the number of parallel jobs and memory corresponding to your available resources.

.. caution::
    ``gatk`` settings are specific to version ``4.5.0.0-gcc-13.2.0`` and may not work for other versions.


.. raw:: html

   <details>
   <summary><a>Optional hard filtering.</a></summary>

.. warning::
    These filters may be too stringent for this use case. Only use if you are seeing an excess of false positives. Note that it is possible to later filter out false positives.

.. code-block:: shell

    # Function to apply GATK hard filters
    apply_stringent_filters() {
        local sample=$1
        local infile_dir="${RESULTS_DIR}/02_chromosome_filtered"
        local outfile_dir="${RESULTS_DIR}/03_stringent_filtered"
        
        echo "Applying hard filters to: $sample"
        
        # Apply GATK hard filters
        gatk --java-options "${GATK_JAVA_OPTS}" VariantFiltration \
            -R ${REFERENCE_FA} \
            -V "${infile_dir}/${sample}.vcf" \
            -O "${outfile_dir}/${sample}_filtered.vcf" \
            --filter-expression "QD < 2.0" --filter-name "QD2" \
            --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
            --filter-expression "SOR > 3.0" --filter-name "SOR3" \
            --filter-expression "FS > 60.0" --filter-name "FS60" \
            --filter-expression "MQ < 40.0" --filter-name "MQ40" \
            --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
        
        # select only PASS variants
        gatk --java-options "${GATK_JAVA_OPTS}" SelectVariants \
            -R ${REFERENCE_FA} \
            -V "${outfile_dir}/${sample}_filtered.vcf" \
            -O "${outfile_dir}/${sample}.vcf" \
            --exclude-filtered
        
        # clean up intermediate files
        rm -f "${outfile_dir}/${sample}_filtered.vcf"
        
        echo "Completed stringent filtering for: $sample"
    }
    
    # store hard-filtered variants
    setup
    mkdir -p "${RESULTS_DIR}/03_hard_filtering"
    
    # Process all samples
    for sample in "${ALL_SAMPLES[@]}"; do
        apply_hard_filters $sample
    done

.. raw:: html

    </details>


Rename chromosomes in vcf with map file, and discard non-chromosomes
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

``data/chrom_table.tsv`` contains the chromosome name mappings. Use this to remap the chromosomes.

At the same time, we are only interested in the chromosomes, which start with ``chr`` and go from ``chr1`` to ``chr25``. So we discard all other contigs (including from the header).

.. code-block:: shell

    filter_and_rename_chromosomes() {
        local sample=$1
        local infile_dir="${RESULTS_DIR}/01_variant_called"
        local outfile_dir="${RESULTS_DIR}/02_chromosome_filtered"
        
        echo "Filtering chromosomes for: $sample"
        # 33550336 use extreme range to cover all
        for i in {1..25}; do 
            echo -e "chr${i}\t1\t1111111111111000000000000"
        done > ../data/chromosomes_only.txt

        # we are going to paste text, so leave it decompressed
        # rename all contigs with mapping
        bcftools annotate --threads ${THREADS} \
            --rename-chr ${CHROM_MAP} \
            ${infile_dir}/${sample}.vcf.gz | \
        grep -v '##contig=<ID=FR' > "${outfile_dir}/${sample}_reannot.vcf"
        
        # filter to main chromosomes only and remove non-chromosomal variants from header
        bcftools view -t $(cut -f1 ${DATA_DIR}/chromosomes.txt | paste -sd,) \
            "${outfile_dir}/${sample}_reannot.vcf" \
            --threads ${THREADS} \
            -o "${outfile_dir}/${sample}_chr_temp.vcf"
        
        # clean up header to remove non-chromosomal contigs
        bcftools reheader -f ${DATA_DIR}/chromosomes.txt \
            "${outfile_dir}/${sample}_chr_temp.vcf" | \
        bgzip -@ ${THREADS} -c > "${outfile_dir}/${sample}.vcf.gz"
        
        bcftools index --threads ${THREADS} \
            "${outfile_dir}/${sample}.vcf.gz"

        # clean up intermediate files
        rm -f "${outfile_dir}/${sample}_reannot.vcf" \
            "${outfile_dir}/${sample}_chr_temp.vcf"
        
        echo "Completed chromosome filtering for: $sample"
    }

    # process all samples sequentially
    setup

    ALL_SAMPLES=("${MUT_SAMPLES[@]}" "${REF_SAMPLES[@]}")
    for sample in "${ALL_SAMPLES[@]}"; do
        filter_and_rename_chromosomes $sample
    done

Decompose multiallelic variants into individual instances
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The aim of the normalisation is to preserve differential events that may occur across alleles. For example, if a `reference control` is heterozygous for ``A/T`` point mutation, and the `mutant sample` is homozygous for a ``G`` point mutation, the mutant is considered to have a SNP event.

.. code-block:: shell

    normalise_vcf() {
        local sample=$1
        local infile_dir="${PROCESSED_VCF_DIR}"
        local outfile_dir="${RESULTS_DIR}/04_normalised/"
        local input_vcf="${infile_dir}/${sample}.vcf.gz"
        local output_vcf="${outfile_dir}/${sample}.vcf.gz"
        
        echo "Normalising VCF for: $sample"
        
        # normalize variants: decompose complex variants and left-align indels
        bcftools norm -m-both -f "${REF_FIXED_FA}" "${input_vcf}" \
            --threads ${THREADS} --write-index -Ob -o "${output_vcf}"
        
        echo "Completed normalisation for: $sample"
        return 0
    }
    
    echo "Normalising VCFs before merging..."

    ALL_SAMPLES=("${MUT_SAMPLES[@]}" "${REF_SAMPLES[@]}")
    for sample in "${ALL_SAMPLES[@]}"; do
        normalise_vcf $sample
    done

Merge files with refs
+++++++++++++++++++++

Samples are then merged for side-by-side comparison between `mutant samples` and `reference controls`.

For each sample, merge the four references:
- TL2209397-ENUref-Female.vcf.gz
- TL2209398-ENUref-Male.vcf.gz
- TL2209399-TUref-Female.vcf.gz
- TL2209400-TUref-Male.vcf.gz

Each sample occupies the first slot in the vcf sample columns, followed by the references.

.. caution::
    The order of samples and references is important.

.. code-block:: shell

    merge_vcf() {
        local sample=$1
        local infile_dir="${RESULTS_DIR}/04_normalised/"
        local outfile_dir="${RESULTS_DIR}/05_samples_merged/"

        # create VCF list for this mutant + all references (using normalised files)
        VCF_LIST="${infile_dir}/${sample}.vcf.gz"
        for ref in "${REF_SAMPLES[@]}"; do
            VCF_LIST="${VCF_LIST} ${infile_dir}/${ref}.vcf.gz"
        done

        # merge this mutant with all references
        echo "Merging ${sample} with references..."
        bcftools merge ${VCF_LIST} --threads ${THREADS} --write-index \
            -Ob -o "${outfile_dir}/${sample}.vcf.gz"
    }
    
    for sample in "${MUT_SAMPLES[@]}"; do
        merge_vcf $sample
    done


Filter for SNP only events
++++++++++++++++++++++++++

Other mutations are not of interest since `N`-ethyl-`N`-nitrosourea (ENU) used in the forward genetic screen mainly induces point mutations.

.. code-block:: shell

    snps_vcf() {
        local sample=$1
        local infile_dir="${RESULTS_DIR}/05_samples_merged/"
        local outfile_dir="${RESULTS_DIR}/06_snps_filtered/"

        # filter for SNPs only
        echo "Filtering SNPs for ${sample}..."
        bcftools view -i 'TYPE="snp"' --threads ${THREADS} \
            "${infile_dir}/${sample}.vcf.gz" \
            --write-index -Ob -o "${outfile_dir}/${sample}.vcf.gz"
    }
    
    for sample in "${MUT_SAMPLES[@]}"; do 
        snps_vcf $sample
    done


Identify candidate SNPs
+++++++++++++++++++++++

Obtain SNPs which occur in the sample but not the references. There are several criteria in our filter, summarised below:
1. `mutant sample` must be covered to a depth of >= 1 read.
2. At least one `reference control` must be covered to a depth of >= 1 read.
3. If the `reference control` is multiallelic for a SNP, it should not match the `mutant sample`.

Beyond this, a light mapping quality filter was implemented during the variant calling stage.

The `reference control` may be for example a grandparent or a wild-type sibling. The pipeline is agnostic to input and will run regardless. An arbitrary number of `reference controls` can be used.

We explain the criteria above with some examples.

Showing the concept of coverage in (1) and (2)::

    MUT SAMPLE           PASS                FAIL                PASS

        MAPPED      --------------                          --------------
        GENOME      ==============      ==============      ==============
    
    REF CONTROL 1   

        MAPPED      --------------                          --------------
        GENOME      ==============      ==============      ==============
    
    REF CONTROL 2

        MAPPED      --------------      
        GENOME      ==============      ==============      ==============


Showing the concept of multiallelic SNP events in (3):
- ORI is the `true original` 
- REF is the `reference control`
- MUT is the `mutant sample`

Let us assume, the true canonical nucleotide is ``C`` at a given position.::

    ==================
    PLATONIC IDEAL REF
    ==================

    ORI ------C-------
        ------C-------

    =====EXAMPLE1=====
    REF homozygous
    MUT homozygous
    SNP in ORI identical to REF    
    SNP in REF different to MUT
    This is a SNP event

    ORI ------C-------
        ------C-------

    REF ------C-------
        ------C-------

    MUT ------G-------
        ------G-------

    =====EXAMPLE2=====
    REF homozygous
    MUT homozygous
    SNP in ORI different to REF
    SNP in REF identical to MUT
    This is not a SNP event

    ORI ------C-------
        ------C-------

    REF ------G-------
        ------G-------

    MUT ------G-------
        ------G-------

    =====EXAMPLE3=====
    REF homozygous
    MUT heterozygous
    SNP in ORI different to REF
    SNP in REF different to MUT
    This is a SNP event

    ORI ------C-------
        ------C-------

    REF ------A-------
        ------T-------

    MUT ------G-------
        ------G-------

    =====EXAMPLE4=====
    REF heterozygous
    MUT homozygous
    SNP in ORI different to REF
    SNP in REF partially different to MUT
    This is not a SNP event

    ORI ------C-------
        ------C-------

    REF ------A-------
        ------G-------

    MUT ------G-------
        ------G-------

    =====EXAMPLE5=====
    REF homozygous
    MUT heterozygous
    SNP in ORI different to REF
    SNP in REF partially different to MUT
    This is not a SNP event

    ORI ------C-------
        ------C-------

    REF ------G-------
        ------G-------

    MUT ------A-------
        ------G-------

The three criteria are implemented together.

1. `mutant sample` must be covered to a depth of >= 1 read.
2. At least one `reference control` must be covered to a depth of >= 1 read.
3. If the `reference control` is multiallelic for a SNP, it should not match the `mutant sample`.

.. code-block:: shell

    find_candidates_vcf() {
        local sample=$1
        local infile_dir="${RESULTS_DIR}/06_snps_filtered/"
        local outfile_dir="${RESULTS_DIR}/07_mutant_candidates/"

        NUM_REFS=${#REF_SAMPLES[@]}

        # mutant AND at least one ref has >=1 read depth
        DP_FILTER="FORMAT/DP[0]>=1"
        for ((i=1; i<=NUM_REFS; i++)); do
            if [ $i -eq 1 ]; then
                DP_FILTER="${DP_FILTER} && (FORMAT/DP[${i}]>=1"
            else
                DP_FILTER="${DP_FILTER} || FORMAT/DP[${i}]>=1"
            fi
        done
        # the trailing bracket isnt a typo btw
        DP_FILTER="${DP_FILTER})"

        # instead of a strict genotyping, we relax thresholds
        # this more effectively approximates genotype instead
        # mutant needs read depth of >= 0.85 and refs <= 0.15
        AF_FILTER='FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85'
        for ((i=1; i<=NUM_REFS; i++)); do
            # Missing data in references = assumed reference (good)
            # Only filter if there IS alt allele data and it's too high
            AF_FILTER="${AF_FILTER} && (FORMAT/AD[${i}:1]=\".\" || FORMAT/AD[${i}:0]=\".\" || FORMAT/AD[${i}:1]/(FORMAT/AD[${i}:0]+FORMAT/AD[${i}:1]) <= 0.15)"
        done        

        # combine filters
        STRICT_FILTER="(${DP_FILTER}) && (${AF_FILTER})"

        bcftools filter -i "${STRICT_FILTER}" --threads ${THREADS} --write_index \
            "${infile_dir}/${sample}.vcf.gz" -Ob -o "${outfile_dir}/${sample}.vcf.gz"
    }
    
    for sample in "${MUT_SAMPLES[@]}"; do 
        snps_vcf $sample
    done


Screen impact of SNP
++++++++++++++++++++

Either ``VEP`` or ``snpEff`` will produce similarly formatted results. We run both separately and combine their annotations at the end.

VEP
***

We run ENSEMBL's variant effect predictor ``VEP`` on the data. Install instructions are provided separately on our install page.

.. important::
    Since we are using an ancient zebrafish genome, a specific combination of settings must be used. A prerequisite is installing an older cached version of the genome, which is covered in the install section. If this was not done, it is unlikely that ``VEP`` will work as intended.

.. code-block:: shell

    annot_vep() {
        local sample=$1
        local infile_dir="${RESULTS_DIR}/07_mutant_candidates/"
        local outfile_dir="${RESULTS_DIR}/08_annot_vep/"

        echo "Predicting variant impact for ${sample}..."
        vep \
            --cache \
            --dir_cache ${VEP_CACHE} \
            --species danio_rerio \
            --assembly ${VEP_ASSEMBLY} \
            --cache_version ${VEP_VERSION} \
            --offline \
            --regulatory \
            --vcf \
            --variant_class \
            --fork ${THREADS} \
            --buffer_size ${VEP_BUFFER} \
            --stats_html \
            --stats_text \
            --force_overwrite \
            --compress_output bgzip \
            --input_file "${infile_dir}/${sample}.vcf.gz" \
            --output_file "${outfile_dir}/${sample}.vcf.gz"
        bcftools index --threads ${THREADS} \
            "${outfile_dir}/${sample}.vcf.gz"            
    }
    
    for sample in "${MUT_SAMPLES[@]}"; do 
        annot_vep $sample
    done

snpEff
******

We run ``snpEff`` to screen for SNPs which are predicted to be impactful. Annotations are added to the data in the process.

.. warning:: 
    ``snpEff`` runs on ``java``, and as a result may run into memory and other issues. You may need to set ``java`` memory flags accordingly.

    .. code-block:: shell

        # adjust max values accordingly
        export JAVA_OPTIONS="-Xms512m -Xmx16g"

.. code-block:: shell

    # SnpEff version SnpEff 5.2 (build 2023-09-29 06:17)
    snpEff download Zv9.75

    annot_eff() { 
        local sample=$1
        local infile_dir="${RESULTS_DIR}/07_mutant_candidates/"
        local outfile_dir="${RESULTS_DIR}/08_annot_eff/"

        export JAVA_OPTIONS="-Xms512m -Xmx16g"

        csv_stat=${outfile_dir}/${sample}.csv
        fas_prot=${outfile_dir}/${sample}.fa
        web_stat=${outfile_dir}/${sample}.html
        out_path=${outfile_dir}/${sample}.vcf

        snpEff \
            -c "${SNPEFF_CONFIG}" \
            -csvStats "${csv_stat}" \
            -i vcf \
            -o vcf \
            -fastaProt "${fas_prot}" \
            -s "${web_stat}" \
            "${SNPEFF_GENOME}" \
            "${infile_dir}/${sample}.vcf.gz" > "${out_path}"
        bgzip -@ ${THREADS} ${out_path}
        bcftools index --threads ${THREADS} ${out_path}.gz
    }

    for sample in "${MUT_SAMPLES[@]}"; do 
        annot_eff $sample
    done

Combine SNP annotations
***********************

Both ``VEP`` and ``snpEff`` annotations are now available. These are recombined so that the aggregated information is available in one file. Note that we split the ``VEP`` and ``snpEff`` steps for efficiency, but running either method sequentially will obtain the same result.

.. code-block:: shell

    annot_all() {
        local sample=$1
        local annot_eff="${RESULTS_DIR}/08_annot_eff/"
        local annot_vep="${RESULTS_DIR}/08_annot_vep/"
        local outfile_dir="${RESULTS_DIR}/08_annot_all/"

        bcftools merge --threads ${THREADS} \
            --force-samples \
            ${annot_vep}/${sample}.vcf.gz \
            ${annot_eff}/${sample}.vcf.gz \
            --write-index -Ob -o ${outfile_dir}/${sample}.vcf.gz
    }
    
    for sample in "${MUT_SAMPLES[@]}"; do 
        annot_all $sample
    done

Visualising data
++++++++++++++++

Preparing data for visualisation
********************************

The aim is to visualise the tracks containing all reads at each SNP position. Visualising the entire genome is ideal but the size of the alignment data makes this step impractical. Therefore, we subset the ``bam`` alignment files [-500,500] around each SNP for a portable but still informative ``bam`` file. The full ``bam`` file can also be loaded and mounted on a separate volume if needed.

In this step, we subset the ``bam`` files for portability.

.. code-block:: shell

    subset_bam() {
        local sample=$1
        local infile_dir="${RESULTS_DIR}/08_annot_all/"
        local outfile_dir="${RESULTS_DIR}/10_visualisations/"

        python parse_vcf.py \
            ${infile_dir}/${sample}.vcf.gz \
            ${SAMPLESHEET} \
            ${outfile_dir}/${sample}.csv
    }
    
    for sample in "${MUT_SAMPLES[@]}"; do 
        subset_bam $sample
    done

Here we parse the final ``vcf`` files into a format suitable for visualisation.

.. code-block:: shell

    prep_vis() {
        local sample=$1
        local infile_dir="${RESULTS_DIR}/08_annot_all/"
        local outfile_dir="${RESULTS_DIR}/10_visualisations/"

        python parse_vcf.py \
            "${infile_dir}/${sample}.vcf.gz" \
            "${SAMPLESHEET}" \
            "${outfile_dir}/${sample}.csv" \
            --subset_path "${outfile_dir}/${sample}.bam" \
            --subset_workers ${THREADS} \
            --coverage_report "${outfile_dir}/${sample}_coverage.csv" \
            --chrom_mapping ${CHROM_MAP}
    }
    
    for sample in "${MUT_SAMPLES[@]}"; do 
        prep_vis $sample
    done




