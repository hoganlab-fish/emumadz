BAM to SNP
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

.. note::
    All steps are run in directory ``source``.

Data setup
++++++++++

Original data (reference only)
##############################

.. raw:: html

   <details>
   <summary><a>Original data run preserved for the record.</a></summary>

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

For the purposes of rerunning the postprocessing, we only want the ``bam`` alignment files and original ``vcf`` files, since we will be recreating everything after this step. Running ``3.setup_dataset.sh`` will create the corresponding input and output directories::

   ../data/Alignment_File
   ../data/VCF_Original
   ../results/Sample_Identity
   ../results/Fastq_File
   ../results/Alignment_File
   ../results/VCF_Merged
   ../results/VCF_ChrFixed
   ../results/VCF_Annotated
   ../results/VCF_Candidates
   ../results/Snzl_NoGaps_NBases
   ../results/Snzl_NoGaps_NSnps
   ../results/Snzl_WithGaps_NBases
   ../results/Snzl_WithGaps_NSnps
   ../results/Json
   
While the ``bam`` and ``vcf`` files are not directly included in this repository due to size constraints, ``md5`` sums are preserved in ``data/Alignment_File/bam.md5`` and ``data/VCF_Original/vcf.md5``.

Four samplesheet files are generated:
- ``samplesheet.220701_A01221_0125_BHCCHNDMXY.tsv``
- ``samplesheet.220930_A00692_316_AHVMJFDSX3.tsv``
- ``samplesheet.221014_A00692_0319_BHGJ7CDMXY.tsv``
- ``samplesheet.231201_A00692_0394_231208_A00692_0396.tsv``

The ``samplesheet.220701_A01221_0125_BHCCHNDMXY.tsv`` contains ``F0`` grandparent reference generation metadata. All other samplesheets contain ``F2`` generation metadata. ``samplesheet_original.tsv`` contains these legacy file paths.

``3.setup_dataset.sh`` copies ``bam`` alignment files and ``vcf`` variant calling files over from their specified locations in ``samplesheet_original.tsv`` to ``data/Alignment_File`` and ``data/VCF_Original`` respectively. We generate ``samplesheet_postprocess.tsv`` for use in this postprocessing analysis.

.. raw:: html

   </details>

Data setup
++++++++++

Software requirements:
- bcftools
- gatk
- samtools
- snpEff

.. note::
    Here we assume that ``bcftools,gatk,samtools,snpEff`` are all installed. Refer back to the installation instructions for how to install these.

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
        cd source/

        # various configuration variables
        THREADS=16
        GATK_MEM="64g"
        GATK_JAVA_OPTS="-Xmx${GATK_MEM} -XX:+UseParallelGC"

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
        mkdir -p ${RESULTS_DIR}/{01_variant_calling,02_chromosome_filtering,03_hard_filtering,04_mutant_analysis,05_annotation,06_final_results}

        # set downstream directories and suffixes
        # if step 3 is actioned, edit these as needed
        PROCESSED_VCF_DIR="${RESULTS_DIR}/02_chromosome_filtering"
        VCF_SUFFIX="_chr.vcf"
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
        head -n 25 ${REF_FIXED_FA}.fai > ${DATA_DIR}/chromosomes.txt
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
            -O "${RESULTS_DIR}/01_variant_calling/${sample}_raw.vcf.gz" \
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
    apply_hard_filters() {
        local sample=$1
        local input_dir="${RESULTS_DIR}/02_chromosome_filtering"
        local output_dir="${RESULTS_DIR}/03_hard_filtering"
        
        echo "Applying hard filters to: $sample"
        
        # Apply GATK hard filters
        gatk --java-options "${GATK_JAVA_OPTS}" VariantFiltration \
            -R ${REFERENCE_FA} \
            -V "${input_dir}/${sample}_chr.vcf" \
            -O "${output_dir}/${sample}_filtered.vcf" \
            --filter-expression "QD < 2.0" --filter-name "QD2" \
            --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
            --filter-expression "SOR > 3.0" --filter-name "SOR3" \
            --filter-expression "FS > 60.0" --filter-name "FS60" \
            --filter-expression "MQ < 40.0" --filter-name "MQ40" \
            --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
        
        # Select only PASS variants
        gatk --java-options "${GATK_JAVA_OPTS}" SelectVariants \
            -R ${REFERENCE_FA} \
            -V "${output_dir}/${sample}_filtered.vcf" \
            -O "${output_dir}/${sample}_pass.vcf" \
            --exclude-filtered
        
        # Normalize with bcftools
        bcftools norm -m-both -f ${REFERENCE_FA} "${output_dir}/${sample}_pass.vcf" \
            --threads ${THREADS} -o "${output_dir}/${sample}_chr.vcf"
        
        # Clean up intermediate files
        rm -f "${output_dir}/${sample}_filtered.vcf" "${output_dir}/${sample}_pass.vcf"
        
        echo "Completed hard filtering for: $sample"
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
        local input_dir="${RESULTS_DIR}/01_variant_calling"
        local output_dir="${RESULTS_DIR}/02_chromosome_filtering"
        
        echo "Filtering chromosomes for: $sample"
        # 33550336 use extreme range to cover all
        for i in {1..25}; do 
            echo -e "chr${i}\t1\t1111111111111000000000000"
        done > ../data/chromosomes_only.txt

        # we are going to paste text, so leave it decompressed
        # rename all contigs with mapping
        bcftools annotate --threads ${THREADS} \
            --rename-chr ${CHROM_MAP} \
            ${input_dir}/${sample}_raw.vcf.gz | \
        grep -v '##contig=<ID=FR' > "${output_dir}/${sample}_reannot.vcf"
        
        # filter to main chromosomes only and remove non-chromosomal variants from header
        bcftools view -t $(cut -f1 ${DATA_DIR}/chromosomes.txt | paste -sd,) \
            "${output_dir}/${sample}_reannot.vcf" \
            --threads ${THREADS} \
            -o "${output_dir}/${sample}_chr_temp.vcf"
        
        # clean up header to remove non-chromosomal contigs
        bcftools reheader -f ${DATA_DIR}/chromosomes.txt \
            "${output_dir}/${sample}_chr_temp.vcf" | \
        bgzip -@ ${THREADS} -c > "${output_dir}/${sample}.vcf.gz"
        
        bcftools index --threads ${THREADS} \
            "${output_dir}/${sample}.vcf.gz"

        # clean up intermediate files
        rm -f "${output_dir}/${sample}_reannot.vcf" \
            "${output_dir}/${sample}_chr_temp.vcf"
        
        echo "Completed chromosome filtering for: $sample"
    }

    # process all samples sequentially
    setup

    ALL_SAMPLES=("${MUT_SAMPLES[@]}" "${REF_SAMPLES[@]}")
    for sample in "${ALL_SAMPLES[@]}"; do
        filter_and_rename_chromosomes $sample
    done

Strict candidate pipeline F0-F2
+++++++++++++++++++++++++++++++

The following steps are carried out for a single mutant:
1. To retain multiallelic snps, normalise all samples
2. To compare `mutants` with `references`, merge into same file
3. 


.. code-block:: shell

    # Function to normalize VCF files before merging
    normalize_vcf() {
        local sample=$1
        local input_dir="${PROCESSED_VCF_DIR}"
        local output_dir="${RESULTS_DIR}/04_mutant_analysis"
        local input_vcf="${input_dir}/${sample}${VCF_SUFFIX}"
        local output_vcf="${output_dir}/${sample}_normalized.vcf"
        
        echo "Normalizing VCF for: $sample"
        
        # Normalize variants: decompose complex variants and left-align indels
        bcftools norm -m-both -f ${REFERENCE_FA} "${input_vcf}" \
            --threads ${THREADS} \
            -o "${output_vcf}"
        
        echo "Completed normalization for: $sample"
        return 0
    }

    # Function to process a single mutant
    process_mutant() {
        local mutant=$1
        local output_dir="${RESULTS_DIR}/04_mutant_analysis"
        
        echo "================================================"
        echo "Processing mutant: $mutant against all references"
        echo "================================================"
        
        # First normalize all VCFs that will be merged (mutant + all references)
        echo "Normalizing VCFs before merging..."
        
        # Normalize mutant VCF
        normalize_vcf $mutant
        
        # Normalize all reference VCFs
        for ref in "${REFERENCE_SAMPLES[@]}"; do
            normalize_vcf $ref
        done
        
        # Create VCF list for this mutant + all references (using normalized files)
        VCF_LIST="${output_dir}/${mutant}_normalized.vcf"
        for ref in "${REFERENCE_SAMPLES[@]}"; do
            VCF_LIST="${VCF_LIST} ${output_dir}/${ref}_normalized.vcf"
        done
        
        # Merge this mutant with all references
        echo "Merging ${mutant} with references..."
        bcftools merge ${VCF_LIST} --threads ${THREADS} -o "${output_dir}/${mutant}_merged.vcf"
        
        # Filter for SNPs only
        echo "Filtering SNPs for ${mutant}..."
        bcftools view -v snps "${output_dir}/${mutant}_merged.vcf" \
            --threads ${THREADS} -o "${output_dir}/${mutant}_snps.vcf"
        
        # Generate dynamic filter expressions
        NUM_REFS=${#REFERENCE_SAMPLES[@]}
        
        # Build genotype filter (mutant=1/1, all refs don't contain "1")
        GT_FILTER='GT[0]="1/1"'
        for ((i=1; i<=NUM_REFS; i++)); do
            GT_FILTER="${GT_FILTER} && GT[${i}]!~\"1\""
        done
        
        # Build depth filter (at least one ref has DP>=1)
        DP_FILTER=""
        for ((i=1; i<=NUM_REFS; i++)); do
            if [ $i -eq 1 ]; then
                DP_FILTER="FORMAT/DP[${i}]>=1"
            else
                DP_FILTER="${DP_FILTER} || FORMAT/DP[${i}]>=1"
            fi
        done
        
        # Combine filters
        STRICT_FILTER="${GT_FILTER} && (${DP_FILTER})"
        
        echo "Applying strict candidate filter for ${mutant}..."
        bcftools filter -i "${STRICT_FILTER}" \
            "${output_dir}/${mutant}_snps.vcf" \
            --threads ${THREADS} \
            -o "${output_dir}/${mutant}_strict_candidates.vcf"
        
        # Generate AF filter for mutant (high AF) and refs (low AF)
        AF_FILTER='FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.8'
        
        for ((i=1; i<=NUM_REFS; i++)); do
            AF_FILTER="${AF_FILTER} && FORMAT/AD[${i}:1]/(FORMAT/AD[${i}:0]+FORMAT/AD[${i}:1]) <= 0.15"
        done
        
        echo "Applying allele frequency filter for ${mutant}..."
        bcftools filter -i "${AF_FILTER}" \
            "${output_dir}/${mutant}_strict_candidates.vcf" \
            --threads ${THREADS} \
            -o "${output_dir}/${mutant}_af_filtered.vcf"
        
        # Clean up intermediate files for this mutant
        rm -f "${output_dir}/${mutant}_normalized.vcf" "${output_dir}/${mutant}_merged.vcf" "${output_dir}/${mutant}_snps.vcf" "${output_dir}/${mutant}_strict_candidates.vcf"
        
        # Clean up normalized reference files (they'll be recreated for the next mutant)
        for ref in "${REFERENCE_SAMPLES[@]}"; do
            rm -f "${output_dir}/${ref}_normalized.vcf"
        done
        
        echo "Completed analysis for mutant: $mutant"
    }

    # Process mutants sequentially
    for mutant in "${MUTANT_SAMPLES[@]}"; do
        process_mutant $mutant
    done

    echo "Per-mutant analysis completed"


Merge files with refs
+++++++++++++++++++++

For each sample, merge the four references:
- TL2209397-ENUref-Female.vcf.gz
- TL2209398-ENUref-Male.vcf.gz
- TL2209399-TUref-Female.vcf.gz
- TL2209400-TUref-Male.vcf.gz

For one example:

.. code-block:: shell

    THREADS=16
    REF_PATHS="../results/VCF_ChrFixed/TL2209397-ENUref-Female.vcf.gz
    ../results/VCF_ChrFixed/TL2209398-ENUref-Male.vcf.gz
    ../results/VCF_ChrFixed/TL2209399-TUref-Female.vcf.gz
    ../results/VCF_ChrFixed/TL2209400-TUref-Male.vcf.gz"

    bcftools merge --threads ${THREADS} \
        ../results/VCF_ChrFixed/TL2312073-163-4L.vcf.gz ${REF_PATHS} \
        -Ob -o ../results/VCF_Merged/TL2312073-163-4L.vcf.gz
    bcftools index --threads ${THREADS} -t \
        ../results/VCF_Merged/TL2312073-163-4L.merged.vcf.gz

For an example on all samples:

.. code-block:: shell

    THREADS=16
    INFILE_DIR="../results/VCF_ChrFixed/*gz"
    OUTFILE_DIR="../results/VCF_Merged/"
    REF_PATHS="../results/VCF_ChrFixed/TL2209397-ENUref-Female.vcf.gz
    ../results/VCF_ChrFixed/TL2209398-ENUref-Male.vcf.gz
    ../results/VCF_ChrFixed/TL2209399-TUref-Female.vcf.gz
    ../results/VCF_ChrFixed/TL2209400-TUref-Male.vcf.gz"

    find ${INFILE_DIR} | grep -v 'ref' | \
        while IFS="\t" read -r line; do
            in=$(echo $line)
            out=${OUTFILE_DIR}$(basename ${in})
            bcftools merge --threads ${THREADS} \
                ${in} ${REF_PATHS} -Ob -o ${out}
            bcftools index --threads ${THREADS} -t ${out}
        done

Filter out non-SNPs
+++++++++++++++++++

For one example:

.. code-block:: shell

    THREADS=16
    bcftools view --threads ${THREADS} -v snps \
        ../results/VCF_Merged/TL2312073-163-4L.vcf.gz \
        -Ob -o ../results/VCF_Snps/TL2312073-163-4L.vcf.gz
    bcftools index --threads ${THREADS} -t \
        ../results/VCF_SNPs/TL2312073-163-4L.vcf.gz

For an example on all samples:

.. code-block:: shell

    THREADS=16
    INFILE_DIR="../results/VCF_Merged/*gz"
    OUTFILE_DIR="../results/VCF_Snps/"

    find ${INFILE_DIR} | sort | \
        while IFS="\t" read -r line; do
            in=$(echo $line)
            out=${OUTFILE_DIR}$(basename ${in})
            bcftools view --threads ${THREADS} -v snps \
                ${in} -Ob -o ${out}
            bcftools index --threads ${THREADS} -t ${out}
        done

Obtain SNPs which occur in the sample but not the references
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. raw:: html

   <details>
   <summary><a>This unrefined implementation is preserved for the record.</a></summary>

Unrefined implementation (works only if grandparents have the reference allele)
*******************************************************************************

Ideally we only should include samples homozygous for the SNPs. However, false heterozygotes appear in the data due to reads with low mapping quality introducing false heterozygoisty. Therefore, this filter is intentionally relaxed to accommodate heterozygotes. The information can be verified by inspecting the ``bam`` file directly.

For one example:

.. code-block:: shell

    THREADS=16
    bcftools view --threads ${THREADS} \
        -i 'GT[0]="hom" && F_MISSING=0.8 || GT[0]="het" && F_MISSING=0.8' \
        ../results/VCF_Snps/TL2312073-163-4L.vcf.gz \
        -Ob -o ../results/VCF_Candidates/TL2312073-163-4L.vcf.gz

To verify that the conditional statements are formulated correctly and that this returns expected results:

.. code-block:: shell

    # both return the same line count of 1481032
    THREADS=16
    bcftools view --threads ${THREADS} \
        -i 'GT[0]="hom" && F_MISSING=0.8 || GT[0]="het" && F_MISSING=0.8' \
        ../results/VCF_Snps/TL2312073-163-4L.vcf.gz | \
        grep -v '#' | wc -l
    bcftools view --threads ${THREADS}\
        -i 'GT[0]="hom" && F_MISSING=0.8 || GT[0]="het" && F_MISSING=0.8' \
        ../results/VCF_Snps/TL2312073-163-4L.vcf.gz | \
        grep -Pc './.:.:.:.:.\t./.:.:.:.:.\t./.:.:.:.:.\t./.:.:.:.:.$'

For an example on all samples:

.. code-block:: shell

    THREADS=16
    INFILE_DIR="../results/VCF_Snps/*gz"
    OUTFILE_DIR="../results/VCF_Candidates/"

    find ${INFILE_DIR} | sort | \
        while IFS="\t" read -r line; do
            in=$(echo $line)
            out=${OUTFILE_DIR}$(basename ${in})
            bcftools view --threads ${THREADS} \
                -i 'GT[0]="hom" && F_MISSING=0.8 || GT[0]="het" && F_MISSING=0.8' \
                ${in} -Ob -o ${out}
            bcftools index --threads ${THREADS} -t ${out}
        done

.. note::
    To be refined. Cannot just obtain all missing snps in refs. e.g:
        - Reference allele may be A
        - Grandparents may be C
        - Mutant may be G

.. raw:: html

   </details>

Refined implementation (works if grandparents have the reference allele)
************************************************************************

``vcf`` file specifications:
- We have a multisample vcf file
- The first sample is the mutant
- The remaining 4 samples are the grandparent `F0` samples used as a reference
- Each sample is the result of variant calling against the primary reference genome

.. raw:: html

   <details>
   <summary><a>Original criteria is preserved for the record.</a></summary>

Original claimed criteria
^^^^^^^^^^^^^^^^^^^^^^^^^

Candidate criteria:
    1. The mutant must be covered to a depth of >= 1 read
    2. At least one of the unaffected sibling or other samples must be covered to a depth of >= 1 read
    3. If the mutant and sibling have a single allele, it must not be the same in both
    4. The mutant must have a majority allele that is:
        a. not the majority allele in any of the 'other' samples'''

Strict candidate criteria:
    1. Mutant is homozygous
    2. All of the 'other' samples are called
    3. Mutant allele is not present in any of the 'other' samples

Reworded criteria
^^^^^^^^^^^^^^^^^

.. raw:: html

   </details>

Candidate criteria:
    1. The mutant must be covered to a depth of >= 1 read
    2. At least one of the other samples must be covered to a depth of >= 1 read
    3. If the mutant and grandparent reference have a single allele, it must not be the same in both
    4. The mutant must have a majority allele that is not the majority allele in any of the F0 references

Strict candidate criteria:
    1. Find positions where F2 has a homozygous mutation
    2. None of the F0s have that mutation
    3. At least one F0 has coverage

Modular ``bcftools``-native implementation of *strict candidate criteria* filter for one sample:

.. code-block:: shell

    THREADS=16
    INFILE_DIR="../results/VCF_Snps/*gz"
    OUTFILE_DIR="../results/VCF_Candidates_Strict/"

    # NOTE: for debugging purposes only
    # are there 1/1 genotypes in the mutant
    bcftools query -f '%CHROM:%POS[\t%GT]\n' SAMPLE.vcf | \
        awk '$2=="1/1"' | wc -l

    # check the range of genotypes of the mutant
    bcftools query -f '%CHROM:%POS[\t%GT]\n' SAMPLE.vcf | \
        cut -f2 | sort | uniq -c

    # GT[0]="1/1" == mutant is homozygous
    #   position 1 is the sample, therefore mutation
    # GT[1+]!~"1" == mutant allele is absent in any references
    #   !~ == doesnt contain 1, 0/0 or ./ OK
    # FORMAT/DP[1+]>=1 == at least one read hits one reference
    #   at least 1 reference has read(s) mapped to the region
    bcftools filter --threads ${THREADS} \
        -i 'GT[0]="1/1" && GT[1]!~"1" && GT[2]!~"1" && GT[3]!~"1" && GT[4]!~"1" && (FORMAT/DP[1]>=1 || FORMAT/DP[2]>=1 || FORMAT/DP[3]>=1 || FORMAT/DP[4]>=1)' \
        ../results/VCF_Snps/TL2312073-163-4L.vcf.gz \
        -Ob -o ../results/VCF_Candidates_Strict/TL2312073-163-4L.vcf.gz
    bcftools index --threads ${THREADS} -t \
        ../results/VCF_Candidates_Strict/TL2312073-163-4L.vcf.gz

For an example on all samples:

.. code-block:: shell

    THREADS=16
    INFILE_DIR="../results/VCF_Snps/*gz"
    OUTFILE_DIR="../results/VCF_Candidates_Strict/"

    find ${INFILE_DIR} | sort | \
        while IFS="\t" read -r line; do
            in=$(echo $line)
            out=${OUTFILE_DIR}$(basename ${in})
            bcftools filter --threads ${THREADS} \
                -i 'GT[0]="1/1" && GT[1]!~"1" && GT[2]!~"1" && GT[3]!~"1" && GT[4]!~"1" && (FORMAT/DP[1]>=1 || FORMAT/DP[2]>=1 || FORMAT/DP[3]>=1 || FORMAT/DP[4]>=1)' \
                ${in} -Ob -o ${out}
            bcftools index --threads ${THREADS} -t ${out}
        done

Screen impact of SNP
++++++++++++++++++++

We run ``snpEff`` to screen for SNPs which are predicted to be impactful. Annotations are added to the data in the process.

.. warning:: 
    ``snpEff`` runs on ``java``, and as a result may run into memory and other issues. You may need to set ``java`` memory flags accordingly.

    .. code-block:: shell

        # adjust max values accordingly
        export JAVA_OPTIONS="-Xms512m -Xmx16g"

.. code-block:: shell

    snpEff -i vcf -v Zv9.75 TL2312073-163-4L.vcf.gz > bar.vcf

Filter by SNP impact
++++++++++++++++++++


TEMP
++++

Quick, unrefined way to pull out strict candidates for reference only. Not to be considered as final results!

.. code-block:: shell

    outdir="../results/tmp/"
    pattern='#|./.:.:.:.:.\t./.:.:.:.:.\t./.:.:.:.:.\t./.:.:.:.:.$'
    metadata="../data/samplesheet_original.tsv"
    mkdir -p ${outdir}
    cut -f6 ${metadata} | tail -n +6 | \
        while IFS="\t" read -r line; do
            in=$(echo $line | cut -d ' ' -f3)
            out="${outdir}$(basename ${in} | cut -d '-' -f1-3).vcf"
            bgzip -cd ${in} | grep -P ${pattern} > ${out}
        done