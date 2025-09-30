#!/usr/bin/env nextflow

/*
 * EMUMADZ Pipeline - BAM to SNP for F0-F2 comparisons
 * 
 * A Nextflow pipeline for variant calling and SNP identification in zebrafish
 * mutagenesis experiments using ENU-induced mutations.
 * 
 * Based on the shell scripts from case_study_1.rst documentation.
 * 
 * Authors: Tyrone Chen, Richard Lupat, Michelle Meier, Maia Zethoven, 
 *          Greg Baillie, Scott Paterson, Oguzhan Baltaci, Cas Simons, 
 *          Jason Li, Benjamin Hogan
 * 
 * License: MIT
 */

nextflow.enable.dsl = 2

// Pipeline parameters
params.samplesheet = null
params.reference_fa = null
params.ref_fixed_fa = null
params.chrom_map = null
params.outdir = './results'
params.threads = 16
params.gatk_mem = '128g'
params.java_options = '-Xms512m -Xmx16g'
params.vep_assembly = 'Zv9'
params.vep_version = '79'
params.vep_cache = "${HOME}/.vep"
params.vep_buffer = '8192'
params.snpeff_config = null
params.snpeff_genome = 'Zv9.75'
params.min_mapping_quality = 20
params.min_base_quality = 20
params.stand_call_conf = 30
params.max_reads_per_alignment_start = 50
params.max_mnp_distance = 1
params.apply_custom_filters = false
params.gatk_filters_config = null
params.help = false

// Print help message
if (params.help) {
    log.info """
    EMUMADZ Pipeline - BAM to SNP for reference-mutant comparisons
    
    Usage:
        nextflow run emumadz_pipeline.nf --samplesheet samplesheet.tsv --reference_fa reference.fna --ref_fixed_fa reference_fixed.fna --chrom_map chrom_table.tsv
    
    Required parameters:
        --samplesheet          Path to samplesheet with sample information
        --reference_fa         Path to reference genome FASTA file
        --ref_fixed_fa         Path to fixed reference genome FASTA file
        --chrom_map            Path to chromosome mapping file
    
    Optional parameters:
        --outdir               Output directory (default: ./results)
        --threads              Number of threads (default: 16)
        --gatk_mem             GATK memory allocation (default: 64g)
        --apply_custom_filters Apply stringent GATK filters (default: false)
        --gatk_filters_config  Path to GATK filters configuration file (optional)
        --snpeff_config        Path to snpEff configuration file (optional)
        --help                 Show this help message
    """
    exit 0
}

// GATK filter parameters
params.qd_threshold = 2.0
params.qual_threshold = 30.0
params.sor_threshold = 3.0
params.fs_threshold = 60.0
params.mq_threshold = 40.0
params.mqranksum_threshold = -12.5
params.readposranksum_threshold = -8.0

// Load custom GATK filters config if provided
if (params.gatk_filters_config) {
    includeConfig params.gatk_filters_config
}

// Process definitions

process setup_metadata {
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/metadata", mode: 'copy'
    
    input:
    path ref_fa
    path ref_fixed_fa
    
    output:
    path ref_fa, emit: ref_fa_indexed
    path "${ref_fa}.fai", emit: ref_fai
    path "${ref_fa.baseName}.dict", emit: ref_dict
    path ref_fixed_fa, emit: ref_fixed_indexed
    path "${ref_fixed_fa}.fai", emit: ref_fixed_fai
    path "chromosomes.txt", emit: chromosomes
    
    script:
    """
    # Create FASTA index for reference
    samtools faidx ${ref_fa}
    
    # Create sequence dictionary
    gatk CreateSequenceDictionary -R ${ref_fa}
    
    # Create FASTA index for fixed reference
    samtools faidx ${ref_fixed_fa}
    
    # Create chromosome list
    grep -P "^chr[0-9]+\t" ${ref_fixed_fa}.fai > chromosomes.txt
    """
}

process call_variants {
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/variants", mode: 'copy'
    cpus params.threads
    memory params.gatk_mem
    
    input:
    tuple val(sample_id), path(bam_file), path(bam_index), val(sample_type)
    path ref_fa
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), val(sample_type), emit: vcf_with_index
    
    script:
    gatk_java_opts = "-Xmx${params.gatk_mem} -XX:+UseParallelGC"
    """
    # Create FASTA index for reference
    samtools faidx ${ref_fa}
    
    # Create sequence dictionary
    gatk CreateSequenceDictionary -R ${ref_fa}

    gatk --java-options "${gatk_java_opts}" HaplotypeCaller \\
        -R ${ref_fa} \\
        -I ${bam_file} \\
        -O ${sample_id}.vcf.gz \\
        --minimum-mapping-quality ${params.min_mapping_quality} \\
        --min-base-quality-score ${params.min_base_quality} \\
        --stand-call-conf ${params.stand_call_conf} \\
        --max-reads-per-alignment-start ${params.max_reads_per_alignment_start} \\
        --max-mnp-distance ${params.max_mnp_distance} \\
        --native-pair-hmm-threads ${params.threads} \\
        --verbosity INFO
    """
}

process filter_chromosomes {
    container 'staphb/bcftools:1.22'
    publishDir "${params.outdir}/filtered", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index), val(sample_type)
    path chrom_map
    path chromosomes
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), val(sample_type), emit: vcf_with_index
    
    script:
    """
    # Create chromosomes_only.txt
    for i in {1..25}; do 
        echo -e "chr\${i}\t1\t1111111111111000000000000"
    done > chromosomes_only.txt
    
    # Rename chromosomes and filter
    bcftools annotate --threads ${params.threads} \\
        --rename-chr ${chrom_map} \\
        ${vcf_file} | \\
    grep -v '##contig=<ID=FR' > ${sample_id}_reannot.vcf
    
    # Filter to main chromosomes only
    bcftools view -t \$(cut -f1 ${chromosomes} | paste -sd,) \\
        ${sample_id}_reannot.vcf \\
        --threads ${params.threads} \\
        -o ${sample_id}_chr_temp.vcf
    
    # Clean up header and compress
    bcftools reheader -f ${chromosomes} \\
        ${sample_id}_chr_temp.vcf | \\
    bgzip -@ ${params.threads} -c > ${sample_id}.vcf.gz
    
    bcftools index --threads ${params.threads} ${sample_id}.vcf.gz
    
    # Clean up intermediate files
    rm -f ${sample_id}_reannot.vcf ${sample_id}_chr_temp.vcf
    """
}

process apply_custom_filters {
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/custom_filtered", mode: 'copy'
    cpus params.threads
    memory params.gatk_mem
    
    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index), val(sample_type)
    path ref_fa
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), val(sample_type), emit: vcf_with_index
    
    script:
    gatk_java_opts = "-Xmx${params.gatk_mem} -XX:+UseParallelGC"
    """
    # Apply GATK custom filters with configurable thresholds
    gatk --java-options "${gatk_java_opts}" VariantFiltration \\
        -R ${ref_fa} \\
        -V ${vcf_file} \\
        -O ${sample_id}_filtered.vcf \\
        --filter-expression "QD < ${params.qd_threshold}" --filter-name "QD${params.qd_threshold}" \\
        --filter-expression "QUAL < ${params.qual_threshold}" --filter-name "QUAL${params.qual_threshold}" \\
        --filter-expression "SOR > ${params.sor_threshold}" --filter-name "SOR${params.sor_threshold}" \\
        --filter-expression "FS > ${params.fs_threshold}" --filter-name "FS${params.fs_threshold}" \\
        --filter-expression "MQ < ${params.mq_threshold}" --filter-name "MQ${params.mq_threshold}" \\
        --filter-expression "MQRankSum < ${params.mqranksum_threshold}" --filter-name "MQRankSum${params.mqranksum_threshold}" \\
        --filter-expression "ReadPosRankSum < ${params.readposranksum_threshold}" --filter-name "ReadPosRankSum${params.readposranksum_threshold}"
    
    # Select only PASS variants
    gatk --java-options "${gatk_java_opts}" SelectVariants \\
        -R ${ref_fa} \\
        -V ${sample_id}_filtered.vcf \\
        -O ${sample_id}.vcf \\
        --exclude-filtered
    
    # Compress and index
    bgzip ${sample_id}.vcf
    bcftools index ${sample_id}.vcf.gz
    
    # Clean up intermediate files
    rm -f ${sample_id}_filtered.vcf
    """
}

process normalise_vcf {
    container 'staphb/bcftools:1.22'
    publishDir "${params.outdir}/normalized", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index), val(sample_type)
    path ref_fixed_fa
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), val(sample_type), emit: vcf_with_index
    
    script:
    """
    # Normalize variants: decompose complex variants and left-align indels
    bcftools norm -m-both -f ${ref_fixed_fa} ${vcf_file} \\
        --threads ${params.threads} --write-index -Ob -o ${sample_id}.vcf.gz
    """
}

process merge_vcf {
    container 'staphb/bcftools:1.22'
    publishDir "${params.outdir}/merged", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index), val(sample_type), path(reference_vcfs)
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf_with_index
    
    script:
    """
    # Create VCF list for this mutant + all references
    echo "${vcf_file}" > vcf_list.txt
    for ref_vcf in ${reference_vcfs}; do
        echo "\${ref_vcf}" >> vcf_list.txt
    done
    
    # Merge this mutant with all references
    bcftools merge -l vcf_list.txt --threads ${params.threads} --write-index \\
        -Ob -o ${sample_id}.vcf.gz
    """
}

process filter_snps {
    container 'staphb/bcftools:1.22'
    publishDir "${params.outdir}/snps", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index)
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf_with_index
    
    script:
    """
    # Filter for SNPs only
    bcftools view -i 'TYPE="snp"' --threads ${params.threads} \\
        ${vcf_file} \\
        --write-index -Ob -o ${sample_id}.vcf.gz
    """
}

process find_candidates {
    container 'staphb/bcftools:1.22'
    publishDir "${params.outdir}/candidates", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index), val(num_refs)
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf_with_index
    
    script:
    """
    # Build depth filter: mutant AND at least one ref has >=1 read depth
    DP_FILTER="FORMAT/DP[0]>=1"
    for ((i=1; i<=${num_refs}; i++)); do
        if [ \$i -eq 1 ]; then
            DP_FILTER="\${DP_FILTER} && (FORMAT/DP[\${i}]>=1"
        else
            DP_FILTER="\${DP_FILTER} || FORMAT/DP[\${i}]>=1"
        fi
    done
    DP_FILTER="\${DP_FILTER})"
    
    # Build allele frequency filter
    AF_FILTER='FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85'
    for ((i=1; i<=${num_refs}; i++)); do
        AF_FILTER="\${AF_FILTER} && (FORMAT/AD[\${i}:1]=\".\" || FORMAT/AD[\${i}:0]=\".\" || FORMAT/AD[\${i}:1]/(FORMAT/AD[\${i}:0]+FORMAT/AD[\${i}:1]) <= 0.15)"
    done
    
    # Combine filters
    STRICT_FILTER="(\${DP_FILTER}) && (\${AF_FILTER})"
    
    # Apply filters
    bcftools filter -i "\${STRICT_FILTER}" --threads ${params.threads} --write-index \\
        ${vcf_file} -Ob -o ${sample_id}.vcf.gz
    """
}

process annotate_vep {
    container 'ensemblorg/ensembl-vep'
    publishDir "${params.outdir}/vep", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index)
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf_with_index
    path "*.html", emit: stats_html
    path "*.txt", emit: stats_text
    
    script:
    """
    # Run VEP
    vep \\
        --cache \\
        --dir_cache ${params.vep_cache} \\
        --species danio_rerio \\
        --assembly ${params.vep_assembly} \\
        --cache_version ${params.vep_version} \\
        --offline \\
        --regulatory \\
        --vcf \\
        --variant_class \\
        --fork ${params.threads} \\
        --buffer_size ${params.vep_buffer} \\
        --stats_html \\
        --stats_text \\
        --force_overwrite \\
        --compress_output bgzip \\
        --input_file ${vcf_file} \\
        --output_file ${sample_id}.vcf.gz
    
    bcftools index --threads ${params.threads} ${sample_id}.vcf.gz
    """
}

process annotate_snpeff {
    container 'staphb/snpeff:5.2f'
    publishDir "${params.outdir}/snpeff", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index)
    path snpeff_config
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf_with_index
    path "*.csv", emit: stats_csv
    path "*.fa", emit: fasta_prot
    path "*.html", emit: stats_html
    
    script:
    config_opt = snpeff_config ? "-c ${snpeff_config}" : ""
    """
    export JAVA_OPTIONS="${params.java_options}"
    
    # Run snpEff
    snpEff \\
        ${config_opt} \\
        -csvStats ${sample_id}.csv \\
        -i vcf \\
        -o vcf \\
        -fastaProt ${sample_id}.fa \\
        -s ${sample_id}.html \\
        ${params.snpeff_genome} \\
        ${vcf_file} > ${sample_id}.vcf
    
    bgzip -@ ${params.threads} ${sample_id}.vcf
    bcftools index --threads ${params.threads} ${sample_id}.vcf.gz
    """
}

process combine_annotations {
    container 'staphb/bcftools:1.22'
    publishDir "${params.outdir}/annotated", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(vep_vcf), path(vep_index)
    tuple val(sample_id), path(snpeff_vcf), path(snpeff_index)
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf_with_index
    
    script:
    """
    # Merge VEP and snpEff annotations
    bcftools merge --threads ${params.threads} \\
        --force-samples \\
        ${vep_vcf} \\
        ${snpeff_vcf} \\
        --write-index -Ob -o ${sample_id}.vcf.gz
    """
}

process prepare_visualization {
    container 'tyronechen/pysam_pandas'
    publishDir "${params.outdir}/visualization", mode: 'copy'
    cpus params.threads
    
    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index)
    path samplesheet
    path chrom_map
    path parse_script
    
    output:
    path "${sample_id}.csv", emit: csv
    path "${sample_id}.json", emit: json
    path "${sample_id}_coverage.csv", emit: coverage
    path "${sample_id}/", emit: bam_subsets
    
    script:
    """
    # Create output directory for BAM subsets
    mkdir -p ${sample_id}
    
    # Run the parse_vcf.py script
    python ${parse_script} \\
        ${vcf_file} \\
        ${samplesheet} \\
        ${sample_id}.csv \\
        --subset_path ${sample_id} \\
        --subset_workers ${params.threads} \\
        --coverage_report ${sample_id}_coverage.csv \\
        --chrom_mapping ${chrom_map} \\
        --force_overwrite
    """
}

// Main workflow
workflow {
    // Validate required parameters
    if (!params.samplesheet || !params.reference_fa || !params.ref_fixed_fa || !params.chrom_map) {
        log.error "Missing required parameters. Use --help for usage information."
        exit 1
    }

    // Create file objects after validation
    ref_fa = file(params.reference_fa)
    ref_fixed_fa = file(params.ref_fixed_fa)
    chrom_map = file(params.chrom_map)
    samplesheet = file(params.samplesheet)
    parse_script = file("${projectDir}/bin/parse_vcf.py")
    
    // Parse samplesheet and create channels
    all_samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row -> [row.sample_identity, row.alignment_file, row.sample_type] }
        .map { sample_id, bam_path, sample_type -> 
            [sample_id, file(bam_path), file("${bam_path}.bai"), sample_type]
        }
    
    // Setup metadata
    setup_metadata(ref_fa, ref_fixed_fa)
    
    // Variant calling for all samples
    call_variants(all_samples, ref_fa)
    
    // Chromosome filtering for all samples
    filter_chromosomes(call_variants.out.vcf_with_index, chrom_map, setup_metadata.out.chromosomes)
    
    // Apply conditional filtering
    if (params.apply_custom_filters) {
        filtered_vcf = apply_custom_filters(filter_chromosomes.out.vcf_with_index, ref_fa)
    } else {
        filtered_vcf = filter_chromosomes.out.vcf_with_index
    }
    
    // Normalize VCFs
    normalise_vcf(filtered_vcf, ref_fixed_fa)
    
    // Separate mutant and reference samples
    mutant_vcfs = normalise_vcf.out.vcf_with_index.filter { sample_id, vcf, index, sample_type -> sample_type == 'mutant' }
    reference_vcfs = normalise_vcf.out.vcf_with_index
        .filter { sample_id, vcf, index, sample_type -> sample_type == 'reference' }
        .map { sample_id, vcf, index, sample_type -> vcf }
        .collect()
    
    // Combine each mutant with collected references
    mutant_with_refs = mutant_vcfs.combine(reference_vcfs)
    
    // Merge mutant samples with references  
    merge_vcf(mutant_with_refs)
    
    // Filter for SNPs only
    filter_snps(merge_vcf.out.vcf_with_index)
    
    // Count references for filtering
    ref_count = reference_vcfs.map { refs -> refs.size() }
    
    // Find candidate SNPs
    find_candidates(filter_snps.out.vcf_with_index.combine(ref_count))
    
    // Annotate with VEP
    annotate_vep(find_candidates.out.vcf_with_index)
    
    // Annotate with snpEff
    snpeff_config_file = params.snpeff_config ? file(params.snpeff_config) : []
    annotate_snpeff(find_candidates.out.vcf_with_index, snpeff_config_file)
    
    // Combine annotations
    combine_annotations(annotate_vep.out.vcf_with_index, annotate_snpeff.out.vcf_with_index)
    
    // Prepare visualization data
    prepare_visualization(combine_annotations.out.vcf_with_index, samplesheet, chrom_map, parse_script)
}

// Workflow completion
workflow.onComplete {
    log.info """
    ================================================
    EMUMADZ Pipeline completed successfully!
    ================================================
    
    Results are available in: ${params.outdir}
    
    Pipeline steps completed:
    - Variant calling
    - Chromosome filtering
    ${params.apply_custom_filters ? '- Custom filtering' : ''}
    - VCF normalization
    - Sample merging
    - SNP filtering
    - Candidate identification
    - VEP annotation
    - snpEff annotation
    - Annotation combination
    - Visualization preparation
    """
}

workflow.onError {
    log.error """
    ================================================
    EMUMADZ Pipeline failed!
    ================================================
    
    Error: ${workflow.errorMessage}
    """
}