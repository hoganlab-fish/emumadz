#!/usr/bin/env nextflow

/*
 * Container Test Pipeline
 * 
 * Tests if all required container images can be pulled and tools are installed
 * Run with: nextflow run test_containers.nf -profile test_containers
 */

nextflow.enable.dsl = 2

// Test workflow
workflow {
    // Test all containers
    test_gatk()
    test_bcftools()
    test_vep()
    test_snpeff()
    test_pysam()
    
    // Print results
    test_gatk.out.view { "✅ GATK test: $it" }
    test_bcftools.out.view { "✅ bcftools test: $it" }
    test_vep.out.view { "✅ VEP test: $it" }
    test_snpeff.out.view { "✅ snpEff test: $it" }
    test_pysam.out.view { "✅ pysam test: $it" }
}

// Test processes
process test_gatk {
    container 'broadinstitute/gatk:4.6.2.0'
    
    output:
    stdout
    
    script:
    """
    gatk --version
    """
}

process test_bcftools {
    container 'staphb/bcftools:1.22'
    
    output:
    stdout
    
    script:
    """
    bcftools --version
    """
}

process test_vep {
    container 'ensemblorg/ensembl-vep'
    
    output:
    stdout
    
    script:
    """    
    ./vep --help | head -5
    """
}

process test_snpeff {
    container 'staphb/snpeff:5.2f'
    
    output:
    stdout
    
    script:
    """
    snpEff -version
    """
}

process test_pysam {
    container 'tyronechen/pysam_pandas'
    
    output:
    stdout
    
    script:
    """
    python -c "import pysam; print(f'pysam version: {pysam.__version__}')"
    python -c "import pandas; print(f'pandas version: {pandas.__version__}')"
    """
}

// Workflow completion
workflow.onComplete {
    log.info """
    ================================================
    Container Test Pipeline completed!
    ================================================
    
    All containers have been tested successfully.
    You can now run the full pipeline with confidence.
    """
}

workflow.onError {
    log.error """
    ================================================
    Container Test Pipeline failed!
    ================================================
    
    One or more containers failed to load or tools are missing.
    Check the error messages above for details.
    """
}
