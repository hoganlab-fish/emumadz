#!/usr/bin/env python3
# take a vcf, extract seqs around each SNP, parse annotations for visualisation
import argparse
import json
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
import pandas as pd
import pysam

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_csq_field(csq_string: str) -> list:
    """Parse VEP CSQ field"""
    if not csq_string:
        return list()
    
    annotations = list()
    for ann in csq_string.split(','):
        fields = ann.split('|')
        if len(fields) >= 4:
            annotations.append({
                'consequence': fields[1] if len(fields) > 1 else '',
                'impact': fields[2] if len(fields) > 2 else '',
                'symbol': fields[3] if len(fields) > 3 else '',
                'gene': fields[4] if len(fields) > 4 else '',
                'transcript': fields[6] if len(fields) > 6 else '',
                'protein_pos': fields[9] if len(fields) > 9 else '',
                'amino_acids': fields[10] if len(fields) > 10 else '',
            })
    return annotations

def parse_ann_field(ann_string: str) -> list:
    """Parse SnpEff ANN field"""
    if not ann_string:
        return list()
    
    annotations = list()
    for ann in ann_string.split(','):
        fields = ann.split('|')
        if len(fields) >= 4:
            annotations.append({
                'effect': fields[1] if len(fields) > 1 else '',
                'impact': fields[2] if len(fields) > 2 else '',
                'gene': fields[3] if len(fields) > 3 else '',
                'transcript': fields[6] if len(fields) > 6 else '',
                'protein_pos': fields[9] if len(fields) > 9 else '',
                'amino_acids': fields[10] if len(fields) > 10 else '',
            })
    return annotations

def calculate_allele_freq(genotype: str) -> float:
    """Calculate allele frequency from genotype"""
    if not genotype or genotype == './.':
        return None, None
    
    gt_calls = str(genotype).split('/') if '/' in str(genotype) else str(genotype).split('|')
    ref_count = sum(1 for call in gt_calls if call == '0')
    alt_count = sum(1 for call in gt_calls if call != '0' and call != '.')
    total = ref_count + alt_count
    
    if total == 0:
        return None, None
    
    return alt_count / total, f"{alt_count}/{total}"

def get_impact_color(impact: str) -> str:
    """Get color coding for impact"""
    colors = {
        'HIGH': '#FF0000',
        'MODERATE': '#FF8C00', 
        'LOW': '#32CD32',
        'MODIFIER': '#87CEEB'
    }
    return colors.get(impact.upper(), '#CCCCCC')

def create_bam_subset(bam_path: str, variant_positions: list, output_path: str, padding: int=50, max_workers: int=4, chrom_mapping: dict=None):
    """
    Create a subset BAM file containing reads around variant positions
    
    Args:
        bam_path: Path to original BAM file
        variant_positions: List of tuples (chromosome, position)
        output_path: Path for subset BAM output
        padding: Base pairs to extend around each variant (default: 50bp)
        chrom_mapping: Dict mapping VCF chroms to BAM chroms
    """
    try:
        # ensure output directory exists
        output_dir = os.path.dirname(output_path)
        if output_dir:  # Only create directory if path contains one
            os.makedirs(output_dir, exist_ok=True)
        
        # open input BAM
        with pysam.AlignmentFile(bam_path, "rb", threads=max_workers) as inbam:
            
            # Create modified header if chromosome mapping provided
            header = inbam.header.to_dict()
            if chrom_mapping:
                # Create reverse mapping (BAM -> VCF)
                reverse_mapping = {v: k for k, v in chrom_mapping.items()}
                
                # Update sequence names in header
                if 'SQ' in header:
                    for sq in header['SQ']:
                        old_name = sq['SN']
                        if old_name in reverse_mapping:
                            sq['SN'] = reverse_mapping[old_name]
            
            # create output BAM with modified header
            with pysam.AlignmentFile(output_path, "wb", header=pysam.AlignmentHeader.from_dict(header), threads=max_workers) as outbam:
                
                # track processed regions to avoid duplicates
                processed_regions = set()
                
                for chrom, pos in variant_positions:
                    # Map chromosome name if mapping provided
                    bam_chrom = chrom_mapping.get(chrom, chrom) if chrom_mapping else chrom
                    
                    # define region with padding
                    start = max(0, pos - padding)
                    end = pos + padding
                    region_key = (bam_chrom, start, end)
                    
                    # skip if we've already processed this region
                    if region_key in processed_regions:
                        continue
                    processed_regions.add(region_key)
                    
                    try:
                        # fetch reads in region using mapped chromosome name
                        for read in inbam.fetch(bam_chrom, start, end):
                            # Write read directly - chromosome renaming handled by header
                            outbam.write(read)
                            
                    except Exception as e:
                        logger.warning(f"Could not fetch reads for {bam_chrom}:{start}-{end}: {e}")
                        continue
        
        # Index the subset BAM
        try:
            # Sort the BAM file first
            sorted_path = output_path.replace('.bam', '_sorted.bam')
            pysam.sort("-o", sorted_path, output_path, threads=max_workers)
            
            # Replace original with sorted version
            os.rename(sorted_path, output_path)
            
            # Index the sorted BAM
            pysam.index(output_path)
            logger.info(f"Created, sorted and indexed subset BAM: {output_path}")
            return True
        except Exception as e:
            logger.error(f"Failed to sort/index {output_path}: {e}")
            return False
            
    except Exception as e:
        logger.error(f"Failed to create subset BAM {output_path}: {e}")
        return False

def process_sample_variants(sample_data: dict, padding: int=50, max_workers: int=4, chrom_mapping: dict=None) -> bool:
    """
    Process all variants for a single sample and create subset BAM
    
    Args:
        sample_data: Dictionary containing sample info and variants with output_path as full BAM path
        padding: Base pairs around variants
    """
    sample_id = sample_data['sample_id']
    bam_path = sample_data['bam_path']
    output_path = sample_data['output_path']  # Now expects full BAM path
    variants = sample_data['variants']
    
    # check if BAM file exists
    if not os.path.exists(bam_path):
        logger.error(f"BAM file not found: {bam_path}")
        return False
    
    # check if BAM is indexed
    bam_index = bam_path + ".bai"
    if not os.path.exists(bam_index):
        try:
            logger.info(f"Indexing BAM file: {bam_path}")
            pysam.index(bam_path, threads=max_workers)
        except Exception as e:
            logger.error(f"Failed to index BAM {bam_path}: {e}")
            return False
    
    # extract variant positions
    variant_positions = [(row['chromosome'], row['position']) for _, row in variants.iterrows()]
    logger.info(f"Processing {len(variant_positions)} variants for sample {sample_id}")    
    success = create_bam_subset(bam_path, variant_positions, output_path, padding, max_workers, chrom_mapping)
    
    if success:
        logger.info(f"Successfully created subset BAM for {sample_id}")
    else:
        logger.error(f"Failed to create subset BAM for {sample_id}")
    
    return success

def validate_bam_subset(subset_path: str, expected_variants: list, chrom_mapping: dict=None, max_workers: int=4) -> dict:
    """
    Validate that a BAM subset contains reads covering expected variant positions
    
    Args:
        subset_path: Path to subset BAM file
        expected_variants: List of (chromosome, position) tuples
    """
    if not os.path.exists(subset_path):
        return False, "Subset BAM file not found"
    
    try:
        with pysam.AlignmentFile(subset_path, "rb", threads=max_workers) as bam:
            coverage_info = {}
            
            for chrom, pos in expected_variants:
                # Map chromosome name if mapping provided
                # print(chrom, pos)
                # bam_chrom = chrom_mapping.get(chrom, chrom) if chrom_mapping else chrom
                # print(bam_chrom)
                read_count = sum(1 for _ in bam.fetch(chrom, pos-1, pos+1))
                coverage_info[f"{chrom}:{pos}"] = read_count
            
            total_reads = sum(coverage_info.values())
            covered_variants = sum(1 for count in coverage_info.values() if count > 0)
            
            return True, {
                'total_reads': total_reads,
                'covered_variants': covered_variants,
                'total_variants': len(expected_variants),
                'coverage_detail': coverage_info
            }
            
    except Exception as e:
        return False, f"Error validating BAM subset: {e}"

def generate_coverage_report(variants_df: pd.DataFrame, subset_path: str, chrom_mapping: dict=None, max_workers: int=4) -> pd.DataFrame:
    """
    Generate coverage report for all BAM subsets
    """
    report = list()
    
    sample_groups = variants_df.groupby('sample_id')
    
    for sample_id, group in sample_groups:
        variant_positions = [(row['chromosome'], row['position']) for _, row in group.iterrows()]
        
        is_valid, result = validate_bam_subset(subset_path, variant_positions, chrom_mapping, max_workers)
        
        if is_valid:
            report.append({
                'sample_id': sample_id,
                'subset_exists': True,
                'total_variants': result['total_variants'],
                'covered_variants': result['covered_variants'],
                'coverage_rate': result['covered_variants'] / result['total_variants'] if result['total_variants'] > 0 else 0,
                'total_reads': result['total_reads']
            })
        else:
            report.append({
                'sample_id': sample_id,
                'subset_exists': False,
                'error': result,
                'total_variants': len(variant_positions),
                'covered_variants': 0,
                'coverage_rate': 0,
                'total_reads': 0
            })
    
    return pd.DataFrame(report)

def process_vcf(vcf_file, samplesheet_file, output_file, subset_path=None, subset_padding=50, max_workers: int=4, chrom_mapping_file=None):
    """Process VCF and create analysis table"""
    
    # load chromosome mapping if provided
    chrom_mapping = None
    if chrom_mapping_file:
        mapping_df = pd.read_csv(chrom_mapping_file, sep='\t', header=None, names=['bam_chrom', 'vcf_chrom'])
        chrom_mapping = dict(zip(mapping_df['vcf_chrom'], mapping_df['bam_chrom']))
        logger.info(f"Loaded chromosome mapping for {len(chrom_mapping)} chromosomes")
    
    # load samplesheet
    samplesheet = pd.read_csv(samplesheet_file, sep='\t')
    sample_info = dict(zip(samplesheet['sample_id'], 
                       zip(samplesheet['bam_path'],
                       samplesheet['sample_type'])))
    # parse VCF
    vcf_reader = pysam.VariantFile(vcf_file, threads=max_workers)
    variants = list()
    
    for record in vcf_reader:
        # Get basic variant info
        variant_id = f"{record.chrom}:{record.pos}:{record.ref}:{','.join([str(alt) for alt in record.alts])}"
        
        # Parse VEP annotations
        vep_annotations = list()
        if 'CSQ' in record.info:
            csq_data = record.info['CSQ'][0] if isinstance(record.info['CSQ'], tuple) else record.info['CSQ']
            vep_annotations = parse_csq_field(csq_data)
        
        # Parse SnpEff annotations  
        snpeff_annotations = list()
        if 'ANN' in record.info:
            ann_data = record.info['ANN'][0] if isinstance(record.info['ANN'], tuple) else record.info['ANN']
            snpeff_annotations = parse_ann_field(ann_data)
        
        # Find first mutant sample (skip controls for SNPs)
        mutant_sample = None
        for sample_name in record.samples:
            # Try exact match first
            if sample_name in sample_info:
                bam_path, sample_type = sample_info[sample_name]
            else:
                # Try fuzzy match - find samplesheet entry that's a prefix of VCF sample
                matched = False
                for sheet_sample in sample_info:
                    if sample_name.startswith(sheet_sample):
                        bam_path, sample_type = sample_info[sheet_sample]
                        matched = True
                        break
                
                if not matched:
                    continue
            
            # Process only the first mutant sample
            if sample_type.lower() == 'mutant':
                mutant_sample = sample_name
                break
        
        if not mutant_sample:
            continue
            
        # Process only the mutant sample
        for sample_name in [mutant_sample]:
            # Get sample info (already validated above)
            if sample_name in sample_info:
                bam_path, sample_type = sample_info[sample_name]
            else:
                for sheet_sample in sample_info:
                    if sample_name.startswith(sheet_sample):
                        bam_path, sample_type = sample_info[sheet_sample]
                        break
                
            sample = record.samples[sample_name]
            
            # Get genotype info
            gt = sample.get('GT', (None, None))
            dp = sample.get('DP', 0)
            gq = sample.get('GQ', 0)
            
            # Format genotype
            gt_str = '/'.join([str(g) if g is not None else '.' for g in gt]) if gt else './.'
            
            # Calculate allele frequency
            af_freq, af_count = calculate_allele_freq(gt_str)
            
            # Get top consequence from each tool
            vep_top = vep_annotations[0] if vep_annotations else {}
            snpeff_top = snpeff_annotations[0] if snpeff_annotations else {}
            
            # Get full annotations as JSON strings
            vep_full = json.dumps(vep_annotations) if vep_annotations else ''
            snpeff_full = json.dumps(snpeff_annotations) if snpeff_annotations else ''
            
            variant_data = {
                'variant_id': variant_id,
                'chromosome': record.chrom,
                'position': record.pos,
                'ref': record.ref,
                'alt': ','.join([str(alt) for alt in record.alts]),
                'sample_id': sample_name,
                'sample_type': sample_type,
                'bam_path': bam_path,
                'bam_subset_path': subset_path if subset_path else "",
                'genotype': gt_str,
                'depth': dp if dp else 0,
                'quality': gq if gq else 0,
                'allele_freq': af_freq,
                'allele_count': af_count,
                'vep_consequence': vep_top.get('consequence', ''),
                'vep_impact': vep_top.get('impact', ''),
                'vep_gene': vep_top.get('symbol', ''),
                'vep_protein_change': vep_top.get('amino_acids', ''),
                'snpeff_effect': snpeff_top.get('effect', ''),
                'snpeff_impact': snpeff_top.get('impact', ''),
                'snpeff_gene': snpeff_top.get('gene', ''),
                'snpeff_protein_change': snpeff_top.get('amino_acids', ''),
                'vep_full_annotations': vep_full,
                'snpeff_full_annotations': snpeff_full,
                'custom_notes': '',
                'vep_impact_color': get_impact_color(vep_top.get('impact', '')),
                'snpeff_impact_color': get_impact_color(snpeff_top.get('impact', ''))
            }
            
            variants.append(variant_data)
    
    # Create DataFrame and save
    df = pd.DataFrame(variants)
    df.to_csv(output_file, index=False)
    
    # Save as JSON for web interface
    json_output = output_file.replace('.csv', '.json')
    df.to_json(json_output, orient='records', indent=2)
    
    # Create BAM subsets if requested
    if subset_path:
        logger.info("Creating BAM subsets...")
        # For single mutant sample, pass the full path directly
        if len(df) > 0:
            sample_data = {
                'sample_id': df.iloc[0]['sample_id'],
                'bam_path': df.iloc[0]['bam_path'], 
                'output_path': subset_path,
                'variants': df
            }
            success = process_sample_variants(sample_data, subset_padding, max_workers,chrom_mapping)
        else:
            print("No variants to process")
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Parse annotated VCF for variant analysis with optional BAM subsetting'
        )
    parser.add_argument('vcf', type=str,
                        help='Merged VCF file with VEP and SnpEff annotations')
    parser.add_argument('samplesheet', type=str,
                        help='Sample metadata file')
    parser.add_argument('outfile_path', type=str,
                        help='Output CSV file')
    parser.add_argument('--subset_path', type=str, default=None, 
                        help='Path to subset BAM')
    parser.add_argument('--subset_padding', type=int, default=25, 
                        help='Base pairs around variants for BAM subsets')
    parser.add_argument('--subset_workers', type=int, default=4, 
                        help='Number of parallel workers for BAM subsetting')
    parser.add_argument('--chrom_mapping', type=str, default=None,
                        help='Two-column TSV mapping VCF::BAM chromosome names')
    parser.add_argument('--coverage_report', type=str, default=None,
                        help='Generate coverage report CSV file')
    
    args = parser.parse_args()
    
    # parse vcf into a human readable table
    # if there is 
    df = process_vcf(
        args.vcf, 
        args.samplesheet, 
        args.outfile_path, 
        subset_path=args.subset_path,
        subset_padding=args.subset_padding,
        max_workers=args.subset_workers,
        chrom_mapping_file=args.chrom_mapping
    )
    
    # Generate coverage report if requested and subsets were created
    if args.coverage_report and args.subset_path:
        logger.info("Generating coverage report")
        # Load chromosome mapping for coverage report
        chrom_mapping = None
        if args.chrom_mapping:
            mapping_df = pd.read_csv(
                args.chrom_mapping, sep='\t', header=None, 
                names=['bam_chrom', 'vcf_chrom']
                )
            chrom_mapping = dict(zip(
                mapping_df['vcf_chrom'], mapping_df['bam_chrom']
                ))
        
        report_df = generate_coverage_report(
            df, args.subset_path, chrom_mapping, args.subset_workers
            )
        report_df.to_csv(args.coverage_report, index=False)
        logger.info(f"Coverage report saved to {args.coverage_report}")
        
        # Print summary
        print(f"\nSUMMARY:")
        print(f"Samples processed: {len(report_df)}")
        print(f"Successful subsets: {report_df['subset_exists'].sum()}")
        print(f"Average coverage rate: {report_df['coverage_rate'].mean():.2%}")
    
    print(f"Processed {len(df)} variant-sample combinations")
    print(f"Results saved to {args.outfile_path}")