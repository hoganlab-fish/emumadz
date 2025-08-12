#!/usr/bin/env python3
import argparse
import json
import os
import pandas as pd
import pysam
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


class VCFParser:
    def __init__(self, vcf_file, samplesheet_file, chrom_mapping_file=None):
        self.vcf_file = vcf_file
        self.samplesheet = pd.read_csv(samplesheet_file, sep='\t')
        self.sample_info = dict(zip(self.samplesheet['sample_id'], 
                                  zip(self.samplesheet['bam_path'], 
                                      self.samplesheet['sample_type'])))
        
        self.chrom_mapping = None
        if chrom_mapping_file:
            mapping_df = pd.read_csv(chrom_mapping_file, sep='\t', header=None, 
                                   names=['bam_chrom', 'renamed_chrom'])
            self.chrom_mapping = dict(zip(mapping_df['bam_chrom'], mapping_df['renamed_chrom']))
    
    def parse_csq_field(self, csq_string):
        """Parse VEP CSQ field"""
        if not csq_string:
            return []
        
        annotations = []
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
    
    def parse_ann_field(self, ann_string):
        """Parse SnpEff ANN field"""
        if not ann_string:
            return []
        
        annotations = []
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
    
    def calculate_allele_freq(self, genotype):
        """Calculate allele frequency from genotype"""
        if not genotype or genotype == './.':
            return None, None
        
        gt_calls = str(genotype).split('/') if '/' in str(genotype) else str(genotype).split('|')
        ref_count = sum(1 for call in gt_calls if call == '0')
        alt_count = sum(1 for call in gt_calls if call not in ['0', '.'])
        total = ref_count + alt_count
        
        if total == 0:
            return None, None
        
        return alt_count / total, f"{alt_count}/{total}"
    
    def get_impact_color(self, impact):
        """Get color coding for impact"""
        colors = {
            'HIGH': '#FF0000',
            'MODERATE': '#FF8C00', 
            'LOW': '#32CD32',
            'MODIFIER': '#87CEEB'
        }
        return colors.get(impact.upper(), '#CCCCCC')
    
    def match_sample(self, vcf_sample):
        """Match VCF sample to samplesheet entry"""
        # Exact match first
        if vcf_sample in self.sample_info:
            return self.sample_info[vcf_sample]
        
        # Fuzzy match - find samplesheet entry that's a prefix of VCF sample
        for sheet_sample in self.sample_info:
            if vcf_sample.startswith(sheet_sample):
                return self.sample_info[sheet_sample]
        
        return None, None
    
    def create_bam_subset(self, bam_path, variant_positions, output_path, padding=50):
        """Create subset BAM containing reads around variant positions"""
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Index BAM if needed
        if not os.path.exists(bam_path + ".bai"):
            pysam.index(bam_path)
        
        with pysam.AlignmentFile(bam_path, "rb") as inbam:
            # Modify header for chromosome renaming
            header = inbam.header.to_dict()
            if self.chrom_mapping and 'SQ' in header:
                for sq in header['SQ']:
                    old_name = sq['SN']
                    if old_name in self.chrom_mapping:
                        sq['SN'] = self.chrom_mapping[old_name]
            
            with pysam.AlignmentFile(output_path, "wb", header=pysam.AlignmentHeader.from_dict(header)) as outbam:
                processed_regions = set()
                
                for chrom, pos in variant_positions:
                    # Get original BAM chromosome name for fetching
                    bam_chrom = chrom
                    if self.chrom_mapping:
                        # Find reverse mapping (new -> old)
                        for old_name, new_name in self.chrom_mapping.items():
                            if new_name == chrom:
                                bam_chrom = old_name
                                break
                    
                    start = max(0, pos - padding)
                    end = pos + padding
                    region_key = (bam_chrom, start, end)
                    
                    if region_key in processed_regions:
                        continue
                    processed_regions.add(region_key)
                    
                    for read in inbam.fetch(bam_chrom, start, end):
                        if self.chrom_mapping and read.reference_name in self.chrom_mapping:
                            # Create new read with renamed reference
                            new_read = pysam.AlignedSegment(outbam.header)
                            new_read.query_name = read.query_name
                            new_read.query_sequence = read.query_sequence
                            new_read.flag = read.flag
                            new_read.reference_id = outbam.header.get_tid(self.chrom_mapping[read.reference_name])
                            new_read.reference_start = read.reference_start
                            new_read.mapping_quality = read.mapping_quality
                            new_read.cigar = read.cigar
                            new_read.next_reference_id = read.next_reference_id
                            new_read.next_reference_start = read.next_reference_start
                            new_read.template_length = read.template_length
                            new_read.query_qualities = read.query_qualities
                            new_read.tags = read.tags
                            outbam.write(new_read)
                        else:
                            outbam.write(read)
        
        # Sort and index
        sorted_path = output_path.replace('.bam', '_sorted.bam')
        pysam.sort("-o", sorted_path, output_path)
        os.rename(sorted_path, output_path)
        pysam.index(output_path)
        
        return True
    
    def process_sample_variants(self, sample_data, padding=50, force_overwrite=False):
        """Process all variants for a single sample and create subset BAM"""
        sample_id = sample_data['sample_id']
        bam_path = sample_data['bam_path']
        output_path = sample_data['output_path']
        variants = sample_data['variants']
        
        if os.path.exists(output_path) and not force_overwrite:
            return True
        
        if not os.path.exists(bam_path):
            return False
        
        variant_positions = [(row['chromosome'], row['position']) for _, row in variants.iterrows()]
        return self.create_bam_subset(bam_path, variant_positions, output_path, padding)
    
    def process_vcf(self, output_file, subset_path=None, subset_padding=50, max_workers=4, force_overwrite=False, coverage_report=None):
        """Process VCF and create analysis table"""
        vcf_reader = pysam.VariantFile(self.vcf_file)
        variants = []
        
        for record in vcf_reader:
            variant_id = f"{record.chrom}:{record.pos}:{record.ref}:{','.join([str(alt) for alt in record.alts])}"
            
            # Parse annotations
            vep_annotations = []
            if 'CSQ' in record.info:
                csq_data = record.info['CSQ'][0] if isinstance(record.info['CSQ'], tuple) else record.info['CSQ']
                vep_annotations = self.parse_csq_field(csq_data)
            
            snpeff_annotations = []
            if 'ANN' in record.info:
                ann_data = record.info['ANN'][0] if isinstance(record.info['ANN'], tuple) else record.info['ANN']
                snpeff_annotations = self.parse_ann_field(ann_data)
            
            # Find mutant and reference samples
            mutant_sample = None
            reference_samples = []
            
            for sample_name in record.samples:
                bam_path, sample_type = self.match_sample(sample_name)
                if not bam_path:
                    continue
                
                if sample_type.lower() == 'mutant' and mutant_sample is None:
                    mutant_sample = (sample_name, bam_path, sample_type)
                elif sample_type.lower() in ['control', 'reference']:
                    reference_samples.append((sample_name, bam_path, sample_type))
            
            if not mutant_sample:
                continue
            
            # Process mutant sample
            sample_name, bam_path, sample_type = mutant_sample
            sample = record.samples[sample_name]
            
            # Get genotype info
            gt = sample.get('GT', (None, None))
            dp = sample.get('DP', 0)
            gq = sample.get('GQ', 0)
            
            gt_str = '/'.join([str(g) if g is not None else '.' for g in gt]) if gt else './.'
            af_freq, af_count = self.calculate_allele_freq(gt_str)
            
            # Get top consequences
            vep_top = vep_annotations[0] if vep_annotations else {}
            snpeff_top = snpeff_annotations[0] if snpeff_annotations else {}
            
            # Create reference BAM paths
            reference_bam_paths = []
            if subset_path and reference_samples:
                reference_bam_paths = [f"{ref_name}_reference.bam" for ref_name, _, _ in reference_samples]
            
            variant_data = {
                'variant_id': variant_id,
                'chromosome': record.chrom,
                'position': record.pos,
                'ref': record.ref,
                'alt': ','.join([str(alt) for alt in record.alts]),
                'sample_id': sample_name,
                'sample_type': sample_type,
                'bam_path': bam_path,
                'mutant_bam_subset_path': f"{sample_name}_mutant.bam" if subset_path else "",
                'reference_bam_paths': '|'.join(reference_bam_paths),
                'genotype': gt_str,
                'depth': dp or 0,
                'quality': gq or 0,
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
                'vep_full_annotations': json.dumps(vep_annotations) if vep_annotations else '',
                'snpeff_full_annotations': json.dumps(snpeff_annotations) if snpeff_annotations else '',
                'custom_notes': '',
                'vep_impact_color': self.get_impact_color(vep_top.get('impact', '')),
                'snpeff_impact_color': self.get_impact_color(snpeff_top.get('impact', ''))
            }
            
            variants.append(variant_data)
        
        # Create DataFrame and save
        df = pd.DataFrame(variants)
        df.to_csv(output_file, index=False)
        
        # Save JSON
        json_output = output_file.replace('.csv', '.json')
        df.to_json(json_output, orient='records', indent=2)
        
        # Create BAM subsets if requested
        if subset_path and len(df) > 0:
            self._create_bam_subsets(df, subset_path, subset_padding, max_workers, force_overwrite)
            
            # Generate coverage report if requested
            if coverage_report:
                report_df = self.generate_coverage_report(df, subset_path)
                report_df.to_csv(coverage_report, index=False)
                print(f"Coverage report saved to {coverage_report}")
                print(f"Regions processed: {len(report_df)}")
                print(f"Successful subsets: {report_df['subset_exists'].sum()}")
                print(f"Average coverage rate: {report_df['coverage_rate'].mean():.2%}")
        
        return df
    
    def validate_bam_subset(self, subset_path, expected_variants):
        """Validate that a BAM subset contains reads covering expected variant positions"""
        if not os.path.exists(subset_path):
            return False, "Subset BAM file not found"
        
        with pysam.AlignmentFile(subset_path, "rb") as bam:
            coverage_info = {}
            
            for chrom, pos in expected_variants:
                # Use renamed chromosome if mapping provided
                bam_chrom = self.chrom_mapping.get(chrom, chrom) if self.chrom_mapping else chrom
                
                # Expand region for better read detection (VCF is 1-based, BAM is 0-based)
                read_count = sum(1 for _ in bam.fetch(bam_chrom, pos-10, pos+10))
                coverage_info[f"{chrom}:{pos}"] = read_count
            
            total_reads = sum(coverage_info.values())
            covered_variants = sum(1 for count in coverage_info.values() if count > 0)
            
            return True, {
                'total_reads': total_reads,
                'covered_variants': covered_variants,
                'total_variants': len(expected_variants),
                'coverage_detail': coverage_info
            }
    
    def generate_coverage_report(self, df, subset_dir):
        """Generate coverage report for all BAM subsets"""
        report = []
        
        # Group by mutant sample only
        mutant_variants = df[df['sample_type'].str.lower() == 'mutant']
        
        for _, row in mutant_variants.iterrows():
            variant_positions = [(row['chromosome'], row['position'])]
            
            # Check mutant BAM
            mutant_subset_path = os.path.join(subset_dir, f"{row['sample_id']}_mutant.bam")
            is_valid, result = self.validate_bam_subset(mutant_subset_path, variant_positions)
            
            if is_valid:
                report.append({
                    'variant_id': row['variant_id'],
                    'sample_id': row['sample_id'],
                    'sample_type': 'mutant',
                    'subset_exists': True,
                    'total_variants': result['total_variants'],
                    'covered_variants': result['covered_variants'],
                    'coverage_rate': result['covered_variants'] / result['total_variants'] if result['total_variants'] > 0 else 0,
                    'total_reads': result['total_reads']
                })
            else:
                report.append({
                    'variant_id': row['variant_id'],
                    'sample_id': row['sample_id'],
                    'sample_type': 'mutant',
                    'subset_exists': False,
                    'error': result,
                    'total_variants': 1,
                    'covered_variants': 0,
                    'coverage_rate': 0,
                    'total_reads': 0
                })
            
            # Check reference BAMs for this variant
            if row['reference_bam_paths']:
                for ref_path in row['reference_bam_paths'].split('|'):
                    if not ref_path:
                        continue
                    ref_sample_id = ref_path.replace('_reference.bam', '')
                    ref_subset_path = os.path.join(subset_dir, ref_path)
                    is_valid, result = self.validate_bam_subset(ref_subset_path, variant_positions)
                    
                    if is_valid:
                        report.append({
                            'variant_id': row['variant_id'],
                            'sample_id': ref_sample_id,
                            'sample_type': 'reference',
                            'subset_exists': True,
                            'total_variants': result['total_variants'],
                            'covered_variants': result['covered_variants'],
                            'coverage_rate': result['covered_variants'] / result['total_variants'] if result['total_variants'] > 0 else 0,
                            'total_reads': result['total_reads']
                        })
                    else:
                        report.append({
                            'variant_id': row['variant_id'],
                            'sample_id': ref_sample_id,
                            'sample_type': 'reference',
                            'subset_exists': False,
                            'error': result,
                            'total_variants': 1,
                            'covered_variants': 0,
                            'coverage_rate': 0,
                            'total_reads': 0
                        })
        
        return pd.DataFrame(report)

    def _create_bam_subsets(self, df, subset_path, subset_padding, max_workers, force_overwrite):
        """Create BAM subsets for all samples"""
        sample_data_list = []
        
        for _, row in df.iterrows():
            # Add mutant sample
            sample_data_list.append({
                'sample_id': row['sample_id'],
                'bam_path': row['bam_path'],
                'output_path': os.path.join(subset_path, f"{row['sample_id']}_mutant.bam"),
                'variants': pd.DataFrame([row])
            })
            
            # Add reference samples
            if row['reference_bam_paths']:
                ref_paths = row['reference_bam_paths'].split('|')
                for ref_path in ref_paths:
                    if not ref_path:
                        continue
                    ref_sample_id = ref_path.replace('_reference.bam', '')
                    if ref_sample_id in self.sample_info:
                        ref_bam_path = self.sample_info[ref_sample_id][0]
                        sample_data_list.append({
                            'sample_id': ref_sample_id,
                            'bam_path': ref_bam_path,
                            'output_path': os.path.join(subset_path, ref_path),
                            'variants': pd.DataFrame([row])
                        })
        
        # Remove duplicates and merge variants
        seen_samples = set()
        unique_sample_data = []
        for sample_data in sample_data_list:
            key = (sample_data['sample_id'], sample_data['output_path'])
            if key not in seen_samples:
                seen_samples.add(key)
                unique_sample_data.append(sample_data)
            else:
                for existing in unique_sample_data:
                    if (existing['sample_id'], existing['output_path']) == key:
                        existing['variants'] = pd.concat([existing['variants'], sample_data['variants']], ignore_index=True)
                        break
        
        # Process in parallel
        success_count = 0
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_sample = {
                executor.submit(self.process_sample_variants, sample_data, subset_padding, force_overwrite): sample_data['sample_id']
                for sample_data in unique_sample_data
            }
            
            for future in as_completed(future_to_sample):
                if future.result():
                    success_count += 1
        
        print(f"BAM subsets created for {success_count} samples")


def main():
    parser = argparse.ArgumentParser(description='Parse annotated VCF for variant analysis with optional BAM subsetting')
    parser.add_argument('vcf', help='Merged VCF file with VEP and SnpEff annotations')
    parser.add_argument('samplesheet', help='Sample metadata file')
    parser.add_argument('outfile_path', help='Output CSV file')
    parser.add_argument('--subset_path', default=None, help='Directory path for subset BAM files')
    parser.add_argument('--subset_padding', type=int, default=10, help='Base pairs around variants for BAM subsets')
    parser.add_argument('--subset_workers', type=int, default=4, help='Number of parallel workers for BAM subsetting')
    parser.add_argument('--chrom_mapping', default=None, help='Two-column TSV mapping BAM chromosome names')
    parser.add_argument('--coverage_report', default=None, help='Generate coverage report CSV file')
    parser.add_argument('--force_overwrite', action='store_true', help='Overwrite existing subset BAM files')
    
    args = parser.parse_args()
    
    # Create parser and process VCF
    parser = VCFParser(args.vcf, args.samplesheet, args.chrom_mapping)
    df = parser.process_vcf(
        args.outfile_path,
        subset_path=args.subset_path,
        subset_padding=args.subset_padding,
        max_workers=args.subset_workers,
        force_overwrite=args.force_overwrite,
        coverage_report=args.coverage_report
    )
    
    print(f"Processed {len(df)} mutant variant entries")
    print(f"Results saved to {args.outfile_path}")

if __name__ == "__main__":
    main()