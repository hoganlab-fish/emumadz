#!/usr/bin/env python3
# loads a vcf with annotations, parses into a table, subset and rename alignments if needed
import argparse
import json
import os
import pandas as pd
import pysam
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union

class VCFParser:
    def __init__(
            self, 
            vcf_file: str, 
            samplesheet_file: str, 
            chrom_mapping_file: Optional[str] = None
        ):
        """Parse annotated VCF files and create BAM subsets for variant analysis.
        
        This class processes VCF files with VEP/SnpEff annotations, matches samples 
        to BAM files, and optionally creates chromosome-renamed BAM subsets around 
        variant positions for visualization.
        
        :param vcf_file: Path to annotated VCF file
        :type vcf_file: str
        :param samplesheet_file: Tab-separated file mapping sample IDs to BAM paths and types
        :type samplesheet_file: str
        :param chrom_mapping_file: Optional TSV file mapping old to new chromosome names
        :type chrom_mapping_file: Optional[str]
        
        Example:
            >>> parser = VCFParser('variants.vcf', 'samples.tsv', 'chr_mapping.tsv')
            >>> df = parser.process_vcf('output.csv', subset_path='bam_subsets/')
        """
        self.vcf_file = vcf_file
        self.samplesheet = pd.read_csv(samplesheet_file, sep='\t')
        self.sample_info = dict(zip(self.samplesheet['sample_identity'], 
                                zip(self.samplesheet['alignment_file'], 
                                    self.samplesheet['sample_type'])))
        
        self.chrom_mapping = None
        if chrom_mapping_file:
            mapping_data = pd.read_csv(
                chrom_mapping_file, sep='\t', header=None, 
                names=['bam_chrom', 'renamed_chrom']
                )
            self.chrom_mapping = dict(zip(
                mapping_data['bam_chrom'], mapping_data['renamed_chrom']
                ))
        
        # Get reference samples for shared processing
        self.reference_samples = {
            sample_id: bam_file for sample_id, (bam_file, sample_type) 
            in self.sample_info.items() if sample_type.lower() in ['control', 'reference']
            }
    
    def get_relative_bam_file(
            self, output_file: str, bam_filename: str, subset_path: str
        ) -> str:
        """Get the correct relative path from output JSON to BAM file.
        
        :param output_file: Path to the output JSON/CSV file
        :type output_file: str
        :param bam_filename: BAM filename (e.g., "sample_mutant.bam")
        :type bam_filename: str
        :param subset_path: Directory where BAM subsets are created
        :type subset_path: str
        :returns: Relative path from JSON to BAM file
        :rtype: str
        """
        if not subset_path or not bam_filename:
            return bam_filename
        
        # Get absolute paths
        output_dir = os.path.dirname(os.path.abspath(output_file))
        bam_abs_path = os.path.abspath(os.path.join(subset_path, bam_filename))
        
        # Calculate relative path from output directory to BAM file
        try:
            relative_path = os.path.relpath(bam_abs_path, output_dir)
            return relative_path
        except ValueError:
            # If on different drives (Windows), return absolute path
            return bam_abs_path
    
    def parse_csq_field(self, csq_string: str) -> List[Dict[str, str]]:
        """Parse VEP CSQ field.
        
        :param csq_string: VEP CSQ annotation string
        :type csq_string: str
        :returns: List of parsed annotation dictionaries
        :rtype: List[Dict[str, str]]
        """
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
    
    def parse_ann_field(self, ann_string: str) -> List[Dict[str, str]]:
        """Parse SnpEff ANN field.
        
        :param ann_string: SnpEff ANN annotation string
        :type ann_string: str
        :returns: List of parsed annotation dictionaries
        :rtype: List[Dict[str, str]]
        """
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
    
    def calculate_allele_freq(self, genotype: str) -> Tuple[Optional[float], Optional[str]]:
        """Calculate allele frequency from genotype.
        
        :param genotype: Genotype string (e.g., "0/1", "1/1")
        :type genotype: str
        :returns: Tuple of (frequency_float, count_string)
        :rtype: Tuple[Optional[float], Optional[str]]
        """
        if not genotype or genotype == './.':
            return None, None
        
        gt_calls = str(genotype).split('/') if '/' in str(genotype) else str(genotype).split('|')
        ref_count = sum(1 for call in gt_calls if call == '0')
        alt_count = sum(1 for call in gt_calls if call not in ['0', '.'])
        total = ref_count + alt_count
        
        if total == 0:
            return None, None
        
        return alt_count / total, f"{alt_count}/{total}"
    
    def get_impact_color(self, impact: str) -> str:
        """Get color coding for impact level.
        
        :param impact: Impact level string (HIGH, MODERATE, LOW, MODIFIER)
        :type impact: str
        :returns: Hex color code
        :rtype: str
        """
        colors = {
            'HIGH': '#FF0000',
            'MODERATE': '#FF8C00', 
            'LOW': '#32CD32',
            'MODIFIER': '#87CEEB'
        }
        return colors.get(impact.upper(), '#CCCCCC')
    
    def match_sample(self, vcf_sample: str) -> Tuple[Optional[str], Optional[str]]:
        """Match VCF sample to samplesheet entry.
        
        :param vcf_sample: Sample name from VCF
        :type vcf_sample: str
        :returns: Tuple of (bam_file, sample_type)
        :rtype: Tuple[Optional[str], Optional[str]]
        """
        # Exact match first
        if vcf_sample in self.sample_info:
            return self.sample_info[vcf_sample]
        
        # Fuzzy match - find samplesheet entry that's a prefix of VCF sample
        for sheet_sample in self.sample_info:
            if vcf_sample.startswith(sheet_sample):
                return self.sample_info[sheet_sample]
        
        return None, None
    
    def create_reference_bam_subsets(
            self, 
            subset_path: str, 
            variant_positions: List[Tuple[str, int]], 
            force_overwrite: bool = False
        ) -> Dict[str, str]:
        """Create chromosome-renamed reference BAM subsets for all reference samples.
        
        :param subset_path: Output directory for subset BAMs
        :type subset_path: str
        :param variant_positions: List of (chromosome, position) tuples
        :type variant_positions: List[Tuple[str, int]]
        :param force_overwrite: Overwrites bam files
        :type force_overwrite: bool
        :returns: Dictionary mapping sample_id to output BAM path
        :rtype: Dict[str, str]
        """
        reference_bam_files = {}
        
        for ref_sample_id, ref_bam_file in self.reference_samples.items():
            if not os.path.exists(ref_bam_file):
                continue
                
            output_path = os.path.join(subset_path, f"{ref_sample_id}_reference.bam")
            
            # Skip if already exists
            if os.path.exists(output_path) and not force_overwrite:
                reference_bam_files[ref_sample_id] = f"{ref_sample_id}_reference.bam"
                continue
                
            # Pass ALL variant positions to create_bam_subset, not just the first one
            success = self.create_bam_subset(ref_bam_file, variant_positions, output_path)
            if success:
                reference_bam_files[ref_sample_id] = f"{ref_sample_id}_reference.bam"
        
        return reference_bam_files
    
    def create_bam_subset(
            self, 
            bam_file: str, 
            variant_positions: List[Tuple[str, int]], 
            output_path: str 
        ) -> bool:
        """Create subset BAM containing reads around variant positions.
        
        :param bam_file: Path to source BAM file
        :type bam_file: str
        :param variant_positions: List of (chromosome, position) tuples
        :type variant_positions: List[Tuple[str, int]]
        :param output_path: Output BAM file path
        :type output_path: str
        :returns: Success status
        :rtype: bool
        """
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Index BAM if needed
        if not os.path.exists(bam_file + ".bai") or \
            not os.path.exists(bam_file.replace(".bam", ".bai")):
            pysam.index(bam_file)
        
        with pysam.AlignmentFile(bam_file, "rb") as inbam:
            # Modify header for chromosome renaming
            header = inbam.header.to_dict()
            if self.chrom_mapping and 'SQ' in header:
                for sq in header['SQ']:
                    old_name = sq['SN']
                    if old_name in self.chrom_mapping:
                        sq['SN'] = self.chrom_mapping[old_name]
            
            with pysam.AlignmentFile(
                    output_path, 
                    "wb", 
                    header=pysam.AlignmentHeader.from_dict(header)
                ) as outbam:
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
                    
                    start = max(0, pos - 1)
                    end = pos + 1
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
    
    def process_sample_variants(
            self, 
            sample_data: Dict[str, Union[str, pd.DataFrame]], 
            force_overwrite: bool = False
        ) -> bool:
        """Process all variants for a single sample and create subset BAM.
        
        :param sample_data: Dictionary containing sample info and variants
        :type sample_data: Dict[str, Union[str, pd.DataFrame]]
        :param force_overwrite: Whether to overwrite existing files
        :type force_overwrite: bool
        :returns: Success status
        :rtype: bool
        """
        sample_id = sample_data['sample_identity']
        bam_file = sample_data['alignment_file']
        output_path = sample_data['output_path']
        variants = sample_data['variants']
        
        if os.path.exists(output_path) and not force_overwrite:
            return True
        
        if not os.path.exists(bam_file):
            return False
        
        variant_positions = [
            (row['chromosome'], row['position']) for _, row in variants.iterrows()
        ]
        return self.create_bam_subset(
            bam_file, variant_positions, output_path
        )
    
    def validate_bam_subset(
            self, 
            subset_path: str, 
            expected_variants: List[Tuple[str, int]]
        ) -> Tuple[bool, Union[str, Dict]]:
        """Validate that a BAM subset contains reads covering expected variant positions.
        
        :param subset_path: Path to subset BAM file
        :type subset_path: str
        :param expected_variants: List of (chromosome, position) tuples
        :type expected_variants: List[Tuple[str, int]]
        :returns: Tuple of (success_bool, result_dict)
        :rtype: Tuple[bool, Union[str, Dict]]
        """
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
    
    def generate_coverage_report(
            self, data: pd.DataFrame, subset_dir: str
        ) -> pd.DataFrame:
        """Generate coverage report for all BAM subsets.
        
        :param data: DataFrame with variant data
        :type data: pd.DataFrame
        :param subset_dir: Directory containing subset BAMs
        :type subset_dir: str
        :returns: Coverage report DataFrame
        :rtype: pd.DataFrame
        """
        report = []
        
        # Group by mutant sample only
        mutant_variants = data[data['sample_type'].str.lower() == 'mutant']
        
        for _, row in mutant_variants.iterrows():
            variant_positions = [(row['chromosome'], row['position'])]
            
            # Check mutant BAM
            mutant_subset_path = os.path.join(subset_dir, f"{row['sample_id']}_mutant.bam")
            is_valid, result = self.validate_bam_subset(mutant_subset_path, variant_positions)
            
            if is_valid:
                report.append({
                    'variant_id': row['variant_id'],
                    'sample_id': row['sample_identity'],
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
                    'sample_id': row['sample_identity'],
                    'sample_type': 'mutant',
                    'subset_exists': False,
                    'error': result,
                    'total_variants': 1,
                    'covered_variants': 0,
                    'coverage_rate': 0,
                    'total_reads': 0
                })
            
            # Check reference BAMs for this variant
            if row['reference_bam_files']:
                for ref_path in row['reference_bam_files'].split('|'):
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

    def process_vcf(
            self,
            output_file: str,
            subset_path: Optional[str] = None, 
            max_workers: int = 4, 
            force_overwrite: bool = False, 
            coverage_report: Optional[str] = None
        ) -> pd.DataFrame:
        """Process VCF and create analysis table.
        
        :param output_file: Output CSV file path
        :type output_file: str
        :param subset_path: Directory for subset BAM files
        :type subset_path: Optional[str]
        :param max_workers: Number of parallel workers for BAM subsetting
        :type max_workers: int
        :param force_overwrite: Whether to overwrite existing files
        :type force_overwrite: bool
        :param coverage_report: Path for coverage report CSV
        :type coverage_report: Optional[str]
        :returns: DataFrame with processed variants
        :rtype: pd.DataFrame
        """
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
                bam_file, sample_type = self.match_sample(sample_name)
                if not bam_file:
                    continue
                
                if sample_type.lower() == 'mutant' and mutant_sample is None:
                    mutant_sample = (sample_name, bam_file, sample_type)
                elif sample_type.lower() in ['control', 'reference']:
                    reference_samples.append((sample_name, bam_file, sample_type))
            
            if not mutant_sample:
                continue
            
            # Process mutant sample
            sample_name, bam_file, sample_type = mutant_sample
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
            
            # Create reference BAM paths list - process once for all variants
            reference_bam_files = []
            if subset_path and reference_samples:
                # Create shared reference BAMs if chromosome mapping exists
                reference_bam_files = [
                    f"{ref_sample_id}_reference.bam" 
                    for ref_sample_id in self.reference_samples.keys()
                ]
            
            # Get full annotations as JSON strings
            vep_full = json.dumps(vep_annotations) if vep_annotations else ''
            snpeff_full = json.dumps(snpeff_annotations) if snpeff_annotations else ''
            
            # Generate BAM subset filenames (will be converted to relative paths later)
            mutant_bam_filename = f"{sample_name}_mutant.bam" if subset_path else ""
            
            variant_data = {
                'variant_id': variant_id,
                'chromosome': record.chrom,
                'position': record.pos,
                'ref': record.ref,
                'alt': ','.join([str(alt) for alt in record.alts]),
                'sample_identity': sample_name,
                'sample_type': sample_type,
                'alignment_file': bam_file,
                'mutant_bam_subset_path': mutant_bam_filename,
                'reference_bam_files': '|'.join(reference_bam_files),
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
                'vep_full_annotations': vep_full,
                'snpeff_full_annotations': snpeff_full,
                'custom_notes': '',
                'vep_impact_color': self.get_impact_color(vep_top.get('impact', '')),
                'snpeff_impact_color': self.get_impact_color(snpeff_top.get('impact', ''))
            }
            
            variants.append(variant_data)
        
        # Create DataFrame
        data = pd.DataFrame(variants)
        
        # Create BAM subsets if requested
        if subset_path and len(data) > 0:
            # Create shared reference BAMs first if chromosome mapping is provided
            if self.chrom_mapping and self.reference_samples:
                all_variants = [(row['chromosome'], row['position']) for _, row in data.iterrows()]
                self.create_reference_bam_subsets(subset_path, all_variants)
            
            self._create_bam_subsets(data, subset_path, max_workers, force_overwrite)
            
            # Update BAM paths to relative paths after subset creation
            data['mutant_bam_subset_path'] = data['mutant_bam_subset_path'].apply(
                lambda x: self.get_relative_bam_file(output_file, x, subset_path) if x else ''
            )
            
            # Update reference BAM paths to relative paths
            data['reference_bam_files'] = data['reference_bam_files'].apply(
                lambda x: '|'.join([
                    self.get_relative_bam_file(output_file, ref_path, subset_path) 
                    for ref_path in x.split('|') if ref_path
                ]) if x else ''
            )
            
            # Generate coverage report if requested
            if coverage_report:
                report_data = self.generate_coverage_report(data, subset_path)
                report_data.to_csv(coverage_report, index=False)
                print(f"Coverage report saved to {coverage_report}")
                print(f"Regions processed: {len(report_data)}")
                print(f"Successful subsets: {report_data['subset_exists'].sum()}")
                print(f"Average coverage rate: {report_data['coverage_rate'].mean():.2%}")
        
        # Save files
        data.to_csv(output_file, index=False)
        
        # Save JSON with correct relative paths
        json_output = output_file.replace('.csv', '.json')
        data.to_json(json_output, orient='records', indent=2)
        
        return data
    
    def _create_bam_subsets(
            self, 
            data: pd.DataFrame, 
            subset_path: str, 
            max_workers: int, 
            force_overwrite: bool
        ) -> None:
        """Create BAM subsets for mutant samples only (references handled separately).
        
        :param data: DataFrame with variant data
        :type data: pd.DataFrame
        :param subset_path: Output directory path
        :type subset_path: str
        :param max_workers: Number of parallel workers
        :type max_workers: int
        :param force_overwrite: Whether to overwrite existing files
        :type force_overwrite: bool
        """
        sample_data_list = []
        
        for _, row in data.iterrows():
            # Add only mutant samples
            sample_data_list.append({
                'sample_id': row['sample_identity'],
                'alignment_file': row['alignment_file'],
                'output_path': os.path.join(subset_path, f"{row['sample_identity']}_mutant.bam"),
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
                executor.submit(
                    self.process_sample_variants, 
                    sample_data, 
                    force_overwrite
                    ): sample_data['sample_id']
                for sample_data in unique_sample_data
            }
            
            for future in as_completed(future_to_sample):
                if future.result():
                    success_count += 1
        
        print(f"BAM subsets created for {success_count} samples")


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Parse annotated VCF into tables with optional BAM subsetting'
    )
    parser.add_argument('vcf', type=str,
                        help='Merged VCF file with VEP and SnpEff annotations')
    parser.add_argument('samplesheet', type=str,
                        help='Sample metadata file')
    parser.add_argument('outfile_path', type=str,
                        help='Output CSV file')
    parser.add_argument('--subset_path', type=str, default=None, 
                        help='Directory path for subset BAM files')
    parser.add_argument('--subset_workers', type=int, default=4, 
                        help='Number of parallel workers for BAM subsetting')
    parser.add_argument('--chrom_mapping', type=str, default=None, 
                        help='Two-column TSV mapping BAM chromosome names')
    parser.add_argument('--coverage_report', type=str, default=None, 
                        help='Generate coverage report CSV file')
    parser.add_argument('--force_overwrite', action='store_true', 
                        help='Overwrite existing subset BAM files')
    
    args = parser.parse_args()
    
    # Create parser and process VCF
    vcf_parser = VCFParser(args.vcf, args.samplesheet, args.chrom_mapping)
    data = vcf_parser.process_vcf(
        args.outfile_path,
        subset_path=args.subset_path,
        max_workers=args.subset_workers,
        force_overwrite=args.force_overwrite,
        coverage_report=args.coverage_report
    )
    
    print(f"Processed {len(data)} mutant variant entries")
    print(f"Results saved to {args.outfile_path}")

if __name__ == "__main__":
    main()
