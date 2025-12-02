#!/usr/bin/env python3
# Identifies regions where mutant is homozygous while references are heterozygous.
# Uses: https://doi.org/10.1016/j.ymeth.2013.05.015
import argparse
import os
import json
import subprocess
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import partial
from multiprocessing import Pool
from typing import Dict, List, Tuple
import pandas as pd
import pysam
from tqdm import tqdm

@dataclass
class VariantData:
    """Container for variant information"""
    chrom: str
    pos: int
    is_snp: bool
    
    def get_sample_data(self, record, sample_name: str) -> Dict:
        """Get sample data on demand"""
        return dict(record.samples[sample_name])

class HomozygosityError(Exception):
    """Custom exception for homozygosity analysis errors"""
    pass

def add_homozygosity_to_variant_json(
    variant_json_path: str,
    tdf_file: str,
    bedgraph_file: str,
    window_size: int,
    display_name: str,
    overwrite: bool = True
    ) -> None:
    """
    Add or update homozygosity track information in variant JSON file.
    
    Args:
        variant_json_path: Path to variant JSON file from VCF parser
        tdf_file: Path to TDF file (or None if not generated)
        bedgraph_file: Path to bedgraph file
        window_size: Analysis window size
        display_name: Sample display name (from samplesheet)
        overwrite: If True, overwrite existing homozygosity data; if False, append new track
    """
    # Load existing variant data
    with open(variant_json_path, 'r') as f:
        variant_data = json.load(f)
    
    # Get directory for relative path calculation
    json_dir = os.path.dirname(os.path.abspath(variant_json_path))
    
    # Calculate relative paths
    if tdf_file and os.path.exists(tdf_file):
        tdf_rel = os.path.relpath(tdf_file, json_dir)
        track_url = tdf_rel
        track_format = "tdf"
    else:
        bedgraph_rel = os.path.relpath(bedgraph_file, json_dir)
        track_url = bedgraph_rel
        track_format = "bedgraph"
    
    # Create homozygosity track config
    homozygosity_track = {
        "name": f"Homozygosity Mapping - {display_name}",
        "type": "wig",
        "format": track_format,
        "url": track_url,
        "height": 100,
        "color": "rgb(0, 150, 200)",
        "altColor": "rgb(0, 100, 150)",
        "min": 0,
        "max": 50,
        "autoscale": False,
        "displayMode": "COLLAPSED",
        "visibilityWindow": -1
    }
    
    # Update each variant record
    updated_count = 0
    added_count = 0
    
    for variant in variant_data:
        if overwrite:
            # Remove old homozygosity keys and replace
            old_keys = [k for k in variant.keys() if k.startswith('homozygosity_')]
            for key in old_keys:
                del variant[key]
            
            # Set as single track
            variant['homozygosity_track'] = homozygosity_track
            variant['homozygosity_window_size'] = window_size
            variant['homozygosity_display_name'] = display_name
            variant['homozygosity_updated'] = pd.Timestamp.now().isoformat()
            updated_count += 1
        else:
            # Append mode: add to list of tracks
            if 'homozygosity_tracks' not in variant:
                variant['homozygosity_tracks'] = []
            
            # Check if this track already exists (by URL)
            existing_urls = [t.get('url') for t in variant.get('homozygosity_tracks', [])]
            if track_url not in existing_urls:
                variant['homozygosity_tracks'].append(homozygosity_track)
                
                # Store metadata in a list as well
                if 'homozygosity_metadata' not in variant:
                    variant['homozygosity_metadata'] = []
                
                variant['homozygosity_metadata'].append({
                    'window_size': window_size,
                    'display_name': display_name,
                    'updated': pd.Timestamp.now().isoformat()
                })
                added_count += 1
    
    # Save updated JSON
    with open(variant_json_path, 'w') as f:
        json.dump(variant_data, f, indent=2)
    
    if overwrite:
        print(f"✓ Replaced homozygosity track info in: {variant_json_path}")
        print(f"  Variants updated: {updated_count}")
    else:
        print(f"✓ Added additional homozygosity track to: {variant_json_path}")
        print(f"  Variants updated: {added_count}")
    
    print(f"  Track URL: {track_url}")
    print(f"  Display Name: {display_name}")
    print(f"  Format: {track_format}")

def load_samplesheet(samplesheet_path: str) -> Dict[str, str]:
    """
    Load samplesheet to map sample names to VCF sample IDs.
    
    Args:
        samplesheet_path: Path to tab-separated samplesheet
    
    Returns:
        Dictionary mapping sample_identity to VCF sample names
    
    Example samplesheet format:
        sample_identity    vcf_sample_name    alignment_file    sample_type
        mutant_1           TL2312073-163-4L   /path/to.bam      mutant
    """
    if not samplesheet_path or not os.path.exists(samplesheet_path):
        return {}
    
    samplesheet = pd.read_csv(samplesheet_path, sep='\t')
    
    # Create mapping: sample_identity -> vcf_sample_name
    # This allows consistent naming even if VCF sample names change
    sample_mapping = {}
    for _, row in samplesheet.iterrows():
        sample_identity = row['sample_identity']
        vcf_sample = row.get('vcf_sample_name', sample_identity)
        sample_mapping[sample_identity] = vcf_sample
    
    return sample_mapping

def parse_allele_depths(ad_field) -> List[int]:
    """Parse AD field with comprehensive error handling"""
    if not ad_field or ad_field == '.' or ad_field is None:
        raise HomozygosityError("Missing AD field")
    try:
        if isinstance(ad_field, (list, tuple)):
            return [int(x) if x is not None else 0 for x in ad_field]
        elif isinstance(ad_field, str):
            return [int(x) for x in ad_field.split(',')]
        else:
            return [int(x) for x in str(ad_field).split(',')]
    except (ValueError, AttributeError) as e:
        raise HomozygosityError(f"Invalid AD format: {ad_field}") from e

def bedgraph_to_tdf(bedgraph_file: str, fasta_index: str, output_tdf: str) -> bool:
    """
    Convert bedgraph to TDF (Tiled Data Format) for efficient IGV visualization.
    
    TDF is a binary format that loads much faster than text bedgraph for large
    genomic datasets. Requires IGVTools to be installed and in PATH.
    
    Args:
        bedgraph_file: Input bedgraph file path
        fasta_index: Reference genome .fai index file (e.g., hg38.fa.fai)
        output_tdf: Output TDF file path
    
    Returns:
        True if conversion succeeded, False if IGVTools not available or failed
    
    Example:
        >>> success = bedgraph_to_tdf("scores.bedgraph", "genome.fa.fai", "scores.tdf")
        >>> if success:
        ...     print("Load scores.tdf in IGV for visualization")
    """
    try:
        subprocess.run(['igvtools'], capture_output=True, check=False)
        cmd = ['igvtools', 'toTDF', bedgraph_file, output_tdf, fasta_index]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            print(f"Successfully created TDF: {output_tdf}")
            return True
        else:
            print(f"IGV tools error: {result.stderr}")
            return False
    except FileNotFoundError:
        print("IGV tools not found - keeping bedgraph format only")
        return False
class HomozygosityScorer(ABC):
    """Abstract base class for homozygosity scoring methods"""
    
    def __init__(self):
        self._validate_config()
    
    @abstractmethod
    def _validate_config(self) -> None:
        pass
    
    @abstractmethod
    def requires_references(self) -> bool:
        pass
    
    @abstractmethod
    def calculate_score(self, record, variant: VariantData, sample_name: str, 
                       ref_samples: List[str]) -> Dict[str, int]:
        pass
    
    @abstractmethod
    def get_name(self) -> str:
        pass

class HenkeScorer(HomozygosityScorer):
    """
    Scorer: (hom / het) × (missing_ref_alleles / present_ref_alleles)
    
    Implements the Henke et al. (2013) mapping score formula, which identifies
    regions where mutants are homozygous for alleles that differ from reference
    strains. High scores indicate candidate linked regions.
    
    The score is calculated as:
        (homozygous_SNPs / heterozygous_SNPs) × (ref_alleles_absent / ref_alleles_present)
    
    Attributes:
        min_coverage: Minimum read depth required for a sample to be considered
        homozygosity_threshold: Allele fraction threshold for "soft" homozygosity
            (0.85 means ≥85% of reads support one allele counts as homozygous)
    
    Example:
        scorer = HenkeScorer(min_coverage=10, homozygosity_threshold=0.85)
    """
    
    def __init__(self, min_coverage: int = 1, homozygosity_threshold: float = 0.85):
        self.min_coverage = min_coverage
        self.homozygosity_threshold = homozygosity_threshold
        super().__init__()
    
    def _validate_config(self) -> None:
        if self.min_coverage < 1:
            raise ValueError(f"Min coverage must be >= 1, got {self.min_coverage}")
        if not 0.5 <= self.homozygosity_threshold <= 1.0:
            raise ValueError(f"Homozygosity threshold must be 0.5-1.0, got {self.homozygosity_threshold}")
    
    def requires_references(self) -> bool:
        return True
    
    def _is_homozygous(self, covered: set, depths: List[int]) -> bool:
        """
        Determine if a genotype is homozygous using strict or soft thresholding.
    
        Strict homozygosity: Only one allele has any coverage
        Soft homozygosity: One allele has ≥ homozygosity_threshold fraction of reads
        
        Args:
            covered: Set of allele indices with non-zero coverage
            depths: List of read depths per allele (e.g., [100, 5, 0])
        
        Returns:
            True if genotype meets homozygosity criteria, False otherwise
        
        Example:
            >>> scorer._is_homozygous({0, 1}, [95, 5, 0])  # 95% ref allele
            True  # if threshold is 0.85
        """
        if len(covered) == 1:
            return True
        total = sum(depths)
        return (max(depths) / total) >= self.homozygosity_threshold if total > 0 else False
    
    def calculate_score(self, record, variant: VariantData, sample_name: str, 
                       ref_samples: List[str]) -> Dict[str, int]:
        """
        Calculate per-variant homozygosity counts for window aggregation.
        
        Examines one VCF record and returns count flags that are later summed
        across genomic windows to compute the Henke score.
        
        Args:
            record: pysam VCF record object
            variant: VariantData container with position info
            sample_name: Name of mutant sample
            ref_samples: List of reference sample names
        
        Returns:
            Dictionary with binary flags (0 or 1):
                - mut_is_hom: 1 if mutant is homozygous
                - mut_is_het: 1 if mutant is heterozygous
                - mut_has_map: 1 if all reference alleles present in mutant
                - mut_miss_map: 1 if any reference alleles missing from mutant
        
        Example:
            >>> counts = scorer.calculate_score(record, variant, "mutant1", ["ref1", "ref2"])
            >>> counts
            {'mut_is_hom': 1, 'mut_is_het': 0, 'mut_has_map': 0, 'mut_miss_map': 1}
        """
        counts = {
            'mut_is_hom': 0,
            'mut_is_het': 0,
            'mut_has_map': 0,
            'mut_miss_map': 0
        }
        
        try:
            # Get mutant data
            mutant_call = variant.get_sample_data(record, sample_name)
            mutant_depths = parse_allele_depths(mutant_call.get('AD'))
            if sum(mutant_depths) < self.min_coverage:
                return counts
            
            mutant_covered = set(i for i, d in enumerate(mutant_depths) if d > 0)
            mutant_is_hom = self._is_homozygous(mutant_covered, mutant_depths)
            
            counts['mut_is_hom' if mutant_is_hom else 'mut_is_het'] = 1
            
            # Aggregate reference samples
            ref_combined_covered = set()
            for ref_sample in ref_samples:
                try:
                    ref_call = variant.get_sample_data(record, ref_sample)
                    ref_depths = parse_allele_depths(ref_call.get('AD'))
                    if sum(ref_depths) >= self.min_coverage:
                        ref_covered = set(i for i, d in enumerate(ref_depths) if d > 0)
                        ref_combined_covered.update(ref_covered)
                except HomozygosityError:
                    continue
            
            if not ref_combined_covered:
                return counts
            
            # Track mapping alleles (ref alleles not in mutant)
            mapping_alleles = ref_combined_covered - mutant_covered
            counts['mut_miss_map' if len(mapping_alleles) > 0 else 'mut_has_map'] = 1
            
            return counts
        except HomozygosityError:
            return counts
    
    def get_name(self) -> str:
        return f"henke_cov{self.min_coverage}_hom{self.homozygosity_threshold}"

class HomozygosityAnalyser:
    """
    Processes VCF files to generate homozygosity scores across genomic windows.
    
    Supports two windowing strategies:
        - BP-based: Fixed genomic coordinates (e.g., 10kb windows)
        - SNP-based: Fixed number of variants (e.g., 200 SNP windows)
    
    BP windows are useful for uniform genome coverage, while SNP windows
    normalize for variable mutation density across the genome.
    
    Attributes:
        scorer: HomozygosityScorer instance (e.g., HenkeScorer)
        window_size: Window size in bp (BP mode) or SNP count (SNP mode)
        step_size: Step size in bp (BP mode) or SNP count (SNP mode)
        use_snp_windows: If True, use SNP-count windows; if False, use BP windows
    
    Example:
        >>> scorer = HenkeScorer(min_coverage=10)
        >>> analyser = HomozygosityAnalyser(scorer, window_size=50000, step_size=5000)
        >>> variants_df, windows_df, sample = analyser.process_vcf("data.vcf", False, 8)
    """
    
    def __init__(self, scorer: HomozygosityScorer, 
                 window_size: int = 10000, step_size: int = 1000,
                 use_snp_windows: bool = False, snp_window_size: int = 200, 
                 snp_step_size: int = 20, samplesheet_path: str = None):
        self.scorer = scorer
        self.window_size = window_size
        self.step_size = step_size
        self.use_snp_windows = use_snp_windows
        self.snp_window_size = snp_window_size
        self.snp_step_size = snp_step_size
        
        # Load samplesheet for name mapping
        self.sample_mapping = load_samplesheet(samplesheet_path) if samplesheet_path else {}
        self._validate_config()
    
    def _infer_samples(self, vcf_path: str) -> Tuple[str, List[str], str]:
        """
        Infer sample and reference names from VCF, return display name.
        
        Returns:
            Tuple of (vcf_sample_name, ref_samples_list, display_name)
        """
        with pysam.VariantFile(vcf_path) as vcf:
            samples = [s for s in vcf.header.samples if ':' not in s]
            if len(samples) < 2:
                raise ValueError("VCF must contain at least 2 samples")
            
            vcf_sample_name = samples[0]
            ref_samples = samples[1:]
            
            # Find display name from samplesheet
            display_name = vcf_sample_name
            if self.sample_mapping:
                # Reverse lookup: find sample_identity for this VCF sample
                for identity, vcf_name in self.sample_mapping.items():
                    if vcf_name == vcf_sample_name:
                        display_name = identity
                        break
            
            return vcf_sample_name, ref_samples, display_name
    
    def _get_chromosomes(self, vcf_path: str) -> List[str]:
        """Extract chromosome list from VCF header"""
        with pysam.VariantFile(vcf_path) as vcf:
            return sorted(list(vcf.header.contigs.keys())) if vcf.header.contigs else []
    
    def _process_chromosome(self, vcf_path: str, chrom: str, all_variants: bool) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Process a single chromosome from VCF file.
    
        This is called in parallel for each chromosome when using multiprocessing.
        Extracts variants, calculates per-variant scores, and generates windows.
        
        Args:
            vcf_path: Path to VCF file
            chrom: Chromosome name to process (e.g., "chr1", "scaffold_10")
            all_variants: If True, include indels; if False, SNPs only
        
        Returns:
            Tuple of (per-variant DataFrame, windowed DataFrame)
            Per-variant DF contains individual variant counts
            Windowed DF contains aggregated scores per genomic window
        """
        sample_name, ref_samples = self._infer_samples(vcf_path)
        variant_data = []
        
        with pysam.VariantFile(vcf_path) as vcf:
            for record in vcf:
                if record.chrom != chrom:
                    continue
                
                variant = VariantData(
                    chrom=record.chrom,
                    pos=record.pos,
                    is_snp=(len(record.ref) == 1 and 
                           all(len(str(alt)) == 1 for alt in record.alts if alt))
                )
                
                if not all_variants and not variant.is_snp:
                    continue
                
                counts = self.scorer.calculate_score(record, variant, sample_name, ref_samples)
                if any(counts.values()):
                    variant_data.append({'chrom': variant.chrom, 'pos': variant.pos, **counts})
        
        if not variant_data:
            return pd.DataFrame(), pd.DataFrame()
        
        variants_df = pd.DataFrame(variant_data)
        windowed_df = self._generate_windows(variants_df)
        return variants_df, windowed_df
    
    def _generate_windows(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Generate windows based on mode"""
        if variants_df.empty:
            return pd.DataFrame()
        
        chrom = variants_df['chrom'].iloc[0]
        chrom_data = variants_df.sort_values('pos')
        windowed_data = []
        
        if self.use_snp_windows:
            # SNP-count windows
            for i in range(0, len(chrom_data) - self.snp_window_size + 1, self.snp_step_size):
                window_vars = chrom_data.iloc[i:i + self.snp_window_size]
                score = self._calculate_henke_score(window_vars)
                windowed_data.append({
                    'chrom': chrom,
                    'start': int(window_vars['pos'].min()),
                    'end': int(window_vars['pos'].max()) + 1,
                    'score': score,
                    'variant_count': len(window_vars)
                })
        else:
            # BP windows
            min_pos, max_pos = int(chrom_data['pos'].min()), int(chrom_data['pos'].max())
            window_start = min_pos
            while window_start < max_pos:
                window_end = min(window_start + self.window_size, max_pos + 1)
                window_vars = chrom_data[(chrom_data['pos'] >= window_start) & 
                                        (chrom_data['pos'] < window_end)]
                if len(window_vars) > 0:
                    score = self._calculate_henke_score(window_vars)
                    windowed_data.append({
                        'chrom': chrom,
                        'start': window_start,
                        'end': window_end,
                        'score': score,
                        'variant_count': len(window_vars)
                    })
                window_start += self.step_size
        
        return pd.DataFrame(windowed_data)
    
    def _calculate_henke_score(self, window_vars: pd.DataFrame) -> float:
        """
        Calculate Henke score from aggregated variant counts in a window.
        
        Formula: (homozygous_positions / heterozygous_positions) × 
                 (missing_ref_alleles  / present_ref_alleles)
        
        Pseudocounts (+1) prevent division by zero and provide conservative scoring
        when windows have sparse data.
        
        Args:
            window_vars: DataFrame of variants in window with count columns:
                mut_is_hom, mut_is_het, mut_has_map, mut_miss_map
        
        Returns:
            Henke score as float. Higher scores indicate stronger evidence of linkage.
            Typical range: 0.1-100+, with scores >10 often indicating linked regions.
        
        Example:
            Window with 80 hom, 20 het, 60 missing_map, 30 has_map:
            Score = (80+1)/(20+1) × (60+1)/(30+1) = 3.86 × 1.97 = 7.6
        """
        mut_is_hom = window_vars['mut_is_hom'].sum() + 1
        mut_is_het = window_vars['mut_is_het'].sum() + 1
        mut_has_map = window_vars['mut_has_map'].sum() + 1
        mut_miss_map = window_vars['mut_miss_map'].sum() + 1
        return (float(mut_is_hom) / float(mut_is_het)) * (float(mut_miss_map) / float(mut_has_map))
    
    def process_vcf(self, vcf_path: str, all_variants: bool, ncpu: int) -> Tuple[pd.DataFrame, pd.DataFrame, str]:
        """
        Process entire VCF file with parallel chromosome processing.
        
        Automatically infers sample names from VCF header (first sample = mutant,
        remaining samples = references). Processes each chromosome in parallel
        and combines results.
        
        Args:
            vcf_path: Path to VCF file (must have AD format field)
            all_variants: If True, process all variant types; if False, SNPs only
            ncpu: Number of CPU cores for parallel processing (one per chromosome optimal)
        
        Returns:
            Tuple of (combined_variants_df, combined_windows_df, sample_name, display_name)
            - combined_variants_df: All per-variant counts across genome
            - combined_windows_df: All window scores across genome
            - sample_name: Name of mutant sample
            - display_name: Name shown on display
        
        Raises:
            ValueError: If VCF has <2 samples or samples not found in header
        
        Example:
            >>> analyser = HomozygosityAnalyser(scorer, window_size=10000)
            >>> vars_df, wins_df, sample = analyser.process_vcf("input.vcf", False, 8)
            >>> print(f"Analyzed {len(wins_df)} windows for sample {sample}")
        """
        vcf_sample_name, ref_samples, display_name = self._infer_samples(vcf_path)
        chroms = self._get_chromosomes(vcf_path)
        self._validate_samples(vcf_path, vcf_sample_name, ref_samples)
        
        print(f"Processing {len(chroms)} chromosomes with {ncpu} cores")
        print(f"VCF Sample: {vcf_sample_name}")
        print(f"Display Name: {display_name}")
        print(f"References: {ref_samples}")
        print(f"Scorer: {self.scorer.get_name()}")
        
        # Process chromosomes in parallel
        with Pool(ncpu) as pool:
            process_func = partial(self._process_chromosome, vcf_path, all_variants=all_variants)
            results = list(tqdm(pool.imap(process_func, chroms), total=len(chroms), 
                              desc="Processing chromosomes"))
        
        all_variants_dfs = [v for v, w in results if not v.empty]
        all_windows_dfs = [w for v, w in results if not w.empty]
        
        combined_variants = pd.concat(all_variants_dfs, ignore_index=True) if all_variants_dfs else pd.DataFrame()
        combined_windows = pd.concat(all_windows_dfs, ignore_index=True) if all_windows_dfs else pd.DataFrame()
        
        print(f"Found {len(combined_variants)} variants across all chromosomes")
        return combined_variants, combined_windows, vcf_sample_name, display_name

def write_bedgraph(df: pd.DataFrame, output_file: str, track_name: str, description: str):
    """
    Write windowed scores to UCSC BedGraph format for genome browser visualization.
    
    BedGraph format: chrom <tab> start <tab> end <tab> score
    Output can be loaded directly into IGV, UCSC Genome Browser, or converted to TDF.
    
    Args:
        df: DataFrame with columns: chrom, start, end, score
        output_file: Path for output .bedgraph file
        track_name: Track name for genome browser display
        description: Track description for genome browser
    
    Example:
        >>> write_bedgraph(windows_df, "output.bedgraph", 
        ...                "mutant1_henke", "Henke homozygosity scores")
    """
    with open(output_file, 'w') as f:
        f.write(f'track type=bedGraph name="{track_name}" description="{description}"\n')
        for _, row in tqdm(df.iterrows(), total=len(df), desc=f"Writing {output_file}"):
            f.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['score']:.3f}\n")

def main():
    parser = argparse.ArgumentParser(
        description="""
        Homozygosity mapping with Henke-style scoring.  
        Identifies regions where mutant is homozygous 
        while references are heterozygous. 
        Python implementation of the scorer published in Henke et al, 2013: 
        https://doi.org/10.1016/j.ymeth.2013.05.015. 
        Please also cite the manuscript if using this software.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('vcf_path', help='Input VCF file')
    parser.add_argument('output_prefix', help='Output file prefix')
    parser.add_argument('samplesheet', help='Samplesheet TSV for sample name mapping')    
    parser.add_argument('-c', '--min_coverage', type=int, default=1,
                       help='Minimum coverage threshold')
    parser.add_argument('-t', '--homozygosity_threshold', type=float, default=0.85,
                       help='Allele fraction for soft homozygosity calling (0.5-1.0)')
    parser.add_argument('-w', '--window_size', type=int, default=10000,
                       help='Window size (bp or SNP count)')
    parser.add_argument('--step_size', type=int, default=1000,
                       help='Step size (bp or SNP count)')
    parser.add_argument('--use_snp_windows', action='store_true',
                       help='Use SNP-count windows instead of bp windows')
    parser.add_argument('--snp_window_size', type=int, default=200,
                       help='SNP window size')
    parser.add_argument('--snp_step_size', type=int, default=20,
                       help='SNP step size')
    parser.add_argument('--all_variants', action='store_true',
                       help='Include all variants (default: SNPs only)')
    parser.add_argument('-i', '--fasta_index', type=str, default=None,
                       help='Genome .fai file for TDF generation')
    parser.add_argument('-g', '--genome', type=str, default="danRer7",
                       help='Genome assembly (e.g., hg38, mm10, danRer11)')
    parser.add_argument('-n', '--ncpu', type=int, default=8,
                       help='Number of CPUs')
    parser.add_argument('-j', '--variant_json', type=str, default=None,
                       help='Variant JSON file to add homozygosity tracks to')
    parser.add_argument('--overwrite_homozygosity', action='store_true', default=True,
                       help='Overwrite existing homozygosity data in JSON')    
    
    args = parser.parse_args()
    
    try:
        scorer = HenkeScorer(
            min_coverage=args.min_coverage,
            homozygosity_threshold=args.homozygosity_threshold
        )
        
        analyser = HomozygosityAnalyser(
            scorer=scorer,
            window_size=args.window_size,
            step_size=args.step_size,
            use_snp_windows=args.use_snp_windows,
            snp_window_size=args.snp_window_size,
            snp_step_size=args.snp_step_size,
            samplesheet_path=args.samplesheet,
        )
        
        variants_df, windowed_df, vcf_sample_name, display_name = analyser.process_vcf(
            args.vcf_path, args.all_variants, args.ncpu
        )
        
        if windowed_df.empty:
            print("No results to write")
            return
        
        window_file = f"{args.output_prefix}_windows.bedgraph"
        window_type = "SNP" if args.use_snp_windows else "BP"
        write_bedgraph(windowed_df, window_file,
                      f"{display_name}_{scorer.get_name()}",
                      f"Homozygosity score ({window_type} windows)")
        
        # Generate TDF if fasta index provided
        tdf_file = None
        if args.fasta_index:
            tdf_file = window_file.replace(".bedgraph", ".tdf")
            if not bedgraph_to_tdf(window_file, args.fasta_index, tdf_file):
                tdf_file = None  # Failed to create TDF
        
        # Add homozygosity track to variant JSON
        if args.variant_json and os.path.exists(args.variant_json):
            add_homozygosity_to_variant_json(
                args.variant_json,
                tdf_file,
                window_file,
                args.window_size,
                display_name,
                overwrite=args.overwrite_homozygosity
            )
            print(f"\n✓ Homozygosity track added to {args.variant_json}")
            print("  The variant viewer will now show homozygosity mapping alongside BAM tracks")
        
        print(f"\nGenerated files:")
        print(f"  - {window_file} (bedgraph)")
        if tdf_file:
            print(f"  - {tdf_file} (TDF for IGV)")
        
    except Exception as e:
        print(f"Error: {e}")
        raise

if __name__ == "__main__":
    main()