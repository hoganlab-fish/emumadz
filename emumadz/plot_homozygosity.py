#!/usr/bin/env python3
# loads a vcf with annotations, create a tdf file of mutation hotspots
import argparse
import logging
import subprocess
import threading
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import partial
from multiprocessing import Pool
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm

@dataclass
class VariantData:
    """Container for variant information"""
    chrom: str
    pos: int
    ref: str
    alts: List[str]
    is_snp: bool
    
    def get_sample_data(self, record, sample_name: str) -> Dict:
        """Get sample data on demand"""
        return dict(record.samples[sample_name])


class HomozygosityError(Exception):
    """Custom exception for homozygosity analysis errors"""
    pass


def parse_allele_depths(ad_field) -> List[int]:
    """Parse AD field with comprehensive error handling"""
    if not ad_field or ad_field == '.' or ad_field is None:
        raise HomozygosityError("Missing AD field")
    # handle cases of (1,2,None)
    try:
        if isinstance(ad_field, (list, tuple)):
            depths = []
            for x in ad_field:
                if x is None:
                    depths.append(0)  # Treat None as zero depth
                else:
                    depths.append(int(x))
            return depths
        elif isinstance(ad_field, str):
            return [int(x) for x in ad_field.split(',')]
        else:
            return [int(x) for x in str(ad_field).split(',')]
    except (ValueError, AttributeError) as e:
        raise HomozygosityError(f"Invalid AD format: {ad_field}") from e

def bedgraph_to_tdf(bedgraph_file: str, genome_file: str, output_tdf: str) -> bool:
    """
    Convert bedgraph to TDF using IGV tools.

    :param bedgraph_file: Input bedgraph file path
    :type bedgraph_file: str
    :param genome_file: Genome .fai file path
    :type genome_file: str
    :param output_tdf: Output TDF file path
    :type output_tdf: str
    :return: True if successful, False otherwise
    :rtype: bool
    """
    try:
        # Check if igvtools is available
        subprocess.run(['igvtools'], capture_output=True, check=False)
        
        # Convert to TDF
        cmd = [
            'igvtools', 'toTDF',
            bedgraph_file,
            output_tdf,
            genome_file
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            print(f"Successfully created TDF: {output_tdf}")
            return True
        else:
            print(f"IGV tools error: {result.stderr}")
            return False
            
    except FileNotFoundError:
        print("IGV tools not found - keeping bedgraph format")
        return False

def sigmoid_transform(score, k: int=8, x0: float=0.3, normalise: bool=True):
    """
    Sigmoid transformation to enhance contrast between signal and noise

        k=5, x0=0.3 - Gentle transition
        k=8, x0=0.3 - Moderate transition
        k=12, x0=0.2 - Sharp transition
        
    :param score: input score (0.0 to 1.0)
    :type score: float
    :param k: steepness of the curve (higher = sharper transition)
    :type k: int
    :param x0: midpoint where transition occurs (0.0 to 1.0)
    :type x0: float
    :param normalise:
    :type bool: True
    :return: transformed score (0.0 to 1.0)
    :rtype: float
    """
    if score <= 0:
        return 0.0
    if score >= 1:
        return 1.0
    
    # Apply sigmoid transformation
    if normalise is True:
        normalised = (score - x0) * 2
        return 1 / (1 + np.exp(-k * normalised))
    else:
        return 1 / (1 + np.exp(-k * (score - x0)))

class HomozygosityScorer(ABC):
    """Abstract base class for homozygosity scoring methods"""
    
    def __init__(self):
        """Initialize scorer and validate configuration"""
        self._validate_config()
    
    @abstractmethod
    def _validate_config(self) -> None:
        """Validate scorer configuration"""
        pass
    
    @abstractmethod
    def requires_references(self) -> bool:
        """Whether this scorer requires reference samples"""
        pass
    
    @abstractmethod
    def calculate_score(self, record, variant: VariantData, sample_name: str, 
                       ref_samples: List[str]) -> float:
        """Calculate homozygosity score for a variant"""
        pass
    
    @abstractmethod
    def get_name(self) -> str:
        """Return scorer name for output files"""
        pass


class AlleleDepthScorer(HomozygosityScorer):
    """Simple allele depth ratio scorer (original script 1 method)"""
    
    def __init__(self, threshold: float = 0.9):
        self.threshold = threshold
        super().__init__()
    
    def _validate_config(self) -> None:
        if not 0.0 <= self.threshold <= 1.0:
            raise ValueError(f"Threshold must be between 0.0 and 1.0, got {self.threshold}")
    
    def requires_references(self) -> bool:
        return False
    
    def calculate_score(self, record, variant: VariantData, sample_name: str, 
                       ref_samples: List[str]) -> float:
        sample_call = variant.get_sample_data(record, sample_name)
        ad_values = sample_call.get('AD')
        
        try:
            depths = parse_allele_depths(ad_values)
            total_depth = sum(depths)
            if total_depth == 0:
                raise HomozygosityError("Zero total depth")
            
            max_allele_fraction = max(depths) / total_depth
            return max_allele_fraction if max_allele_fraction >= self.threshold else 0.0
            
        except HomozygosityError:
            return 0.0
    
    def get_name(self) -> str:
        return f"allele_depth_t{self.threshold}"


class MutantReferenceScorer(HomozygosityScorer):
    """Mutant vs reference pool comparison scorer"""
    
    def __init__(self, min_coverage: int = 1, require_all_coverage: bool = False):
        self.min_coverage = min_coverage
        self.require_all_coverage = require_all_coverage
        super().__init__()
    
    def _validate_config(self) -> None:
        if self.min_coverage < 1:
            raise ValueError(f"Minimum coverage must be >= 1, got {self.min_coverage}")
    
    def requires_references(self) -> bool:
        return True
    
    def calculate_score(self, record, variant: VariantData, sample_name: str, 
                       ref_samples: List[str]) -> float:
        if not ref_samples:
            raise ValueError("MutantReferenceScorer requires reference samples")
        
        try:
            # Get mutant sample data
            mutant_call = variant.get_sample_data(record, sample_name)
            mutant_depths = parse_allele_depths(mutant_call.get('AD'))
            mutant_total = sum(mutant_depths)
            
            if mutant_total < self.min_coverage:
                return 0.0
            
            # Process reference samples
            ref_alt_fractions = []
            
            for ref_sample in ref_samples:
                try:
                    ref_call = variant.get_sample_data(record, ref_sample)
                    ref_depths = parse_allele_depths(ref_call.get('AD'))
                    ref_total = sum(ref_depths)
                    
                    if ref_total < self.min_coverage:
                        if self.require_all_coverage:
                            return 0.0
                        continue
                    
                    # Calculate alt allele fraction
                    ref_alt_fraction = sum(ref_depths[1:]) / ref_total if len(ref_depths) > 1 else 0.0
                    ref_alt_fractions.append(ref_alt_fraction)
                    
                except HomozygosityError:
                    if self.require_all_coverage:
                        return 0.0
                    continue
            
            if not ref_alt_fractions:
                return 0.0
            
            # Calculate score
            mutant_alt_fraction = sum(mutant_depths[1:]) / mutant_total if len(mutant_depths) > 1 else 0.0
            avg_ref_alt_fraction = sum(ref_alt_fractions) / len(ref_alt_fractions)
            score = mutant_alt_fraction - avg_ref_alt_fraction
            
            return max(0.0, min(1.0, score))
            
        except HomozygosityError:
            return 0.0
    
    def get_name(self) -> str:
        return f"mut_vs_ref_cov{self.min_coverage}"


class HomozygosityAnalyser:
    """Main homozygosity analysis class with streaming processing"""
    
    def __init__(self, scorer: HomozygosityScorer, 
                 window_size: int = 10000, step_size: int = 1000,
                 use_snp_windows: bool = False, snp_window_size: int = 200, 
                 snp_step_size: int = 20, sigmoid: dict = {"k": 8, "x0": 0.3},
                 normalise: bool = False):
        self.scorer = scorer
        self.window_size = window_size
        self.step_size = step_size
        self.use_snp_windows = use_snp_windows
        self.snp_window_size = snp_window_size
        self.snp_step_size = snp_step_size
        self.k = sigmoid["k"]
        self.x0 = sigmoid["x0"]
        self.normalise = normalise
        self.total = None
        
        self._validate_config()
    
    def _validate_config(self) -> None:
        """Validate analyser configuration"""
        if self.window_size <= 0 or self.step_size <= 0:
            raise ValueError("Window and step sizes must be positive")
        if self.snp_window_size <= 0 or self.snp_step_size <= 0:
            raise ValueError("SNP window and step sizes must be positive")
        if self.step_size > self.window_size:
            logging.warning("Step size larger than window size - windows won't overlap")

    def _infer_samples(self, vcf_path: str) -> Tuple[str, List[str]]:
        """Infer sample and reference names from VCF"""
        with pysam.VariantFile(vcf_path) as vcf:
            samples = list(vcf.header.samples)
            if len(samples) < 2:
                raise ValueError("VCF must contain at least 2 samples")
            return samples[0], samples[1:]

    def _count_vcf_records(self, vcf_path: str) -> int:
        """Count entries in VCF file in an efficient way"""
        try:
            with tqdm(desc=f"Counting variants in {vcf_path}", bar_format="{desc}: {elapsed}") as pbar:            
                def _update_timer():
                    while not stop_event.is_set():
                        pbar.update(0)
                        time.sleep(1)
                
                stop_event = threading.Event()
                timer_thread = threading.Thread(target=_update_timer)
                timer_thread.start()
                
                total = subprocess.run(
                    ['bcftools', 'view', '-H', vcf_path], 
                    capture_output=True, text=True, check=True
                )
                
                stop_event.set()
                timer_thread.join()
                
                self.total = len([
                    line for line in total.stdout.strip().split('\n') if line
                    ])
            
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.total = None

    def _get_chromosomes(self, vcf_path: str) -> List[str]:
        """Extract chromosome list from VCF header"""
        with pysam.VariantFile(vcf_path) as vcf:
            if vcf.header.contigs:
                return sorted(list(vcf.header.contigs.keys()))
            else:
                chroms = set()
                for record in tqdm(vcf, desc="Counting chromosmomes"):
                    chroms.add(record.chrom)
                return sorted(list(chroms))

    def _process_chromosome(self, vcf_path: str, chrom: str, all_variants: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Process single chromosome and return variants + windows"""
        sample_name, ref_samples = self._infer_samples(vcf_path)
        
        variant_scores = []
        with pysam.VariantFile(vcf_path) as vcf:
            for record in vcf:
                if record.chrom != chrom:
                    continue
                    
                variant = VariantData(
                    chrom=record.chrom,
                    pos=record.pos,
                    ref=record.ref,
                    alts=[str(alt) for alt in record.alts] if record.alts else [],
                    is_snp=(len(record.ref) == 1 and 
                        all(len(str(alt)) == 1 for alt in record.alts if alt is not None))
                )
                
                if not all_variants and not variant.is_snp:
                    continue
                
                score = self.scorer.calculate_score(record, variant, sample_name, ref_samples)

                if self.k or self.x0:
                    score = sigmoid_transform(score, self.k, self.x0, self.normalise)
                
                if score > 0:
                    variant_scores.append({
                        'chrom': variant.chrom,
                        'pos': variant.pos,
                        'score': score,
                        'is_snp': variant.is_snp
                    })
        
        if not variant_scores:
            return pd.DataFrame(), pd.DataFrame()
        
        # Convert to DataFrame and deduplicate
        variants_df = pd.DataFrame(variant_scores)
        variants_df = variants_df.groupby(['chrom', 'pos'], as_index=False)['score'].max()
        
        # Generate windows for this chromosome
        windowed_df = self._generate_windows_single_chrom(variants_df)
        
        return variants_df, windowed_df

    def process_vcf(self, vcf_path: str, all_variants: bool = True, ncpu: int = 1) -> Tuple[pd.DataFrame, pd.DataFrame, str]:
        """Process VCF with parallel chromosome processing"""
        sample_name, ref_samples = self._infer_samples(vcf_path)
        chroms = self._get_chromosomes(vcf_path)
        self._validate_samples(vcf_path, sample_name, ref_samples)
                
        print(f"Processing {len(chroms)} chromosomes with {ncpu} cores")
        print(f"Sample: {sample_name}")
        print(f"References: {ref_samples}")
        print(f"Scorer: {self.scorer.get_name()}")
        
        # Process chromosomes in parallel
        with Pool(ncpu) as pool:
            process_func = partial(self._process_chromosome, vcf_path, all_variants=all_variants)
            results = list(tqdm(
                pool.imap(process_func, chroms), 
                total=len(chroms), 
                desc="Processing chromosomes"
            ))
        
        # Combine results
        all_variants_dfs = []
        all_windows_dfs = []
        
        for variants_df, windowed_df in results:
            if not variants_df.empty:
                all_variants_dfs.append(variants_df)
            if not windowed_df.empty:
                all_windows_dfs.append(windowed_df)
        
        combined_variants = pd.concat(all_variants_dfs, ignore_index=True) if all_variants_dfs else pd.DataFrame()
        combined_windows = pd.concat(all_windows_dfs, ignore_index=True) if all_windows_dfs else pd.DataFrame()
        
        print(f"Found {len(combined_variants)} variants across all chromosomes")
        
        return combined_variants, combined_windows, sample_name

    def _validate_samples(self, vcf_path: str, sample_name: str, ref_samples: List[str]) -> None:
        """Validate sample names exist in VCF"""
        with pysam.VariantFile(vcf_path) as vcf:
            samples = list(vcf.header.samples)
            
            if sample_name not in samples:
                raise ValueError(f"Sample {sample_name} not found in VCF. Available: {samples}")
            
            if self.scorer.requires_references() and not ref_samples:
                raise ValueError(f"Scorer {self.scorer.get_name()} requires reference samples")
            
            for ref in ref_samples:
                if ref not in samples:
                    raise ValueError(f"Reference sample {ref} not found in VCF. Available: {samples}")

    def _generate_windows_single_chrom(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Generate windows for single chromosome"""
        if self.use_snp_windows:
            return self._generate_snp_windows_single_chrom(variants_df)
        else:
            return self._generate_bp_windows_single_chrom(variants_df)

    def _generate_bp_windows_single_chrom(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Generate BP windows for single chromosome"""
        if variants_df.empty:
            return pd.DataFrame()
        
        windowed_data = []
        chrom = variants_df['chrom'].iloc[0]  # Single chromosome
        chrom_data = variants_df.sort_values('pos')
        
        min_pos = int(chrom_data['pos'].min())
        max_pos = int(chrom_data['pos'].max())
        
        window_start = min_pos
        while window_start < max_pos:
            window_end = min(window_start + self.window_size, max_pos + 1)
            
            window_variants = chrom_data[
                (chrom_data['pos'] >= window_start) & 
                (chrom_data['pos'] < window_end)
            ]
            
            if len(window_variants) > 0:
                windowed_data.append({
                    'chrom': chrom,
                    'start': window_start,
                    'end': window_end,
                    'score': window_variants['score'].mean(),
                    'variant_count': len(window_variants)
                })
            
            window_start += self.step_size
        
        return pd.DataFrame(windowed_data)

    def _generate_snp_windows_single_chrom(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """Generate SNP windows for single chromosome"""
        if variants_df.empty or len(variants_df) < self.snp_window_size:
            return pd.DataFrame()
        
        windowed_data = []
        chrom = variants_df['chrom'].iloc[0]
        chrom_data = variants_df.sort_values('pos')
        
        for i in range(0, len(chrom_data) - self.snp_window_size + 1, self.snp_step_size):
            window_variants = chrom_data.iloc[i:i + self.snp_window_size]
            
            windowed_data.append({
                'chrom': chrom,
                'start': int(window_variants['pos'].min()),
                'end': int(window_variants['pos'].max()) + 1,
                'score': window_variants['score'].mean(),
                'variant_count': len(window_variants)
            })
        
        return pd.DataFrame(windowed_data)



def write_bedgraph(df: pd.DataFrame, output_file: str, track_name: str, 
                  description: str, score_col: str = 'score') -> None:
    """Write DataFrame to bedgraph format"""
    with open(output_file, 'w') as f:
        f.write(f'track type=bedGraph name="{track_name}" description="{description}"\n')
        
        for _, row in tqdm(df.iterrows(), total=len(df), desc=f"Writing {output_file}"):
            if score_col in row:
                if 'start' in row:  # Windowed data
                    f.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row[score_col]:.3f}\n")
                else:  # Per-variant data
                    f.write(f"{row['chrom']}\t{row['pos']-1}\t{row['pos']}\t{row[score_col]:.3f}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Enhanced homozygosity analysis with multiple scoring methods",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('vcf_path', help='Input VCF file')
    parser.add_argument('output_prefix', help='Output file prefix')
    
    # Scorer options
    parser.add_argument('--scorer', choices=['allele_depth', 'mut_vs_ref'], 
                       default='mut_vs_ref', help='Homozygosity scoring method')
    parser.add_argument('-t', '--threshold', type=float, default=0.9,
                       help='Threshold for allele_depth scorer')
    parser.add_argument('-c', '--min-coverage', type=int, default=1,
                       help='Minimum coverage for mut_vs_ref scorer')
    parser.add_argument('--require-all-coverage', action='store_true',
                       help='Require all reference samples meet min coverage')
    
    # Window options
    parser.add_argument('-w', '--window-size', type=int, default=10000,
                       help='Window size (bp or SNP count)')
    parser.add_argument('--step-size', type=int, default=1000,
                       help='Step size (bp or SNP count)')
    parser.add_argument('--use-snp-windows', action='store_true',
                       help='Use SNP-count windows instead of bp windows')
    parser.add_argument('--snp-window-size', type=int, default=200,
                       help='SNP window size')
    parser.add_argument('--snp-step-size', type=int, default=20,
                       help='SNP step size')
    
    # Transform options
    parser.add_argument('-k', type=int, default=5,
                        help='Sigmoid function steepness')
    parser.add_argument('-x', type=float, default=0.3,
                        help='Sigmoid function steepness')
    parser.add_argument('--normalise', action='store_true',
                        help='Normalise sigmoid function (default: False)')
    
    # General options
    parser.add_argument('--all-variants', action='store_true',
                       help='Include all variants (default: SNPs only)')
    parser.add_argument('-i', '--fasta-index', type=str, default=None,
                       help='Path to fasta.fai index (if absent, skip tdf generation)')
    parser.add_argument('-n', '--ncpu', type=int, default=8,
                        help="Number of cpus (ideally one per chr)")
    parser.add_argument('-v', '--verbose', action='store_true')
    
    args = parser.parse_args()
    
    # Set up logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    
    try:
        # Create scorer
        if args.scorer == 'allele_depth':
            scorer = AlleleDepthScorer(threshold=args.threshold)
        elif args.scorer == 'mut_vs_ref':
            scorer = MutantReferenceScorer(min_coverage=args.min_coverage,
                                       require_all_coverage=args.require_all_coverage)
        else:
            raise ValueError(f"Unknown scorer: {args.scorer}")
        
        # Create analyser
        analyser = HomozygosityAnalyser(
            scorer=scorer,
            window_size=args.window_size,
            step_size=args.step_size,
            use_snp_windows=args.use_snp_windows,
            snp_window_size=args.snp_window_size,
            snp_step_size=args.snp_step_size,
            sigmoid={"k": args.k, "x0": args.x0},
            normalise=args.normalise
        )
        
        # Process VCF
        variants_df, windowed_df, sample_name = analyser.process_vcf(
            args.vcf_path, args.all_variants, args.ncpu
        )
        
        if variants_df.empty:
            print("No results to write")
            return
        
        # Write outputs
        variant_file = f"{args.output_prefix}_variants.bedgraph"
        window_file = f"{args.output_prefix}_windows.bedgraph"
        
        write_bedgraph(variants_df, variant_file, 
                      f"{sample_name}_{scorer.get_name()}", 
                      "Per-variant homozygosity scores")
        
        if args.fasta_index:
            variant_tdf = variant_file.replace(".bedgraph", ".tdf")
            bedgraph_to_tdf(variant_file, args.fasta_index, variant_tdf)
        
        if not windowed_df.empty:
            window_type = "SNP" if args.use_snp_windows else "BP"
            write_bedgraph(windowed_df, window_file,
                          f"{sample_name}_{scorer.get_name()}_windowed",
                          f"Windowed homozygosity ({window_type} windows)")
            if args.fasta_index:
                window_tdf = window_file.replace(".bedgraph", ".tdf")
                bedgraph_to_tdf(window_file, args.fasta_index, window_tdf)            
        
        print(f"\nGenerated files:")
        print(f"  - {variant_file}")
        if not windowed_df.empty:
            print(f"  - {window_file}")
            
    except Exception as e:
        print(f"Error: {e}")
        if args.verbose:
            raise

if __name__ == "__main__":
    main()