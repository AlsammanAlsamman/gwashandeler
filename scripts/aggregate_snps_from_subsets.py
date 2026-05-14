#!/usr/bin/env python3
"""
aggregate_snps_from_subsets.py
Extracts SNPs (rsids) from all subset GWAS files and creates unique sorted lists per chromosome.
"""

import sys
import argparse
import pandas as pd
import os
import glob
import re


def _find_column(columns, patterns):
    for pattern in patterns:
        for col in columns:
            if re.search(pattern, col, re.IGNORECASE):
                return col
    return None

def aggregate_snps(subset_dir, output_dir):
    """
    Scan all chromosome-split subset GWAS files, extract unique SNPs per chromosome.
    Output: results/snp_lists/chr{N}_snps.txt
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all subset files
    subset_files = glob.glob(os.path.join(subset_dir, "*", "*", "*chr*_subset.tsv"))
    print(f"Found {len(subset_files)} subset GWAS files")
    
    if len(subset_files) == 0:
        raise ValueError(f"No subset files found in {subset_dir}")
    
    # Dictionary to store SNPs per chromosome
    snps_by_chr = {}
    all_chromosomes = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
    
    for filepath in subset_files:
        print(f"Processing: {filepath}")
        try:
            df = pd.read_csv(filepath, sep='\t', dtype=str, low_memory=False)

            chrom_col = _find_column(df.columns, [r'^CHROM$', r'^CHR$', r'^chrom$', r'^chr$'])
            snp_col = _find_column(df.columns, [r'^ID$', r'^rsid$', r'^snpid$', r'^markername$', r'^varid$'])

            if chrom_col is None or snp_col is None:
                raise ValueError(f"Could not detect chromosome/SNP columns in {filepath}")

            for _, row in df[[chrom_col, snp_col]].dropna().iterrows():
                chrom = str(row[chrom_col]).strip()
                snp_id = str(row[snp_col]).strip()

                if not chrom or not snp_id or snp_id.upper() == 'NA':
                    continue

                if chrom not in snps_by_chr:
                    snps_by_chr[chrom] = set()
                snps_by_chr[chrom].add(snp_id)
        except Exception as e:
            print(f"Warning: Error processing {filepath}: {e}")
            continue
    
    # Write unique sorted SNP lists per chromosome, including empty files for missing chromosomes
    output_files = []
    for chrom in all_chromosomes:
        snps = sorted(list(snps_by_chr.get(chrom, set())))
        output_file = os.path.join(output_dir, f"chr{chrom}_snps.txt")
        with open(output_file, 'w') as f:
            for snp in snps:
                f.write(f"{snp}\n")
        print(f"Wrote {len(snps)} unique SNPs to {output_file}")
        output_files.append(output_file)
    
    print(f"SNP aggregation complete: {len(output_files)} chromosome files created")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate unique SNPs from subset GWAS files by chromosome")
    parser.add_argument("--subset-dir", required=True, help="Base directory containing subset GWAS files")
    parser.add_argument("--output-dir", required=True, help="Output directory for SNP lists")
    
    args = parser.parse_args()
    
    try:
        aggregate_snps(args.subset_dir, args.output_dir)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
