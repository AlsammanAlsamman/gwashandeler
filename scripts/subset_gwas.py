#!/usr/bin/env python3
"""
subset_gwas.py
Filters GWAS data by p-value threshold, keeping only variants with p <= threshold
"""

import sys
import argparse
import pandas as pd
import re

def find_pvalue_column(df, column_names):
    """
    Identify p-value column from DataFrame columns.
    Tries common patterns: p_*, pval*, P, p, etc.
    """
    patterns = [
        r'^p_',
        r'^pval',
        r'^P$',
        r'^p$'
    ]
    
    for pattern in patterns:
        matching = [col for col in column_names if re.search(pattern, col, re.IGNORECASE)]
        if matching:
            return matching[0]
    
    raise ValueError(f"Could not find p-value column in: {column_names}")

def subset_gwas(input_file, output_dir, threshold, dataset, table):
    """
    Read GWAS file, filter by p-value threshold, split by chromosome, write per-chromosome files.
    """
    print(f"Reading GWAS data from: {input_file}")
    df = pd.read_csv(input_file, sep='\t', dtype=str, low_memory=False)
    
    print(f"Data loaded: {len(df)} rows, {len(df.columns)} columns")
    print(f"Available columns: {', '.join(df.columns)}")
    
    # Find p-value column
    pval_col = find_pvalue_column(df, df.columns)
    print(f"Using p-value column: {pval_col}")
    
    # Find chromosome column
    chr_patterns = [r'^chrom', r'^chr', r'^CHROM', r'^CHR']
    chr_col = None
    for pattern in chr_patterns:
        matching = [col for col in df.columns if re.search(pattern, col, re.IGNORECASE)]
        if matching:
            chr_col = matching[0]
            break
    
    if chr_col is None:
        raise ValueError(f"Could not find chromosome column in: {', '.join(df.columns)}")
    print(f"Using chromosome column: {chr_col}")
    
    # Convert p-value to numeric
    df[pval_col] = pd.to_numeric(df[pval_col], errors='coerce')
    
    # Filter by threshold
    df_subset = df[df[pval_col] <= threshold].copy()
    
    print(f"After filtering (p-value <= {threshold}): {len(df_subset)} rows")
    print(f"Removed: {len(df) - len(df_subset)} rows")
    
    # Split by chromosome and write separate files
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    chromosomes = sorted(df_subset[chr_col].unique())
    print(f"Chromosomes in subset: {chromosomes}")
    
    output_files = []
    for chrom in chromosomes:
        df_chr = df_subset[df_subset[chr_col] == chrom].copy()
        output_file = os.path.join(output_dir, f"{dataset}_{table}_chr{chrom}_subset.tsv")
        print(f"Writing {len(df_chr)} rows to: {output_file}")
        df_chr.to_csv(output_file, sep='\t', index=False, na_rep='NA')
        output_files.append(output_file)
    
    print(f"Subsetting complete for {dataset}/{table}: {len(output_files)} chromosome files created")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subset GWAS by p-value threshold, split by chromosome")
    parser.add_argument("--input", required=True, help="Input extracted GWAS TSV file")
    parser.add_argument("--output-dir", required=True, help="Output directory for chromosome-split TSV files")
    parser.add_argument("--threshold", type=float, default=5e-3, help="P-value threshold (default 5e-3)")
    parser.add_argument("--dataset", required=True, help="Dataset name for logging")
    parser.add_argument("--table", required=True, help="Table name for logging")
    
    args = parser.parse_args()
    
    try:
        subset_gwas(args.input, args.output_dir, args.threshold, args.dataset, args.table)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
