#!/usr/bin/env python3
"""
extract_gwas.py - Extract specified columns from GWAS files using fast file reading.
Uses data.table::fread-equivalent (pandas with optimized reading).
"""

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    print("Error: pandas not installed. Install with: pip install pandas", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Extract specified columns from GWAS TSV file"
    )
    parser.add_argument("--input", required=True, help="Input GWAS TSV file")
    parser.add_argument("--output", required=True, help="Output extracted TSV file")
    parser.add_argument("--columns", required=True, help="Comma-separated column names to extract")
    parser.add_argument("--dataset", required=True, help="Dataset name (for logging)")
    parser.add_argument("--table", required=True, help="Table name (for logging)")
    
    args = parser.parse_args()
    
    # Parse columns
    columns_to_extract = [col.strip() for col in args.columns.split(',')]
    
    print(f"Extracting table '{args.table}' from dataset '{args.dataset}'")
    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}")
    print(f"Columns to extract: {', '.join(columns_to_extract)}")
    
    # Check input file exists
    if not Path(args.input).exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    try:
        # Read file using fast pandas reading (equivalent to R's fread)
        print(f"Reading {args.input}...")
        df = pd.read_csv(
            args.input,
            sep='\t',
            dtype=str,  # Read all as string first to handle flexible formats
            low_memory=False  # Prevent dtype inference warnings
        )
        
        print(f"File loaded: {len(df)} rows, {len(df.columns)} columns")
        print(f"Available columns: {', '.join(df.columns.tolist())}")
        
        # Check if all requested columns exist
        missing_cols = set(columns_to_extract) - set(df.columns)
        if missing_cols:
            print(f"Error: Missing columns in {args.dataset}/{args.table}: {', '.join(missing_cols)}", 
                  file=sys.stderr)
            print(f"Available columns: {', '.join(df.columns)}", file=sys.stderr)
            sys.exit(1)
        
        # Extract and convert numeric columns where appropriate
        df_extracted = df[columns_to_extract].copy()
        
        # Try to convert numeric columns
        for col in columns_to_extract:
            if col.lower() in ['p', 'p_value', 'pval', 'pos', 'position', 'bp', 
                               'beta', 'se', 'or', 'z_stat', 'n', 'eaf']:
                try:
                    df_extracted[col] = pd.to_numeric(df_extracted[col], errors='coerce')
                except:
                    pass  # Keep as is if conversion fails
        
        # Write output
        print(f"Writing {len(df_extracted)} rows to {args.output}...")
        df_extracted.to_csv(args.output, sep='\t', index=False, na_rep='NA')
        
        print(f"Extraction complete. Output file: {args.output}")
        print(f"Output size: {len(df_extracted)} rows x {len(df_extracted.columns)} columns")
        
    except Exception as e:
        print(f"Error during extraction: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
