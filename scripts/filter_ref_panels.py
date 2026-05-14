#!/usr/bin/env python3
"""
filter_ref_panels.py
Filters reference panels using PLINK to keep only SNPs present in subset GWAS data.
"""

import sys
import argparse
import os
import subprocess

def filter_ref_panel(snp_list_file, ref_panel_path, chrom, output_prefix):
    """
    Use PLINK to filter reference panel by SNP list.
    Input: ref_panel_path (without extension, e.g. /path/to/panel)
    Output: {output_prefix}.bed/bim/fam
    """
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    
    print(f"Filtering reference panel for chromosome {chrom}")
    print(f"SNP list: {snp_list_file}")
    print(f"Ref panel: {ref_panel_path}")
    print(f"Output: {output_prefix}")

    # Verify input files exist
    if not os.path.exists(snp_list_file):
        raise FileNotFoundError(f"SNP list not found: {snp_list_file}")
    
    if not os.path.exists(f"{ref_panel_path}.bed"):
        raise FileNotFoundError(f"Ref panel .bed not found: {ref_panel_path}.bed")
    if not os.path.exists(f"{ref_panel_path}.bim"):
        raise FileNotFoundError(f"Ref panel .bim not found: {ref_panel_path}.bim")
    if not os.path.exists(f"{ref_panel_path}.fam"):
        raise FileNotFoundError(f"Ref panel .fam not found: {ref_panel_path}.fam")
    
    # Count SNPs in list
    with open(snp_list_file, 'r') as f:
        snp_count = sum(1 for _ in f)
    print(f"Extracted {snp_count} SNPs from list")
    
    # Run PLINK to filter; prefer plink2 and fall back to plink.
    plink_candidates = ["plink2", "plink"]
    run_ok = False
    last_error = None

    for plink_bin in plink_candidates:
        plink_cmd = [
            plink_bin,
            '--bfile', ref_panel_path,
            '--extract', snp_list_file,
            '--make-bed',
            '--out', output_prefix,
            '--threads', '4'
        ]

        print(f"Running: {' '.join(plink_cmd)}")
        try:
            result = subprocess.run(plink_cmd, capture_output=True, text=True)
        except FileNotFoundError as e:
            last_error = e
            print(f"{plink_bin} not found in PATH, trying next candidate...")
            continue

        if result.returncode == 0:
            print("PLINK stdout:", result.stdout)
            run_ok = True
            break

        last_error = RuntimeError(f"{plink_bin} failed with return code {result.returncode}: {result.stderr}")
        print(f"{plink_bin} failed, trying next candidate...")

    if not run_ok:
        raise RuntimeError(f"PLINK filtering failed for all candidates: {last_error}")
    
    # Verify output files exist
    for ext in ['.bed', '.bim', '.fam']:
        if not os.path.exists(f"{output_prefix}{ext}"):
            raise FileNotFoundError(f"Output file not created: {output_prefix}{ext}")
    
    print(f"Filtered reference panel created: {output_prefix}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter reference panel by SNP list using PLINK")
    parser.add_argument("--snp-list", required=True, help="SNP list file (one per line)")
    parser.add_argument("--ref-panel", required=True, help="Reference panel path (without extension)")
    parser.add_argument("--chrom", required=True, help="Chromosome number")
    parser.add_argument("--output", required=True, help="Output prefix for PLINK files")
    
    args = parser.parse_args()
    
    try:
        filter_ref_panel(args.snp_list, args.ref_panel, args.chrom, args.output)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
