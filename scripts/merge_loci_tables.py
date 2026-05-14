#!/usr/bin/env python3
"""
merge_loci_tables.py
Merge chromosome-level loci files into one final table per dataset/table.
"""

import argparse
import glob
import os
import re
import sys
import pandas as pd


def chr_key(val):
    s = str(val)
    if s.isdigit():
        return (0, int(s))
    return (1, s)


def merge_loci(dataset, table, input_dir, output_file):
    pattern = os.path.join(input_dir, "chr*_loci.tsv")
    files = sorted(glob.glob(pattern))

    if not files:
        out = pd.DataFrame(columns=[
            "genomic_locus", "dataset", "table", "reference_panel", "chromosome",
            "start", "end", "size_bp", "top_snp", "top_p", "n_sig_snps", "n_ref_matched_snps"
        ])
        out.to_csv(output_file, sep="\t", index=False)
        return

    frames = []
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t")
            if not df.empty:
                frames.append(df)
        except Exception:
            continue

    if not frames:
        out = pd.DataFrame(columns=[
            "genomic_locus", "dataset", "table", "reference_panel", "chromosome",
            "start", "end", "size_bp", "top_snp", "top_p", "n_sig_snps", "n_ref_matched_snps"
        ])
        out.to_csv(output_file, sep="\t", index=False)
        return

    merged = pd.concat(frames, ignore_index=True)
    merged = merged.sort_values(
        by=["chromosome", "start", "top_p"],
        key=lambda col: col.map(chr_key) if col.name == "chromosome" else col,
    ).reset_index(drop=True)

    merged.insert(0, "genomic_locus", range(1, len(merged) + 1))
    merged.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge chromosome loci tables")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--table", required=True)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    try:
        merge_loci(args.dataset, args.table, args.input_dir, args.output)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
