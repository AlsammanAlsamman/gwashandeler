#!/usr/bin/env python3
"""
merge_dataset_loci.py
Merge all per-table loci files into one table per dataset.
"""

import argparse
import os
import re
import sys
import pandas as pd


def chr_key(val):
    s = str(val)
    if s.isdigit():
        return (0, int(s))
    return (1, s)


def infer_table_name(path):
    name = os.path.basename(path)
    m = re.match(r"(.+)_loci\.tsv$", name)
    return m.group(1) if m else name


def merge_dataset(dataset, inputs, output_file):
    frames = []

    for path in inputs:
        if not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path, sep="\t", low_memory=False)
        except Exception:
            continue

        if df.empty:
            continue

        if "table" not in df.columns:
            df["table"] = infer_table_name(path)
        if "dataset" not in df.columns:
            df["dataset"] = dataset

        frames.append(df)

    if not frames:
        out = pd.DataFrame(columns=[
            "dataset_locus", "dataset", "table", "reference_panel", "chromosome",
            "start", "end", "size_bp", "top_snp", "top_p", "n_sig_snps", "n_ref_matched_snps"
        ])
        out.to_csv(output_file, sep="\t", index=False)
        return

    merged = pd.concat(frames, ignore_index=True, sort=False)

    if "genomic_locus" in merged.columns:
        merged = merged.drop(columns=["genomic_locus"])

    sort_cols = [c for c in ["chromosome", "start", "top_p", "table"] if c in merged.columns]
    if sort_cols:
        merged = merged.sort_values(
            by=sort_cols,
            key=lambda col: col.map(chr_key) if col.name == "chromosome" else col,
        ).reset_index(drop=True)

    merged.insert(0, "dataset_locus", range(1, len(merged) + 1))
    merged.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge per-table loci files into one dataset-level loci table")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    try:
        merge_dataset(args.dataset, args.inputs, args.output)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
