#!/usr/bin/env python3
"""
Combine per-chromosome loci details TSV files into a single audit TSV.
"""

import argparse
import os
import re
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Combine per-chromosome loci details files")
    parser.add_argument("--inputs", nargs="+", required=True, help="Per-chromosome details TSV files")
    parser.add_argument("--output-tsv", required=True, help="Combined details TSV output")
    return parser.parse_args()


def norm_chr(value):
    return re.sub(r"^chr", "", str(value).strip(), flags=re.IGNORECASE)


def chr_sort_key(value):
    text = norm_chr(value).upper()
    if text == "X":
        return (23, "")
    if text == "Y":
        return (24, "")
    if text in {"M", "MT"}:
        return (25, "")
    if text.isdigit():
        return (int(text), "")
    return (999, text)


def main():
    args = parse_args()

    frames = []
    for path in args.inputs:
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path, sep="\t", low_memory=False)
        if not df.empty:
            frames.append(df)

    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    if not out.empty and "chr" in out.columns:
        out["_chr_order"] = out["chr"].map(lambda x: chr_sort_key(x)[0] if pd.notna(x) else 999)

        numeric_sort_col = None
        for c in ["start", "left_start", "source_start"]:
            if c in out.columns:
                out["_start_num"] = pd.to_numeric(out[c], errors="coerce")
                numeric_sort_col = "_start_num"
                break

        sort_cols = ["_chr_order"]
        if "record_type" in out.columns:
            sort_cols.append("record_type")
        if numeric_sort_col:
            sort_cols.append(numeric_sort_col)

        out = out.sort_values(sort_cols, na_position="last").reset_index(drop=True)

        drop_cols = [c for c in ["_chr_order", "_start_num"] if c in out.columns]
        if drop_cols:
            out = out.drop(columns=drop_cols)

    Path(os.path.dirname(args.output_tsv)).mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
