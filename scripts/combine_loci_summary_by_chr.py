#!/usr/bin/env python3
"""
Combine per-chromosome loci summary TSV files into a single TSV/XLSX.
"""

import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import PatternFill

HIGHLIGHT_FILL = PatternFill(start_color="FFF2CC", end_color="FFF2CC", fill_type="solid")


def parse_args():
    parser = argparse.ArgumentParser(description="Combine per-chromosome loci summary files")
    parser.add_argument("--inputs", nargs="+", required=True, help="Per-chromosome summary TSV files")
    parser.add_argument("--output-tsv", required=True, help="Combined TSV output")
    parser.add_argument("--output-xlsx", required=True, help="Combined XLSX output")
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


def write_excel_with_highlight(df, xlsx_path, highlight_cols, threshold=5e-8):
    Path(os.path.dirname(xlsx_path)).mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="dataset_loci_summary")

    wb = load_workbook(xlsx_path)
    ws = wb["dataset_loci_summary"]

    col_idx = {cell.value: cell.column for cell in ws[1]}
    for col_name in highlight_cols:
        idx = col_idx.get(col_name)
        if idx is None:
            continue
        for row in range(2, ws.max_row + 1):
            value = ws.cell(row=row, column=idx).value
            try:
                if value is not None and float(value) < threshold:
                    ws.cell(row=row, column=idx).fill = HIGHLIGHT_FILL
            except (ValueError, TypeError):
                continue

    wb.save(xlsx_path)


def main():
    args = parse_args()

    frames = []
    for path in args.inputs:
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path, sep="\t", low_memory=False)
        if not df.empty:
            frames.append(df)

    if frames:
        out = pd.concat(frames, ignore_index=True)
    else:
        out = pd.DataFrame()

    if not out.empty and "chr" in out.columns and "start" in out.columns:
        out["_chr_order"] = out["chr"].map(lambda x: chr_sort_key(x)[0])
        out["_start_num"] = pd.to_numeric(out["start"], errors="coerce")
        out = out.sort_values(["_chr_order", "_start_num"], na_position="last").reset_index(drop=True)
        out = out.drop(columns=["_chr_order", "_start_num"])

    if not out.empty and "locus_code" in out.columns:
        out["locus_code"] = [f"Locus_{i:04d}" for i in range(1, len(out) + 1)]

    highlight_cols = [
        c for c in out.columns
        if str(c).endswith("__top_p") or str(c).endswith("__top_snp_p")
    ]

    Path(os.path.dirname(args.output_tsv)).mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output_tsv, sep="\t", index=False)
    write_excel_with_highlight(out, args.output_xlsx, highlight_cols, threshold=5e-8)


if __name__ == "__main__":
    main()
