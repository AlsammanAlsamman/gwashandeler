#!/usr/bin/env python3
"""
Build merged loci from FUMA region files and summarize minimum p-values per locus.
"""

import argparse
import glob
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from openpyxl import load_workbook
from openpyxl.styles import PatternFill


RED_FILL = PatternFill(start_color="FFF2CC", end_color="FFF2CC", fill_type="solid")


def parse_args():
    parser = argparse.ArgumentParser(description="Merge FUMA loci and summarize min p-values")
    parser.add_argument("--config", required=True, help="Path to analysis.yml")
    parser.add_argument("--output-xlsx", required=True, help="Output Excel file")
    parser.add_argument("--output-loci", required=True, help="Output merged loci TSV")
    parser.add_argument("--distance", type=int, default=250000, help="Merge gap distance in bp")
    return parser.parse_args()


def read_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def norm_chr(value):
    text = str(value).strip()
    text = re.sub(r"^chr", "", text, flags=re.IGNORECASE)
    return text


def chr_sort_key(value):
    text = norm_chr(value)
    upper = text.upper()
    if upper == "X":
        return (23, "")
    if upper == "Y":
        return (24, "")
    if upper in {"M", "MT"}:
        return (25, "")
    if upper.isdigit():
        return (int(upper), "")
    return (999, upper)


def first_matching_column(columns, candidates):
    lower_map = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def load_fuma_loci(fuma_folder):
    txt_files = sorted(glob.glob(os.path.join(fuma_folder, "*.txt")))
    if not txt_files:
        raise FileNotFoundError(f"No .txt FUMA files found in: {fuma_folder}")

    records = []
    for path in txt_files:
        df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        chr_col = first_matching_column(df.columns, ["chr", "chrom", "chromosome"])
        start_col = first_matching_column(df.columns, ["start", "locus_start", "region_start"])
        end_col = first_matching_column(df.columns, ["end", "locus_end", "region_end"])

        if not chr_col or not start_col or not end_col:
            raise ValueError(
                f"File {path} is missing chr/start/end columns. "
                f"Found columns: {', '.join(df.columns)}"
            )

        sub = df[[chr_col, start_col, end_col]].copy()
        sub.columns = ["chr", "start", "end"]
        sub["chr"] = sub["chr"].map(norm_chr)
        sub["start"] = pd.to_numeric(sub["start"], errors="coerce")
        sub["end"] = pd.to_numeric(sub["end"], errors="coerce")
        sub = sub.dropna(subset=["chr", "start", "end"])

        if sub.empty:
            continue

        # Ensure start <= end for every row.
        sub["start_fixed"] = np.minimum(sub["start"], sub["end"]).astype(int)
        sub["end_fixed"] = np.maximum(sub["start"], sub["end"]).astype(int)
        sub["source_file"] = os.path.basename(path)

        records.append(sub[["chr", "start_fixed", "end_fixed", "source_file"]])

    if not records:
        raise ValueError("No valid loci found across FUMA files.")

    loci = pd.concat(records, ignore_index=True)
    loci = loci.rename(columns={"start_fixed": "start", "end_fixed": "end"})
    return loci


def merge_loci(loci_df, distance_bp):
    merged = []

    loci_df = loci_df.copy()
    loci_df["chr_order"] = loci_df["chr"].map(lambda x: chr_sort_key(x)[0])
    loci_df["chr_label"] = loci_df["chr"].map(lambda x: chr_sort_key(x)[1])
    loci_df = loci_df.sort_values(["chr_order", "chr_label", "start", "end"]).reset_index(drop=True)

    for chr_value, sub in loci_df.groupby("chr", sort=False):
        sub = sub.sort_values(["start", "end"]).reset_index(drop=True)

        cur_start = int(sub.loc[0, "start"])
        cur_end = int(sub.loc[0, "end"])
        sources = {sub.loc[0, "source_file"]}

        for i in range(1, len(sub)):
            s = int(sub.loc[i, "start"])
            e = int(sub.loc[i, "end"])
            src = sub.loc[i, "source_file"]

            if s <= cur_end + distance_bp:
                cur_end = max(cur_end, e)
                sources.add(src)
            else:
                merged.append((chr_value, cur_start, cur_end, ";".join(sorted(sources))))
                cur_start, cur_end = s, e
                sources = {src}

        merged.append((chr_value, cur_start, cur_end, ";".join(sorted(sources))))

    out = pd.DataFrame(merged, columns=["chr", "start", "end", "source_files"])
    out["chr_order"] = out["chr"].map(lambda x: chr_sort_key(x)[0])
    out["chr_label"] = out["chr"].map(lambda x: chr_sort_key(x)[1])
    out = out.sort_values(["chr_order", "chr_label", "start", "end"]).reset_index(drop=True)
    out["locus_code"] = [f"Locus_{i:03d}" for i in range(1, len(out) + 1)]
    out["length_bp"] = out["end"] - out["start"] + 1
    return out[["locus_code", "chr", "start", "end", "length_bp", "source_files"]]


def detect_pvalue_columns(columns):
    pcols = []
    for col in columns:
        low = col.lower()
        if low == "p" or low.startswith("p_") or "pvalue" in low or low.startswith("pval"):
            pcols.append(col)
    return pcols


def summarize_table_min_p(table_path, loci_df):
    head = pd.read_csv(table_path, sep="\t", nrows=0)
    cols = list(head.columns)

    chr_col = first_matching_column(cols, ["chr", "chrom", "chromosome", "chrom_col", "chromosome_name"])
    pos_col = first_matching_column(cols, ["pos", "position", "bp"])
    p_cols = detect_pvalue_columns(cols)

    if not chr_col or not pos_col or not p_cols:
        return {}, []

    usecols = [chr_col, pos_col] + p_cols
    df = pd.read_csv(table_path, sep="\t", usecols=usecols, dtype=str, low_memory=False)
    df[chr_col] = df[chr_col].map(norm_chr)
    df[pos_col] = pd.to_numeric(df[pos_col], errors="coerce")
    df = df.dropna(subset=[chr_col, pos_col]).copy()
    df[pos_col] = df[pos_col].astype(int)

    for pcol in p_cols:
        df[pcol] = pd.to_numeric(df[pcol], errors="coerce")

    indexed = {}
    for chr_value, sub in df.groupby(chr_col):
        sub = sub.sort_values(pos_col)
        entry = {"pos": sub[pos_col].to_numpy()}
        for pcol in p_cols:
            entry[pcol] = sub[pcol].to_numpy(dtype=float)
        indexed[chr_value] = entry

    summary = {pcol: [] for pcol in p_cols}
    for _, locus in loci_df.iterrows():
        chr_value = norm_chr(locus["chr"])
        start = int(locus["start"])
        end = int(locus["end"])

        chr_data = indexed.get(chr_value)
        if chr_data is None:
            for pcol in p_cols:
                summary[pcol].append(np.nan)
            continue

        pos = chr_data["pos"]
        left = np.searchsorted(pos, start, side="left")
        right = np.searchsorted(pos, end, side="right")

        if right <= left:
            for pcol in p_cols:
                summary[pcol].append(np.nan)
            continue

        for pcol in p_cols:
            vals = chr_data[pcol][left:right]
            with np.errstate(all="ignore"):
                minv = np.nanmin(vals) if vals.size else np.nan
            if np.isinf(minv):
                minv = np.nan
            summary[pcol].append(minv)

    return summary, p_cols


def write_excel_with_highlight(df, xlsx_path, highlight_cols, threshold=5e-8):
    Path(os.path.dirname(xlsx_path)).mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="merged_loci_min_p")

    wb = load_workbook(xlsx_path)
    ws = wb["merged_loci_min_p"]

    col_idx = {cell.value: cell.column for cell in ws[1]}
    for col_name in highlight_cols:
        idx = col_idx.get(col_name)
        if idx is None:
            continue
        for row in range(2, ws.max_row + 1):
            value = ws.cell(row=row, column=idx).value
            try:
                if value is not None and float(value) < threshold:
                    ws.cell(row=row, column=idx).fill = RED_FILL
            except (ValueError, TypeError):
                continue

    wb.save(xlsx_path)


def main():
    args = parse_args()
    cfg = read_yaml(args.config)

    fuma_folder = cfg.get("fuma_folder", "input/fumas")
    results_dir = cfg.get("results_dir", "results")
    gwastables = cfg.get("gwastables", [])

    loci_raw = load_fuma_loci(fuma_folder)
    loci_merged = merge_loci(loci_raw, args.distance)

    summary_df = loci_merged.copy()
    highlight_cols = []

    for dataset_entry in gwastables:
        dataset = dataset_entry["name"]
        tables = dataset_entry.get("tables", [])
        for table_entry in tables:
            table = table_entry["table_name"]
            table_path = os.path.join(results_dir, "extracted", dataset, f"{table}.tsv")

            if not os.path.exists(table_path):
                continue

            p_summary, p_cols = summarize_table_min_p(table_path, loci_merged)
            for pcol in p_cols:
                out_col = f"{dataset}__{table}__{pcol}_minp"
                summary_df[out_col] = p_summary[pcol]
                highlight_cols.append(out_col)

    Path(os.path.dirname(args.output_loci)).mkdir(parents=True, exist_ok=True)
    loci_merged.to_csv(args.output_loci, sep="\t", index=False)
    write_excel_with_highlight(summary_df, args.output_xlsx, highlight_cols, threshold=5e-8)

    print(f"Merged loci: {len(loci_merged)}")
    print(f"Loci TSV: {args.output_loci}")
    print(f"Excel summary: {args.output_xlsx}")


if __name__ == "__main__":
    main()
