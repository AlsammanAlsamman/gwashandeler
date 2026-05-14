#!/usr/bin/env python3
"""
identify_loci_per_chromosome.py
Identify loci on a single chromosome from subset GWAS data,
optionally restricted to SNPs present in the filtered reference panel BIM.
"""

import argparse
import os
import re
import sys
import pandas as pd


def find_col(columns, patterns):
    for pattern in patterns:
        for col in columns:
            if re.search(pattern, col, re.IGNORECASE):
                return col
    return None


def load_ref_snp_ids(ref_bim_path):
    if not os.path.exists(ref_bim_path):
        raise FileNotFoundError(f"Reference BIM file not found: {ref_bim_path}")
    bim = pd.read_csv(ref_bim_path, sep=r"\s+", header=None, dtype=str)
    if bim.shape[1] < 2:
        raise ValueError(f"Unexpected BIM format in {ref_bim_path}")
    return set(bim.iloc[:, 1].astype(str).tolist())


def identify_loci(
    subset_file,
    ref_bim,
    output_file,
    dataset,
    table,
    ref_panel,
    lead_p_threshold,
    merge_distance_kb,
    min_snps_per_locus,
    require_refpanel_match,
):
    if not os.path.exists(subset_file):
        raise FileNotFoundError(f"Subset file not found: {subset_file}")

    df = pd.read_csv(subset_file, sep="\t", dtype=str, low_memory=False)
    if df.empty:
        out = pd.DataFrame(columns=[
            "dataset", "table", "reference_panel", "chromosome", "locus_chr_index",
            "start", "end", "size_bp", "top_snp", "top_p", "n_sig_snps", "n_ref_matched_snps"
        ])
        out.to_csv(output_file, sep="\t", index=False)
        return

    chr_col = find_col(df.columns, [r"^CHROM$", r"^CHR$", r"^chrom$", r"^chr$"])
    pos_col = find_col(df.columns, [r"^POS$", r"^BP$", r"^pos$", r"^bp$"])
    p_col = find_col(df.columns, [r"^P$", r"^p$", r"^p_", r"^pval", r"pvalue"])
    id_col = find_col(df.columns, [r"^ID$", r"^rsid$", r"^snpid$", r"^markername$", r"^varid$"])

    missing = [
        name for name, col in [
            ("chromosome", chr_col),
            ("position", pos_col),
            ("pvalue", p_col),
            ("snp_id", id_col),
        ] if col is None
    ]
    if missing:
        raise ValueError(f"Missing required columns in {subset_file}: {', '.join(missing)}")

    df[pos_col] = pd.to_numeric(df[pos_col], errors="coerce")
    df[p_col] = pd.to_numeric(df[p_col], errors="coerce")
    df = df.dropna(subset=[pos_col, p_col, id_col, chr_col]).copy()

    df = df[(df[p_col] > 0) & (df[p_col] <= lead_p_threshold)].copy()

    ref_ids = load_ref_snp_ids(ref_bim)
    if require_refpanel_match:
        df = df[df[id_col].astype(str).isin(ref_ids)].copy()

    if df.empty:
        out = pd.DataFrame(columns=[
            "dataset", "table", "reference_panel", "chromosome", "locus_chr_index",
            "start", "end", "size_bp", "top_snp", "top_p", "n_sig_snps", "n_ref_matched_snps"
        ])
        out.to_csv(output_file, sep="\t", index=False)
        return

    df = df.sort_values(by=[pos_col, p_col]).reset_index(drop=True)
    merge_distance_bp = int(merge_distance_kb * 1000)

    loci = []
    locus_idx = 0
    cur_start = int(df.loc[0, pos_col])
    cur_end = int(df.loc[0, pos_col])
    cur_rows = [0]

    for i in range(1, len(df)):
        pos = int(df.loc[i, pos_col])
        if pos - cur_end <= merge_distance_bp:
            cur_end = max(cur_end, pos)
            cur_rows.append(i)
        else:
            if len(cur_rows) >= min_snps_per_locus:
                locus_idx += 1
                rows = df.loc[cur_rows]
                top = rows.sort_values(by=p_col).iloc[0]
                loci.append({
                    "dataset": dataset,
                    "table": table,
                    "reference_panel": ref_panel,
                    "chromosome": str(top[chr_col]),
                    "locus_chr_index": locus_idx,
                    "start": cur_start,
                    "end": cur_end,
                    "size_bp": cur_end - cur_start + 1,
                    "top_snp": str(top[id_col]),
                    "top_p": float(top[p_col]),
                    "n_sig_snps": len(rows),
                    "n_ref_matched_snps": len(rows),
                })
            cur_start = pos
            cur_end = pos
            cur_rows = [i]

    if len(cur_rows) >= min_snps_per_locus:
        locus_idx += 1
        rows = df.loc[cur_rows]
        top = rows.sort_values(by=p_col).iloc[0]
        loci.append({
            "dataset": dataset,
            "table": table,
            "reference_panel": ref_panel,
            "chromosome": str(top[chr_col]),
            "locus_chr_index": locus_idx,
            "start": cur_start,
            "end": cur_end,
            "size_bp": cur_end - cur_start + 1,
            "top_snp": str(top[id_col]),
            "top_p": float(top[p_col]),
            "n_sig_snps": len(rows),
            "n_ref_matched_snps": len(rows),
        })

    out = pd.DataFrame(loci)
    if out.empty:
        out = pd.DataFrame(columns=[
            "dataset", "table", "reference_panel", "chromosome", "locus_chr_index",
            "start", "end", "size_bp", "top_snp", "top_p", "n_sig_snps", "n_ref_matched_snps"
        ])
    out.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify loci on one chromosome from subset GWAS file")
    parser.add_argument("--subset", required=True)
    parser.add_argument("--ref-bim", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--table", required=True)
    parser.add_argument("--ref-panel", required=True)
    parser.add_argument("--lead-p-threshold", type=float, required=True)
    parser.add_argument("--merge-distance-kb", type=float, required=True)
    parser.add_argument("--min-snps-per-locus", type=int, required=True)
    parser.add_argument("--require-refpanel-match", type=int, choices=[0, 1], required=True)

    args = parser.parse_args()

    try:
        identify_loci(
            subset_file=args.subset,
            ref_bim=args.ref_bim,
            output_file=args.output,
            dataset=args.dataset,
            table=args.table,
            ref_panel=args.ref_panel,
            lead_p_threshold=args.lead_p_threshold,
            merge_distance_kb=args.merge_distance_kb,
            min_snps_per_locus=args.min_snps_per_locus,
            require_refpanel_match=bool(args.require_refpanel_match),
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
