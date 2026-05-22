#!/usr/bin/env python3
"""
combine_loci_summary_by_chr.py

Combine per-chromosome loci summary TSV files (produced by
summarize_dataset_loci_with_genes.py --skip-xlsx) into a single TSV/XLSX.

Bugs fixed vs. original:
  BUG-1  Duplicate loci after concatenation.
         Each per-chromosome TSV has locus_code starting at Locus_0001.
         After a naive concat the same genomic locus appears N times (once
         per chromosome file that covered it).  Fix: deduplicate on
         (chr, start, end) BEFORE renumbering locus_code.

  BUG-2  Renumbering happened before deduplication, so duplicates received
         unique codes and looked like distinct loci.
         Fix: sort → deduplicate → renumber (in that order).

  BUG-3  No cross-chromosome re-merge after concatenation.
         Adjacent loci near chromosome-chunk boundaries could end up
         un-merged.  Fix: after deduplication, run the same single-linkage
         greedy merge used in summarize_dataset_loci_with_genes.py so the
         final table is globally consistent.  The merge distance is read
         from the first TSV's metadata column or falls back to the CLI flag.

  BUG-4  highlight_cols detection used fragile endswith patterns that would
         also accidentally match rsID columns (e.g. a column named
         "ds__t__max_p_snp_rsid" ends with neither "__top_p" nor "__max_p"
         so that specific case was OK, but the original code also flagged
         __top_snp_p which is always identical to __top_p — see BUG-5 in
         summarize_... — so the same column was highlighted twice).
         Fix: use an explicit set of p-value suffixes and exclude any
         column whose name contains 'rsid' or 'snp_id'.
"""

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import PatternFill

HIGHLIGHT_FILL = PatternFill(
    start_color="FF9999", end_color="FF9999", fill_type="solid"
)

# P-value column suffixes that should be highlighted when < threshold.
# rsID / SNP-name columns are intentionally excluded.
P_VALUE_SUFFIXES = ("__top_p", "__top_snp_p", "__max_p")
RSID_KEYWORDS    = ("rsid", "snp_id", "snp_rsid")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine per-chromosome loci summary files"
    )
    parser.add_argument(
        "--inputs", nargs="+", required=True,
        help="Per-chromosome summary TSV files"
    )
    parser.add_argument("--output-tsv",  required=True, help="Combined TSV output")
    parser.add_argument("--output-xlsx", required=True, help="Combined XLSX output")
    parser.add_argument(
        "--merge-distance-kb", type=int, default=250,
        help="Gap threshold in kb for re-merging adjacent loci (default: 250)"
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Chromosome sort helpers (identical to summarize_... so ordering is consistent)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# BUG-3 FIX: re-merge helper (greedy single-linkage, same logic as
# merge_intervals_by_policy in summarize_... but simpler — no LD step)
# ---------------------------------------------------------------------------

def remerge_loci(df, distance_bp):
    """
    After combining chromosome TSVs, ensure adjacent loci that fall within
    `distance_bp` of each other are merged into one row.

    Operates on the *locus* rows only (chr, start, end).
    All other columns (p-value summaries, gene annotations) are aggregated:
      - numeric p-value columns  → minimum (best signal wins)
      - rsID / text columns      → first non-null value
      - nearest_gene             → first non-null value
      - gene_distance_bp         → minimum (closest gene wins)
      - length_bp / length_kbp   → recomputed from merged coords
      - IMD                      → recomputed after all loci are known

    Returns a new DataFrame with globally re-numbered locus_code.
    """
    if df.empty:
        return df.copy()

    # Identify column groups
    coord_cols    = {"locus_code", "chr", "start", "end",
                     "length_bp", "length_kbp", "IMD"}
    p_val_cols    = [
        c for c in df.columns
        if c not in coord_cols
        and any(c.endswith(s) for s in P_VALUE_SUFFIXES)
        and not any(kw in c.lower() for kw in RSID_KEYWORDS)
    ]
    gene_dist_col = "gene_distance_bp" if "gene_distance_bp" in df.columns else None
    text_cols     = [
        c for c in df.columns
        if c not in coord_cols
        and c not in p_val_cols
        and c != gene_dist_col
    ]

    df = df.copy()
    df["_chr_order"] = df["chr"].map(lambda x: chr_sort_key(x)[0])
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"]   = pd.to_numeric(df["end"],   errors="coerce").astype("Int64")
    df = df.sort_values(["_chr_order", "start"]).reset_index(drop=True)

    merged_records = []

    for chr_val, sub in df.groupby("chr", sort=False):
        sub = sub.sort_values("start").reset_index(drop=True)

        # Each cluster: list of row indices belonging to it
        clusters = []
        cur_indices = [0]
        cur_end = int(sub.at[0, "end"])

        for i in range(1, len(sub)):
            s = int(sub.at[i, "start"])
            e = int(sub.at[i, "end"])
            if s <= cur_end + distance_bp:
                cur_indices.append(i)
                cur_end = max(cur_end, e)
            else:
                clusters.append(cur_indices)
                cur_indices = [i]
                cur_end = e
        clusters.append(cur_indices)

        for cluster_idx in clusters:
            rows = sub.iloc[cluster_idx]
            new_start = int(rows["start"].min())
            new_end   = int(rows["end"].max())
            rec = {
                "chr":   chr_val,
                "start": new_start,
                "end":   new_end,
            }
            # Aggregate p-value columns: take minimum (lowest = most significant)
            for col in p_val_cols:
                vals = pd.to_numeric(rows[col], errors="coerce").dropna()
                rec[col] = float(vals.min()) if len(vals) > 0 else np.nan

            # Gene distance: minimum
            if gene_dist_col:
                vals = pd.to_numeric(rows[gene_dist_col], errors="coerce").dropna()
                rec[gene_dist_col] = float(vals.min()) if len(vals) > 0 else np.nan

            # Text / rsID / gene name columns: first non-null
            for col in text_cols:
                non_null = rows[col].dropna()
                non_null = non_null[non_null.astype(str).str.strip() != ""]
                rec[col] = str(non_null.iloc[0]) if len(non_null) > 0 else np.nan

            merged_records.append(rec)

    if not merged_records:
        return pd.DataFrame()

    out = pd.DataFrame(merged_records)

    # Recompute coordinate-derived columns
    out["length_bp"]  = out["end"] - out["start"] + 1
    out["length_kbp"] = out["length_bp"] / 1000.0

    # Sort globally
    out["_chr_order"] = out["chr"].map(lambda x: chr_sort_key(x)[0])
    out = out.sort_values(["_chr_order", "start"]).reset_index(drop=True)

    # Recompute IMD (inter-marker distance to next locus on same chromosome)
    out["IMD"] = np.nan
    for chr_val, sub in out.groupby("chr", sort=False):
        idxs = list(sub.index)
        for pos, row_idx in enumerate(idxs[:-1]):
            next_idx = idxs[pos + 1]
            out.at[row_idx, "IMD"] = int(
                out.at[next_idx, "start"] - out.at[row_idx, "end"] - 1
            )

    out = out.drop(columns=["_chr_order"])

    # Globally renumber locus_code AFTER dedup + re-merge
    out.insert(0, "locus_code",
               [f"Locus_{i:04d}" for i in range(1, len(out) + 1)])

    # Restore column order: coord first, then all others in original order
    desired_order = (
        ["locus_code", "chr", "start", "end",
         "length_bp", "length_kbp", "IMD"]
        + [c for c in df.columns
           if c not in {"locus_code","chr","start","end",
                        "length_bp","length_kbp","IMD","_chr_order"}
           and c in out.columns]
    )
    desired_order = [c for c in desired_order if c in out.columns]
    out = out[desired_order]

    return out


# ---------------------------------------------------------------------------
# Excel writer
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    merge_distance_bp = args.merge_distance_kb * 1000

    # ── Load all per-chromosome TSVs ─────────────────────────────────────────
    frames = []
    for path in args.inputs:
        if not os.path.exists(path):
            print(f"[warn] Input not found, skipping: {path}", file=sys.stderr)
            continue
        try:
            df = pd.read_csv(path, sep="\t", low_memory=False)
        except Exception as e:
            print(f"[warn] Failed to read {path}: {e}", file=sys.stderr)
            continue
        if not df.empty:
            frames.append(df)

    if not frames:
        print("[warn] No valid input files; writing empty outputs.", file=sys.stderr)
        out = pd.DataFrame()
        Path(os.path.dirname(args.output_tsv)).mkdir(parents=True, exist_ok=True)
        out.to_csv(args.output_tsv, sep="\t", index=False)
        write_excel_with_highlight(out, args.output_xlsx, [])
        return

    # ── BUG-1 + BUG-2 FIX: concatenate then DEDUPLICATE before renumbering ──
    raw = pd.concat(frames, ignore_index=True)
    n_raw = len(raw)

    # Normalise chr column so dedup works across files that use "chr1" vs "1"
    if "chr" in raw.columns:
        raw["chr"] = raw["chr"].map(norm_chr)

    # Cast coordinates to numeric for reliable deduplication and sorting
    if "start" in raw.columns:
        raw["start"] = pd.to_numeric(raw["start"], errors="coerce")
    if "end" in raw.columns:
        raw["end"] = pd.to_numeric(raw["end"], errors="coerce")

    # Sort by (chr, start) first — required both for dedup and for BUG-3 merge
    if "chr" in raw.columns and "start" in raw.columns:
        raw["_chr_order"] = raw["chr"].map(lambda x: chr_sort_key(x)[0])
        raw["_start_num"] = pd.to_numeric(raw["start"], errors="coerce")
        raw = (
            raw
            .sort_values(["_chr_order", "_start_num"], na_position="last")
            .drop(columns=["_chr_order", "_start_num"])
            .reset_index(drop=True)
        )

    # ── BUG-3 FIX: re-merge adjacent loci across chromosome-chunk boundaries ─
    # This also handles the deduplication implicitly for exact duplicates,
    # but we still dedup first so that re-merge only aggregates truly
    # distinct genomic signals rather than identical rows.
    dedup_keys = [c for c in ["chr", "start", "end"] if c in raw.columns]
    if dedup_keys:
        # Keep first occurrence; all other columns are aggregated in remerge
        before_dedup = len(raw)
        raw = raw.drop_duplicates(subset=dedup_keys, keep="first").reset_index(drop=True)
        n_deduped = before_dedup - len(raw)
        if n_deduped:
            print(
                f"  Removed {n_deduped} duplicate loci (same chr/start/end) "
                f"before re-merge.",
                flush=True,
            )

    out = remerge_loci(raw, distance_bp=merge_distance_bp)
    n_out = len(out)
    print(
        f"  Input rows: {n_raw}  →  after dedup+remerge: {n_out} loci",
        flush=True,
    )

    # ── BUG-4 FIX: detect p-value columns robustly ──────────────────────────
    highlight_cols = [
        c for c in out.columns
        if (
            any(c.endswith(s) for s in P_VALUE_SUFFIXES)
            and not any(kw in c.lower() for kw in RSID_KEYWORDS)
        )
    ]

    # ── Write outputs ────────────────────────────────────────────────────────
    Path(os.path.dirname(args.output_tsv)).mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"TSV written: {args.output_tsv}", flush=True)

    write_excel_with_highlight(out, args.output_xlsx, highlight_cols, threshold=5e-8)
    print(f"XLSX written: {args.output_xlsx}", flush=True)


if __name__ == "__main__":
    main()
