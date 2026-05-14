#!/usr/bin/env python3
"""
summarize_dataset_loci_with_genes.py

Reads per-dataset merged loci TSV files produced by merge_dataset_loci.py
(rule 06), pools every locus interval across all datasets/tables, merges
intervals that overlap or are within a configurable gap (default 250 kb),
annotates each merged locus with the nearest gene from a MAGMA-format
.gene.loc file, and for every dataset x table records the best (minimum)
top_p within the merged locus.

Outputs:
  - a merged-loci TSV  (--output-tsv)
  - a highlighted Excel workbook (--output-xlsx)
    p-value columns are highlighted when value < 5e-8
"""

import argparse
import os
import re
import sys
import subprocess
import tempfile
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from openpyxl import load_workbook
from openpyxl.styles import PatternFill

HIGHLIGHT_FILL = PatternFill(start_color="FFF2CC", end_color="FFF2CC", fill_type="solid")


# ---------------------------------------------------------------------------
# Helper: Find PLINK executable
# ---------------------------------------------------------------------------

def find_plink_executable():
    """
    Try to find PLINK executable in PATH and common locations.
    Returns path to plink or None if not found.
    """
    # Try PATH first
    plink_path = shutil.which("plink")
    if plink_path:
        return plink_path
    
    # Try common locations
    common_paths = [
        "/usr/bin/plink",
        "/usr/local/bin/plink",
        "/opt/plink/plink",
        os.path.expanduser("~/bin/plink"),
        os.path.expanduser("~/plink/plink"),
    ]
    for p in common_paths:
        if os.path.exists(p) and os.access(p, os.X_OK):
            return p
    
    return None


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Merge dataset loci across all datasets, annotate with nearest gene, "
            "and summarise best p-value per dataset x table."
        )
    )
    parser.add_argument("--config", required=True, help="Path to analysis.yml")
    parser.add_argument(
        "--inputs", nargs="+", required=True,
        help="Per-dataset merged loci TSV files (one per dataset, from merge_dataset_loci.py)"
    )
    parser.add_argument("--output-tsv", required=True, help="Output merged-loci TSV")
    parser.add_argument("--output-xlsx", required=False, help="Output Excel workbook")
    parser.add_argument(
        "--distance", type=int, default=250000,
        help="Maximum gap in bp between adjacent loci to be merged (default: 250000)"
    )
    parser.add_argument(
        "--use-ld", action="store_true",
        help="Use LD-based criteria for merging (requires SNP data and PLINK)"
    )
    parser.add_argument(
        "--chromosome", required=False,
        help="Optional chromosome filter (e.g. 1..22, X)."
    )
    parser.add_argument(
        "--skip-xlsx", action="store_true",
        help="Skip writing Excel output (useful for per-chromosome parallel jobs)."
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_yaml(path):
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def norm_chr(value):
    """Strip leading 'chr' prefix and whitespace."""
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
# Step 1 - load all per-dataset loci
# ---------------------------------------------------------------------------

def load_all_loci(input_files):
    """Concatenate all per-dataset loci TSV files into one DataFrame."""
    frames = []
    for path in input_files:
        if not os.path.exists(path):
            print(f"[warn] Input not found, skipping: {path}", file=sys.stderr)
            continue
        df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        if df.empty:
            continue
        frames.append(df)

    if not frames:
        raise ValueError(
            "No valid loci found across all input files. "
            "Check that rule 06 completed successfully."
        )

    all_loci = pd.concat(frames, ignore_index=True)

    all_loci["chromosome"] = all_loci["chromosome"].map(norm_chr)
    all_loci["start"] = pd.to_numeric(all_loci["start"], errors="coerce")
    all_loci["end"] = pd.to_numeric(all_loci["end"], errors="coerce")
    all_loci["top_p"] = pd.to_numeric(all_loci["top_p"], errors="coerce")
    all_loci = all_loci.dropna(subset=["chromosome", "start", "end"])
    all_loci["start"] = all_loci["start"].astype(int)
    all_loci["end"] = all_loci["end"].astype(int)

    return all_loci


# ---------------------------------------------------------------------------
# Step 1b - load SNPs for each locus (for LD-based merging)
# ---------------------------------------------------------------------------

def load_snps_for_loci():
    """
    Load SNP data from all GWAS files in input/ directory.
    Handles multiple column name formats (CHR/chrom, POS/pos, ID/rsid/snpid/varid).
    Files with .gwas.tsv, .tsv suffixes are loaded.
    Returns dict: {dataset: DataFrame with CHR, POS, ID columns}
    """
    loci_snps = {}
    input_dir = "input"
    
    if not os.path.exists(input_dir):
        print(f"[warn] Input directory not found: {input_dir}", file=sys.stderr)
        return loci_snps
    
    # Column name mappings for different file formats
    chr_variants = ["CHR", "CHROM", "chrom", "Chr", "chromosome"]
    pos_variants = ["POS", "POS_BUILD37", "pos", "Pos", "position"]
    id_variants = ["ID", "SNP", "SNPID", "snpid", "rsid", "varid", "markername", "varid"]
    
    # Scan for GWAS files
    gwas_files = []
    for fname in os.listdir(input_dir):
        if fname.endswith('.gwas.tsv') or (fname.endswith('.tsv') and not fname.startswith('.')):
            full_path = os.path.join(input_dir, fname)
            if os.path.isfile(full_path):
                dataset_name = fname.replace('.gwas.tsv', '').replace('.tsv', '')
                gwas_files.append((dataset_name, full_path))
    
    for dataset, gwas_file in sorted(gwas_files):
        try:
            gwas_df = pd.read_csv(gwas_file, sep="\t", dtype=str, low_memory=False, nrows=1000000)
        except Exception as e:
            print(f"[warn] Failed to read GWAS file {gwas_file}: {e}", file=sys.stderr)
            continue

        # Find the correct column names for CHR, POS, ID
        chr_col = None
        pos_col = None
        id_col = None
        
        for col in gwas_df.columns:
            if col in chr_variants and chr_col is None:
                chr_col = col
            if col in pos_variants and pos_col is None:
                pos_col = col
            if col in id_variants and id_col is None:
                id_col = col
        
        # Skip files missing required columns
        if not chr_col or not pos_col:
            print(f"[warn] GWAS file {gwas_file} missing chromosome/position columns", file=sys.stderr)
            continue
        
        if not id_col:
            print(f"[warn] GWAS file {gwas_file} missing SNP ID column (tried: {', '.join(id_variants)})", file=sys.stderr)
            continue
        
        # Extract and normalize
        try:
            snps_df = gwas_df[[chr_col, pos_col, id_col]].copy()
            snps_df.columns = ["CHR", "POS", "ID"]
            snps_df["CHR"] = snps_df["CHR"].astype(str).map(norm_chr)
            snps_df["POS"] = pd.to_numeric(snps_df["POS"], errors="coerce")
            snps_df = snps_df.dropna(subset=["CHR", "POS", "ID"])
            snps_df["POS"] = snps_df["POS"].astype(int)
            snps_df = snps_df[snps_df["ID"].astype(str).str.strip() != ""]
            
            if snps_df.empty:
                print(f"[warn] No valid SNPs found in {gwas_file}", file=sys.stderr)
                continue
            
            loci_snps[dataset] = snps_df[["CHR", "POS", "ID"]].copy()
            print(f"  Loaded {len(snps_df)} SNPs from {dataset}", flush=True)
        except Exception as e:
            print(f"[warn] Error processing SNPs from {gwas_file}: {e}", file=sys.stderr)
            continue

    return loci_snps


def collect_boundary_snps(loci_snps, chr_val, start, end, n_snps=10, at_end=True):
    """
    Collect boundary SNP IDs across all loaded dataset SNP frames.
    Returns up to n_snps unique SNP IDs from the requested boundary.
    """
    frames = []
    chr_norm = norm_chr(chr_val)

    for snps_df in loci_snps.values():
        if snps_df is None or snps_df.empty:
            continue
        sub = snps_df[snps_df["CHR"] == chr_norm]
        sub = sub[(sub["POS"] >= start) & (sub["POS"] <= end)]
        if not sub.empty:
            frames.append(sub[["POS", "ID"]].copy())

    if not frames:
        return []

    merged = pd.concat(frames, ignore_index=True)
    merged = merged.dropna(subset=["ID"]).drop_duplicates(subset=["ID"]) 
    merged = merged.sort_values("POS").reset_index(drop=True)

    if at_end:
        return merged["ID"].tail(n_snps).tolist()
    return merged["ID"].head(n_snps).tolist()


# ---------------------------------------------------------------------------
# Step 2 - merge overlapping / near loci
# ---------------------------------------------------------------------------

def merge_loci(all_loci_df, distance_bp):
    """
    Pool every (chr, start, end) interval and merge those that overlap or
    are separated by at most distance_bp.  Returns a DataFrame with columns:
      locus_code, chr, start, end, length_bp, length_kbp, IMD
    """
    df = all_loci_df[["chromosome", "start", "end"]].copy()
    df["chr_order"] = df["chromosome"].map(lambda x: chr_sort_key(x)[0])
    df = df.sort_values(["chr_order", "chromosome", "start", "end"]).reset_index(drop=True)

    records = []
    for chr_value, sub in df.groupby("chromosome", sort=False):
        sub = sub.sort_values("start").reset_index(drop=True)
        cur_start = int(sub.loc[0, "start"])
        cur_end = int(sub.loc[0, "end"])

        for i in range(1, len(sub)):
            s = int(sub.loc[i, "start"])
            e = int(sub.loc[i, "end"])
            if s <= cur_end + distance_bp:
                cur_end = max(cur_end, e)
            else:
                records.append((chr_value, cur_start, cur_end))
                cur_start, cur_end = s, e
        records.append((chr_value, cur_start, cur_end))

    merged = pd.DataFrame(records, columns=["chr", "start", "end"])
    merged["chr_order"] = merged["chr"].map(lambda x: chr_sort_key(x)[0])
    merged = merged.sort_values(["chr_order", "start"]).reset_index(drop=True)
    merged["locus_code"] = [f"Locus_{i:04d}" for i in range(1, len(merged) + 1)]
    merged["length_bp"] = merged["end"] - merged["start"] + 1
    merged["length_kbp"] = merged["length_bp"] / 1000.0

    merged["IMD"] = np.nan
    for chr_value, sub in merged.groupby("chr", sort=False):
        idxs = list(sub.index)
        for pos, row_idx in enumerate(idxs[:-1]):
            next_idx = idxs[pos + 1]
            merged.at[row_idx, "IMD"] = int(merged.at[next_idx, "start"] - merged.at[row_idx, "end"] - 1)

    merged = merged.drop(columns=["chr_order"])
    return merged[["locus_code", "chr", "start", "end", "length_bp", "length_kbp", "IMD"]]


def merge_and_filter_loci(all_loci_df, distance_bp, min_width_bp):
    """Merge loci, drop loci narrower than min_width_bp, then re-merge survivors."""
    merged = merge_loci(all_loci_df, distance_bp)
    filtered = merged[merged["length_bp"] >= min_width_bp].copy()
    if filtered.empty:
        return filtered

    if len(filtered) != len(merged):
        filtered_input = filtered.rename(columns={"chr": "chromosome"})[["chromosome", "start", "end"]].copy()
        filtered = merge_loci(filtered_input, distance_bp)

    return filtered


def compute_ld_with_plink(snp_ids_1, snp_ids_2, ref_panel_prefix, temp_dir=None):
    """
    Compute cross-set LD between SNP set 1 and SNP set 2 using PLINK.
    Returns max R2 across cross-set pairs, or None if unavailable.
    """
    if not snp_ids_1 or not snp_ids_2:
        return None

    plink_exe = find_plink_executable()
    if not plink_exe:
        print(f"[warn] PLINK executable not found in PATH or common locations", file=sys.stderr)
        return None

    set1 = list(dict.fromkeys([str(x) for x in snp_ids_1 if str(x).strip()]))
    set2 = list(dict.fromkeys([str(x) for x in snp_ids_2 if str(x).strip()]))
    if not set1 or not set2:
        return None

    tmp_root = temp_dir if temp_dir else None

    try:
        with tempfile.TemporaryDirectory(dir=tmp_root, prefix="plink_ld_") as tmpdir:
            set1_file = os.path.join(tmpdir, "set1.snps")
            set2_file = os.path.join(tmpdir, "set2.snps")
            all_file = os.path.join(tmpdir, "all.snps")
            out_prefix = os.path.join(tmpdir, "ld_out")

            with open(set1_file, "w", encoding="utf-8") as f:
                f.write("\n".join(set1) + "\n")
            with open(set2_file, "w", encoding="utf-8") as f:
                f.write("\n".join(set2) + "\n")
            with open(all_file, "w", encoding="utf-8") as f:
                f.write("\n".join(list(dict.fromkeys(set1 + set2))) + "\n")

            cmd = [
                plink_exe,
                "--bfile", ref_panel_prefix,
                "--extract", all_file,
                "--r2",
                "--ld-window-kb", "10000",
                "--ld-window", "99999",
                "--ld-window-r2", "0",
                "--allow-no-sex",
                "--out", out_prefix,
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
            if result.returncode != 0:
                print(f"[warn] PLINK error: {result.stderr[:200]}", file=sys.stderr)
                return None

            ld_file = f"{out_prefix}.ld"
            if not os.path.exists(ld_file):
                return None

            ld_df = pd.read_csv(ld_file, sep=r"\s+", dtype=str)
            if ld_df.empty or "R2" not in ld_df.columns or "SNP_A" not in ld_df.columns or "SNP_B" not in ld_df.columns:
                return None

            set1_lookup = set(set1)
            set2_lookup = set(set2)

            cross = ld_df[
                ((ld_df["SNP_A"].isin(set1_lookup)) & (ld_df["SNP_B"].isin(set2_lookup))) |
                ((ld_df["SNP_A"].isin(set2_lookup)) & (ld_df["SNP_B"].isin(set1_lookup)))
            ].copy()

            if cross.empty:
                return None

            cross["R2"] = pd.to_numeric(cross["R2"], errors="coerce")
            cross = cross.dropna(subset=["R2"])
            if cross.empty:
                return None

            return float(cross["R2"].max())

    except subprocess.TimeoutExpired:
        print("[warn] PLINK LD computation timed out", file=sys.stderr)
    except FileNotFoundError:
        print("[warn] PLINK executable not found in PATH", file=sys.stderr)
    except Exception as e:
        print(f"[warn] LD computation failed: {e}", file=sys.stderr)

    return None


def pick_reference_panel(ref_panels):
    """Pick the first reference panel prefix that has .bed/.bim/.fam files."""
    for prefix in ref_panels.values():
        if all(os.path.exists(f"{prefix}.{ext}") for ext in ("bed", "bim", "fam")):
            return prefix
    return None


def merge_intervals_by_policy(intervals_df, auto_merge_gap_bp, ld_test_min_gap_bp, ld_test_max_gap_bp,
                              ld_merge_threshold, loci_snps, ref_panel_prefix, n_boundary_snps):
    """
    Merge intervals with policy:
      1) overlap/touch => merge
      2) gap <= auto_merge_gap_bp => merge
      3) ld_test_min_gap_bp <= gap <= ld_test_max_gap_bp => run PLINK LD and merge if R2 >= ld_merge_threshold
      4) otherwise => do not merge
    """
    records = []

    for chr_value, sub in intervals_df.groupby("chromosome", sort=False):
        sub = sub.sort_values(["start", "end"]).reset_index(drop=True)
        cur_start = int(sub.loc[0, "start"])
        cur_end = int(sub.loc[0, "end"])

        for i in range(1, len(sub)):
            s = int(sub.loc[i, "start"])
            e = int(sub.loc[i, "end"])
            gap = s - cur_end - 1

            should_merge = False
            reason = ""

            if gap <= 0:
                should_merge = True
                reason = "overlap"
            elif gap <= auto_merge_gap_bp:
                should_merge = True
                reason = f"gap<=auto({auto_merge_gap_bp})"
            elif ld_test_min_gap_bp <= gap <= ld_test_max_gap_bp and loci_snps and ref_panel_prefix:
                left_snps = collect_boundary_snps(
                    loci_snps, chr_value, cur_start, cur_end, n_snps=n_boundary_snps, at_end=True
                )
                right_snps = collect_boundary_snps(
                    loci_snps, chr_value, s, e, n_snps=n_boundary_snps, at_end=False
                )
                ld_val = compute_ld_with_plink(left_snps, right_snps, ref_panel_prefix)

                if ld_val is not None and ld_val >= ld_merge_threshold:
                    should_merge = True
                    reason = f"LD={ld_val:.4f}>=thr({ld_merge_threshold})"
                else:
                    reason = f"LD={(ld_val if ld_val is not None else 'NA')}<thr({ld_merge_threshold})"

                print(
                    f"  LD check chr{chr_value}: [{cur_start}-{cur_end}] vs [{s}-{e}], gap={gap}, {reason}",
                    flush=True,
                )

            if should_merge:
                cur_end = max(cur_end, e)
            else:
                records.append((chr_value, cur_start, cur_end))
                cur_start, cur_end = s, e

        records.append((chr_value, cur_start, cur_end))

    merged = pd.DataFrame(records, columns=["chr", "start", "end"])
    merged["chr_order"] = merged["chr"].map(lambda x: chr_sort_key(x)[0])
    merged = merged.sort_values(["chr_order", "start", "end"]).reset_index(drop=True)
    merged["locus_code"] = [f"Locus_{i:04d}" for i in range(1, len(merged) + 1)]
    merged["length_bp"] = merged["end"] - merged["start"] + 1
    merged["length_kbp"] = merged["length_bp"] / 1000.0
    merged["IMD"] = np.nan
    for chr_value, sub in merged.groupby("chr", sort=False):
        idxs = list(sub.index)
        for pos, row_idx in enumerate(idxs[:-1]):
            next_idx = idxs[pos + 1]
            merged.at[row_idx, "IMD"] = int(merged.at[next_idx, "start"] - merged.at[row_idx, "end"] - 1)
    merged = merged.drop(columns=["chr_order"])
    return merged[["locus_code", "chr", "start", "end", "length_bp", "length_kbp", "IMD"]]


# ---------------------------------------------------------------------------
# Step 3 - assign each source-locus row to its merged locus
# ---------------------------------------------------------------------------

def assign_loci(all_loci_df, merged_df):
    """
    For every row in all_loci_df return the locus_code of the merged locus
    that contains it (by chromosome + overlapping coordinates).
    """
    merged_by_chr = {}
    for _, row in merged_df.iterrows():
        c = norm_chr(row["chr"])
        merged_by_chr.setdefault(c, []).append(
            (int(row["start"]), int(row["end"]), row["locus_code"])
        )

    codes = []
    for _, row in all_loci_df.iterrows():
        c = norm_chr(row["chromosome"])
        s = int(row["start"])
        e = int(row["end"])
        code = None
        for (ms, me, mc) in merged_by_chr.get(c, []):
            if s <= me and e >= ms:
                code = mc
                break
        codes.append(code)
    return codes


# ---------------------------------------------------------------------------
# Step 4 - gene annotation
# ---------------------------------------------------------------------------

def load_gene_loc(path):
    """
    Load a MAGMA-format .gene.loc file (whitespace-separated, no header).
    Expected layout: GeneID  Chr  Start  Stop  [Strand]  [Symbol]
    Returns a DataFrame with columns: gene_id, gene_name, chr, start, end, center.
    """
    if not path or not os.path.exists(path):
        return None

    df = pd.read_csv(path, sep=r"\s+", header=None, dtype=str, comment="#")
    n = df.shape[1]
    if n < 4:
        return None

    col_names = ["gene_id", "chr", "start", "end"]
    if n >= 5:
        col_names.append("strand")
    if n >= 6:
        col_names.append("gene_name")
    for i in range(len(col_names), n):
        col_names.append(f"_col{i}")
    df.columns = col_names

    if "gene_name" not in df.columns:
        df["gene_name"] = df["gene_id"]

    df["chr"] = df["chr"].map(norm_chr)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chr", "start", "end"])
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["center"] = (df["start"] + df["end"]) // 2
    return df[["gene_id", "gene_name", "chr", "start", "end", "center"]]


def find_nearest_gene(genes_df, locus_chr, locus_start, locus_end):
    """Return (gene_name, distance_bp) for the gene nearest to the locus."""
    if genes_df is None:
        return "", np.nan

    sub = genes_df[genes_df["chr"] == norm_chr(locus_chr)]
    if sub.empty:
        return "", np.nan

    locus_center = (locus_start + locus_end) // 2

    # Genes that overlap the locus are at distance 0
    overlap = sub[(sub["start"] <= locus_end) & (sub["end"] >= locus_start)]
    if not overlap.empty:
        best = overlap.iloc[
            (overlap["center"] - locus_center).abs().to_numpy().argmin()
        ]
        return str(best["gene_name"]), 0

    dist = np.minimum(
        np.abs(sub["start"].values - locus_end),
        np.abs(sub["end"].values - locus_start),
    )
    idx = int(np.argmin(dist))
    return str(sub.iloc[idx]["gene_name"]), int(dist[idx])


# ---------------------------------------------------------------------------
# Step 5 - build summary columns (best p-value per dataset x table per locus)
# ---------------------------------------------------------------------------

def build_summary(merged_df, all_loci_df, assigned_codes, gwastables):
    """
        Attach one column per dataset x table to merged_df:
            {dataset}__{table}__top_p  - minimum p-value in the locus

    Returns (summary_df, highlight_col_names).
    """
    tagged = all_loci_df.copy()
    tagged["_merged_locus"] = assigned_codes
    tagged = tagged[tagged["_merged_locus"].notna()].copy()

    summary = merged_df.copy()
    highlight_cols = []

    for dataset_entry in gwastables:
        dataset = str(dataset_entry["name"])
        for table_entry in dataset_entry.get("tables", []):
            table = str(table_entry["table_name"])
            p_col = f"{dataset}__{table}__top_p"

            subset = tagged[
                (tagged["dataset"].astype(str) == dataset) &
                (tagged["table"].astype(str) == table)
            ].copy()

            if subset.empty:
                summary[p_col] = np.nan
            else:
                # Sort ascending by top_p so .first() picks the best SNP per locus
                best = (
                    subset
                    .sort_values("top_p")
                    .groupby("_merged_locus", as_index=False)
                    .first()
                    .set_index("_merged_locus")
                )
                summary[p_col] = summary["locus_code"].map(
                    lambda lc, b=best: b.loc[lc, "top_p"] if lc in b.index else np.nan
                )

            highlight_cols.append(p_col)

    return summary, highlight_cols


# ---------------------------------------------------------------------------
# Step 6 - write Excel with p-value highlight
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
    cfg = read_yaml(args.config)

    if not args.skip_xlsx and not args.output_xlsx:
        raise ValueError("--output-xlsx is required unless --skip-xlsx is used")

    gwastables = cfg.get("gwastables", [])
    gene_loc_path = cfg.get("gene_annotation", {}).get("gencode_gtf", "")
    locus_cfg = cfg.get("loci_identiffication", {})
    if "min_locus_width_kbp" in locus_cfg:
        min_locus_width_bp = int(round(float(locus_cfg["min_locus_width_kbp"]) * 1000))
    else:
        min_locus_width_bp = int(locus_cfg.get("min_locus_width_bp", 1000))

    use_ld = args.use_ld
    ld_threshold = float(locus_cfg.get("ld_threshold", 0.2))
    n_boundary_snps = int(locus_cfg.get("n_boundary_snps", 10))
    ref_panels = cfg.get("ref_panels", {})
    auto_merge_gap_bp = int(float(locus_cfg.get("auto_merge_gap_kb", 100)) * 1000)
    ld_test_min_gap_bp = int(float(locus_cfg.get("ld_test_gap_min_kb", 150)) * 1000)
    ld_test_max_gap_bp = int(float(locus_cfg.get("ld_test_gap_max_kb", args.distance / 1000.0)) * 1000)

    if auto_merge_gap_bp < 0:
        raise ValueError("auto_merge_gap_kb must be >= 0")
    if not (0 <= ld_threshold <= 1):
        raise ValueError("ld_threshold must be between 0 and 1")
    if n_boundary_snps < 1:
        raise ValueError("n_boundary_snps must be >= 1")
    if ld_test_min_gap_bp > ld_test_max_gap_bp:
        raise ValueError("ld_test_gap_min_kb cannot be larger than ld_test_gap_max_kb")

    # 1. Load all per-dataset loci
    print("Loading per-dataset loci ...", flush=True)
    all_loci = load_all_loci(args.inputs)

    if args.chromosome:
        chr_filter = norm_chr(args.chromosome)
        all_loci = all_loci[all_loci["chromosome"] == chr_filter].copy()
        print(f"  Chromosome filter applied: chr{chr_filter}", flush=True)

    if all_loci.empty:
        print("  No loci after filtering; writing empty outputs.", flush=True)

        merged = pd.DataFrame(columns=[
            "locus_code", "chr", "start", "end", "length_bp", "nearest_gene", "gene_distance_bp", "length_kbp", "IMD"
        ])
        summary, highlight_cols = build_summary(merged, all_loci, [], gwastables)

        Path(os.path.dirname(args.output_tsv)).mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.output_tsv, sep="\t", index=False)
        print(f"TSV written:  {args.output_tsv}", flush=True)

        if not args.skip_xlsx:
            write_excel_with_highlight(summary, args.output_xlsx, highlight_cols, threshold=5e-8)
            print(f"XLSX written: {args.output_xlsx}", flush=True)
        return

    print(f"  Total source loci loaded: {len(all_loci)}", flush=True)

    loci_snps = {}
    if use_ld:
        print("Loading SNP data for LD-based merging ...", flush=True)
        loci_snps = load_snps_for_loci()
        print(f"  Loaded SNP data for {len(loci_snps)} datasets", flush=True)

    # 2. Merge loci by configured policy
    print(f"Merging loci with policy: overlap or gap<={auto_merge_gap_bp}bp merge directly; "
          f"{ld_test_min_gap_bp}-{ld_test_max_gap_bp}bp uses LD (if enabled)", flush=True)

    if use_ld:
        ref_panel_prefix = pick_reference_panel(ref_panels)
        if ref_panel_prefix:
            print(f"Using PLINK reference panel: {ref_panel_prefix}", flush=True)
        else:
            print("[warn] No usable reference panel found (.bed/.bim/.fam). LD checks will be skipped.", file=sys.stderr)
    else:
        ref_panel_prefix = None

    merged_initial = merge_intervals_by_policy(
        all_loci[["chromosome", "start", "end"]].copy(),
        auto_merge_gap_bp=auto_merge_gap_bp,
        ld_test_min_gap_bp=ld_test_min_gap_bp,
        ld_test_max_gap_bp=ld_test_max_gap_bp,
        ld_merge_threshold=ld_threshold,
        loci_snps=loci_snps if use_ld else {},
        ref_panel_prefix=ref_panel_prefix,
        n_boundary_snps=n_boundary_snps,
    )
    print(f"  Merged loci before width filter: {len(merged_initial)}", flush=True)

    merged = merged_initial[merged_initial["length_bp"] >= min_locus_width_bp].copy()
    merged.reset_index(drop=True, inplace=True)
    merged["locus_code"] = [f"Locus_{i:04d}" for i in range(1, len(merged) + 1)]

    # Recompute IMD after width filtering.
    merged["IMD"] = np.nan
    for chr_value, sub in merged.groupby("chr", sort=False):
        idxs = list(sub.index)
        for pos, row_idx in enumerate(idxs[:-1]):
            next_idx = idxs[pos + 1]
            merged.at[row_idx, "IMD"] = int(merged.at[next_idx, "start"] - merged.at[row_idx, "end"] - 1)

    print(f"  Merged loci after width filter (>= {min_locus_width_bp} bp): {len(merged)}", flush=True)

    # 3. Assign each source locus to a merged locus
    assigned = assign_loci(all_loci, merged)
    n_unassigned = sum(c is None for c in assigned)
    if n_unassigned:
        print(f"  [warn] {n_unassigned} source loci could not be assigned "
              "to a merged locus", file=sys.stderr)

    # 4. Annotate with nearest gene
    print("Annotating with nearest gene ...", flush=True)
    genes = load_gene_loc(gene_loc_path)
    if genes is None:
        print(f"  [warn] Gene annotation file not found or empty: {gene_loc_path}",
              file=sys.stderr)

    gene_names, gene_dists = [], []
    for _, locus in merged.iterrows():
        name, dist = find_nearest_gene(
            genes, locus["chr"], int(locus["start"]), int(locus["end"])
        )
        gene_names.append(name)
        gene_dists.append(dist)

    insert_at = merged.columns.get_loc("length_bp") + 1
    merged.insert(insert_at,     "nearest_gene",     gene_names)
    merged.insert(insert_at + 1, "gene_distance_bp", gene_dists)

    # 5. Build per-dataset x table summary columns
    print("Building per-dataset x table summary ...", flush=True)
    summary, highlight_cols = build_summary(merged, all_loci, assigned, gwastables)

    # 6. Write outputs
    Path(os.path.dirname(args.output_tsv)).mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"TSV written:  {args.output_tsv}", flush=True)

    if not args.skip_xlsx:
        write_excel_with_highlight(summary, args.output_xlsx, highlight_cols, threshold=5e-8)
        print(f"XLSX written: {args.output_xlsx}", flush=True)


if __name__ == "__main__":
    main()

