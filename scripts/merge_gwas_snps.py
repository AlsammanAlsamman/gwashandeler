#!/usr/bin/env python3
"""Merge all GWAS datasets into a comprehensive SNP-level table.

Inputs are Step 00 extracted GWAS tables and Step 7 loci summaries.
Outputs SNPs that are significant (p < threshold) in any dataset with full
p-value and effect size data across all GWAS sources.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side


P_THRESHOLD = 5e-5
GWS = 5e-8


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge all GWAS SNP data across datasets")
    parser.add_argument("--extracted", nargs="+", required=True, help="Step 00 extracted GWAS TSV files")
    parser.add_argument("--loci-summary", required=True, help="Step 7 loci summary TSV file")
    parser.add_argument("--datasets", required=True, help="Comma-separated dataset names")
    parser.add_argument("--p-threshold", type=float, default=5e-5, help="P-value threshold for filtering")
    parser.add_argument("--output-xlsx", required=True, help="Formatted Excel workbook output")
    parser.add_argument("--output-tsv", required=True, help="TSV export of merged SNPs")
    parser.add_argument("--config", default="configs/analysis.yml", help="Pipeline analysis YAML config")
    return parser.parse_args()


def read_yaml(path: str) -> dict:
    config: dict[str, dict[str, str] | str] = {}
    current_section: str | None = None
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip()
            if not line or line.lstrip().startswith("#"):
                continue
            if not raw_line.startswith((" ", "\t")) and ":" in line:
                key, value = line.split(":", 1)
                key = key.strip()
                value = value.strip().strip('"').strip("'")
                if value:
                    config[key] = value
                    current_section = None
                else:
                    config[key] = {}
                    current_section = key
                continue
            if current_section and ":" in line:
                key, value = line.split(":", 1)
                key = key.strip()
                value = value.strip().strip('"').strip("'")
                section = config.setdefault(current_section, {})
                if isinstance(section, dict):
                    section[key] = value
    return config


def _pick_first(columns: list[str], candidates: list[str]) -> str | None:
    """Find first candidate column name (case-insensitive)."""
    lower_to_original = {str(col).lower(): str(col) for col in columns}
    for candidate in candidates:
        match = lower_to_original.get(candidate.lower())
        if match is not None:
            return match
    return None


def _normalize_chr(value) -> int | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text:
        return None
    text = re.sub(r"^chr", "", text, flags=re.IGNORECASE)
    if text.isdigit():
        return int(text)
    return None


def load_loci_map(loci_summary_path: str) -> pd.DataFrame:
    """Load loci summary file to get chr, start, end, locus_id, gene_name mapping."""
    if not Path(loci_summary_path).exists():
        print(f"[warn] Loci summary file not found: {loci_summary_path}")
        return pd.DataFrame(columns=["chr", "start", "end", "locus_id", "gene_name"])
    
    loci_map = pd.read_csv(loci_summary_path, sep="\t", low_memory=False)
    
    # Normalize column names
    loci_map.columns = [col.lower() for col in loci_map.columns]
    
    # Ensure required columns exist
    required = ["chr", "start", "end"]
    for col in required:
        if col not in loci_map.columns:
            # Try to find similar column
            cols_lower = [c.lower() for c in loci_map.columns]
            if "chromosome" in cols_lower and col == "chr":
                idx = cols_lower.index("chromosome")
                loci_map.rename(columns={loci_map.columns[idx]: "chr"}, inplace=True)
            elif "pos_start" in cols_lower and col == "start":
                idx = cols_lower.index("pos_start")
                loci_map.rename(columns={loci_map.columns[idx]: "start"}, inplace=True)
            elif "pos_end" in cols_lower and col == "end":
                idx = cols_lower.index("pos_end")
                loci_map.rename(columns={loci_map.columns[idx]: "end"}, inplace=True)
    
    # Numeric conversion
    loci_map["chr"] = pd.to_numeric(loci_map.get("chr", []), errors="coerce")
    loci_map["start"] = pd.to_numeric(loci_map.get("start", []), errors="coerce")
    loci_map["end"] = pd.to_numeric(loci_map.get("end", []), errors="coerce")
    
    # Default columns if missing
    if "locus_id" not in loci_map.columns:
        loci_map["locus_id"] = ""
    if "gene_name" not in loci_map.columns:
        loci_map["gene_name"] = ""
    
    return loci_map[["chr", "start", "end", "locus_id", "gene_name"]].dropna(subset=["chr", "start", "end"])


def snp_in_locus(snp_chr: int, snp_pos: int, locus_chr: int, locus_start: int, locus_end: int) -> bool:
    """Check if SNP is within locus bounds."""
    return snp_chr == locus_chr and locus_start <= snp_pos <= locus_end


def load_and_merge_gwas(extracted_paths: list[str], datasets: list[str], p_threshold: float) -> pd.DataFrame:
    """Load all extracted GWAS tables and merge by SNP."""
    
    dataset_snp_data = {}
    
    for path in extracted_paths:
        if not Path(path).exists():
            continue
        
        # Infer dataset from path: results/extracted/{dataset}/{table}.tsv
        path_str = str(path)
        dataset = None
        for d in datasets:
            if d in path_str:
                dataset = d
                break
        
        if not dataset:
            print(f"[warn] Could not infer dataset from path: {path}")
            continue
        
        # Read header to discover columns
        header = pd.read_csv(path, sep="\t", nrows=0)
        columns = list(header.columns)
        
        chr_col = _pick_first(columns, ["chromosome", "chrom", "chr", "CHROM", "CHR"])
        pos_col = _pick_first(columns, ["position", "pos", "bp", "POS"])
        p_col = _pick_first(columns, ["p", "pvalue", "p_value", "pval", "P", "P_VALUE"])
        snp_col = _pick_first(columns, ["snpid", "rsid", "id", "markername", "varid", "SNP", "ID"])
        or_col = _pick_first(columns, ["or", "OR", "odds_ratio"])
        beta_col = _pick_first(columns, ["beta", "BETA", "effect"])
        
        if chr_col is None or pos_col is None or p_col is None:
            print(f"[warn] Could not detect chr/pos/p columns in {path}; skipping")
            continue
        
        usecols = [chr_col, pos_col, p_col]
        if snp_col and snp_col not in usecols:
            usecols.append(snp_col)
        if or_col and or_col not in usecols:
            usecols.append(or_col)
        if beta_col and beta_col not in usecols:
            usecols.append(beta_col)
        
        gwas = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
        gwas["_chr"] = gwas[chr_col].map(_normalize_chr)
        gwas["_pos"] = pd.to_numeric(gwas[pos_col], errors="coerce")
        gwas["_p"] = pd.to_numeric(gwas[p_col], errors="coerce")
        
        if snp_col:
            gwas["_snp"] = gwas[snp_col].astype(str)
        else:
            gwas["_snp"] = gwas["_chr"].astype("Int64").astype(str) + ":" + gwas["_pos"].astype("Int64").astype(str)
        
        if or_col:
            gwas["_or"] = pd.to_numeric(gwas[or_col], errors="coerce")
        elif beta_col:
            gwas["_or"] = np.exp(pd.to_numeric(gwas[beta_col], errors="coerce"))
        else:
            gwas["_or"] = np.nan
        
        gwas = gwas.dropna(subset=["_chr", "_pos", "_p"]).copy()
        
        # Create unique SNP key
        gwas["_key"] = gwas["_chr"].astype(str) + ":" + gwas["_pos"].astype(str)
        
        if dataset not in dataset_snp_data:
            dataset_snp_data[dataset] = {}
        
        for _, row in gwas.iterrows():
            key = row["_key"]
            if key not in dataset_snp_data[dataset]:
                dataset_snp_data[dataset][key] = {
                    "chr": row["_chr"],
                    "pos": row["_pos"],
                    "snp": row["_snp"],
                    "p": row["_p"],
                    "or": row["_or"],
                }
            else:
                # Keep best p-value if multiple entries for same SNP in same dataset
                if row["_p"] < dataset_snp_data[dataset][key]["p"]:
                    dataset_snp_data[dataset][key]["p"] = row["_p"]
                    dataset_snp_data[dataset][key]["snp"] = row["_snp"]
                    dataset_snp_data[dataset][key]["or"] = row["_or"]
    
    # Find all SNPs with p < threshold in any dataset
    sig_snps = set()
    for dataset, snps in dataset_snp_data.items():
        for key, data in snps.items():
            if pd.notna(data["p"]) and data["p"] < p_threshold:
                sig_snps.add(key)
    
    # Build output dataframe
    rows = []
    for key in sorted(sig_snps):
        chr_pos = key.split(":")
        chr_val = int(chr_pos[0])
        pos_val = int(chr_pos[1])
        
        row = {"chr": chr_val, "pos": pos_val}
        
        # Add per-dataset columns
        for dataset in datasets:
            if dataset in dataset_snp_data and key in dataset_snp_data[dataset]:
                data = dataset_snp_data[dataset][key]
                row[f"{dataset}_snp"] = data["snp"]
                row[f"{dataset}_p"] = data["p"]
                row[f"{dataset}_or"] = data["or"]
            else:
                row[f"{dataset}_snp"] = ""
                row[f"{dataset}_p"] = np.nan
                row[f"{dataset}_or"] = np.nan
        
        rows.append(row)
    
    return pd.DataFrame(rows)


def assign_loci_and_genes(snp_df: pd.DataFrame, loci_map: pd.DataFrame) -> pd.DataFrame:
    """Assign locus ID and gene name to each SNP."""
    
    snp_df["locus_id"] = ""
    snp_df["gene_name"] = ""
    
    for idx, snp_row in snp_df.iterrows():
        snp_chr = int(snp_row["chr"])
        snp_pos = int(snp_row["pos"])
        
        # Find locus containing this SNP
        matches = loci_map[
            (loci_map["chr"] == snp_chr) & 
            (loci_map["start"] <= snp_pos) & 
            (loci_map["end"] >= snp_pos)
        ]
        
        if not matches.empty:
            best_match = matches.iloc[0]
            snp_df.at[idx, "locus_id"] = str(best_match.get("locus_id", ""))
            snp_df.at[idx, "gene_name"] = str(best_match.get("gene_name", ""))
    
    # Reorder columns: chr, pos, snp, locus_id, gene_name, then dataset columns
    base_cols = ["chr", "pos", "locus_id", "gene_name"]
    dataset_cols = [col for col in snp_df.columns if col not in base_cols]
    
    return snp_df[base_cols + dataset_cols]


def write_workbook(df: pd.DataFrame, output_xlsx: str, datasets: list[str]) -> None:
    """Write formatted Excel workbook."""
    wb = Workbook()
    ws = wb.active
    ws.title = "Significant SNPs"
    ws.sheet_view.showGridLines = False
    ws.freeze_panes = "A2"
    
    # Headers
    headers = ["Chr", "Pos", "Locus ID", "Gene"]
    col_idx = {}
    current_col = 1
    
    for h in headers:
        cell = ws.cell(row=1, column=current_col, value=h)
        cell.fill = PatternFill("solid", start_color="1F3864", end_color="1F3864")
        cell.font = Font(bold=True, color="FFFFFF", name="Calibri", size=10)
        cell.alignment = Alignment(horizontal="center", vertical="center")
        col_idx[h.lower().replace(" ", "_")] = current_col
        current_col += 1
    
    # Add per-dataset headers
    for dataset in datasets:
        for suffix in ["_snp", "_p", "_or"]:
            col_name = f"{dataset}{suffix}"
            header_text = {
                "_snp": f"{dataset} SNP",
                "_p": f"{dataset} p-value",
                "_or": f"{dataset} OR"
            }[suffix]
            cell = ws.cell(row=1, column=current_col, value=header_text)
            cell.fill = PatternFill("solid", start_color="2C3E50", end_color="2C3E50")
            cell.font = Font(bold=True, color="FFFFFF", name="Calibri", size=8.5)
            cell.alignment = Alignment(horizontal="center", vertical="center")
            col_idx[col_name] = current_col
            current_col += 1
    
    ws.row_dimensions[1].height = 20
    
    # Data rows
    for row_index, (_, row) in enumerate(df.iterrows(), start=2):
        base_fill = ["F8F9FA", "FFFFFF"][row_index % 2]
        
        # Base columns
        ws.cell(row=row_index, column=col_idx["chr"], value=int(row["chr"]))
        ws.cell(row=row_index, column=col_idx["pos"], value=int(row["pos"]))
        ws.cell(row=row_index, column=col_idx["locus_id"], value=row.get("locus_id", ""))
        ws.cell(row=row_index, column=col_idx["gene_name"], value=row.get("gene_name", ""))
        
        # Dataset columns
        for dataset in datasets:
            snp_col = col_idx.get(f"{dataset}_snp")
            p_col = col_idx.get(f"{dataset}_p")
            or_col = col_idx.get(f"{dataset}_or")
            
            snp_val = row.get(f"{dataset}_snp", "")
            p_val = row.get(f"{dataset}_p", np.nan)
            or_val = row.get(f"{dataset}_or", np.nan)
            
            if snp_col:
                ws.cell(row=row_index, column=snp_col, value=snp_val if pd.notna(snp_val) and snp_val != "" else None)
            
            if p_col:
                if pd.notna(p_val):
                    p_str = f"{float(p_val):.2e}"
                    cell = ws.cell(row=row_index, column=p_col, value=p_str)
                    if p_val < GWS:
                        cell.fill = PatternFill("solid", start_color="C6EFCE", end_color="C6EFCE")
                        cell.font = Font(bold=True, name="Calibri", size=8.5, color="276221")
                    elif p_val < 1e-5:
                        cell.fill = PatternFill("solid", start_color="FFEB9C", end_color="FFEB9C")
                        cell.font = Font(name="Calibri", size=8.5, color="9C6500")
                    else:
                        cell.fill = PatternFill("solid", start_color=base_fill, end_color=base_fill)
                        cell.font = Font(name="Calibri", size=8.5)
                else:
                    ws.cell(row=row_index, column=p_col, value=None)
            
            if or_col:
                if pd.notna(or_val) and or_val > 0:
                    ws.cell(row=row_index, column=or_col, value=f"{float(or_val):.4f}")
                else:
                    ws.cell(row=row_index, column=or_col, value=None)
    
    # Adjust column widths
    ws.column_dimensions['A'].width = 5
    ws.column_dimensions['B'].width = 12
    ws.column_dimensions['C'].width = 12
    ws.column_dimensions['D'].width = 18
    for col_num in range(5, current_col):
        ws.column_dimensions[chr(64 + col_num)].width = 14
    
    Path(output_xlsx).parent.mkdir(parents=True, exist_ok=True)
    wb.save(output_xlsx)


def main() -> None:
    args = parse_args()
    
    datasets = [d.strip() for d in args.datasets.split(",")]
    
    print(f"Loading loci summary from: {args.loci_summary}")
    loci_map = load_loci_map(args.loci_summary)
    print(f"  Loaded {len(loci_map)} locus regions")
    
    print(f"Loading and merging {len(args.extracted)} GWAS extracted files for {len(datasets)} datasets...")
    snp_df = load_and_merge_gwas(args.extracted, datasets, args.p_threshold)
    print(f"  Found {len(snp_df)} SNPs with p < {args.p_threshold} in any dataset")
    
    print(f"Assigning loci and genes...")
    snp_df = assign_loci_and_genes(snp_df, loci_map)
    
    print(f"Writing outputs...")
    output_tsv = Path(args.output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    snp_df.to_csv(output_tsv, sep="\t", index=False)
    print(f"  Saved TSV: {args.output_tsv}")
    
    write_workbook(snp_df, args.output_xlsx, datasets)
    print(f"  Saved XLSX: {args.output_xlsx}")
    
    print(f"\nSummary:")
    print(f"  Total significant SNPs: {len(snp_df)}")
    for dataset in datasets:
        n_sig = (pd.to_numeric(snp_df[f"{dataset}_p"], errors="coerce") < args.p_threshold).sum()
        print(f"  - {dataset}: {n_sig} SNPs with p < {args.p_threshold}")


if __name__ == "__main__":
    main()
