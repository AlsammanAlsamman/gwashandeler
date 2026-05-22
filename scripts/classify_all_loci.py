#!/usr/bin/env python3
"""Classify merged loci across pairwise Tractor datasets and the Hispanic GWAS.

Inputs are the dataset-level merged loci TSVs produced by Step 6. The script
merges nearby loci across datasets, classifies each merged locus by its driving
ancestry, and writes both TSV and formatted XLSX outputs.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

GWS = 5e-8
MERGE_GAP = 250_000
ANCESTRIES = ["AMR", "AFR", "EUR", "EAS"]
DEFAULT_SOURCES = ["Hisp", "AFR_AMR", "AFR_EAS", "AFR_EUR", "AMR_EAS", "AMR_EUR", "EUR_EAS"]

# Global list of allowed sources - will be populated from config in main()
ALLOWED_SOURCES = DEFAULT_SOURCES.copy()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Classify merged loci across Step 6 dataset outputs")
    parser.add_argument("--inputs", nargs="+", required=True, help="Dataset-level merged loci TSV files from Step 6")
    parser.add_argument("--extracted", nargs="*", default=[], help="Step 00 extracted GWAS TSV files for primary datasets")
    parser.add_argument("--output-xlsx", required=True, help="Formatted Excel workbook output")
    parser.add_argument("--output-tsv", required=True, help="TSV export of the main classified table")
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


def norm_chr(value) -> str:
    return re.sub(r"^chr", "", str(value).strip(), flags=re.IGNORECASE)


def load_gene_loc(path: str | None) -> pd.DataFrame | None:
    if not path or not Path(path).exists():
        return None

    df = pd.read_csv(path, sep=r"\s+", header=None, dtype=str, comment="#")
    if df.shape[1] < 4:
        return None

    col_names = ["gene_id", "chr", "start", "end"]
    if df.shape[1] >= 5:
        col_names.append("strand")
    if df.shape[1] >= 6:
        col_names.append("gene_name")
    for index in range(len(col_names), df.shape[1]):
        col_names.append(f"_col{index}")
    df.columns = col_names

    if "gene_name" not in df.columns:
        df["gene_name"] = df["gene_id"]

    df["chr"] = df["chr"].map(norm_chr)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["chr", "start", "end"]).copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["center"] = (df["start"] + df["end"]) // 2
    return df[["gene_id", "gene_name", "chr", "start", "end", "center"]]


def find_nearest_gene(genes_df: pd.DataFrame | None, locus_chr, locus_start: int, locus_end: int) -> tuple[str, float]:
    if genes_df is None:
        return "", np.nan

    sub = genes_df[genes_df["chr"] == norm_chr(locus_chr)]
    if sub.empty:
        return "", np.nan

    locus_center = (locus_start + locus_end) // 2
    overlap = sub[(sub["start"] <= locus_end) & (sub["end"] >= locus_start)]
    if not overlap.empty:
        best = overlap.iloc[(overlap["center"] - locus_center).abs().to_numpy().argmin()]
        return str(best["gene_name"]), 0

    dist = np.minimum(np.abs(sub["start"].values - locus_end), np.abs(sub["end"].values - locus_start))
    idx = int(np.argmin(dist))
    return str(sub.iloc[idx]["gene_name"]), int(dist[idx])


def infer_source_from_path(path: str) -> str:
    return Path(path).stem.replace("_loci", "")


def infer_source_from_extracted_path(path: str) -> str:
    return Path(path).parent.name


def _pick_first(columns: list[str], candidates: list[str]) -> str | None:
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


def _discover_gwas_columns(path: str) -> tuple[str | None, str | None, str | None, str | None]:
    header = pd.read_csv(path, sep="\t", nrows=0)
    columns = list(header.columns)
    chr_col = _pick_first(columns, ["chromosome", "chrom", "chr", "CHROM", "CHR"])
    pos_col = _pick_first(columns, ["position", "pos", "bp", "POS"])
    p_col = _pick_first(columns, ["p", "pvalue", "p_value", "pval", "P"])
    snp_col = _pick_first(columns, ["snpid", "rsid", "id", "markername", "varid", "SNP"])
    return chr_col, pos_col, snp_col, p_col


def append_primary_region_details(summary_df: pd.DataFrame, extracted_paths: list[str], region_sources: list[str]) -> pd.DataFrame:
    if summary_df.empty or not extracted_paths or not region_sources:
        return summary_df

    out = summary_df.copy()
    source_set = set(region_sources)

    for source in region_sources:
        out[f"{source}_region_best_p"] = np.nan
        out[f"{source}_region_best_snp"] = ""
        out[f"{source}_region_best_table"] = ""

    for path in extracted_paths:
        source = infer_source_from_extracted_path(path)
        if source not in source_set:
            continue
        if not Path(path).exists():
            print(f"[warn] Extracted GWAS file missing for source '{source}': {path}")
            continue

        chr_col, pos_col, snp_col, p_col = _discover_gwas_columns(path)
        if chr_col is None or pos_col is None or p_col is None:
            print(f"[warn] Could not detect chr/pos/p columns in {path}; skipping")
            continue

        usecols = [chr_col, pos_col, p_col]
        if snp_col is not None and snp_col not in usecols:
            usecols.append(snp_col)

        gwas = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
        gwas["_chr"] = gwas[chr_col].map(_normalize_chr)
        gwas["_pos"] = pd.to_numeric(gwas[pos_col], errors="coerce")
        gwas["_p"] = pd.to_numeric(gwas[p_col], errors="coerce")
        if snp_col is None:
            gwas["_snp"] = gwas["_chr"].astype("Int64").astype(str) + ":" + gwas["_pos"].astype("Int64").astype(str)
        else:
            gwas["_snp"] = gwas[snp_col].astype(str)

        gwas = gwas.dropna(subset=["_chr", "_pos", "_p"]).copy()
        if gwas.empty:
            continue

        table_name = Path(path).stem
        for idx, row in out.iterrows():
            rchr = _normalize_chr(row.get("chr"))
            rstart = pd.to_numeric(row.get("start"), errors="coerce")
            rend = pd.to_numeric(row.get("end"), errors="coerce")
            if rchr is None or pd.isna(rstart) or pd.isna(rend):
                continue

            subset = gwas[(gwas["_chr"] == rchr) & (gwas["_pos"] >= rstart) & (gwas["_pos"] <= rend)]
            if subset.empty:
                continue

            best = subset.loc[subset["_p"].idxmin()]
            best_p = float(best["_p"])
            best_snp = str(best["_snp"])
            current_best = pd.to_numeric(out.at[idx, f"{source}_region_best_p"], errors="coerce")
            if pd.isna(current_best) or best_p < current_best:
                out.at[idx, f"{source}_region_best_p"] = best_p
                out.at[idx, f"{source}_region_best_snp"] = best_snp
                out.at[idx, f"{source}_region_best_table"] = table_name

    return out


def get_ancestry(table: str) -> str:
    if table.startswith("local_ancestry_"):
        return "local_ancestry"
    return table.split("_")[0]


def get_sigtype(table: str) -> str:
    return "local_ancestry" if table.startswith("local_ancestry_") else "dosage"


def load_pool(inputs: list[str]) -> pd.DataFrame:
    rows = []
    for path in inputs:
        src = infer_source_from_path(path)
        if src not in ALLOWED_SOURCES:
            continue
        df = pd.read_csv(path, sep="\t", low_memory=False)
        if df.empty:
            continue
        id_col = "dataset_locus" if "dataset_locus" in df.columns else "genomic_locus"
        if id_col not in df.columns or "table" not in df.columns:
            continue
        df = df.rename(columns={id_col: "orig_id"}).copy()
        df["source"] = src
        df["ancestry"] = df["table"].apply(get_ancestry)
        df["sig_type"] = df["table"].apply(get_sigtype)
        rows.append(df[["orig_id", "source", "table", "chromosome", "start", "end", "top_snp", "top_p", "n_sig_snps", "ancestry", "sig_type"]])

    if not rows:
        return pd.DataFrame(columns=["orig_id", "source", "table", "chromosome", "start", "end", "top_snp", "top_p", "n_sig_snps", "ancestry", "sig_type"])

    pool = pd.concat(rows, ignore_index=True)
    pool["chromosome"] = pd.to_numeric(pool["chromosome"], errors="coerce")
    pool["start"] = pd.to_numeric(pool["start"], errors="coerce")
    pool["end"] = pd.to_numeric(pool["end"], errors="coerce")
    pool["top_p"] = pd.to_numeric(pool["top_p"], errors="coerce")
    pool["n_sig_snps"] = pd.to_numeric(pool["n_sig_snps"], errors="coerce").fillna(0).astype(int)
    pool = pool.dropna(subset=["chromosome", "start", "end", "top_p"]).copy()
    pool["chromosome"] = pool["chromosome"].astype(int)
    pool["start"] = pool["start"].astype(int)
    pool["end"] = pool["end"].astype(int)
    return pool

CLASS_COL = {
    "AMR": "FFB347",
    "EUR": "74B2E0",
    "AFR": "81C995",
    "EAS": "98D8C8",
    "AMR + EUR": "C3A6E8",
    "AMR + AFR": "FFAA80",
    "AMR + EAS": "FFE0A0",
    "EUR + EAS": "A0C8F0",
    "Multi-ancestry (HLA)": "F28B82",
    "Multi-ancestry (3)": "E07070",
    "Multi-ancestry (4)": "C04040",
    "Admixture mapping signal": "FFD580",
    "Hisp GWAS (shared / unclear)": "D3D3D3",
    "Not classified": "EEEEEE",
}

# Placeholder sections - will be built dynamically in main()
SECTIONS = []
SEC_HDR_COL = {}
COL_LABELS = {}
COL_WIDTHS = {}
P_COLS = set()

def build_output_schema(primary_sources: list[str]) -> None:
    """Build output schema (SECTIONS, COL_LABELS, COL_WIDTHS, SEC_HDR_COL, P_COLS) dynamically from primary sources."""
    global SECTIONS, SEC_HDR_COL, COL_LABELS, COL_WIDTHS, P_COLS
    
    # Build primary GWAS section
    primary_cols = []
    primary_labels = {}
    primary_widths = {}
    for source in primary_sources:
        p_col = f"{source}_p"
        snp_col = f"{source}_snp"
        region_best_p_col = f"{source}_region_best_p"
        region_best_snp_col = f"{source}_region_best_snp"
        region_best_table_col = f"{source}_region_best_table"
        primary_cols.extend([p_col, snp_col, region_best_p_col, region_best_snp_col, region_best_table_col])
        primary_labels[p_col] = f"{source} p-value"
        primary_labels[snp_col] = f"{source} lead SNP"
        primary_labels[region_best_p_col] = f"{source} best region p"
        primary_labels[region_best_snp_col] = f"{source} best region SNP"
        primary_labels[region_best_table_col] = f"{source} best-hit table"
        primary_widths[p_col] = 13
        primary_widths[snp_col] = 13
        primary_widths[region_best_p_col] = 13
        primary_widths[region_best_snp_col] = 20
        primary_widths[region_best_table_col] = 16
    
    # Build sections
    SECTIONS = [
        ("Locus", ["merged_id", "chr", "start", "end", "size_kb", "gene_name"]),
        ("Origin", ["n_sources", "sources", "orig_loci"]),
        ("Primary GWAS", primary_cols),
        ("AMR Signal", ["AMR_dosage_p", "AMR_dosage_snp", "AMR_dosage_pair", "AMR_n_sig"]),
        ("AFR Signal", ["AFR_dosage_p", "AFR_dosage_snp", "AFR_dosage_pair", "AFR_n_sig"]),
        ("EUR Signal", ["EUR_dosage_p", "EUR_dosage_snp", "EUR_dosage_pair", "EUR_n_sig"]),
        ("EAS Signal", ["EAS_dosage_p", "EAS_dosage_snp", "EAS_dosage_pair", "EAS_n_sig"]),
        ("Local Ancestry", ["LA_p", "LA_snp", "LA_pair"]),
        ("Classification", ["Class", "Reason"]),
    ]
    
    # Build headers
    SEC_HDR_COL = {
        "Locus": "1F3864",
        "Origin": "17375E",
        "Primary GWAS": "1C4587",
        "AMR Signal": "7F4B00",
        "AFR Signal": "7F1D1D",
        "EUR Signal": "1A237E",
        "EAS Signal": "1B5E20",
        "Local Ancestry": "4A148C",
        "Classification": "212121",
    }
    
    # Build labels
    COL_LABELS = {
        "merged_id": "Locus#",
        "chr": "Chr",
        "start": "Start (bp)",
        "end": "End (bp)",
        "size_kb": "Size (kb)",
        "gene_name": "Gene",
        "n_sources": "# Source files",
        "sources": "Source files",
        "orig_loci": "Original loci (merged)",
        "AMR_dosage_p": "AMR p (dosage)",
        "AMR_dosage_snp": "AMR lead SNP",
        "AMR_dosage_pair": "AMR source pair",
        "AMR_n_sig": "AMR #sig SNPs",
        "AFR_dosage_p": "AFR p (dosage)",
        "AFR_dosage_snp": "AFR lead SNP",
        "AFR_dosage_pair": "AFR source pair",
        "AFR_n_sig": "AFR #sig SNPs",
        "EUR_dosage_p": "EUR p (dosage)",
        "EUR_dosage_snp": "EUR lead SNP",
        "EUR_dosage_pair": "EUR source pair",
        "EUR_n_sig": "EUR #sig SNPs",
        "EAS_dosage_p": "EAS p (dosage)",
        "EAS_dosage_snp": "EAS lead SNP",
        "EAS_dosage_pair": "EAS source pair",
        "EAS_n_sig": "EAS #sig SNPs",
        "LA_p": "Local Ancestry p",
        "LA_snp": "LA lead SNP",
        "LA_pair": "LA source pair",
        "Class": "Class",
        "Reason": "Classification Reason",
        **primary_labels
    }
    
    # Build widths
    COL_WIDTHS = {
        "merged_id": 7, "chr": 5, "start": 13, "end": 13, "size_kb": 9, "gene_name": 18,
        "n_sources": 9, "sources": 16, "orig_loci": 60,
        "AMR_dosage_p": 13, "AMR_dosage_snp": 13, "AMR_dosage_pair": 13, "AMR_n_sig": 9,
        "AFR_dosage_p": 13, "AFR_dosage_snp": 13, "AFR_dosage_pair": 13, "AFR_n_sig": 9,
        "EUR_dosage_p": 13, "EUR_dosage_snp": 13, "EUR_dosage_pair": 13, "EUR_n_sig": 9,
        "EAS_dosage_p": 13, "EAS_dosage_snp": 13, "EAS_dosage_pair": 13, "EAS_n_sig": 9,
        "LA_p": 13, "LA_snp": 13, "LA_pair": 13,
        "Class": 22, "Reason": 75,
        **primary_widths
    }
    
    # P-columns for highlighting
    P_COLS = set()
    for source in primary_sources:
        P_COLS.add(f"{source}_p")
        P_COLS.add(f"{source}_region_best_p")
    P_COLS.update({"AMR_dosage_p", "AFR_dosage_p", "EUR_dosage_p", "EAS_dosage_p", "LA_p"})
ROW_FILLS = ["F8F9FA", "FFFFFF"]


def merge_loci(df: pd.DataFrame, gap: int = MERGE_GAP) -> tuple[pd.DataFrame, int]:
    df = df.copy().sort_values(["chromosome", "start"]).reset_index(drop=True)
    df["merged_id"] = -1
    cluster = 0
    cur_chr = -1
    cur_end = -1
    for i, row in df.iterrows():
        chromosome, start, end = row["chromosome"], row["start"], row["end"]
        if chromosome != cur_chr or start > cur_end + gap:
            cluster += 1
            cur_chr = chromosome
            cur_end = end
        else:
            cur_end = max(cur_end, end)
        df.at[i, "merged_id"] = cluster
    return df, cluster


def fmt_p(p: float) -> str:
    return f"{p:.2e}" if pd.notna(p) else "NS"


def classify(row: pd.Series) -> tuple[str, str]:
    ch = row["chr"]
    start = row["start"]
    amr_p = row["AMR_dosage_p"]
    afr_p = row["AFR_dosage_p"]
    eur_p = row["EUR_dosage_p"]
    eas_p = row["EAS_dosage_p"]
    la_p = row["LA_p"]
    hp = row["Hisp_p"]
    amr = pd.notna(amr_p) and amr_p < GWS
    afr = pd.notna(afr_p) and afr_p < GWS
    eur = pd.notna(eur_p) and eur_p < GWS
    eas = pd.notna(eas_p) and eas_p < GWS
    la = pd.notna(la_p) and la_p < GWS
    inh = pd.notna(hp) and hp < GWS
    gws_ancs = [a for a, g in zip(ANCESTRIES, [amr, afr, eur, eas]) if g]

    if ch == 6 and 24_000_000 <= start <= 47_000_000:
        involved = ", ".join([f"{a} p={fmt_p(row[f'{a}_dosage_p'])}" for a in ANCESTRIES if pd.notna(row[f"{a}_dosage_p"])])
        return ("Multi-ancestry (HLA)", f"HLA region (chr6:24–47 Mb). Complex multi-ancestry autoimmune locus with extensive LD. Involved ancestries: {involved}. Local ancestry p={fmt_p(la_p)}. Hisp GWAS p={fmt_p(hp)}.")
    if len(gws_ancs) == 4:
        return ("Multi-ancestry (4)", f"All 4 ancestries GWS: AMR p={fmt_p(amr_p)}, AFR p={fmt_p(afr_p)}, EUR p={fmt_p(eur_p)}, EAS p={fmt_p(eas_p)}. Likely a deep, ancient causal variant shared across all populations.")
    if len(gws_ancs) == 3:
        return ("Multi-ancestry (3)", f"3 ancestries GWS: {', '.join(gws_ancs)} (AMR p={fmt_p(amr_p)}, AFR p={fmt_p(afr_p)}, EUR p={fmt_p(eur_p)}, EAS p={fmt_p(eas_p)}). Broad multi-ethnic signal. Hisp GWAS p={fmt_p(hp)}.")
    if amr and eur and not afr and not eas:
        dominant = "AMR-dominant" if amr_p < eur_p else "EUR-dominant"
        return ("AMR + EUR", f"Both AMR (p={fmt_p(amr_p)}) and EUR (p={fmt_p(eur_p)}) dosage GWS. {dominant} by p-value. Likely shared causal variant with ancestry-specific effect sizes. AFR p={fmt_p(afr_p)}, EAS p={fmt_p(eas_p)}. {'Local ancestry GWS (p=' + fmt_p(la_p) + ') supports haplotype block effect. ' if la else ''}Hisp GWAS p={fmt_p(hp)}.")
    if amr and afr and not eur and not eas:
        return ("AMR + AFR", f"Both AMR (p={fmt_p(amr_p)}) and AFR (p={fmt_p(afr_p)}) dosage GWS. EUR p={fmt_p(eur_p)}, EAS p={fmt_p(eas_p)}. Signal present on both Native American and African haplotypes. Hisp GWAS p={fmt_p(hp)}.")
    if eur and eas and not amr and not afr:
        return ("EUR + EAS", f"Both EUR (p={fmt_p(eur_p)}) and EAS (p={fmt_p(eas_p)}) dosage GWS. AMR p={fmt_p(amr_p)}, AFR p={fmt_p(afr_p)}. Signal present on Eurasian lineage haplotypes. Hisp GWAS p={fmt_p(hp)}.")
    if amr and eas and not eur and not afr:
        return ("AMR + EAS", f"Both AMR (p={fmt_p(amr_p)}) and EAS (p={fmt_p(eas_p)}) dosage GWS. Consistent with Native American / East Asian shared ancestry. EUR p={fmt_p(eur_p)}, AFR p={fmt_p(afr_p)}. Hisp GWAS p={fmt_p(hp)}.")
    if amr and not eur and not afr and not eas:
        return ("AMR", f"AMR dosage GWS (p={fmt_p(amr_p)}); EUR p={fmt_p(eur_p)}, AFR p={fmt_p(afr_p)}, EAS p={fmt_p(eas_p)}. Signal driven by Native American haplotypes. {'Local ancestry GWS (p=' + fmt_p(la_p) + '). ' if la else ''}Hisp GWAS p={fmt_p(hp)}.")
    if eur and not amr and not afr and not eas:
        return ("EUR", f"EUR dosage GWS (p={fmt_p(eur_p)}); AMR p={fmt_p(amr_p)}, AFR p={fmt_p(afr_p)}, EAS p={fmt_p(eas_p)}. Signal driven by European haplotypes. {'Local ancestry GWS (p=' + fmt_p(la_p) + '). ' if la else ''}Hisp GWAS p={fmt_p(hp)}.")
    if afr and not amr and not eur and not eas:
        return ("AFR", f"AFR dosage GWS (p={fmt_p(afr_p)}); AMR p={fmt_p(amr_p)}, EUR p={fmt_p(eur_p)}, EAS p={fmt_p(eas_p)}. Signal on African haplotypes. {'Local ancestry GWS (p=' + fmt_p(la_p) + '). ' if la else ''}Hisp GWAS p={fmt_p(hp)}.")
    if eas and not amr and not eur and not afr:
        return ("EAS", f"EAS dosage GWS (p={fmt_p(eas_p)}); AMR p={fmt_p(amr_p)}, EUR p={fmt_p(eur_p)}, AFR p={fmt_p(afr_p)}. Signal on East Asian haplotypes. Hisp GWAS p={fmt_p(hp)}.")
    if la and not any([amr, eur, afr, eas]):
        return ("Admixture mapping signal", f"No ancestry-specific dosage reaches GWS (AMR p={fmt_p(amr_p)}, EUR p={fmt_p(eur_p)}, AFR p={fmt_p(afr_p)}, EAS p={fmt_p(eas_p)}), but local ancestry is GWS (p={fmt_p(la_p)}). Hisp GWAS p={fmt_p(hp)}.")
    if inh and not any([amr, eur, afr, eas, la]):
        return ("Hisp GWAS (shared / unclear)", f"Significant in Hispanic GWAS (p={fmt_p(hp)}) but no single ancestry reaches GWS in pairwise Tractor analyses. AMR p={fmt_p(amr_p)}, EUR p={fmt_p(eur_p)}, AFR p={fmt_p(afr_p)}, EAS p={fmt_p(eas_p)}.")
    if not inh and not any([amr, eur, afr, eas, la]):
        return ("Not classified", f"No GWS signal in any test. AMR p={fmt_p(amr_p)}, EUR p={fmt_p(eur_p)}, AFR p={fmt_p(afr_p)}, EAS p={fmt_p(eas_p)}, LA p={fmt_p(la_p)}, Hisp p={fmt_p(hp)}.")
    return (" + ".join(gws_ancs) if gws_ancs else "Suggestive only", f"GWS ancestries: {gws_ancs}. AMR p={fmt_p(amr_p)}, EUR p={fmt_p(eur_p)}, AFR p={fmt_p(afr_p)}, EAS p={fmt_p(eas_p)}, LA p={fmt_p(la_p)}, Hisp p={fmt_p(hp)}.")


def build_summary(pool: pd.DataFrame, n_clusters: int, genes_df: pd.DataFrame | None, primary_sources: list[str] = None) -> pd.DataFrame:
    if primary_sources is None:
        primary_sources = ["Hisp"]
    
    records = []
    for mid in range(1, n_clusters + 1):
        grp = pool[pool["merged_id"] == mid]
        chromosome = grp["chromosome"].iloc[0]
        start = grp["start"].min()
        end = grp["end"].max()
        best_overall = grp.loc[grp["top_p"].idxmin()]
        gene_name, _ = find_nearest_gene(
            genes_df,
            best_overall["chromosome"],
            int(best_overall["start"]),
            int(best_overall["end"]),
        )
        orig_parts = []
        for source in sorted(grp["source"].unique()):
            ids = sorted(grp.loc[grp["source"] == source, "orig_id"].astype(str).unique())
            orig_parts.append(f"{source}: {','.join(ids)}")
        orig_loci = " | ".join(orig_parts)

        anc_data = {}
        for anc in ANCESTRIES:
            dosage = grp[(grp["ancestry"] == anc) & (grp["sig_type"] == "dosage")]
            if dosage.empty:
                anc_data[f"{anc}_dosage_p"] = np.nan
                anc_data[f"{anc}_dosage_snp"] = ""
                anc_data[f"{anc}_dosage_pair"] = ""
                anc_data[f"{anc}_n_sig"] = 0
            else:
                best = dosage.loc[dosage["top_p"].idxmin()]
                anc_data[f"{anc}_dosage_p"] = float(best["top_p"])
                anc_data[f"{anc}_dosage_snp"] = str(best["top_snp"])
                anc_data[f"{anc}_dosage_pair"] = str(best["source"])
                anc_data[f"{anc}_n_sig"] = int(dosage["n_sig_snps"].max())

        local_ancestry = grp[grp["sig_type"] == "local_ancestry"]
        
        # Extract data for each primary source
        primary_data = {}
        for source in primary_sources:
            source_rows = grp[grp["source"] == source]
            if source_rows.empty:
                primary_data[f"{source}_p"] = np.nan
                primary_data[f"{source}_snp"] = ""
            else:
                best_source = source_rows.loc[source_rows["top_p"].idxmin()]
                primary_data[f"{source}_p"] = float(best_source["top_p"])
                primary_data[f"{source}_snp"] = str(best_source["top_snp"])
        
        records.append({
            "merged_id": mid,
            "chr": chromosome,
            "start": start,
            "end": end,
            "size_bp": end - start,
            "size_kb": round((end - start) / 1000, 1),
            "gene_name": gene_name,
            "n_input_rows": len(grp),
            "n_sources": grp["source"].nunique(),
            "sources": "|".join(sorted(grp["source"].unique())),
            "orig_loci": orig_loci,
            **primary_data,
            **anc_data,
            "LA_p": float(local_ancestry["top_p"].min()) if not local_ancestry.empty else np.nan,
            "LA_snp": str(local_ancestry.loc[local_ancestry["top_p"].idxmin(), "top_snp"]) if not local_ancestry.empty else "",
            "LA_pair": str(local_ancestry.loc[local_ancestry["top_p"].idxmin(), "source"]) if not local_ancestry.empty else "",
        })

    summary = pd.DataFrame(records)
    if not summary.empty:
        classified = summary.apply(classify, axis=1)
        summary["Class"] = [item[0] for item in classified]
        summary["Reason"] = [item[1] for item in classified]
    return summary


def thin() -> Border:
    side = Side(style="thin", color="D0D0D0")
    return Border(left=side, right=side, top=side, bottom=side)


def header_fill(hex_color: str) -> PatternFill:
    return PatternFill("solid", start_color=hex_color, end_color=hex_color)


def pfmt(value) -> str:
    if pd.isna(value) or value == "":
        return ""
    try:
        return f"{float(value):.2e}"
    except Exception:
        return str(value)


def write_workbook(df: pd.DataFrame, output_xlsx: str) -> None:
    wb = Workbook()
    ws = wb.active
    ws.title = "All Loci Classified"
    ws.sheet_view.showGridLines = False
    ws.freeze_panes = "A4"

    all_cols = [col for _, cols in SECTIONS for col in cols]
    col_idx = {col: i + 1 for i, col in enumerate(all_cols)}

    current_column = 1
    for section_name, section_cols in SECTIONS:
        start_col = current_column
        end_col = current_column + len(section_cols) - 1
        if start_col != end_col:
            ws.merge_cells(start_row=1, start_column=start_col, end_row=1, end_column=end_col)
        cell = ws.cell(row=1, column=start_col, value=section_name)
        cell.fill = header_fill(SEC_HDR_COL.get(section_name, "1F3864"))
        cell.font = Font(bold=True, color="FFFFFF", name="Calibri", size=10)
        cell.alignment = Alignment(horizontal="center", vertical="center")
        cell.border = thin()
        current_column += len(section_cols)
    ws.row_dimensions[1].height = 20

    for column_name, column_index in col_idx.items():
        cell = ws.cell(row=2, column=column_index, value=COL_LABELS.get(column_name, column_name))
        cell.fill = header_fill("2C3E50")
        cell.font = Font(bold=True, color="FFFFFF", name="Calibri", size=8.5)
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = thin()
    ws.row_dimensions[2].height = 30

    ws.cell(row=3, column=col_idx.get("AMR_dosage_p", 1), value="GWS: p < 5×10⁻⁸  |  Yellow: p < 1×10⁻⁵  |  Gap for merging: 250 kb").font = Font(italic=True, color="888888", name="Calibri", size=7.5)
    ws.row_dimensions[3].height = 12

    for row_index, (_, row) in enumerate(df.iterrows(), start=4):
        class_name = row.get("Class", "")
        class_fill = CLASS_COL.get(class_name, "FFFFFF")
        base_fill = ROW_FILLS[row_index % 2]
        for column_name, column_index in col_idx.items():
            raw = row.get(column_name, None)
            if column_name in P_COLS:
                value = pfmt(raw)
            elif column_name == "in_hisp_gwas":
                value = "Yes" if raw else "-"
            elif column_name in ("AMR_n_sig", "AFR_n_sig", "EUR_n_sig", "EAS_n_sig"):
                value = int(raw) if raw and raw > 0 else ""
            elif column_name == "size_kb":
                value = f"{raw:,.1f}" if pd.notna(raw) else ""
            elif column_name in ("start", "end", "n_input_rows"):
                value = int(raw) if pd.notna(raw) else ""
            else:
                value = raw if pd.notna(raw) and raw != "" else ""
            cell = ws.cell(row=row_index, column=column_index, value=value if value != "" else None)
            cell.font = Font(name="Calibri", size=8.5)
            cell.border = thin()
            if column_name == "Class":
                cell.fill = PatternFill("solid", start_color=class_fill, end_color=class_fill)
                cell.font = Font(bold=True, name="Calibri", size=8.5, color="000000")
                cell.alignment = Alignment(horizontal="center", vertical="center")
            elif column_name == "Reason":
                cell.fill = PatternFill("solid", start_color="FAFAFA", end_color="FAFAFA")
                cell.font = Font(name="Calibri", size=7.5, color="333333")
                cell.alignment = Alignment(wrap_text=True, vertical="top")
            elif column_name in P_COLS:
                # Highlight p-values based on GWS threshold:
                # - p < GWS (5e-8): Green highlight + bold (genome-wide significant)
                # - GWS <= p < 1e-5: Yellow highlight (suggestive significance)
                # - p >= 1e-5: Light gray background (mention but not highlighted)
                # - NS (no value): White background with "NS" text
                if pd.isna(raw) or raw == "":
                    value = "NS"
                    cell.fill = PatternFill("solid", start_color="FFFFFF", end_color="FFFFFF")
                    cell.font = Font(name="Calibri", size=8.5, color="999999")
                elif raw < GWS:
                    cell.fill = PatternFill("solid", start_color="C6EFCE", end_color="C6EFCE")
                    cell.font = Font(bold=True, name="Calibri", size=8.5, color="276221")
                elif raw < 1e-5:
                    cell.fill = PatternFill("solid", start_color="FFEB9C", end_color="FFEB9C")
                    cell.font = Font(name="Calibri", size=8.5, color="9C6500")
                else:
                    # Non-significant p-values: light gray background to make them visible
                    cell.fill = PatternFill("solid", start_color="F0F0F0", end_color="F0F0F0")
                    cell.font = Font(name="Calibri", size=8.5, color="555555")
                cell.alignment = Alignment(horizontal="center", vertical="center")
            elif column_name == "orig_loci":
                cell.fill = PatternFill("solid", start_color=base_fill, end_color=base_fill)
                cell.font = Font(name="Calibri", size=7.5, color="444444")
                cell.alignment = Alignment(wrap_text=True, vertical="top")
            elif column_name.endswith("_region_best_snp") or column_name.endswith("_region_best_table"):
                cell.fill = PatternFill("solid", start_color=base_fill, end_color=base_fill)
                cell.font = Font(name="Calibri", size=7.5, color="444444")
                cell.alignment = Alignment(wrap_text=True, vertical="top")
            else:
                cell.fill = PatternFill("solid", start_color=base_fill, end_color=base_fill)
                cell.alignment = Alignment(horizontal="center", vertical="center")
        ws.row_dimensions[row_index].height = 40

    for column_name, column_index in col_idx.items():
        ws.column_dimensions[get_column_letter(column_index)].width = COL_WIDTHS.get(column_name, 12)

    ws2 = wb.create_sheet("Summary by Class")
    ws2.sheet_view.showGridLines = False
    ws2.cell(row=1, column=1, value="Locus Classification Summary — All Pairwise Sources")
    ws2.cell(row=1, column=1).font = Font(bold=True, size=13, color="1F3864", name="Calibri")
    for ci, header in enumerate(["Class", "# Loci", "Locus numbers"], start=1):
        cell = ws2.cell(row=2, column=ci, value=header)
        cell.fill = header_fill("1F3864")
        cell.font = Font(bold=True, color="FFFFFF", name="Calibri", size=10)
        cell.alignment = Alignment(horizontal="center", vertical="center")
        cell.border = thin()

    summary = df.groupby("Class").agg(N=("merged_id", "count"), IDs=("merged_id", lambda values: ", ".join(f"#{value}" for value in sorted(values)))).reset_index().sort_values("N", ascending=False) if not df.empty else pd.DataFrame(columns=["Class", "N", "IDs"])
    for row_index, (_, row) in enumerate(summary.iterrows(), start=3):
        cell = ws2.cell(row=row_index, column=1, value=row["Class"])
        cell.fill = PatternFill("solid", start_color=CLASS_COL.get(row["Class"], "FFFFFF"), end_color=CLASS_COL.get(row["Class"], "FFFFFF"))
        cell.font = Font(bold=True, name="Calibri", size=10)
        for column_index, value in enumerate([int(row["N"]), row["IDs"]], start=2):
            ws2.cell(row=row_index, column=column_index, value=value)

    Path(output_xlsx).parent.mkdir(parents=True, exist_ok=True)
    wb.save(output_xlsx)


def main() -> None:
    global GWS, MERGE_GAP, ALLOWED_SOURCES
    args = parse_args()
    cfg = read_yaml(args.config)
    
    # Read classification parameters from config
    classification_cfg = cfg.get("classification", {})
    primary_sources = sorted({infer_source_from_extracted_path(path) for path in args.extracted})
    if not primary_sources:
        primary_sources = sorted({infer_source_from_path(path) for path in args.inputs})

    # Always allow all datasets present in Step 8 inputs.
    ALLOWED_SOURCES = sorted({infer_source_from_path(path) for path in args.inputs})

    if isinstance(classification_cfg, dict):
        if "gws_threshold" in classification_cfg:
            GWS = float(classification_cfg["gws_threshold"])
        if "merge_gap_bp" in classification_cfg:
            MERGE_GAP = int(classification_cfg["merge_gap_bp"])
    
    # Build output schema with primary sources
    build_output_schema(primary_sources)
    
    gene_loc_path = cfg.get("gene_annotation", {}).get("gencode_gtf", "")
    genes_df = load_gene_loc(gene_loc_path)
    if genes_df is None:
        print(f"[warn] Gene annotation file not found or empty: {gene_loc_path}")
    pool = load_pool(args.inputs)
    output_tsv = Path(args.output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    if pool.empty:
        df = pd.DataFrame(columns=[column for _, columns in SECTIONS for column in columns])
        df.to_csv(output_tsv, sep="\t", index=False)
        write_workbook(df, args.output_xlsx)
        print("No loci available for classification; wrote empty outputs.")
        return
    pool, n_clusters = merge_loci(pool, MERGE_GAP)
    df = build_summary(pool, n_clusters, genes_df, primary_sources)
    df = append_primary_region_details(df, args.extracted, primary_sources)
    df.to_csv(output_tsv, sep="\t", index=False)
    write_workbook(df, args.output_xlsx)
    print(f"Saved Excel workbook: {args.output_xlsx}")
    print(f"Saved TSV table: {args.output_tsv}")
    print(f"Total merged loci: {n_clusters}")
    if not df.empty:
        primary_p_cols = [f"{source}_p" for source in primary_sources if f"{source}_p" in df.columns]
        if primary_p_cols:
            print("Loci with genome-wide signal in primary GWAS sources (p < 5e-8):")
            for p_col in primary_p_cols:
                source = p_col[:-2]
                n_gws = int((pd.to_numeric(df[p_col], errors="coerce") < GWS).sum())
                print(f"  - {source}: {n_gws}")
        best_p_cols = [f"{source}_region_best_p" for source in primary_sources if f"{source}_region_best_p" in df.columns]
        if best_p_cols:
            print("Loci with at least one Step 00 SNP inside region by source:")
            for p_col in best_p_cols:
                source = p_col.replace("_region_best_p", "")
                n_hits = int(pd.to_numeric(df[p_col], errors="coerce").notna().sum())
                print(f"  - {source}: {n_hits}")
        print("Class breakdown:")
        print(df["Class"].value_counts().to_string())


if __name__ == "__main__":
    main()
