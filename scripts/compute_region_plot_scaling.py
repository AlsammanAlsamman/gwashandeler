#!/usr/bin/env python3
"""Compute global region-plot scaling across all configured GWAS datasets/tables.

Outputs a TSV with fixed y-axis and threshold positions so all region plots are
directly comparable and keep headroom for the most significant SNPs.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import yaml


P_PATTERNS = ("p", "pval", "pvalue")
EXCLUDE_NAMES = {"pos", "position", "chrom", "chr", "bp", "allele"}


def detect_p_column(columns: list[str]) -> str | None:
    for col in columns:
        cl = col.strip().lower()
        if cl in EXCLUDE_NAMES:
            continue
        if any(tok in cl for tok in P_PATTERNS):
            return col
    return None


def min_positive_pvalue_from_tsv(tsv_path: Path, p_col: str) -> float | None:
    local_min = None
    with tsv_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header or p_col not in header:
            return None
        p_idx = header.index(p_col)

        for row in reader:
            if p_idx >= len(row):
                continue
            cell = row[p_idx].strip()
            if not cell:
                continue
            try:
                p_val = float(cell)
            except ValueError:
                continue

            if 0 < p_val <= 1:
                if local_min is None or p_val < local_min:
                    local_min = p_val

    return local_min


def scientific_label_from_y(y_val: float) -> str:
    p = 10 ** (-y_val)
    return f"{p:.1e}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    cfg_path = Path(args.config)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with cfg_path.open("r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    gwastables = cfg.get("gwastables", [])
    if not gwastables:
        raise ValueError("No gwastables found in config")

    global_min_p = None

    for ds in gwastables:
        file_path = Path(ds["file"])
        if not file_path.exists():
            continue

        for table in ds.get("tables", []):
            table_cols = [str(c) for c in table.get("columns", [])]
            p_col = detect_p_column(table_cols)
            if p_col is None:
                continue

            local_min = min_positive_pvalue_from_tsv(file_path, p_col)
            if local_min is not None and (global_min_p is None or local_min < global_min_p):
                global_min_p = local_min

    if global_min_p is None:
        # Safe fallback if no valid p-values were found.
        global_max_log10p = 8.0
    else:
        global_max_log10p = -math.log10(global_min_p)

    # Keep room above the strongest SNP so peaks are never clipped.
    global_ymax = max(8.0, round(global_max_log10p * 1.2 + 0.5, 2))

    # Fixed thresholds in y-space derived from the global significance range.
    threshold_low_y = max(1.5, round(global_max_log10p * 0.35, 2))
    threshold_high_y = max(threshold_low_y + 0.5, round(global_max_log10p * 0.60, 2))

    # Keep thresholds inside visible y-range with headroom.
    threshold_high_y = min(threshold_high_y, round(global_ymax * 0.90, 2))
    threshold_low_y = min(threshold_low_y, round(threshold_high_y - 0.5, 2))

    rows = [
        ("global_min_p", global_min_p if global_min_p is not None else "NA"),
        ("global_max_log10p", global_max_log10p),
        ("global_ymax", global_ymax),
        ("threshold_low_y", threshold_low_y),
        ("threshold_high_y", threshold_high_y),
        ("threshold_low_label", scientific_label_from_y(threshold_low_y)),
        ("threshold_high_label", scientific_label_from_y(threshold_high_y)),
    ]

    with out_path.open("w", encoding="utf-8") as fh:
        fh.write("key\tvalue\n")
        for key, val in rows:
            fh.write(f"{key}\t{val}\n")

    print(f"Wrote global region plot scaling to {out_path}")


if __name__ == "__main__":
    main()
