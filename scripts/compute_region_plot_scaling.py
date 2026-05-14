#!/usr/bin/env python3
"""Compute global region-plot scaling across all configured GWAS datasets/tables.

Outputs a TSV with fixed y-axis and threshold positions so all region plots are
directly comparable and keep headroom for the most significant SNPs.
"""

from __future__ import annotations

import argparse
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


def min_positive_from_file_for_columns(tsv_path: Path, candidate_cols: list[str]) -> float | None:
    """Scan a TSV once and return the minimum positive p-value across candidate columns."""
    local_min = None

    with tsv_path.open("r", encoding="utf-8", newline="") as fh:
        header_line = fh.readline()
        if not header_line:
            return None

        header = [h.strip() for h in header_line.rstrip("\n\r").split("\t")]
        idxs = [header.index(col) for col in candidate_cols if col in header]
        if not idxs:
            return None

        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n\r").split("\t")
            for idx in idxs:
                if idx >= len(parts):
                    continue
                cell = parts[idx].strip()
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
        ds_name = ds.get("name", "unknown")
        file_path = Path(ds["file"])
        if not file_path.exists():
            print(f"[WARN] Missing GWAS file for dataset {ds_name}: {file_path}", flush=True)
            continue

        p_cols = []
        for table in ds.get("tables", []):
            table_cols = [str(c) for c in table.get("columns", [])]
            p_col = detect_p_column(table_cols)
            if p_col:
                p_cols.append(p_col)

        # Unique while preserving order
        seen = set()
        p_cols = [c for c in p_cols if not (c in seen or seen.add(c))]

        if not p_cols:
            print(f"[WARN] No p-value columns found for dataset {ds_name}", flush=True)
            continue

        print(f"[INFO] Scanning {ds_name}: {file_path} with columns {p_cols}", flush=True)
        local_min = min_positive_from_file_for_columns(file_path, p_cols)
        if local_min is not None and (global_min_p is None or local_min < global_min_p):
            global_min_p = local_min
            print(f"[INFO] Updated global minimum p-value to {global_min_p}", flush=True)

    if global_min_p is None:
        # Safe fallback if no valid p-values were found.
        global_max_log10p = 8.0
    else:
        global_max_log10p = -math.log10(global_min_p)

    # Keep room above the strongest SNP so peaks are never clipped.
    global_ymax = max(8.0, round(global_max_log10p * 1.2 + 0.5, 2))

    # Lower fixed thresholds in y-space, per user request.
    threshold_low_y = 2.5
    threshold_high_y = 4.5

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

    print(f"Wrote global region plot scaling to {out_path}", flush=True)


if __name__ == "__main__":
    main()
