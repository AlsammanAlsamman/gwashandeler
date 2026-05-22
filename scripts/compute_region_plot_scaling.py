#!/usr/bin/env python3
"""Compute region Manhattan y-axis scaling from region subset files.

Reads pre-extracted subset TSVs for one region, finds the most significant
p-value across datasets/tables, and sets that region's y-axis max to 10%
above -log10(min_p).
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


P_PATTERNS = ("p_", "pval", "pvalue", "^p$")
EXCLUDE_NAMES = {"pos", "position", "chrom", "chr", "bp", "allele"}


def find_col(header: list[str], patterns: tuple[str, ...], exclude: set[str] | None = None) -> str | None:
    exclude = exclude or set()
    for pattern in patterns:
        for col in header:
            col_norm = col.strip().lower()
            if col_norm in exclude:
                continue
            if re.search(pattern, col_norm, re.IGNORECASE):
                return col
    return None


def min_positive_pvalue_from_subset(tsv_path: Path) -> float | None:
    local_min = None
    with tsv_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header:
            return None

        p_col = find_col(header, P_PATTERNS, EXCLUDE_NAMES)
        if not p_col:
            return None

        p_idx = header.index(p_col)

        for row in reader:
            if p_idx >= len(row):
                continue
            p_cell = row[p_idx].strip()
            if not p_cell:
                continue
            try:
                p_val = float(p_cell)
            except ValueError:
                continue

            if 0 < p_val <= 1 and (local_min is None or p_val < local_min):
                local_min = p_val

    return local_min


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--subsets-root", required=True)
    parser.add_argument("--region-id", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    subsets_root = Path(args.subsets_root)
    region_id = args.region_id
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not subsets_root.exists():
        raise FileNotFoundError(f"Subset root not found: {subsets_root}")

    global_min_p = None

    for subset_file in subsets_root.rglob(f"{region_id}__*.tsv"):
        local_min = min_positive_pvalue_from_subset(subset_file)
        if local_min is not None and (global_min_p is None or local_min < global_min_p):
            global_min_p = local_min

    if global_min_p is None:
        # Safe fallback if no valid p-values were found.
        global_max_log10p = 8.0
    else:
        global_max_log10p = -math.log10(global_min_p)

    # 10% headroom above the most significant SNP in configured regions.
    global_ymax = round(max(8.0, global_max_log10p * 1.10), 3)

    rows = [
        ("global_min_p", global_min_p if global_min_p is not None else "NA"),
        ("global_max_log10p", global_max_log10p),
        ("global_ymax", global_ymax),
        ("threshold_suggestive", 5e-5),
        ("threshold_genomewide", 5e-8),
    ]

    with out_path.open("w", encoding="utf-8") as fh:
        fh.write("key\tvalue\n")
        for key, val in rows:
            fh.write(f"{key}\t{val}\n")

    print(f"Wrote region plot scaling for {region_id} to {out_path}")


if __name__ == "__main__":
    main()
