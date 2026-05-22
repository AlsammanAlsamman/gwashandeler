#!/usr/bin/env python3
"""Extract configured regions from one GWAS dataset into per-region/per-table files.

This script reads a GWAS input file once for a dataset and writes subset TSV files
for all configured regions and tables, so downstream plotting does not re-read full files.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import yaml


CHR_PATTERNS = ("chrom", "^chr$")
POS_PATTERNS = ("^pos$", "^position$", "^bp$")


def norm_chr(value: str) -> str:
    return re.sub(r"^chr", "", str(value).strip(), flags=re.IGNORECASE)


def sanitize_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())


def find_col(header: list[str], patterns: tuple[str, ...]) -> str | None:
    for pattern in patterns:
        for col in header:
            if re.search(pattern, col.strip().lower(), re.IGNORECASE):
                return col
    return None


def load_regions(regions_path: Path) -> list[dict[str, object]]:
    regions: list[dict[str, object]] = []
    with regions_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            locus = str(row.get("locus", "")).strip()
            chr_raw = str(row.get("chr", "")).strip()
            start_raw = str(row.get("start", "")).strip()
            end_raw = str(row.get("end", "")).strip()
            if not locus or not chr_raw or not start_raw or not end_raw:
                continue

            try:
                start = int(float(start_raw))
                end = int(float(end_raw))
            except ValueError:
                continue

            if end < start:
                start, end = end, start

            chr_norm = norm_chr(chr_raw)
            region_id = sanitize_name(f"{locus}_chr{chr_norm}_{start}_{end}")
            regions.append(
                {
                    "locus": locus,
                    "chr": chr_norm,
                    "start": start,
                    "end": end,
                    "region_id": region_id,
                }
            )

    return regions


def load_dataset_config(cfg_path: Path, dataset_name: str) -> dict[str, object]:
    with cfg_path.open("r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    for ds in cfg.get("gwastables", []):
        if ds.get("name") == dataset_name:
            return ds
    raise ValueError(f"Dataset '{dataset_name}' not found in gwastables")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--regions", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    cfg_path = Path(args.config)
    dataset_name = args.dataset
    input_path = Path(args.input)
    regions_path = Path(args.regions)
    output_dir = Path(args.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    dataset_cfg = load_dataset_config(cfg_path, dataset_name)
    tables = dataset_cfg.get("tables", [])
    delimiter = str(dataset_cfg.get("delimiter", "\t"))

    regions = load_regions(regions_path)
    if not regions:
        raise ValueError("No valid regions found")

    regions_by_chr: dict[str, list[dict[str, object]]] = {}
    for region in regions:
        regions_by_chr.setdefault(str(region["chr"]), []).append(region)

    handles: dict[tuple[str, str], object] = {}
    writers: dict[tuple[str, str], csv.writer] = {}

    with input_path.open("r", encoding="utf-8", newline="") as infh:
        reader = csv.DictReader(infh, delimiter=delimiter)
        header = [str(h) for h in (reader.fieldnames or [])]
        if not header:
            raise ValueError(f"Input file has no header: {input_path}")

        chr_col = find_col(header, CHR_PATTERNS)
        pos_col = find_col(header, POS_PATTERNS)
        if not chr_col or not pos_col:
            raise ValueError("Could not detect chromosome/position columns in input file")

        table_cols: dict[str, list[str]] = {}
        for table in tables:
            table_name = str(table["table_name"])
            cols = [str(c) for c in table.get("columns", [])]
            missing = [c for c in cols if c not in header]
            if missing:
                raise ValueError(
                    f"Dataset {dataset_name}, table {table_name}: missing columns in input: {', '.join(missing)}"
                )
            table_cols[table_name] = cols

        for region in regions:
            region_id = str(region["region_id"])
            for table_name, cols in table_cols.items():
                out_path = output_dir / f"{region_id}__{table_name}.tsv"
                fh = out_path.open("w", encoding="utf-8", newline="")
                writer = csv.writer(fh, delimiter="\t")
                writer.writerow(cols)
                key = (region_id, table_name)
                handles[key] = fh
                writers[key] = writer

        for row in reader:
            chr_cell = norm_chr(row.get(chr_col, ""))
            if not chr_cell:
                continue
            chr_regions = regions_by_chr.get(chr_cell, [])
            if not chr_regions:
                continue

            pos_raw = str(row.get(pos_col, "")).strip()
            if not pos_raw:
                continue
            try:
                pos_val = int(float(pos_raw))
            except ValueError:
                continue

            for region in chr_regions:
                start = int(region["start"])
                end = int(region["end"])
                if not (start <= pos_val <= end):
                    continue

                region_id = str(region["region_id"])
                for table_name, cols in table_cols.items():
                    writer = writers[(region_id, table_name)]
                    writer.writerow([row.get(col, "") for col in cols])

    for fh in handles.values():
        fh.close()

    print(f"Extracted region subsets for dataset {dataset_name} into {output_dir}")


if __name__ == "__main__":
    main()
