import re
import csv
import sys

sys.path.append("utils")
from bioconfigme import get_results_dir

configfile: "configs/analysis.yml"

results_dir = get_results_dir()


def norm_chr(value):
    return re.sub(r"^chr", "", str(value).strip(), flags=re.IGNORECASE)


def sanitize_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())


def load_regions():
    regions_path = config.get("regions", "input/regions.tsv")
    required = {"locus", "chr", "start", "end"}
    rows = []

    with open(regions_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"regions file {regions_path} is empty or missing a header")

        missing = [col for col in required if col not in reader.fieldnames]
        if missing:
            raise ValueError(
                f"regions file {regions_path} missing required columns: {', '.join(missing)}"
            )

        for row in reader:
            locus = str(row.get("locus", "")).strip()
            chr_value = norm_chr(row.get("chr", ""))
            try:
                start = int(float(row.get("start", "")))
                end = int(float(row.get("end", "")))
            except ValueError:
                continue

            if not locus or not chr_value:
                continue

            rows.append(
                {
                    "locus": locus,
                    "chr": chr_value,
                    "start": start,
                    "end": end,
                    "region_id": sanitize_name(f"{locus}_chr{chr_value}_{start}_{end}"),
                }
            )

    return rows


REGIONS = load_regions()
REGION_BY_ID = {row["region_id"]: row for row in REGIONS}
GLOBAL_SCALING = f"{results_dir}/plots/regions/_global_scaling.tsv"


def get_plot_targets():
    targets = []
    for region in REGIONS:
        for dataset_entry in config["gwastables"]:
            dataset_name = dataset_entry["name"]
            for table_entry in dataset_entry["tables"]:
                table_name = table_entry["table_name"]
                targets.append(
                    f"{results_dir}/plots/regions/{region['region_id']}/{dataset_name}__{table_name}.done"
                )
    return targets


def get_original_table_path(dataset, table):
    for dataset_entry in config["gwastables"]:
        if dataset_entry["name"] == dataset:
            return dataset_entry["file"]
    raise ValueError(f"Dataset {dataset} not found in config")


def get_table_columns(dataset, table):
    for dataset_entry in config["gwastables"]:
        if dataset_entry["name"] == dataset:
            for table_entry in dataset_entry["tables"]:
                if table_entry["table_name"] == table:
                    return ",".join(table_entry["columns"])
    raise ValueError(f"Table {table} not found in dataset {dataset}")


rule all:
    input:
        get_plot_targets()


rule compute_region_plot_scaling:
    output:
        scaling=GLOBAL_SCALING
    resources:
        mem_mb=16000,
        cores=1,
        time="01:30:00"
    log:
        f"{results_dir}/log/compute_region_plot_scaling.log"
    shell:
        """
        mkdir -p $(dirname {output.scaling}) $(dirname {log})
        python3 scripts/compute_region_plot_scaling.py \
            --config configs/analysis.yml \
            --output {output.scaling} \
            2>&1 | tee {log}
        """


rule plot_region_gwas_manhattan:
    input:
        gwas_file=lambda wc: get_original_table_path(wc.dataset, wc.table),
        scaling=GLOBAL_SCALING
    output:
        png=f"{results_dir}/plots/regions/{{region_id}}/{{dataset}}__{{table}}.png",
        done=f"{results_dir}/plots/regions/{{region_id}}/{{dataset}}__{{table}}.done"
    params:
        region_name=lambda wc: REGION_BY_ID[wc.region_id]["locus"],
        region_chr=lambda wc: REGION_BY_ID[wc.region_id]["chr"],
        region_start=lambda wc: REGION_BY_ID[wc.region_id]["start"],
        region_end=lambda wc: REGION_BY_ID[wc.region_id]["end"],
        columns=lambda wc: get_table_columns(wc.dataset, wc.table)
    resources:
        mem_mb=16000,
        cores=1,
        time="00:30:00"
    log:
        f"{results_dir}/log/plot_region_{{region_id}}_{{dataset}}_{{table}}.log"
    shell:
        """
        mkdir -p $(dirname {output.png}) $(dirname {log})
        bash scripts/plot_region_manhattan_wrapper.sh \
            {input.gwas_file} \
            {output.png} \
            "{params.region_name}" \
            {params.region_chr} \
            {params.region_start} \
            {params.region_end} \
            {wildcards.dataset} \
            {wildcards.table} \
            "{params.columns}" \
            {input.scaling} \
            2>&1 | tee {log}
        touch {output.done}
        """