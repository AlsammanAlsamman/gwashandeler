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
DATASETS = [entry["name"] for entry in config["gwastables"]]


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


def region_scaling_path(region_id):
    return f"{results_dir}/plots/regions/{region_id}/mpvalue.txt"


def get_region_scaling_targets():
    return [region_scaling_path(region["region_id"]) for region in REGIONS]


def region_best_hits_tsv_path():
    return f"{results_dir}/plots/regions/region_best_hits.tsv"


def region_best_hits_xlsx_path():
    return f"{results_dir}/plots/regions/region_best_hits.xlsx"


def get_dataset_tables(dataset):
    for dataset_entry in config["gwastables"]:
        if dataset_entry["name"] == dataset:
            return [table_entry["table_name"] for table_entry in dataset_entry["tables"]]
    raise ValueError(f"Dataset {dataset} not found in config")


def region_subset_path(region_id, dataset, table):
    return f"{results_dir}/gwas_region_subsets/{dataset}/{region_id}__{table}.tsv"


def dataset_subset_done_path(dataset):
    return f"{results_dir}/gwas_region_subsets/{dataset}/.extract_complete"


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
        get_region_scaling_targets(),
        region_best_hits_tsv_path(),
        region_best_hits_xlsx_path(),
        get_plot_targets()


rule extract_dataset_region_subsets:
    input:
        gwas_file=lambda wc: get_original_table_path(wc.dataset, None)
    output:
        done=dataset_subset_done_path("{dataset}")
    params:
        dataset=lambda wc: wc.dataset,
        output_dir=lambda wc: f"{results_dir}/gwas_region_subsets/{wc.dataset}"
    resources:
        mem_mb=20000,
        cores=1,
        time="01:30:00"
    log:
        f"{results_dir}/log/extract_dataset_region_subsets_{{dataset}}.log"
    shell:
        """
        mkdir -p {params.output_dir} $(dirname {log})
        python3 scripts/extract_dataset_regions.py \
            --config configs/analysis.yml \
            --dataset {params.dataset} \
            --input {input.gwas_file} \
            --regions input/regions.tsv \
            --output-dir {params.output_dir} \
            2>&1 | tee {log}
        touch {output.done}
        """


rule compute_region_plot_scaling:
    input:
        subset_done=expand(dataset_subset_done_path("{dataset}"), dataset=DATASETS)
    output:
        scaling_tsv=f"{results_dir}/plots/regions/{{region_id}}/mpvalue.txt"
    log:
        f"{results_dir}/log/compute_region_plot_scaling_{{region_id}}.log"
    shell:
        """
        mkdir -p $(dirname {output.scaling_tsv}) $(dirname {log})
        python3 scripts/compute_region_plot_scaling.py \
            --subsets-root {results_dir}/gwas_region_subsets \
            --region-id {wildcards.region_id} \
            --output {output.scaling_tsv} \
            2>&1 | tee {log}
        """


rule summarize_region_best_hits:
    input:
        subset_done=expand(dataset_subset_done_path("{dataset}"), dataset=DATASETS)
    output:
        tsv=region_best_hits_tsv_path(),
        xlsx=region_best_hits_xlsx_path()
    log:
        f"{results_dir}/log/summarize_region_best_hits.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        python3 scripts/summarize_region_best_hits.py \
            --config configs/analysis.yml \
            --regions input/regions.tsv \
            --subsets-root {results_dir}/gwas_region_subsets \
            --output-tsv {output.tsv} \
            --output-xlsx {output.xlsx} \
            2>&1 | tee {log}
        """


rule plot_region_gwas_manhattan:
    input:
        subset_done=lambda wc: dataset_subset_done_path(wc.dataset),
        scaling_tsv=lambda wc: region_scaling_path(wc.region_id)
    output:
        png=f"{results_dir}/plots/regions/{{region_id}}/{{dataset}}__{{table}}.png",
        done=f"{results_dir}/plots/regions/{{region_id}}/{{dataset}}__{{table}}.done"
    params:
        subset_file=lambda wc: region_subset_path(wc.region_id, wc.dataset, wc.table),
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
            {params.subset_file} \
            {output.png} \
            "{params.region_name}" \
            {params.region_chr} \
            {params.region_start} \
            {params.region_end} \
            {wildcards.dataset} \
            {wildcards.table} \
            "{params.columns}" \
            {input.scaling_tsv} \
            2>&1 | tee {log}
        touch {output.done}
        """