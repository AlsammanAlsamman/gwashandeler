import sys
sys.path.append("utils")
from bioconfigme import get_results_dir

configfile: "configs/analysis.yml"

results_dir = get_results_dir()


def get_dataset_entries():
    return config["gwastables"]


def get_dataset_names():
    return [d["name"] for d in get_dataset_entries()]


def get_table_loci_inputs(dataset):
    for d in get_dataset_entries():
        if d["name"] == dataset:
            return [f"{results_dir}/loci/{dataset}/{t['table_name']}_loci.tsv" for t in d["tables"]]
    raise ValueError(f"Dataset not found: {dataset}")


def get_merged_dataset_targets():
    return [f"{results_dir}/loci/{d}/merged_dataset_loci.done" for d in get_dataset_names()]


rule all:
    input:
        get_merged_dataset_targets()


rule merge_loci_per_dataset:
    input:
        table_loci=lambda wc: get_table_loci_inputs(wc.dataset)
    output:
        loci=f"{results_dir}/loci/{{dataset}}/{{dataset}}_loci.tsv",
        done=f"{results_dir}/loci/{{dataset}}/merged_dataset_loci.done"
    resources:
        mem_mb=32000,
        cores=1,
        time="00:30:00"
    log:
        f"{results_dir}/log/loci_merge_dataset_{{dataset}}.log"
    shell:
        """
        mkdir -p $(dirname {output.loci})
        python3 scripts/merge_dataset_loci.py \
            --dataset {wildcards.dataset} \
            --inputs {input.table_loci} \
            --output {output.loci} \
            2>&1 | tee {log}
        touch {output.done}
        """
