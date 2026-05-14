import sys
sys.path.append("utils")
from bioconfigme import get_results_dir

configfile: "configs/analysis.yml"

results_dir = get_results_dir()

def get_manhattan_targets():
    """Generate all Manhattan plot targets from gwastables config."""
    targets = []
    for dataset_entry in config['gwastables']:
        dataset_name = dataset_entry['name']
        for table_entry in dataset_entry['tables']:
            table_name = table_entry['table_name']
            target = f"{results_dir}/plots/{dataset_name}/{table_name}_manhattan.done"
            targets.append(target)
    return targets

rule all:
    input:
        get_manhattan_targets()

rule plot_gwas_manhattan:
    input:
        gwas_file=lambda wc: f"{results_dir}/extracted/{wc.dataset}/{wc.table}.tsv"
    output:
        png=f"{results_dir}/plots/{{dataset}}/{{table}}_manhattan.png",
        done=f"{results_dir}/plots/{{dataset}}/{{table}}_manhattan.done"
    params:
        dataset="{dataset}",
        table="{table}",
        columns=lambda wc: _get_table_columns(wc.dataset, wc.table)
    resources:
        mem_mb=65536,
        cores=2,
        time="02:00:00"
    log:
        f"{results_dir}/log/plot_manhattan_{{dataset}}_{{table}}.log"
    shell:
        """
        mkdir -p $(dirname {output.png})
        bash scripts/plot_manhattan_wrapper.sh \
            {input.gwas_file} \
            {output.png} \
            {params.dataset} \
            {params.table} \
            2>&1 | tee {log}
        touch {output.done}
        """

def _get_gwas_file(dataset):
    """Get the input file for a dataset."""
    for entry in config['gwastables']:
        if entry['name'] == dataset:
            return entry['file']
    raise ValueError(f"Dataset {dataset} not found in config")

def _get_table_columns(dataset, table):
    """Get column info for a specific table."""
    for entry in config['gwastables']:
        if entry['name'] == dataset:
            for t in entry['tables']:
                if t['table_name'] == table:
                    return ','.join(t['columns'])
    raise ValueError(f"Table {table} not found in dataset {dataset}")

