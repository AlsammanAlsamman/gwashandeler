import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_gwastables

configfile: "configs/analysis.yml"

results_dir = get_results_dir()

def get_extraction_targets():
    """Generate all extraction targets from gwastables config."""
    targets = []
    for dataset_entry in config['gwastables']:
        dataset_name = dataset_entry['name']
        for table_entry in dataset_entry['tables']:
            table_name = table_entry['table_name']
            target = f"{results_dir}/extracted/{dataset_name}/{table_name}.done"
            targets.append(target)
    return targets

rule all:
    input:
        get_extraction_targets()

rule extract_gwas_table:
    input:
        gwas_file=lambda wc: _get_gwas_file(wc.dataset)
    output:
        table=f"{results_dir}/extracted/{{dataset}}/{{table}}.tsv",
        done=f"{results_dir}/extracted/{{dataset}}/{{table}}.done"
    params:
        dataset="{dataset}",
        table="{table}",
        columns=lambda wc: _get_table_columns(wc.dataset, wc.table)
    resources:
        mem_mb=32000,
        cores=2,
        time="01:00:00"
    log:
        f"{results_dir}/log/extract_{{dataset}}_{{table}}.log"
    shell:
        """
        mkdir -p $(dirname {output.table})
        python3 scripts/extract_gwas.py \
            --input {input.gwas_file} \
            --output {output.table} \
            --columns {params.columns} \
            --dataset {params.dataset} \
            --table {params.table} \
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
    """Get column names for a specific table as comma-separated string."""
    for entry in config['gwastables']:
        if entry['name'] == dataset:
            for t in entry['tables']:
                if t['table_name'] == table:
                    return ','.join(t['columns'])
    raise ValueError(f"Table {table} not found in dataset {dataset}")
