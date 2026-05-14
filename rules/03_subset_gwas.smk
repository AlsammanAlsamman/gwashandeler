import sys
sys.path.append("utils")
from bioconfigme import get_results_dir
import pandas as pd
import re

configfile: "configs/analysis.yml"

results_dir = get_results_dir()

def get_chromosomes_in_file(filepath):
    """
    Read input file and return sorted list of unique chromosomes.
    """
    try:
        df = pd.read_csv(filepath, sep='\t', dtype=str, low_memory=False, nrows=10000)
        chr_patterns = [r'^chrom', r'^chr', r'^CHROM', r'^CHR']
        chr_col = None
        for pattern in chr_patterns:
            matching = [col for col in df.columns if re.search(pattern, col, re.IGNORECASE)]
            if matching:
                chr_col = matching[0]
                break
        if chr_col is None:
            return ['unknown']
        return sorted(df[chr_col].unique())
    except:
        return ['unknown']

def get_subset_targets():
    """Generate all subset GWAS targets from gwastables config."""
    targets = []
    for dataset_entry in config['gwastables']:
        dataset_name = dataset_entry['name']
        input_file = f"{results_dir}/extracted/{dataset_name}/"
        for table_entry in dataset_entry['tables']:
            table_name = table_entry['table_name']
            target = f"{results_dir}/gwas_subset/{dataset_name}/{table_name}_subset.done"
            targets.append(target)
    return targets

rule all:
    input:
        get_subset_targets()

rule subset_gwas_by_pvalue:
    input:
        gwas_file=lambda wc: f"{results_dir}/extracted/{wc.dataset}/{wc.table}.tsv"
    output:
        subset_files=directory(f"{results_dir}/gwas_subset/{{dataset}}/{{table}}_subset_dir"),
        done=f"{results_dir}/gwas_subset/{{dataset}}/{{table}}_subset.done"
    params:
        dataset="{dataset}",
        table="{table}",
        threshold=lambda wc: config['subsetting']['pvalue_threshold']
    resources:
        mem_mb=64000,
        cores=2,
        time="01:00:00"
    log:
        f"{results_dir}/log/subset_gwas_{{dataset}}_{{table}}.log"
    shell:
        """
        mkdir -p {output.subset_files}
        python3 scripts/subset_gwas.py \
            --input {input.gwas_file} \
            --output-dir {output.subset_files} \
            --threshold {params.threshold} \
            --dataset {params.dataset} \
            --table {params.table} \
            2>&1 | tee {log}
        touch {output.done}
        """
