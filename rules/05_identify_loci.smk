import glob
import os
import re
import sys
sys.path.append("utils")
from bioconfigme import get_results_dir

configfile: "configs/analysis.yml"

results_dir = get_results_dir()


def get_dataset_entries():
    return config["gwastables"]


def get_dataset_table_pairs():
    pairs = []
    for d in get_dataset_entries():
        for t in d["tables"]:
            pairs.append((d["name"], t["table_name"]))
    return pairs


def get_ref_panel_for_dataset(dataset):
    for d in get_dataset_entries():
        if d["name"] == dataset:
            return d["reference_panel"]
    raise ValueError(f"Dataset not found: {dataset}")


def get_chromosomes_for_dataset_table(dataset, table):
    pattern = f"{results_dir}/gwas_subset/{dataset}/{table}_subset_dir/*_chr*_subset.tsv"
    files = glob.glob(pattern)
    chromosomes = set()
    for f in files:
        m = re.search(r"_chr([^_/]+)_subset\.tsv$", os.path.basename(f))
        if m:
            chromosomes.add(m.group(1))

    chrom_mode = config["loci_identiffication"].get("chromosomes", "autosomes")
    if chrom_mode == "autosomes":
        chromosomes = {c for c in chromosomes if c.isdigit() and 1 <= int(c) <= 22}

    def sort_key(ch):
        return (0, int(ch)) if ch.isdigit() else (1, ch)

    return sorted(chromosomes, key=sort_key)


def get_chr_loci_inputs(dataset, table):
    chroms = get_chromosomes_for_dataset_table(dataset, table)
    return [f"{results_dir}/loci_by_chr/{dataset}/{table}/chr{c}_loci.tsv" for c in chroms]


def get_merged_targets():
    targets = []
    for dataset, table in get_dataset_table_pairs():
        targets.append(f"{results_dir}/loci/{dataset}/{table}_loci.done")
    return targets


rule all:
    input:
        get_merged_targets()


rule identify_loci_per_chromosome:
    input:
        subset=lambda wc: f"{results_dir}/gwas_subset/{wc.dataset}/{wc.table}_subset_dir/{wc.dataset}_{wc.table}_chr{wc.chrom}_subset.tsv",
        ref_bim=lambda wc: f"{results_dir}/ref_panels_filtered/{get_ref_panel_for_dataset(wc.dataset)}/{get_ref_panel_for_dataset(wc.dataset)}_chr{wc.chrom}.bim",
        ref_done=lambda wc: f"{results_dir}/ref_panels_filtered/{get_ref_panel_for_dataset(wc.dataset)}/{get_ref_panel_for_dataset(wc.dataset)}_chr{wc.chrom}.done"
    output:
        loci=f"{results_dir}/loci_by_chr/{{dataset}}/{{table}}/chr{{chrom}}_loci.tsv"
    params:
        ref_panel=lambda wc: get_ref_panel_for_dataset(wc.dataset),
        lead_p_threshold=lambda wc: config["loci_identiffication"]["lead_p_threshold"],
        merge_distance_kb=lambda wc: config["loci_identiffication"]["merge_distance_kb"],
        min_snps_per_locus=lambda wc: config["loci_identiffication"]["min_snps_per_locus"],
        require_refpanel_match=lambda wc: 1 if config["loci_identiffication"].get("require_refpanel_match", True) else 0
    resources:
        mem_mb=64000,
        cores=2,
        time="01:00:00"
    log:
        f"{results_dir}/log/loci_chr_{{dataset}}_{{table}}_chr{{chrom}}.log"
    shell:
        """
        mkdir -p $(dirname {output.loci})
        python3 scripts/identify_loci_per_chromosome.py \
            --subset {input.subset} \
            --ref-bim {input.ref_bim} \
            --output {output.loci} \
            --dataset {wildcards.dataset} \
            --table {wildcards.table} \
            --ref-panel {params.ref_panel} \
            --lead-p-threshold {params.lead_p_threshold} \
            --merge-distance-kb {params.merge_distance_kb} \
            --min-snps-per-locus {params.min_snps_per_locus} \
            --require-refpanel-match {params.require_refpanel_match} \
            2>&1 | tee {log}
        """


rule merge_loci_per_table:
    input:
        chr_loci=lambda wc: get_chr_loci_inputs(wc.dataset, wc.table)
    output:
        loci=f"{results_dir}/loci/{{dataset}}/{{table}}_loci.tsv",
        done=f"{results_dir}/loci/{{dataset}}/{{table}}_loci.done"
    resources:
        mem_mb=32000,
        cores=1,
        time="00:30:00"
    log:
        f"{results_dir}/log/loci_merge_{{dataset}}_{{table}}.log"
    shell:
        """
        mkdir -p $(dirname {output.loci})
        python3 scripts/merge_loci_tables.py \
            --dataset {wildcards.dataset} \
            --table {wildcards.table} \
            --input-dir {results_dir}/loci_by_chr/{wildcards.dataset}/{wildcards.table} \
            --output {output.loci} \
            2>&1 | tee {log}
        touch {output.done}
        """
