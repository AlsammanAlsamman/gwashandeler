import sys
import glob
import os
import re
sys.path.append("utils")
from bioconfigme import get_results_dir

configfile: "configs/analysis.yml"

results_dir = get_results_dir()

def get_all_ref_panels():
    """Extract all reference panels from config."""
    return list(config['ref_panels'].keys())

def get_all_chromosomes():
    """Get chromosomes present in subset outputs; fallback to autosomes if not found."""
    subset_pattern = f"{results_dir}/gwas_subset/*/*_subset_dir/*_chr*_subset.tsv"
    files = glob.glob(subset_pattern)
    chromosomes = set()

    for f in files:
        m = re.search(r'_chr([^_/]+)_subset\.tsv$', os.path.basename(f))
        if m:
            chromosomes.add(m.group(1))

    if not chromosomes:
        return [str(i) for i in range(1, 23)]

    def chr_sort_key(ch):
        return (0, int(ch)) if ch.isdigit() else (1, ch)

    return sorted(chromosomes, key=chr_sort_key)

def get_snp_aggregation_targets():
    """Aggregate SNP files per chromosome."""
    return expand(f"{results_dir}/snp_lists/chr{{chrom}}_snps.txt", chrom=get_all_chromosomes())

def get_subset_done_targets():
    """Return all subset completion markers from rule 03 output structure."""
    targets = []
    for dataset_entry in config['gwastables']:
        dataset_name = dataset_entry['name']
        for table_entry in dataset_entry['tables']:
            table_name = table_entry['table_name']
            targets.append(f"{results_dir}/gwas_subset/{dataset_name}/{table_name}_subset.done")
    return targets

def get_ref_panel_filter_targets():
    """Filter ref panels per chromosome per panel."""
    ref_panels = get_all_ref_panels()
    targets = []
    for panel in ref_panels:
        for chrom in get_all_chromosomes():
            target = f"{results_dir}/ref_panels_filtered/{panel}/{panel}_chr{chrom}.done"
            targets.append(target)
    return targets

rule all:
    input:
        get_snp_aggregation_targets(),
        get_ref_panel_filter_targets()

rule aggregate_snps_from_subsets:
    input:
        subset_done=get_subset_done_targets()
    output:
        snp_lists=expand(f"{results_dir}/snp_lists/chr{{chrom}}_snps.txt", chrom=get_all_chromosomes()),
        done=f"{results_dir}/snp_lists/aggregation.done"
    params:
        subset_dir=f"{results_dir}/gwas_subset",
        output_dir=f"{results_dir}/snp_lists"
    resources:
        mem_mb=128000,
        cores=4,
        time="02:00:00"
    log:
        f"{results_dir}/log/aggregate_snps.log"
    shell:
        """
        mkdir -p {params.output_dir}
        python3 scripts/aggregate_snps_from_subsets.py \
            --subset-dir {params.subset_dir} \
            --output-dir {params.output_dir} \
            2>&1 | tee {log}
        touch {output.done}
        """

rule filter_ref_panels_by_snp:
    input:
        snp_list=f"{results_dir}/snp_lists/chr{{chrom}}_snps.txt",
        aggregation_done=f"{results_dir}/snp_lists/aggregation.done"
    output:
        bed=f"{results_dir}/ref_panels_filtered/{{panel}}/{{panel}}_chr{{chrom}}.bed",
        bim=f"{results_dir}/ref_panels_filtered/{{panel}}/{{panel}}_chr{{chrom}}.bim",
        fam=f"{results_dir}/ref_panels_filtered/{{panel}}/{{panel}}_chr{{chrom}}.fam",
        done=f"{results_dir}/ref_panels_filtered/{{panel}}/{{panel}}_chr{{chrom}}.done"
    params:
        ref_panel_path=lambda wc: config['ref_panels'][wc.panel],
        chrom="{chrom}",
        panel="{panel}",
        output_prefix=lambda wc: f"{results_dir}/ref_panels_filtered/{wc.panel}/{wc.panel}_chr{wc.chrom}"
    resources:
        mem_mb=64000,
        cores=4,
        time="01:00:00"
    log:
        f"{results_dir}/log/filter_ref_{{panel}}_chr{{chrom}}.log"
    shell:
        """
        mkdir -p $(dirname {output.bed})
        bash scripts/filter_ref_panels_wrapper.sh \
            {input.snp_list} \
            {params.ref_panel_path} \
            {params.chrom} \
            {params.output_prefix} \
            2>&1 | tee {log}
        touch {output.done}
        """
