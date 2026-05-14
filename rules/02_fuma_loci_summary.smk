import sys
sys.path.append("utils")
from bioconfigme import get_results_dir

configfile: "configs/analysis.yml"

results_dir = get_results_dir()

rule all:
    input:
        f"{results_dir}/fuma/merged_loci_250kb.tsv",
        f"{results_dir}/fuma/merged_loci_minp_summary.xlsx",
        f"{results_dir}/fuma/merged_loci_minp_summary.done"


rule summarize_fuma_loci:
    output:
        loci=f"{results_dir}/fuma/merged_loci_250kb.tsv",
        xlsx=f"{results_dir}/fuma/merged_loci_minp_summary.xlsx",
        done=f"{results_dir}/fuma/merged_loci_minp_summary.done"
    params:
        config="configs/analysis.yml",
        distance=250000
    resources:
        mem_mb=128000,
        cores=2,
        time="02:00:00"
    log:
        f"{results_dir}/log/summarize_fuma_loci.log"
    shell:
        """
        mkdir -p $(dirname {output.loci})
        python3 scripts/summarize_fuma_loci.py \
            --config {params.config} \
            --output-xlsx {output.xlsx} \
            --output-loci {output.loci} \
            --distance {params.distance} \
            2>&1 | tee {log}
        touch {output.done}
        """
