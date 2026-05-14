import sys
sys.path.append("utils")
from bioconfigme import get_results_dir

configfile: "configs/analysis.yml"

results_dir = get_results_dir()


def get_dataset_names():
    return [d["name"] for d in config["gwastables"]]


def get_dataset_loci_done_targets():
    return [f"{results_dir}/loci/{d}/merged_dataset_loci.done" for d in get_dataset_names()]


def get_dataset_loci_tsv_inputs():
    return [f"{results_dir}/loci/{d}/{d}_loci.tsv" for d in get_dataset_names()]


def get_merge_distance_bp():
    return int(config.get("loci_identiffication", {}).get("merge_distance_kb", 250)) * 1000


def get_use_ld():
    # Support both nested and top-level config overrides.
    nested = config.get("loci_identiffication", {}).get("use_ld", None)
    if nested is not None:
        return bool(nested)
    return bool(config.get("use_ld", False))


def get_chromosomes():
    # Autosomes 1..22 only for the parallel chromosome summary stage.
    return [str(i) for i in range(1, 23)]


def get_chr_summary_tsv_targets():
    return [f"{results_dir}/loci_summary/by_chr/chr{c}.summary.tsv" for c in get_chromosomes()]


def get_chr_summary_done_targets():
    return [f"{results_dir}/loci_summary/by_chr/chr{c}.done" for c in get_chromosomes()]


rule all:
    input:
        f"{results_dir}/loci_summary/dataset_loci_gene_summary.tsv",
        f"{results_dir}/loci_summary/dataset_loci_gene_summary.xlsx",
        f"{results_dir}/loci_summary/dataset_loci_gene_summary.done"


rule summarize_dataset_loci_with_genes_by_chr:
    input:
        dataset_done=get_dataset_loci_done_targets(),
        dataset_loci=get_dataset_loci_tsv_inputs()
    output:
        tsv=f"{results_dir}/loci_summary/by_chr/chr{{chromosome}}.summary.tsv",
        done=f"{results_dir}/loci_summary/by_chr/chr{{chromosome}}.done"
    params:
        config="configs/analysis.yml",
        distance=get_merge_distance_bp(),
        use_ld=get_use_ld(),
        chromosome="{chromosome}"
    wildcard_constraints:
        chromosome="\d+"
    resources:
        mem_mb=16000,
        cores=2,
        time="01:00:00"
    log:
        f"{results_dir}/log/loci_dataset_summary_chr{{chromosome}}.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        LD_FLAG=""
        if [ "{params.use_ld}" == "True" ] || [ "{params.use_ld}" == "true" ]; then
            LD_FLAG="--use-ld"
        fi
        python3 scripts/summarize_dataset_loci_with_genes.py \
            --config {params.config} \
            --inputs {input.dataset_loci} \
            --chromosome {params.chromosome} \
            --output-tsv {output.tsv} \
            --distance {params.distance} \
            --skip-xlsx \
            $LD_FLAG \
            2>&1 | tee {log}
        touch {output.done}
        """


rule summarize_dataset_loci_with_genes:
    input:
        chr_tsv=get_chr_summary_tsv_targets(),
        chr_done=get_chr_summary_done_targets()
    output:
        tsv=f"{results_dir}/loci_summary/dataset_loci_gene_summary.tsv",
        xlsx=f"{results_dir}/loci_summary/dataset_loci_gene_summary.xlsx",
        done=f"{results_dir}/loci_summary/dataset_loci_gene_summary.done"
    resources:
        mem_mb=32000,
        cores=2,
        time="00:40:00"
    log:
        f"{results_dir}/log/loci_dataset_summary_combine.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        python3 scripts/combine_loci_summary_by_chr.py \
            --inputs {input.chr_tsv} \
            --output-tsv {output.tsv} \
            --output-xlsx {output.xlsx} \
            2>&1 | tee {log}
        touch {output.done}
        """
