import sys
sys.path.append("utils")
from bioconfigme import get_results_dir

configfile: "configs/analysis.yml"

results_dir = get_results_dir()

# Build dataset list from config: combine primary and ancestral GWAS files
GWAS_DATASETS = []
if "classification" in config and "gwas_files" in config["classification"]:
    gwas_files = config["classification"]["gwas_files"]
    GWAS_DATASETS = gwas_files.get("primary", []) + gwas_files.get("ancestral", [])
else:
    # Fallback to defaults if config section missing
    GWAS_DATASETS = ["Hisp", "WangEAS", "BJEUR", "AFR_AMR", "AFR_EAS", "AFR_EUR", "AMR_EAS", "AMR_EUR", "EUR_EAS"]


def get_extracted_gwas_inputs():
    """Return Step 00 extracted GWAS TSVs for all classification datasets."""
    paths = []
    gwastables = config.get("gwastables", [])
    
    for dataset in GWAS_DATASETS:
        entry = next((item for item in gwastables if item.get("name") == dataset), None)
        if not entry:
            continue
        for table in entry.get("tables", []):
            table_name = table.get("table_name")
            if table_name and "gwas_primary" in table_name or "dosage" in table_name:
                paths.append(f"{results_dir}/extracted/{dataset}/{table_name}.tsv")
    return paths


rule all:
    input:
        f"{results_dir}/gwas_snps_merged/all_gwas_snps_significant.tsv",
        f"{results_dir}/gwas_snps_merged/all_gwas_snps_significant.xlsx"


rule merge_all_gwas_snps:
    input:
        extracted=get_extracted_gwas_inputs(),
        loci_summary=f"{results_dir}/loci_summary/dataset_loci_gene_summary.tsv"
    output:
        tsv=f"{results_dir}/gwas_snps_merged/all_gwas_snps_significant.tsv",
        xlsx=f"{results_dir}/gwas_snps_merged/all_gwas_snps_significant.xlsx"
    params:
        datasets=",".join(GWAS_DATASETS),
        p_threshold=5e-5
    resources:
        mem_mb=64000,
        cores=4,
        time="01:30:00"
    log:
        f"{results_dir}/log/merge_gwas_snps.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        python3 scripts/merge_gwas_snps.py \
            --extracted {input.extracted} \
            --loci-summary {input.loci_summary} \
            --datasets {params.datasets} \
            --p-threshold {params.p_threshold} \
            --config configs/analysis.yml \
            --output-tsv {output.tsv} \
            --output-xlsx {output.xlsx} \
            2>&1 | tee {log}
        """
