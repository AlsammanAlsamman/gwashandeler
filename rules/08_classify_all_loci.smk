import sys
sys.path.append("utils")
from bioconfigme import get_results_dir

configfile: "configs/analysis.yml"

results_dir = get_results_dir()

# Build dataset list from config: combine primary and ancestral GWAS files
CLASSIFY_DATASETS = []
if "classification" in config and "gwas_files" in config["classification"]:
    gwas_files = config["classification"]["gwas_files"]
    CLASSIFY_DATASETS = gwas_files.get("primary", []) + gwas_files.get("ancestral", [])
else:
    # Fallback to defaults if config section missing
    CLASSIFY_DATASETS = ["Hisp", "AFR_AMR", "AFR_EAS", "AFR_EUR", "AMR_EAS", "AMR_EUR", "EUR_EAS", "AFR", "AMR", "EAS", "EUR"]


def get_step6_loci_inputs():
    return [f"{results_dir}/loci/{dataset}/{dataset}_loci.tsv" for dataset in CLASSIFY_DATASETS]


def get_primary_extracted_inputs():
    """Return Step 00 extracted GWAS TSVs for all datasets used in classification."""
    paths = []
    gwastables = config.get("gwastables", [])

    for dataset in CLASSIFY_DATASETS:
        entry = next((item for item in gwastables if item.get("name") == dataset), None)
        if not entry:
            continue
        for table in entry.get("tables", []):
            table_name = table.get("table_name")
            if table_name:
                paths.append(f"{results_dir}/extracted/{dataset}/{table_name}.tsv")
    return paths


rule all:
    input:
        f"{results_dir}/loci_classification/all_loci_classified.tsv",
        f"{results_dir}/loci_classification/all_loci_classified.xlsx"


rule classify_all_loci:
    input:
        loci=get_step6_loci_inputs(),
        extracted=get_primary_extracted_inputs()
    output:
        tsv=f"{results_dir}/loci_classification/all_loci_classified.tsv",
        xlsx=f"{results_dir}/loci_classification/all_loci_classified.xlsx"
    resources:
        mem_mb=32000,
        cores=1,
        time="00:45:00"
    log:
        f"{results_dir}/log/classify_all_loci.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        python3 scripts/classify_all_loci.py \
            --inputs {input.loci} \
            --extracted {input.extracted} \
            --config configs/analysis.yml \
            --output-tsv {output.tsv} \
            --output-xlsx {output.xlsx} \
            2>&1 | tee {log}
        """