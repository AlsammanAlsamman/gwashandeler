#!/bin/bash
# steps.sh - Pipeline execution commands for GWAS handler

# Step 1: Extract GWAS tables for each dataset
./submit.sh --snakefile rules/00_extract_gwas.smk --jobs 20 --cores 8

# Step 2: Plot GWAS Manhattan plots for all tables
./submit.sh --snakefile rules/01_plot_gwas_manhattan.smk --jobs 20 --cores 8

# Step 2b: Plot region Manhattan plots from the original extracted tables
./submit.sh --snakefile rules/01_plot_regions_manhattan.smk --jobs 40 --cores 8

# Step 3: Subset GWAS by p-value threshold
./submit.sh --snakefile rules/03_subset_gwas.smk --jobs 20 --cores 8

# Step 4: Aggregate SNPs and filter reference panels
./submit.sh --snakefile rules/04_filter_ref_panels.smk --jobs 20 --cores 8

# Step 5: Identify loci per chromosome and merge per table
./submit.sh --snakefile rules/05_identify_loci.smk --jobs 40 --cores 8

# Step 6: Merge all tables into one loci table per dataset
./submit.sh --snakefile rules/06_merge_dataset_loci.smk --jobs 10 --cores 4

# Step 7: Summarize loci by chromosome (22 parallel jobs) then combine
./submit.sh --snakefile rules/07_loci_dataset_summary.smk --jobs 22 --cores 8

# Step 8: Classify merged loci across datasets and Hispanic GWAS
./submit.sh --snakefile rules/08_classify_all_loci.smk --jobs 1 --cores 2

# Step 9: Merge all GWAS SNPs (p < 5E-5) with loci and gene annotations
./submit.sh --snakefile rules/09_merge_gwas_snps.smk --jobs 1 --cores 4

# Step 10: Build merged FUMA loci and Excel min-p summary
./submit.sh --snakefile rules/02_fuma_loci_summary.smk --jobs 1 --cores 2
