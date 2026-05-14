#!/bin/bash
# steps.sh - Pipeline execution commands for GWAS handler

# Step 1: Extract GWAS tables for each dataset
./submit.sh --snakefile rules/00_extract_gwas.smk --jobs 10 --cores 8

# Step 2: Plot GWAS Manhattan plots for all tables
./submit.sh --snakefile rules/01_plot_gwas_manhattan.smk --jobs 20 --cores 8

# Step 3: Subset GWAS by p-value threshold
./submit.sh --snakefile rules/03_subset_gwas.smk --jobs 10 --cores 8

# Step 4: Aggregate SNPs and filter reference panels
./submit.sh --snakefile rules/04_filter_ref_panels.smk --jobs 10 --cores 8

# Step 5: Build merged FUMA loci and Excel min-p summary
./submit.sh --snakefile rules/02_fuma_loci_summary.smk --jobs 1 --cores 2
