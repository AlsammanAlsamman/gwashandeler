#!/usr/bin/env Rscript
# Quick test of column detection logic

library(data.table)
library(dplyr)

# Test with AMR_EUR dataset
input_file <- "input/AMR_EUR.tsv"
table_columns_str <- "CHROM,POS,ID,REF,ALT,n,beta_local_ancestry,se_local_ancestry,p_local_ancestry,ancestry_i,ancestry_j"

cat("Testing column detection with:", input_file, "\n")
cat("Expected table columns:", table_columns_str, "\n\n")

gwas_data <- as.data.frame(fread(input_file, sep = "\t", data.table = FALSE))
cat("Data loaded:", nrow(gwas_data), "rows x", ncol(gwas_data), "columns\n")
cat("Column names:", paste(colnames(gwas_data), collapse = ", "), "\n\n")

table_columns <- unique(trimws(unlist(strsplit(table_columns_str, ",", fixed = TRUE))))

# Test p-value detection
cat("--- Testing P-value column detection ---\n")
p_cols <- grep("p_|pval|^p$|pvalue", colnames(gwas_data), ignore.case = TRUE)
cat("Initial grep matches:", paste(colnames(gwas_data)[p_cols], collapse = ", "), "\n")

p_cols <- p_cols[!tolower(colnames(gwas_data)[p_cols]) %in% c("pos", "position", "chrom", "chr", "bp", "allele")]
cat("After exclusion:", paste(colnames(gwas_data)[p_cols], collapse = ", "), "\n")

if (length(table_columns) > 0) {
    preferred_p <- which(tolower(colnames(gwas_data)) %in% tolower(table_columns) &
                         grepl("p_|pval|^p$|pvalue", colnames(gwas_data), ignore.case = TRUE) &
                         !tolower(colnames(gwas_data)) %in% c("pos", "position", "chrom", "chr", "bp", "allele"))
    if (length(preferred_p) > 0) {
        cat("Preferred matches from table_columns:", paste(colnames(gwas_data)[preferred_p], collapse = ", "), "\n")
        p_cols <- unique(c(preferred_p, p_cols))
    }
}

p_col <- colnames(gwas_data)[p_cols[1]]
cat("SELECTED P-VALUE COLUMN:", p_col, "\n\n")

# Test CHR detection
cat("--- Testing CHR column detection ---\n")
chr_cols <- grep("chrom|^chr$", colnames(gwas_data), ignore.case = TRUE)
cat("Grep matches:", paste(colnames(gwas_data)[chr_cols], collapse = ", "), "\n")
chr_col <- colnames(gwas_data)[chr_cols[1]]
cat("SELECTED CHR COLUMN:", chr_col, "\n\n")

# Test POS detection
cat("--- Testing POS column detection ---\n")
pos_cols <- grep("^pos$|^position$|^bp$", colnames(gwas_data), ignore.case = TRUE)
cat("Grep matches:", paste(colnames(gwas_data)[pos_cols], collapse = ", "), "\n")
pos_col <- colnames(gwas_data)[pos_cols[1]]
cat("SELECTED POS COLUMN:", pos_col, "\n\n")

# Show sample data
cat("--- Sample data (first 3 rows) ---\n")
print(head(gwas_data[, c(chr_col, pos_col, p_col)], 3))

# Test region filtering
region_chr <- "3"
region_start <- 169508915
region_end <- 169518455

norm_chr <- function(x) {
    x <- as.character(x)
    gsub("^chr", "", trimws(x), ignore.case = TRUE)
}

plot_data <- gwas_data %>%
    select(all_of(c(chr_col, pos_col, p_col))) %>%
    mutate(
        CHR = norm_chr(.data[[chr_col]]),
        POS = suppressWarnings(as.numeric(.data[[pos_col]])),
        PVAL = suppressWarnings(as.numeric(.data[[p_col]]))
    ) %>%
    filter(!is.na(CHR), !is.na(POS), !is.na(PVAL), PVAL > 0, PVAL <= 1) %>%
    filter(CHR == norm_chr(region_chr), POS >= region_start, POS <= region_end)

cat("\n--- Region filtering test ---\n")
cat("Looking for chr", region_chr, ":", region_start, "-", region_end, "\n")
cat("Found:", nrow(plot_data), "variants\n")
if (nrow(plot_data) > 0) {
    cat("\nFirst few variants in region:\n")
    print(head(plot_data[, c("CHR", "POS", "PVAL")], 5))
}
