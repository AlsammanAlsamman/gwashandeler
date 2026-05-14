#!/usr/bin/env Rscript
# plot_manhattan.R
# Creates Manhattan plots with significance threshold lines and lambda calculation

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: Rscript plot_manhattan.R <input_file> <output_png> <dataset> <table_name>")
}

input_file <- args[1]
output_png <- args[2]
dataset <- args[3]
table_name <- args[4]

# Plot tuning requested by user
threshold_suggestive <- 5e-4
threshold_genomewide <- 5e-6

# Load required libraries
suppressWarnings({
    library(data.table, quietly = TRUE)
    library(ggplot2, quietly = TRUE)
    library(dplyr, quietly = TRUE)
})

# Read the data using fread (fast reading, equivalent to pandas read_csv)
cat("Reading GWAS data from:", input_file, "\n")
gwas_data <- as.data.frame(fread(input_file, sep = "\t", data.table = FALSE))

cat("Data loaded:", nrow(gwas_data), "rows x", ncol(gwas_data), "columns\n")

# Find p-value column (try common column names)
p_cols <- grep("^p_|^pval|^P$|^p$", colnames(gwas_data), ignore.case = TRUE)
if (length(p_cols) == 0) {
    stop("Could not find p-value column in data. Available columns:", paste(colnames(gwas_data), collapse = ", "))
}
p_col <- colnames(gwas_data)[p_cols[1]]

# Find chromosome column
chr_cols <- grep("^chrom|^chr|^CHROM|^CHR", colnames(gwas_data), ignore.case = TRUE)
if (length(chr_cols) == 0) {
    stop("Could not find chromosome column in data.")
}
chr_col <- colnames(gwas_data)[chr_cols[1]]

# Find position column
pos_cols <- grep("^pos|^position|^bp|^POS|^BP", colnames(gwas_data), ignore.case = TRUE)
if (length(pos_cols) == 0) {
    stop("Could not find position column in data.")
}
pos_col <- colnames(gwas_data)[pos_cols[1]]

# Find Z-score column for color mapping
z_cols <- grep("z_stat|zscore|z_score|^z$", colnames(gwas_data), ignore.case = TRUE)
z_col <- if (length(z_cols) > 0) colnames(gwas_data)[z_cols[1]] else NA_character_

# If explicit Z-score is not available, estimate it from beta/se when possible.
if (is.na(z_col)) {
    beta_cols <- grep("^beta$|beta_", colnames(gwas_data), ignore.case = TRUE)
    se_cols <- grep("^se$|se_", colnames(gwas_data), ignore.case = TRUE)
    if (length(beta_cols) > 0 && length(se_cols) > 0) {
        beta_col <- colnames(gwas_data)[beta_cols[1]]
        se_col <- colnames(gwas_data)[se_cols[1]]
        gwas_data$ZSCORE_EST <- suppressWarnings(as.numeric(gwas_data[[beta_col]]) / as.numeric(gwas_data[[se_col]]))
        z_col <- "ZSCORE_EST"
        cat("Using estimated Z-score from", beta_col, "and", se_col, "\n")
    }
}

# Last-resort fallback: unsigned z from p-value.
if (is.na(z_col)) {
    gwas_data$ZSCORE_EST <- suppressWarnings(qnorm(1 - as.numeric(gwas_data[[p_col]]) / 2))
    z_col <- "ZSCORE_EST"
    cat("Using Z-score estimated from p-value only (unsigned)\n")
}

cat("Using columns - CHR:", chr_col, "POS:", pos_col, "PVAL:", p_col, "ZSCORE:", z_col, "\n")

# Prepare data for plotting
plot_data <- gwas_data %>%
    select(all_of(c(chr_col, pos_col, p_col, z_col))) %>%
    filter(!is.na(!!sym(p_col)), !is.na(!!sym(chr_col)), !is.na(!!sym(pos_col)), !is.na(!!sym(z_col))) %>%
    mutate(
        CHR = as.numeric(!!sym(chr_col)),
        POS = as.numeric(!!sym(pos_col)),
        PVAL = as.numeric(!!sym(p_col)),
        ZSCORE = as.numeric(!!sym(z_col))
    ) %>%
    filter(!is.na(CHR), !is.na(POS), !is.na(PVAL), !is.na(ZSCORE), PVAL > 0, PVAL <= 1) %>%
    select(CHR, POS, PVAL, ZSCORE) %>%
    mutate(
        log10p = -log10(PVAL),
        CHR_factor = factor(CHR, levels = sort(unique(CHR)))
    ) %>%
    arrange(CHR, POS)

cat("Data prepared:", nrow(plot_data), "variants\n")

# Calculate lambda (genomic inflation factor)
# Lambda = median(observed chi-square) / median(expected chi-square)
tryCatch({
    chisq_obs <- qchisq(1 - plot_data$PVAL, 1)
    chisq_exp <- qchisq(seq(1, nrow(plot_data)) / (nrow(plot_data) + 1), 1)
    lambda <- median(chisq_obs) / median(chisq_exp)
}, error = function(e) {
    cat("Warning: Could not calculate lambda:", e$message, "\n")
    lambda <<- NA
})

if (is.na(lambda)) {
    lambda_label <- "Lambda = NA"
} else {
    lambda_label <- paste0("Lambda = ", round(lambda, 3))
}

cat("Lambda:", lambda_label, "\n")

# Piecewise y transform to expand the 1-12 region by 5x.
transform_logp <- function(y, low = 1, high = 12, factor = 5) {
    ifelse(
        y < low,
        y,
        ifelse(y <= high,
               low + (y - low) * factor,
               low + (high - low) * factor + (y - high))
    )
}

plot_data <- plot_data %>%
    mutate(log10p_zoom = transform_logp(log10p, low = 1, high = 12, factor = 5))

y_max_raw <- max(plot_data$log10p, na.rm = TRUE)
base_breaks <- seq(0, min(12, ceiling(y_max_raw)), by = 1)
if (ceiling(y_max_raw) >= 15) {
    high_breaks <- seq(15, ceiling(y_max_raw), by = 5)
    y_breaks_raw <- unique(c(base_breaks, high_breaks))
} else {
    y_breaks_raw <- base_breaks
}
y_breaks_zoom <- transform_logp(y_breaks_raw, low = 1, high = 12, factor = 5)

# Create Manhattan plot
p <- ggplot(plot_data, aes(x = POS, y = log10p_zoom, color = ZSCORE)) +
    geom_point(size = 0.8, alpha = 0.65) +
    facet_wrap(. ~ CHR_factor, nrow = 1, scales = "free_x") +
    geom_hline(yintercept = transform_logp(-log10(threshold_suggestive), low = 1, high = 12, factor = 5), linetype = "dashed", color = "#d6820d", linewidth = 1.2, alpha = 0.9) +
    geom_hline(yintercept = transform_logp(-log10(threshold_genomewide), low = 1, high = 12, factor = 5), linetype = "dashed", color = "#b30000", linewidth = 1.2, alpha = 0.9) +
    scale_color_gradient2(
        low = "#d73027",
        mid = "#f7f7f7",
        high = "#4575b4",
        midpoint = 0,
        name = "Z-score"
    ) +
    scale_y_continuous(breaks = y_breaks_zoom, labels = y_breaks_raw) +
    labs(
        title = paste0("Manhattan Plot: ", dataset, " - ", table_name),
        subtitle = paste0(lambda_label, " | y-axis: 1-12 expanded 5x | thresholds: ", threshold_suggestive, " and ", threshold_genomewide),
        x = NULL,
        y = "-log10(p-value)"
    ) +
    theme_minimal() +
    theme(
        legend.position = "right",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10)
    )

# Save plot
cat("Saving PNG to:", output_png, "\n")
ggsave(output_png, plot = p, width = 24, height = 9, dpi = 300, device = "png")

cat("Manhattan plot creation complete.\n")
