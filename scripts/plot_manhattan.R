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

cat("Using columns - CHR:", chr_col, "POS:", pos_col, "PVAL:", p_col, "\n")

# Prepare data for plotting
plot_data <- gwas_data %>%
    select(all_of(c(chr_col, pos_col, p_col))) %>%
    filter(!is.na(!!sym(p_col)), !is.na(!!sym(chr_col)), !is.na(!!sym(pos_col))) %>%
    mutate(
        CHR = as.numeric(!!sym(chr_col)),
        POS = as.numeric(!!sym(pos_col)),
        PVAL = as.numeric(!!sym(p_col))
    ) %>%
    filter(!is.na(CHR), !is.na(POS), !is.na(PVAL), PVAL > 0, PVAL <= 1) %>%
    select(CHR, POS, PVAL) %>%
    mutate(
        log10p = -log10(PVAL),
        CHR_factor = factor(CHR, levels = sort(unique(CHR)))
    ) %>%
    arrange(CHR, POS)

cat("Data prepared:", nrow(plot_data), "variants\n")

# Calculate lambda (genomic inflation factor)
# Lambda = median(observed χ²) / median(expected χ²)
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

inv_transform_logp <- function(y_zoom, low = 1, high = 12, factor = 5) {
    zoom_high <- low + (high - low) * factor
    ifelse(
        y_zoom < low,
        y_zoom,
        ifelse(y_zoom <= zoom_high,
               low + (y_zoom - low) / factor,
               high + (y_zoom - zoom_high))
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
p <- ggplot(plot_data, aes(x = POS, y = log10p_zoom, color = CHR_factor)) +
    geom_point(size = 1.2, alpha = 0.6) +
    facet_wrap(. ~ CHR_factor, nrow = 1, scales = "free_x") +
    geom_hline(yintercept = transform_logp(-log10(5e-5), low = 1, high = 12, factor = 5), linetype = "dashed", color = "blue", size = 0.8, alpha = 0.7) +
    geom_hline(yintercept = transform_logp(-log10(5e-8), low = 1, high = 12, factor = 5), linetype = "dashed", color = "red", size = 0.8, alpha = 0.7) +
    scale_color_manual(values = rep(c("gray75", "black"), length.out = nlevels(plot_data$CHR_factor))) +
    scale_y_continuous(breaks = y_breaks_zoom, labels = y_breaks_raw) +
    labs(
        title = paste0("Manhattan Plot: ", dataset, " - ", table_name),
        subtitle = paste0(lambda_label, " | y-axis: 1-12 expanded 5x"),
        x = NULL,
        y = "-log10(p-value)"
    ) +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10)
    )

# Save plot
cat("Saving PNG to:", output_png, "\n")
ggsave(output_png, plot = p, width = 16, height = 6, dpi = 300, device = "png")

cat("Manhattan plot creation complete.\n")
