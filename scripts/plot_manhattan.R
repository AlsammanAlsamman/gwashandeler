#!/usr/bin/env Rscript
# plot_manhattan.R
# Creates standard Manhattan plots with significance threshold lines and lambda calculation.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: Rscript plot_manhattan.R <input_file> <output_png> <dataset> <table_name>")
}

input_file <- args[1]
output_png <- args[2]
dataset <- args[3]
table_name <- args[4]

# Plot tuning requested by user
threshold_suggestive <- 5e-5
threshold_genomewide <- 5e-8

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
        CHR_RAW = as.character(!!sym(chr_col)),
        POS = as.numeric(!!sym(pos_col)),
        PVAL = as.numeric(!!sym(p_col))
    ) %>%
    mutate(
        CHR = suppressWarnings(as.numeric(gsub("^chr", "", CHR_RAW, ignore.case = TRUE)))
    ) %>%
    filter(!is.na(CHR), !is.na(POS), !is.na(PVAL), PVAL > 0, PVAL <= 1) %>%
    select(CHR, POS, PVAL) %>%
    mutate(
        log10p = -log10(PVAL),
        CHR_factor = factor(CHR, levels = sort(unique(CHR)))
    ) %>%
    arrange(CHR, POS)

cat("Data prepared:", nrow(plot_data), "variants\n")

lambda <- NA_real_
if (nrow(plot_data) >= 3) {
    # Lambda = median(observed chi-square) / expected median chi-square (0.5 quantile)
    tryCatch({
        chisq_obs <- qchisq(1 - plot_data$PVAL, 1)
        lambda <- median(chisq_obs, na.rm = TRUE) / qchisq(0.5, 1)
    }, error = function(e) {
        cat("Warning: Could not calculate lambda:", e$message, "\n")
        lambda <<- NA_real_
    })
}

if (is.na(lambda)) {
    lambda_label <- "Lambda = NA"
} else {
    lambda_label <- paste0("Lambda = ", round(lambda, 3))
}

cat("Lambda:", lambda_label, "\n")

if (nrow(plot_data) == 0) {
    p <- ggplot() +
        annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = paste0("No variants available after filtering\n", dataset, " - ", table_name),
            size = 5
        ) +
        xlim(0, 1) +
        ylim(0, 1) +
        theme_void() +
        labs(title = paste0("Manhattan Plot: ", dataset, " - ", table_name))
} else {
    chr_info <- plot_data %>%
        group_by(CHR) %>%
        summarise(chr_len = max(POS, na.rm = TRUE), .groups = "drop") %>%
        arrange(CHR) %>%
        mutate(tot = cumsum(chr_len) - chr_len)

    plot_data <- plot_data %>%
        left_join(chr_info, by = "CHR") %>%
        mutate(BPcum = POS + tot,
               chr_index = as.integer(factor(CHR, levels = chr_info$CHR)),
               CHR_COLOR = ifelse(chr_index %% 2 == 0, "gray50", "black"))

    axis_df <- plot_data %>%
        group_by(CHR) %>%
        summarise(center = (min(BPcum) + max(BPcum)) / 2, .groups = "drop") %>%
        arrange(CHR)

    y_max <- max(plot_data$log10p, na.rm = TRUE)

    p <- ggplot(plot_data, aes(x = BPcum, y = log10p, color = CHR_COLOR)) +
        geom_point(size = 0.8) +
        geom_hline(yintercept = -log10(threshold_suggestive), linetype = "dashed", color = "#d6820d", linewidth = 1.1) +
        geom_hline(yintercept = -log10(threshold_genomewide), linetype = "dashed", color = "#b30000", linewidth = 1.1) +
        scale_color_identity() +
        scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +
        coord_cartesian(ylim = c(0, y_max * 1.05), expand = FALSE) +
        labs(
            title = paste0("Manhattan Plot: ", dataset, " - ", table_name),
            subtitle = paste0(lambda_label, " | thresholds: ", threshold_suggestive, " and ", threshold_genomewide),
            x = "Chromosome",
            y = "-log10(p-value)"
        ) +
        theme_minimal() +
        theme(
            legend.position = "none",
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(size = 9),
            plot.title = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 10)
        )
}

# Save plot
cat("Saving PNG to:", output_png, "\n")
ggsave(output_png, plot = p, width = 24, height = 9, dpi = 300, device = "png")

cat("Manhattan plot creation complete.\n")
