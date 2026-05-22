#!/usr/bin/env Rscript
# plot_region_manhattan.R
# Plot a single genomic region using the original extracted GWAS table.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
    stop("Usage: Rscript plot_region_manhattan.R <input_file> <output_png> <region_name> <region_chr> <region_start> <region_end> <dataset> <table_name> <columns_csv> [scaling_tsv]")
}

input_file <- args[1]
output_png <- args[2]
region_name <- args[3]
region_chr <- args[4]
region_start <- as.numeric(args[5])
region_end <- as.numeric(args[6])
dataset <- args[7]
table_name <- args[8]
columns_csv <- if (length(args) >= 9) args[9] else ""
scaling_tsv <- if (length(args) >= 10) args[10] else ""

# Plot tuning requested by user
threshold_suggestive <- 5e-5
threshold_genomewide <- 5e-8
global_ymax <- NA_real_

if (nzchar(scaling_tsv) && file.exists(scaling_tsv)) {
    scaling <- tryCatch(
        data.table::fread(scaling_tsv, sep = "\t", data.table = FALSE),
        error = function(e) NULL
    )
    if (!is.null(scaling) && all(c("key", "value") %in% colnames(scaling))) {
        val <- scaling$value[scaling$key == "global_ymax"][1]
        parsed <- suppressWarnings(as.numeric(val))
        if (!is.na(parsed) && is.finite(parsed) && parsed > 0) {
            global_ymax <- parsed
            cat("Using global y-axis max from scaling file:", global_ymax, "\n")
        }
    }
}

table_columns <- unique(trimws(unlist(strsplit(columns_csv, ",", fixed = TRUE))))

suppressWarnings({
    library(data.table, quietly = TRUE)
    library(ggplot2, quietly = TRUE)
    library(dplyr, quietly = TRUE)
})

norm_chr <- function(x) {
    x <- as.character(x)
    gsub("^chr", "", trimws(x), ignore.case = TRUE)
}

cat("Reading GWAS data from:", input_file, "\n")
gwas_data <- as.data.frame(fread(input_file, sep = "\t", data.table = FALSE))
cat("Data loaded:", nrow(gwas_data), "rows x", ncol(gwas_data), "columns\n")

p_cols <- grep("p_|pval|^p$|pvalue", colnames(gwas_data), ignore.case = TRUE)
# Exclude coordinate/allele columns
p_cols <- p_cols[!tolower(colnames(gwas_data)[p_cols]) %in% c("pos", "position", "chrom", "chr", "bp", "allele")]
if (length(table_columns) > 0) {
    preferred_p <- which(tolower(colnames(gwas_data)) %in% tolower(table_columns) &
                         grepl("p_|pval|^p$|pvalue", colnames(gwas_data), ignore.case = TRUE) &
                         !tolower(colnames(gwas_data)) %in% c("pos", "position", "chrom", "chr", "bp", "allele"))
    if (length(preferred_p) > 0) {
        p_cols <- unique(c(preferred_p, p_cols))
    }
}
if (length(p_cols) == 0) {
    stop(paste0("Could not find p-value column in data. Available columns: ", paste(colnames(gwas_data), collapse = ", ")))
}
p_col <- colnames(gwas_data)[p_cols[1]]

chr_cols <- grep("chrom|^chr$", colnames(gwas_data), ignore.case = TRUE)
if (length(table_columns) > 0) {
    preferred_chr <- which(tolower(colnames(gwas_data)) %in% tolower(table_columns) &
                           grepl("chrom|^chr$", colnames(gwas_data), ignore.case = TRUE))
    if (length(preferred_chr) > 0) {
        chr_cols <- unique(c(preferred_chr, chr_cols))
    }
}
if (length(chr_cols) == 0) {
    stop("Could not find chromosome column in data.")
}
chr_col <- colnames(gwas_data)[chr_cols[1]]

pos_cols <- grep("^pos$|^position$|^bp$", colnames(gwas_data), ignore.case = TRUE)
if (length(table_columns) > 0) {
    preferred_pos <- which(tolower(colnames(gwas_data)) %in% tolower(table_columns) &
                           grepl("^pos$|^position$|^bp$", colnames(gwas_data), ignore.case = TRUE))
    if (length(preferred_pos) > 0) {
        pos_cols <- unique(c(preferred_pos, pos_cols))
    }
}
if (length(pos_cols) == 0) {
    stop("Could not find position column in data.")
}
pos_col <- colnames(gwas_data)[pos_cols[1]]

# Try to find Z-score column
z_cols <- grep("z_stat|zscore|z_score|^z$", colnames(gwas_data), ignore.case = TRUE)
if (length(table_columns) > 0) {
    preferred_z <- which(tolower(colnames(gwas_data)) %in% tolower(table_columns) &
                         grepl("z_stat|zscore|z_score|^z$", colnames(gwas_data), ignore.case = TRUE))
    if (length(preferred_z) > 0) {
        z_cols <- unique(c(preferred_z, z_cols))
    }
}
z_col <- if (length(z_cols) > 0) colnames(gwas_data)[z_cols[1]] else NA_character_
beta_col <- NA_character_
se_col <- NA_character_
z_from_beta_se <- FALSE
z_from_p <- FALSE

# If explicit Z-score is not available, estimate it from beta/se when possible.
if (is.na(z_col)) {
    beta_cols <- grep("^beta$|beta_", colnames(gwas_data), ignore.case = TRUE)
    se_cols <- grep("^se$|se_", colnames(gwas_data), ignore.case = TRUE)
    if (length(beta_cols) > 0 && length(se_cols) > 0) {
        beta_col <- colnames(gwas_data)[beta_cols[1]]
        se_col <- colnames(gwas_data)[se_cols[1]]
        z_from_beta_se <- TRUE
        cat("Using estimated Z-score from", beta_col, "and", se_col, "\n")
    }
}

# Last-resort fallback: unsigned z from p-value.
if (is.na(z_col) && !z_from_beta_se) {
    z_from_p <- TRUE
    cat("Using Z-score estimated from p-value only (unsigned)\n")
}

cat("Using columns - CHR:", chr_col, "POS:", pos_col, "PVAL:", p_col)
if (!is.na(z_col)) {
    cat(" ZSCORE:", z_col)
} else if (z_from_beta_se) {
    cat(" ZSCORE: estimated from", beta_col, "and", se_col)
} else if (z_from_p) {
    cat(" ZSCORE: estimated from p-value")
}
cat("\n")
cat("Data sample (first 3 rows of used columns):\n")
sample_cols <- c(chr_col, pos_col, p_col)
if (!is.na(z_col)) sample_cols <- c(sample_cols, z_col)
if (z_from_beta_se) sample_cols <- c(sample_cols, beta_col, se_col)
print(head(gwas_data[, sample_cols], 3))

select_cols <- c(chr_col, pos_col, p_col)
if (!is.na(z_col)) select_cols <- c(select_cols, z_col)
if (z_from_beta_se) select_cols <- c(select_cols, beta_col, se_col)
select_cols <- unique(select_cols)

plot_data <- gwas_data %>%
    select(all_of(select_cols)) %>%
    mutate(
        CHR = norm_chr(.data[[chr_col]]),
        POS = suppressWarnings(as.numeric(.data[[pos_col]]))
    ) %>%
    filter(CHR == norm_chr(region_chr), POS >= region_start, POS <= region_end) %>%
    mutate(PVAL = suppressWarnings(as.numeric(.data[[p_col]])))

if (!is.na(z_col)) {
    plot_data <- plot_data %>%
        mutate(ZSCORE = suppressWarnings(as.numeric(.data[[z_col]])))
} else if (z_from_beta_se) {
    plot_data <- plot_data %>%
        mutate(ZSCORE = suppressWarnings(as.numeric(.data[[beta_col]]) / as.numeric(.data[[se_col]])))
} else if (z_from_p) {
    plot_data <- plot_data %>%
        mutate(ZSCORE = suppressWarnings(qnorm(1 - PVAL / 2)))
} else {
    plot_data <- plot_data %>%
        mutate(ZSCORE = NA_real_)
}

plot_data <- plot_data %>%
    filter(!is.na(CHR), !is.na(POS), !is.na(PVAL), !is.na(ZSCORE), PVAL > 0, PVAL <= 1) %>%
    arrange(POS) %>%
    mutate(log10p = -log10(PVAL))

cat("After region filter (chr", region_chr, ":", region_start, "-", region_end, "):", nrow(plot_data), "variants\n")

lambda_label <- "Lambda = NA"
if (nrow(plot_data) >= 3) {
    tryCatch({
        chisq_obs <- qchisq(1 - plot_data$PVAL, 1)
        lambda <- median(chisq_obs, na.rm = TRUE) / qchisq(0.5, 1)
        if (!is.na(lambda) && is.finite(lambda)) {
            lambda_label <- paste0("Lambda = ", round(lambda, 3))
        }
    }, error = function(e) {
        cat("Warning: Could not calculate lambda:", e$message, "\n")
    })
}

dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)

if (nrow(plot_data) == 0) {
    p <- ggplot() +
        annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = paste0("No variants found in region\n", region_name, "\n", dataset, " - ", table_name),
            size = 5
        ) +
        xlim(0, 1) +
        ylim(0, 1) +
        theme_void() +
        labs(title = paste0("Region Manhattan: ", region_name, " | ", dataset, " - ", table_name))
} else {
    region_width <- max(1, region_end - region_start)
    x_padding <- region_width * 0.05
    x_start <- max(0, region_start - x_padding)
    x_end <- region_end + x_padding
    y_max <- if (!is.na(global_ymax)) global_ymax else max(plot_data$log10p, na.rm = TRUE) * 1.10

    p <- ggplot(plot_data, aes(x = POS, y = log10p, color = ZSCORE)) +
        geom_point(aes(size = log10p), alpha = 0.7) +
        scale_color_gradient2(
            low = "#4575b4",
            mid = "#fee08b",
            high = "#d73027",
            midpoint = 0,
            name = "Z-score"
        ) +
        geom_hline(yintercept = -log10(threshold_suggestive), linetype = "dashed", color = "#d6820d", linewidth = 1.2, alpha = 0.9) +
        geom_hline(yintercept = -log10(threshold_genomewide), linetype = "dashed", color = "#b30000", linewidth = 1.2, alpha = 0.9) +
        annotate("text", x = Inf, y = -log10(threshold_suggestive), label = as.character(threshold_suggestive),
                 vjust = -0.3, hjust = 1.1, size = 3, color = "#d6820d") +
        annotate("text", x = Inf, y = -log10(threshold_genomewide), label = as.character(threshold_genomewide),
                 vjust = -0.3, hjust = 1.1, size = 3, color = "#b30000") +
        scale_size_continuous(
            name = "-log10(p-value)",
            range = c(0.7, 3.8),
            breaks = seq(0, ceiling(max(plot_data$log10p, na.rm = TRUE)), by = 2)
        ) +
        coord_cartesian(xlim = c(x_start, x_end), ylim = c(0, y_max), expand = FALSE) +
        labs(
            title = paste0("Region Manhattan: ", region_name, " | ", dataset, " - ", table_name),
            subtitle = paste0(norm_chr(region_chr), ":", format(region_start, scientific = FALSE, trim = TRUE),
                              "-", format(region_end, scientific = FALSE, trim = TRUE),
                              " | ", nrow(plot_data), " variants | ", lambda_label),
            x = paste0("Position on chr", norm_chr(region_chr)),
            y = "-log10(p-value)"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(size = 9),
            plot.title = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 10),
            legend.position = "right",
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 9)
        )
}

cat("Saving PNG to:", output_png, "\n")
ggsave(output_png, plot = p, width = 18, height = 8, dpi = 300, device = "png")
cat("Region Manhattan plot creation complete.\n")
