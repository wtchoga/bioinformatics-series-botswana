################################################################################
# Organize significant GSEA Hallmark pathways into biological modules
# Input:  Week24_DESeq2_output/08_GSEA_Hallmark_Rebound_vs_NoRebound_Week24.csv
# Filter: nominal p value < 0.05
# Output: module-level and pathway-level CSV files + dot plots
################################################################################

## =============================================================================
## 0. SETUP
## =============================================================================

cran_packages <- c("dplyr", "ggplot2", "readr", "stringr", "tidyr", "forcats")
missing_pkgs <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
  library(forcats)
})

## =============================================================================
## 1. FILE PATHS
## =============================================================================

outdir <- "Week24_DESeq2_output"
gsea_path <- file.path(outdir, "08_GSEA_Hallmark_Rebound_vs_NoRebound_Week24.csv")
module_outdir <- file.path(outdir, "Pathway_Modules_pvalue_lt_0.05")
dir.create(module_outdir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(gsea_path)) {
  stop(
    "Could not find GSEA results file: ", gsea_path, "\n",
    "Run Week24_DESeq2_Rebound_vs_NoRebound_FINAL.R first, or update gsea_path."
  )
}

## =============================================================================
## 2. READ AND CHECK GSEA RESULTS
## =============================================================================

gsea_df <- read.csv(gsea_path, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("pathway", "pval", "padj", "NES")
missing_cols <- setdiff(required_cols, colnames(gsea_df))
if (length(missing_cols) > 0) {
  stop("The GSEA file is missing required column(s): ", paste(missing_cols, collapse = ", "))
}

# Make sure numeric columns are numeric even if the CSV was edited in Excel.
gsea_df <- gsea_df %>%
  mutate(
    pval = as.numeric(.data$pval),
    padj = as.numeric(.data$padj),
    NES = as.numeric(.data$NES)
  )

## =============================================================================
## 3. MODULE ASSIGNMENT FUNCTION
## =============================================================================

assign_hallmark_module <- function(pathway_name) {
  p <- toupper(pathway_name)

  case_when(
    str_detect(p, "INTERFERON|TNFA|NFKB|IL6|JAK|STAT|IL2|COMPLEMENT|INFLAMMATORY|ALLOGRAFT|COAGULATION") ~
      "Immune / inflammation",

    str_detect(p, "E2F|G2M|MITOTIC|MYC|DNA_REPAIR|SPERMATOGENESIS") ~
      "Cell cycle / proliferation",

    str_detect(p, "GLYCOLYSIS|OXIDATIVE_PHOSPHORYLATION|FATTY_ACID|CHOLESTEROL|BILE_ACID|ADIPOGENESIS|XENOBIOTIC|HEME|PEROXISOME|MTORC1|PI3K|AKT|MTOR") ~
      "Metabolism / growth signaling",

    str_detect(p, "HYPOXIA|APOPTOSIS|P53|UNFOLDED|REACTIVE_OXYGEN|UV_RESPONSE") ~
      "Stress response / cell death",

    str_detect(p, "EPITHELIAL_MESENCHYMAL|TGF_BETA|APICAL_JUNCTION|APICAL_SURFACE|ANGIOGENESIS|MYOGENESIS|KRAS") ~
      "Tissue remodeling / EMT",

    str_detect(p, "ESTROGEN|ANDROGEN|WNT|NOTCH|HEDGEHOG|HORMONE") ~
      "Developmental / hormone signaling",

    TRUE ~ "Other Hallmark programs"
  )
}

clean_pathway_label <- function(pathway_name) {
  pathway_name %>%
    str_remove("^HALLMARK_") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

## =============================================================================
## 4. FILTER SIGNIFICANT PATHWAYS AND ADD MODULES
## =============================================================================

sig_pathways <- gsea_df %>%
  filter(!is.na(.data$pval), .data$pval < 0.05, !is.na(.data$NES)) %>%
  mutate(
    Module = assign_hallmark_module(.data$pathway),
    Direction = case_when(
      .data$NES > 0 ~ "Enriched in Rebound",
      .data$NES < 0 ~ "Enriched in NoRebound",
      TRUE ~ "No directional enrichment"
    ),
    pathway_label = clean_pathway_label(.data$pathway),
    neg_log10_pvalue = -log10(.data$pval),
    abs_NES = abs(.data$NES)
  ) %>%
  arrange(.data$Module, .data$Direction, .data$pval)

if (nrow(sig_pathways) == 0) {
  warning("No pathways passed nominal p value < 0.05. Empty output tables will be written.")
}

write.csv(
  sig_pathways,
  file.path(module_outdir, "Significant_Pathways_pvalue_lt_0.05_with_modules.csv"),
  row.names = FALSE
)

## =============================================================================
## 5. MODULE-LEVEL SUMMARY TABLES
## =============================================================================

module_summary <- sig_pathways %>%
  group_by(.data$Module, .data$Direction) %>%
  summarise(
    n_pathways = n(),
    mean_NES = mean(.data$NES, na.rm = TRUE),
    median_NES = median(.data$NES, na.rm = TRUE),
    mean_abs_NES = mean(abs(.data$NES), na.rm = TRUE),
    best_pvalue = min(.data$pval, na.rm = TRUE),
    best_padj = min(.data$padj, na.rm = TRUE),
    mean_neg_log10_pvalue = mean(.data$neg_log10_pvalue, na.rm = TRUE),
    top_pathways = paste(head(.data$pathway_label[order(.data$pval)], 5), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(.data$Direction, .data$best_pvalue)

write.csv(
  module_summary,
  file.path(module_outdir, "Pathway_Module_Summary_pvalue_lt_0.05.csv"),
  row.names = FALSE
)

# Use explicit dplyr::count() and tidyr::pivot_wider().
# This avoids errors when another package, especially plyr, masks dplyr::count().
# Also, pivot_wider() expects bare column names, not .data$ pronouns.
module_wide <- sig_pathways %>%
  as.data.frame() %>%
  dplyr::count(Module, Direction, name = "n_pathways") %>%
  tidyr::pivot_wider(
    names_from = Direction,
    values_from = n_pathways,
    values_fill = list(n_pathways = 0)
  )

write.csv(
  module_wide,
  file.path(module_outdir, "Pathway_Module_Counts_Wide_pvalue_lt_0.05.csv"),
  row.names = FALSE
)

## =============================================================================
## 6. MODULE DOT PLOT
## =============================================================================

if (nrow(module_summary) > 0) {
  module_summary_plot <- module_summary %>%
    mutate(
      Module = fct_reorder(.data$Module, .data$mean_NES),
      DirectionColor = if_else(.data$mean_NES >= 0, "Positive NES", "Negative NES")
    )

  module_dotplot <- ggplot(
    module_summary_plot,
    aes(
      x = .data$mean_NES,
      y = .data$Module,
      size = .data$n_pathways,
      fill = .data$mean_neg_log10_pvalue,
      color = .data$DirectionColor
    )
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    geom_point(shape = 21, stroke = 1.4, alpha = 0.9) +
    scale_color_manual(
      values = c("Positive NES" = "#C9485B", "Negative NES" = "#2A69AC"),
      name = "NES direction"
    ) +
    scale_fill_gradient(low = "white", high = "#5B2C83", name = "Mean -log10(p value)") +
    scale_size_continuous(name = "Number of pathways", range = c(4, 13)) +
    labs(
      title = "Significant Hallmark pathways grouped into biological modules",
      subtitle = "Includes pathways with nominal GSEA p value < 0.05",
      x = "Mean NES per module",
      y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  ggsave(
    file.path(module_outdir, "Pathway_Module_Dotplot_pvalue_lt_0.05.pdf"),
    module_dotplot,
    width = 9,
    height = 6
  )
  ggsave(
    file.path(module_outdir, "Pathway_Module_Dotplot_pvalue_lt_0.05.png"),
    module_dotplot,
    width = 9,
    height = 6,
    dpi = 300
  )
}

## =============================================================================
## 7. PATHWAY-BY-MODULE HEATMAP-STYLE PLOT
## =============================================================================

if (nrow(sig_pathways) > 0) {
  heatmap_df <- sig_pathways %>%
    mutate(
      pathway_label = fct_reorder(.data$pathway_label, .data$NES),
      Module = factor(.data$Module)
    )

  pathway_heatmap <- ggplot(
    heatmap_df,
    aes(x = .data$Module, y = .data$pathway_label, fill = .data$NES)
  ) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient2(
      low = "#2A69AC",
      mid = "white",
      high = "#C9485B",
      midpoint = 0,
      name = "NES"
    ) +
    labs(
      title = "Significant Hallmark pathways by biological module",
      subtitle = "Red = positive NES, enriched in Rebound; blue = negative NES, enriched in NoRebound",
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid = element_blank()
    )

  plot_height <- max(6, 0.22 * nrow(heatmap_df) + 2)

  ggsave(
    file.path(module_outdir, "Pathway_Module_Heatmap_pvalue_lt_0.05.pdf"),
    pathway_heatmap,
    width = 10,
    height = plot_height,
    limitsize = FALSE
  )
  ggsave(
    file.path(module_outdir, "Pathway_Module_Heatmap_pvalue_lt_0.05.png"),
    pathway_heatmap,
    width = 10,
    height = plot_height,
    dpi = 300,
    limitsize = FALSE
  )
}

## =============================================================================
## 8. CONSOLE SUMMARY
## =============================================================================

cat("\nDONE: significant pathway module analysis complete.\n")
cat("Input file: ", gsea_path, "\n", sep = "")
cat("Pathways with nominal p value < 0.05: ", nrow(sig_pathways), "\n", sep = "")
cat("Output folder: ", module_outdir, "\n", sep = "")

if (nrow(module_summary) > 0) {
  print(module_summary)
}
################################################################################
