################################################################################
# WEEK24 BULK RNA-SEQ ANALYSIS: REBOUND vs NOREBOUND
# Method: DESeq2 differential expression + UMAP + heatmap + scatter plot + GSEA
#
# Required input files in the working directory:
#   1) Expression.csv                 genes x samples raw counts
#   2) Meta_Data.csv                  metadata for all samples
#   3) h.all.v7.0.symbols.gmt.txt     MSigDB Hallmark GMT file
#
# Important design choice:
#   - The full expression and metadata are loaded and kept as expr_full/meta_full.
#   - ALL analyses are then performed only on the Week24 subset.
#   - Status is forced to c("NoRebound", "Rebound").
#   - DESeq2 contrast is Rebound vs NoRebound.
#
# Interpretation:
#   positive log2FoldChange = higher in Rebound
#   negative log2FoldChange = higher in NoRebound
################################################################################

## =============================================================================
## 0. SETUP
## =============================================================================

set.seed(42)

outdir <- "Week24_DESeq2_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

expr_path <- "Expression.csv"
meta_path <- "Meta_Data.csv"
gmt_path  <- "h.all.v7.0.symbols.gmt.txt"

required_files <- c(expr_path, meta_path, gmt_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required file(s): ", paste(missing_files, collapse = ", "))
}

cran_packages <- c(
  "ggplot2", "pheatmap", "dplyr", "tibble", "ggrepel",
  "RColorBrewer", "uwot"
)

bioc_packages <- c("DESeq2", "fgsea")
optional_bioc_packages <- c("apeglm")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# apeglm is useful for LFC shrinkage but the analysis can proceed without it.
if (!requireNamespace("apeglm", quietly = TRUE)) {
  tryCatch(
    BiocManager::install("apeglm", update = FALSE, ask = FALSE),
    error = function(e) {
      message("Optional package apeglm could not be installed. The script will use unshrunken log2FC if needed. Reason: ", conditionMessage(e))
    }
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(tibble)
  library(ggrepel)
  library(RColorBrewer)
  library(uwot)
})

## =============================================================================
## 1. LOAD FULL DATA AND KEEP IT
## =============================================================================

# Expression.csv has genes in the first column and samples in the remaining columns.
expr_raw <- read.csv(
  expr_path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Meta_Data.csv has an index/sample-name column in the first column. Use row.names = 1
# so that the blank or unnamed first header does not become a problematic metadata column.
meta_raw <- read.csv(
  meta_path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Keep complete copies of all input data.
expr_full <- expr_raw
meta_full <- meta_raw

# Clean metadata column names and key character columns.
colnames(meta_full) <- trimws(colnames(meta_full))

required_meta_cols <- c("SampleID", "Time", "Status")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta_full))
if (length(missing_meta_cols) > 0) {
  stop("Meta_Data.csv is missing required column(s): ", paste(missing_meta_cols, collapse = ", "))
}

meta_full$SampleID <- trimws(as.character(meta_full$SampleID))
meta_full$Time     <- trimws(as.character(meta_full$Time))
meta_full$Status   <- trimws(as.character(meta_full$Status))

cat("Loaded full expression matrix: ", nrow(expr_full), " genes x ", ncol(expr_full), " samples\n", sep = "")
cat("Loaded full metadata: ", nrow(meta_full), " samples x ", ncol(meta_full), " columns\n", sep = "")
cat("Time distribution in full metadata:\n")
print(table(meta_full$Time, useNA = "ifany"))
cat("Status distribution in full metadata:\n")
print(table(meta_full$Status, useNA = "ifany"))

## =============================================================================
## 2. CREATE WEEK24-ONLY ANALYSIS OBJECTS
## =============================================================================

meta_wk24 <- meta_full %>%
  dplyr::filter(tolower(trimws(.data$Time)) == "week24") %>%
  dplyr::mutate(
    Status = trimws(as.character(.data$Status)),
    Status = dplyr::case_when(
      tolower(Status) %in% c("rebound", "rb") ~ "Rebound",
      tolower(Status) %in% c("norebound", "no rebound", "no_rebound", "nonrebound", "nr") ~ "NoRebound",
      TRUE ~ Status
    ),
    Status = factor(Status, levels = c("NoRebound", "Rebound"))
  )

if (nrow(meta_wk24) == 0) {
  stop("No Week24 samples were found in Meta_Data.csv. Check the Time column values.")
}

if (any(is.na(meta_wk24$Status))) {
  bad_status <- unique(meta_wk24$Status[is.na(meta_wk24$Status)])
  stop("Some Week24 Status values could not be mapped to NoRebound/Rebound. Check Status column.")
}

if (!all(c("NoRebound", "Rebound") %in% as.character(meta_wk24$Status))) {
  stop("Week24 data must contain both Status groups: NoRebound and Rebound.")
}

missing_expr_samples <- setdiff(meta_wk24$SampleID, colnames(expr_full))
if (length(missing_expr_samples) > 0) {
  stop(
    "These Week24 SampleID values are missing from Expression.csv columns: ",
    paste(missing_expr_samples, collapse = ", ")
  )
}

# Week24 expression matrix only, ordered exactly like Week24 metadata.
expr_wk24 <- expr_full[, meta_wk24$SampleID, drop = FALSE]
stopifnot(identical(colnames(expr_wk24), meta_wk24$SampleID))

rownames(meta_wk24) <- meta_wk24$SampleID

cat("\nWeek24 analysis subset:\n")
cat("  Samples: ", nrow(meta_wk24), "\n", sep = "")
cat("  Genes before filtering: ", nrow(expr_wk24), "\n", sep = "")
cat("  Status distribution at Week24:\n")
print(table(meta_wk24$Status))

## =============================================================================
## 3. PREPARE COUNT MATRIX FOR DESEQ2
## =============================================================================

# Convert counts to a numeric matrix, then integer counts for DESeq2.
count_mat <- as.matrix(expr_wk24)
suppressWarnings(storage.mode(count_mat) <- "numeric")

if (any(is.na(count_mat))) {
  stop("Expression matrix contains NA or non-numeric values after conversion to numeric.")
}

if (any(count_mat < 0)) {
  stop("Expression matrix contains negative counts, which are invalid for DESeq2.")
}

# DESeq2 requires integer-like raw counts.
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"

# Aggregate duplicate gene symbols if present. This preserves gene symbols for GSEA.
gene_ids <- trimws(rownames(count_mat))
valid_gene <- !is.na(gene_ids) & gene_ids != ""
count_mat <- count_mat[valid_gene, , drop = FALSE]
gene_ids <- gene_ids[valid_gene]

if (anyDuplicated(gene_ids) > 0) {
  cat("Duplicate gene symbols detected; aggregating duplicate rows by summing counts.\n")
  count_mat <- rowsum(count_mat, group = gene_ids, reorder = FALSE)
} else {
  rownames(count_mat) <- gene_ids
}

# Remove genes with zero counts across Week24 samples before creating DESeq2 object.
keep_nonzero <- rowSums(count_mat) > 0
count_mat <- count_mat[keep_nonzero, , drop = FALSE]

cat("Genes after removing zero-count genes: ", nrow(count_mat), "\n", sep = "")

## =============================================================================
## 4. DESEQ2 DIFFERENTIAL EXPRESSION: REBOUND vs NOREBOUND AT WEEK24
## =============================================================================

coldata <- meta_wk24
coldata$Status <- droplevels(coldata$Status)

# The design compares Rebound vs NoRebound only at Week24.
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = coldata,
  design    = ~ Status
)

# Filter weakly expressed genes. Keep genes with at least 10 counts in at least
# the size of the smaller Week24 group.
min_group_size <- min(table(coldata$Status))
keep_expr <- rowSums(counts(dds) >= 10) >= min_group_size
dds <- dds[keep_expr, ]

cat("Genes after DESeq2 low-count filtering: ", nrow(dds), "\n", sep = "")

# Run DESeq2.
dds <- DESeq(dds)

# Contrast: Rebound vs NoRebound.
res <- results(
  dds,
  contrast = c("Status", "Rebound", "NoRebound"),
  alpha = 0.05
)

cat("\nDESeq2 results names:\n")
print(resultsNames(dds))

# Shrink log2FC using apeglm if possible. Keep unshrunken Wald stat for GSEA.
coef_name <- "Status_Rebound_vs_NoRebound"

res_shrunk <- tryCatch(
  {
    if (requireNamespace("apeglm", quietly = TRUE) && coef_name %in% resultsNames(dds)) {
      lfcShrink(dds, coef = coef_name, type = "apeglm")
    } else {
      message("apeglm is unavailable or coefficient name was not found; using unshrunken DESeq2 log2FC.")
      res
    }
  },
  error = function(e) {
    message("apeglm shrinkage failed; using unshrunken DESeq2 results for log2FC. Reason: ", conditionMessage(e))
    res
  }
)

## =============================================================================
## 5. BUILD FINAL RESULTS TABLE SAFELY
## =============================================================================

# Never self-rename columns like rename(baseMean = baseMean). That can fail.
# Instead, build with base R and explicitly check the expected DESeq2 columns.
res_unshrunken_raw <- as.data.frame(res)
res_unshrunken_raw$Gene <- rownames(res_unshrunken_raw)

needed_cols <- c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
missing_result_cols <- setdiff(needed_cols, colnames(res_unshrunken_raw))
if (length(missing_result_cols) > 0) {
  stop(
    "The unshrunken DESeq2 result is missing required column(s): ",
    paste(missing_result_cols, collapse = ", "),
    ". Available columns are: ", paste(colnames(res_unshrunken_raw), collapse = ", ")
  )
}

res_unshrunken_df <- res_unshrunken_raw[, needed_cols]
colnames(res_unshrunken_df) <- c(
  "Gene", "baseMean", "log2FoldChange_unshrunken", "lfcSE_unshrunken",
  "stat", "pvalue", "padj"
)

res_shrunk_raw <- as.data.frame(res_shrunk)
res_shrunk_raw$Gene <- rownames(res_shrunk_raw)

if (!"log2FoldChange" %in% colnames(res_shrunk_raw)) {
  res_shrunk_df <- data.frame(
    Gene = res_unshrunken_df$Gene,
    log2FoldChange_shrunk = NA_real_,
    lfcSE_shrunk = NA_real_,
    stringsAsFactors = FALSE
  )
} else {
  res_shrunk_df <- data.frame(
    Gene = res_shrunk_raw$Gene,
    log2FoldChange_shrunk = res_shrunk_raw$log2FoldChange,
    lfcSE_shrunk = if ("lfcSE" %in% colnames(res_shrunk_raw)) res_shrunk_raw$lfcSE else NA_real_,
    stringsAsFactors = FALSE
  )
}

res_df <- res_unshrunken_df %>%
  dplyr::left_join(res_shrunk_df, by = "Gene") %>%
  dplyr::mutate(
    log2FoldChange = dplyr::if_else(
      !is.na(.data$log2FoldChange_shrunk),
      .data$log2FoldChange_shrunk,
      .data$log2FoldChange_unshrunken
    ),
    Direction = dplyr::case_when(
      is.na(.data$log2FoldChange) ~ NA_character_,
      .data$log2FoldChange > 0 ~ "Higher in Rebound",
      .data$log2FoldChange < 0 ~ "Higher in NoRebound",
      TRUE ~ "No change"
    ),
    neg_log10_pvalue = -log10(.data$pvalue),
    neg_log10_padj = -log10(.data$padj)
  ) %>%
  dplyr::arrange(.data$padj, .data$pvalue)

# Replace infinite -log10 values for plotting safety.
finite_p <- res_df$neg_log10_pvalue[is.finite(res_df$neg_log10_pvalue)]
if (length(finite_p) > 0) {
  max_finite_p <- max(finite_p, na.rm = TRUE)
  res_df$neg_log10_pvalue[!is.finite(res_df$neg_log10_pvalue)] <- max_finite_p + 1
}
finite_padj <- res_df$neg_log10_padj[is.finite(res_df$neg_log10_padj)]
if (length(finite_padj) > 0) {
  max_finite_padj <- max(finite_padj, na.rm = TRUE)
  res_df$neg_log10_padj[!is.finite(res_df$neg_log10_padj)] <- max_finite_padj + 1
}

write.csv(
  res_df,
  file.path(outdir, "01_DESeq2_Rebound_vs_NoRebound_Week24_all_genes.csv"),
  row.names = FALSE
)

sig_padj <- res_df %>%
  dplyr::filter(!is.na(.data$padj), .data$padj < 0.05)

sig_pvalue <- res_df %>%
  dplyr::filter(!is.na(.data$pvalue), .data$pvalue < 0.05)

write.csv(sig_padj, file.path(outdir, "02_DESeq2_sig_genes_padj_0.05.csv"), row.names = FALSE)
write.csv(sig_pvalue, file.path(outdir, "03_DESeq2_nominal_genes_pvalue_0.05.csv"), row.names = FALSE)

rebound_genes <- sig_pvalue %>%
  dplyr::filter(.data$log2FoldChange > 0) %>%
  dplyr::arrange(dplyr::desc(.data$log2FoldChange), .data$pvalue)

norebound_genes <- sig_pvalue %>%
  dplyr::filter(.data$log2FoldChange < 0) %>%
  dplyr::arrange(.data$log2FoldChange, .data$pvalue)

write.csv(rebound_genes, file.path(outdir, "04_Genes_Higher_In_Rebound_nominal_p_0.05.csv"), row.names = FALSE)
write.csv(norebound_genes, file.path(outdir, "05_Genes_Higher_In_NoRebound_nominal_p_0.05.csv"), row.names = FALSE)

cat("\nDESeq2 summary:\n")
print(summary(res))
cat("Nominal p < 0.05 genes: ", nrow(sig_pvalue), "\n", sep = "")
cat("Adjusted padj < 0.05 genes: ", nrow(sig_padj), "\n", sep = "")

## =============================================================================
## 6. NORMALIZED EXPRESSION FOR UMAP AND HEATMAP
## =============================================================================

vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

## =============================================================================
## 7. FIGURE 1: UMAP AT WEEK24
## =============================================================================

# Use the most variable genes for UMAP to reduce noise.
gene_var <- apply(vsd_mat, 1, var)
top_var_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(5000, length(gene_var)))]

umap_input <- t(vsd_mat[top_var_genes, , drop = FALSE])
num_samples <- nrow(umap_input)
num_neighbors <- min(10, max(2, num_samples - 1))

set.seed(42)
umap_coords <- uwot::umap(
  umap_input,
  n_neighbors = num_neighbors,
  min_dist = 0.3,
  metric = "cosine",
  verbose = FALSE
)

umap_df <- data.frame(
  SampleID = rownames(umap_input),
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    meta_wk24 %>% dplyr::select(SampleID, Status, Group, AnimalID),
    by = "SampleID"
  )

fig1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Status, shape = Status)) +
  geom_point(size = 4, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = AnimalID), size = 3, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = c("NoRebound" = "#3B76B5", "Rebound" = "#D1495B")) +
  labs(
    title = "UMAP of Week24 bulk RNA-seq samples",
    subtitle = "Variance-stabilized expression; Status forced to NoRebound vs Rebound",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme_bw(base_size = 13)

pdf(file.path(outdir, "Fig1_UMAP_Week24.pdf"), width = 7, height = 6)
print(fig1)
dev.off()

ggsave(file.path(outdir, "Fig1_UMAP_Week24.png"), fig1, width = 7, height = 6, dpi = 300)

## =============================================================================
## 8. FIGURE 2: HEATMAP OF TOP 50 GENES BETWEEN REBOUND AND NOREBOUND
## =============================================================================

top50 <- res_df %>%
  dplyr::filter(!is.na(.data$pvalue)) %>%
  dplyr::arrange(.data$pvalue) %>%
  dplyr::slice_head(n = 50)

write.csv(top50, file.path(outdir, "06_Top50_DE_genes_by_pvalue_Week24.csv"), row.names = FALSE)

if (nrow(top50) >= 2) {
  heatmap_genes <- intersect(top50$Gene, rownames(vsd_mat))
  hm_mat <- vsd_mat[heatmap_genes, , drop = FALSE]
  hm_mat_z <- t(scale(t(hm_mat)))
  hm_mat_z[is.na(hm_mat_z)] <- 0

  ann_col <- meta_wk24 %>%
    dplyr::select(Status, Group, AnimalID) %>%
    as.data.frame()
  rownames(ann_col) <- meta_wk24$SampleID

  ann_colors <- list(
    Status = c("NoRebound" = "#3B76B5", "Rebound" = "#D1495B")
  )

  pdf(file.path(outdir, "Fig2_Heatmap_Top50_DE_Genes_Week24.pdf"), width = 9, height = 10)
  pheatmap::pheatmap(
    hm_mat_z,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    show_colnames = FALSE,
    fontsize_row = 7,
    color = colorRampPalette(c("#3B76B5", "white", "#D1495B"))(100),
    main = "Top 50 genes: Rebound vs NoRebound at Week24"
  )
  dev.off()

  png(file.path(outdir, "Fig2_Heatmap_Top50_DE_Genes_Week24.png"), width = 2700, height = 3000, res = 300)
  pheatmap::pheatmap(
    hm_mat_z,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    show_colnames = FALSE,
    fontsize_row = 7,
    color = colorRampPalette(c("#3B76B5", "white", "#D1495B"))(100),
    main = "Top 50 genes: Rebound vs NoRebound at Week24"
  )
  dev.off()
} else {
  warning("Fewer than 2 genes available for top-50 heatmap; heatmap skipped.")
}

## =============================================================================
## 9. FIGURE 3: SCATTER PLOT OF GENES WITH P < 0.05
## =============================================================================

scatter_df <- res_df %>%
  dplyr::filter(!is.na(.data$pvalue), .data$pvalue < 0.05, !is.na(.data$log2FoldChange)) %>%
  dplyr::mutate(
    Regulation = dplyr::case_when(
      .data$log2FoldChange > 0 ~ "Higher in Rebound",
      .data$log2FoldChange < 0 ~ "Higher in NoRebound",
      TRUE ~ "No change"
    )
  )

top10_up <- scatter_df %>%
  dplyr::filter(.data$log2FoldChange > 0) %>%
  dplyr::arrange(dplyr::desc(.data$log2FoldChange), .data$pvalue) %>%
  dplyr::slice_head(n = 10)

top10_down <- scatter_df %>%
  dplyr::filter(.data$log2FoldChange < 0) %>%
  dplyr::arrange(.data$log2FoldChange, .data$pvalue) %>%
  dplyr::slice_head(n = 10)

label_df <- dplyr::bind_rows(top10_up, top10_down) %>%
  dplyr::distinct(.data$Gene, .keep_all = TRUE)

write.csv(label_df, file.path(outdir, "07_Scatter_Top10_Up_and_Top10_Down_Labels.csv"), row.names = FALSE)

if (nrow(scatter_df) > 0) {
  fig3 <- ggplot(scatter_df, aes(x = log2FoldChange, y = neg_log10_pvalue)) +
    geom_point(aes(color = Regulation), alpha = 0.75, size = 2) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = Gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4,
      show.legend = FALSE
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c(
      "Higher in Rebound" = "#D1495B",
      "Higher in NoRebound" = "#3B76B5",
      "No change" = "grey60"
    )) +
    labs(
      title = "Genes with nominal P < 0.05 at Week24",
      subtitle = "Top 10 upregulated and top 10 downregulated genes in Rebound are labeled",
      x = "log2 fold change: Rebound vs NoRebound",
      y = expression(-log[10](P~value)),
      color = NULL
    ) +
    theme_bw(base_size = 13)

  pdf(file.path(outdir, "Fig3_Scatter_AllGenes_Pvalue_lt_0.05_Top10_Labels.pdf"), width = 8, height = 6.5)
  print(fig3)
  dev.off()

  ggsave(file.path(outdir, "Fig3_Scatter_AllGenes_Pvalue_lt_0.05_Top10_Labels.png"), fig3, width = 8, height = 6.5, dpi = 300)
} else {
  warning("No genes with nominal p < 0.05 were found; scatter plot skipped.")
}

## =============================================================================
## 10. GSEA USING ATTACHED MSIGDB HALLMARK GMT FILE
## =============================================================================

pathways <- fgsea::gmtPathways(gmt_path)

# Use DESeq2 Wald statistic as the ranked list.
# Positive ranks represent genes higher in Rebound.
# Negative ranks represent genes higher in NoRebound.
if (!"stat" %in% colnames(res_df)) {
  stop("Column stat is missing from res_df. This should not happen because stat is taken from unshrunken DESeq2 results.")
}

rank_df <- res_df %>%
  dplyr::filter(!is.na(.data$stat), !is.na(.data$Gene), .data$Gene != "") %>%
  dplyr::group_by(.data$Gene) %>%
  dplyr::slice_max(order_by = abs(.data$stat), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

rank_vector <- rank_df$stat
names(rank_vector) <- rank_df$Gene
rank_vector <- sort(rank_vector, decreasing = TRUE)

if (length(rank_vector) < 100) {
  stop("Too few ranked genes for GSEA. Check that gene symbols and DESeq2 statistics were created correctly.")
}

gsea_res <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = rank_vector,
  minSize = 10,
  maxSize = 500
)

gsea_df <- as.data.frame(gsea_res) %>%
  dplyr::mutate(
    leadingEdge = vapply(.data$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)),
    Direction = dplyr::case_when(
      .data$NES > 0 ~ "Upregulated in Rebound",
      .data$NES < 0 ~ "Downregulated in Rebound / Higher in NoRebound",
      TRUE ~ "No direction"
    ),
    neg_log10_pvalue = -log10(.data$pval),
    pathway_label = gsub("^HALLMARK_", "", .data$pathway),
    pathway_label = gsub("_", " ", .data$pathway_label)
  ) %>%
  dplyr::arrange(.data$padj, .data$pval)

finite_gp <- gsea_df$neg_log10_pvalue[is.finite(gsea_df$neg_log10_pvalue)]
if (length(finite_gp) > 0) {
  max_finite_gp <- max(finite_gp, na.rm = TRUE)
  gsea_df$neg_log10_pvalue[!is.finite(gsea_df$neg_log10_pvalue)] <- max_finite_gp + 1
}

write.csv(gsea_df, file.path(outdir, "08_GSEA_Hallmark_Rebound_vs_NoRebound_Week24.csv"), row.names = FALSE)

top20_up_pathways <- gsea_df %>%
  dplyr::filter(.data$NES > 0) %>%
  dplyr::arrange(.data$pval, dplyr::desc(.data$NES)) %>%
  dplyr::slice_head(n = 20)

top20_down_pathways <- gsea_df %>%
  dplyr::filter(.data$NES < 0) %>%
  dplyr::arrange(.data$pval, .data$NES) %>%
  dplyr::slice_head(n = 20)

write.csv(top20_up_pathways, file.path(outdir, "09_GSEA_Top20_Upregulated_Pathways.csv"), row.names = FALSE)
write.csv(top20_down_pathways, file.path(outdir, "10_GSEA_Top20_Downregulated_Pathways.csv"), row.names = FALSE)

## =============================================================================
## 11. FIGURE 4: GSEA TOP 20 UPREGULATED PATHWAYS
## =============================================================================

if (nrow(top20_up_pathways) > 0) {
  fig4 <- ggplot(
    top20_up_pathways,
    aes(
      x = NES,
      y = reorder(pathway_label, NES),
      size = abs(NES),
      fill = neg_log10_pvalue
    )
  ) +
    geom_point(shape = 21, color = "black", alpha = 0.9) +
    scale_fill_gradient(low = "#FDE0DD", high = "#D1495B", name = expression(-log[10](P~value))) +
    scale_size_continuous(name = "|NES|") +
    labs(
      title = "Top 20 Hallmark pathways upregulated in Rebound",
      subtitle = "Positive NES = enriched among genes higher in Rebound",
      x = "Normalized enrichment score (NES)",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 9))

  pdf(file.path(outdir, "Fig4_GSEA_Top20_Upregulated_Pathways.pdf"), width = 9, height = 7)
  print(fig4)
  dev.off()

  ggsave(file.path(outdir, "Fig4_GSEA_Top20_Upregulated_Pathways.png"), fig4, width = 9, height = 7, dpi = 300)
} else {
  warning("No pathways with positive NES were found; upregulated GSEA dot plot skipped.")
}

## =============================================================================
## 12. FIGURE 5: GSEA TOP 20 DOWNREGULATED PATHWAYS
## =============================================================================

if (nrow(top20_down_pathways) > 0) {
  fig5 <- ggplot(
    top20_down_pathways,
    aes(
      x = NES,
      y = reorder(pathway_label, NES),
      size = abs(NES),
      fill = neg_log10_pvalue
    )
  ) +
    geom_point(shape = 21, color = "black", alpha = 0.9) +
    scale_fill_gradient(low = "#DEEBF7", high = "#3B76B5", name = expression(-log[10](P~value))) +
    scale_size_continuous(name = "|NES|") +
    labs(
      title = "Top 20 Hallmark pathways downregulated in Rebound",
      subtitle = "Negative NES = enriched among genes higher in NoRebound",
      x = "Normalized enrichment score (NES)",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 9))

  pdf(file.path(outdir, "Fig5_GSEA_Top20_Downregulated_Pathways.pdf"), width = 9, height = 7)
  print(fig5)
  dev.off()

  ggsave(file.path(outdir, "Fig5_GSEA_Top20_Downregulated_Pathways.png"), fig5, width = 9, height = 7, dpi = 300)
} else {
  warning("No pathways with negative NES were found; downregulated GSEA dot plot skipped.")
}

## =============================================================================
## 13. SAVE OBJECTS AND SESSION INFO
## =============================================================================

saveRDS(dds, file.path(outdir, "DESeq2_dds_Week24.rds"))
saveRDS(vsd, file.path(outdir, "DESeq2_vst_Week24.rds"))
saveRDS(expr_full, file.path(outdir, "expr_full_all_samples.rds"))
saveRDS(meta_full, file.path(outdir, "meta_full_all_samples.rds"))
saveRDS(expr_wk24, file.path(outdir, "expr_wk24_analysis_subset.rds"))
saveRDS(meta_wk24, file.path(outdir, "meta_wk24_analysis_subset.rds"))

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

cat("\nDONE. All Week24 DESeq2, UMAP, heatmap, scatter, and GSEA outputs were written to: ", outdir, "\n", sep = "")
################################################################################
