################################################################################
# RANDOM FOREST ANALYSIS OF TRANSCRIPTOMIC CORRELATES OF HIV/SIV REBOUND
# AFTER ANALYTICAL TREATMENT INTERRUPTION (ATI)
#
# Input:
#   Expression.csv  - genes (rows) x samples (columns), raw RNA-seq counts
#   Meta_Data.csv   - one row per sample; must contain SampleID, Status
#                      (Rebound / NoRebound), Group, SetVL, and longitudinal
#                      plasma viral load columns d10.p.t.i ... d114.p.t.i
#
# Output (written to ./RF_output/):
#   - Biomarker_Table.csv           ranked candidate biomarkers + scores/p-values
#   - Fig1_PCA_QC.pdf
#   - Fig2_LOOCV_ROC.pdf
#   - Fig3_Boruta_Importance.pdf
#   - Fig4_RF_Importance_TopGenes.pdf
#   - Fig5_Heatmap_TopBiomarkers.pdf
#   - Fig6_Boxplots_TopGenes.pdf
#   - Genes_Associated_With_Rebound.csv
#   - Genes_Associated_With_NoRebound.csv
################################################################################


## =============================================================================
## 0. SETUP
## =============================================================================

packages_cran  <- c("randomForest","Boruta","pROC","ggplot2","pheatmap",
                     "reshape2","dplyr","tibble","RColorBrewer","ggrepel")
packages_bioc  <- c("edgeR","limma")
new_cran <- packages_cran[!(packages_cran %in% installed.packages()[,"Package"])]

# Boruta >= 10.0.0 pulls in a new dependency ('fru') that needs a Rust
# toolchain to compile from source. Most machines don't have Rust installed,
# so we install Boruta separately, forcing a binary (pre-compiled) install
# first, and falling back to the older 9.0.0 source release (no Rust needed)
# if no binary is available for this platform/R version.
if ("Boruta" %in% new_cran) {
  ok <- tryCatch({
    install.packages("Boruta", type = "binary")
    requireNamespace("Boruta", quietly = TRUE)
  }, error = function(e) FALSE, warning = function(w) FALSE)

  if (!isTRUE(ok)) {
    message("Binary install of Boruta failed/unavailable — falling back to Boruta 9.0.0 (pre-Rust-dependency release).")
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_version("Boruta", version = "9.0.0", repos = "https://cran.r-project.org")
  }
  new_cran <- setdiff(new_cran, "Boruta")
}

# Install remaining missing CRAN packages ONE AT A TIME, defensively.
# Installing in a single batch call means one package's failure (or an
# attempt to silently upgrade a package that's already loaded elsewhere in
# your R session, e.g. 'survival' via 'fitdistrplus') can abort or corrupt
# the whole install step. IMPORTANT: run this script in a FRESH R session
# (Session > Restart R in RStudio) so no packages are pre-loaded before we
# start installing/upgrading anything.
for (pkg in new_cran) {
  tryCatch({
    install.packages(pkg, dependencies = TRUE)
  }, error = function(e) {
    message("Could not install '", pkg, "': ", conditionMessage(e),
            " -- if this package was already installed, this is likely fine; ",
            "otherwise restart R and try again.")
  })
}

new_bioc <- packages_bioc[!(packages_bioc %in% installed.packages()[,"Package"])]
if (length(new_bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for (pkg in new_bioc) {
    tryCatch({
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }, error = function(e) {
      message("Could not install '", pkg, "': ", conditionMessage(e))
    })
  }
}

invisible(lapply(c(packages_cran, packages_bioc), library, character.only = TRUE))

set.seed(42)

dir.create("RF_output", showWarnings = FALSE)

## Point these at your files
expr_path <- "Expression.csv"
meta_path <- "Meta_Data.csv"


## =============================================================================
## 1. LOAD & ALIGN DATA
## =============================================================================

expr_raw <- read.csv(expr_path, row.names = 1, check.names = FALSE)
meta     <- read.csv(meta_path, row.names = 1, stringsAsFactors = FALSE)

# Make sure every metadata sample has a matching expression column, then
# force identical order (this is the single most common source of silent
# errors in RF/DE pipelines, so we check it explicitly)
stopifnot(all(meta$SampleID %in% colnames(expr_raw)))
expr_raw <- expr_raw[, meta$SampleID]
stopifnot(identical(colnames(expr_raw), meta$SampleID))

# Outcome variable. NoRebound is the reference level so that positive
# importance/ROC direction corresponds to "predicts Rebound"
meta$Status <- factor(meta$Status, levels = c("NoRebound", "Rebound"))
meta$Group  <- factor(meta$Group)

cat("Samples:", nrow(meta), " | Genes (raw):", nrow(expr_raw), "\n")
print(table(meta$Status))


## =============================================================================
## 2. FILTERING & NORMALIZATION (edgeR TMM + voom-style log2 CPM)
## =============================================================================
# RNA-seq counts must not be fed raw into RF: library size differences and the
# extreme right skew of count data will dominate tree splits. We (a) drop
# genes with negligible expression, (b) normalize library composition (TMM),
# (c) log2-transform (CPM) so that Euclidean-style variance filtering and
# correlation-based diagnostics behave sensibly downstream.

# 2a. Hard filter: drop any gene with fewer than 100 total reads summed
# across ALL samples (this is a simple, absolute floor applied before the
# more nuanced group-aware filter below).
min_total_reads <- 100
keep_reads <- rowSums(expr_raw) >= min_total_reads
cat("Genes with total counts >=", min_total_reads, "across all samples:",
    sum(keep_reads), "of", nrow(expr_raw), "\n")
expr_filt <- expr_raw[keep_reads, ]

# 2b. Design-aware filter: also require adequate expression within at least
# one group (Rebound/NoRebound), which additionally protects against genes
# that clear the total-count floor but are only expressed in 1-2 samples.
dge <- DGEList(counts = expr_filt, group = meta$Status)
keep <- filterByExpr(dge, group = meta$Status)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
logCPM <- cpm(dge, log = TRUE, prior.count = 1)   # genes x samples, filtered + normalized

cat("Genes retained after expression filtering:", nrow(logCPM), "\n")

# Sanitize gene names into valid R identifiers. Formula interfaces (used by
# randomForest, Boruta, and rfsrc below) parse "-" as arithmetic subtraction,
# so a gene symbol like "MAMU-DOA" is read as "the variable MAMU minus the
# variable DOA" and fails with "object 'MAMU-DOA' not found". We rename rows
# to syntactically-safe names and keep a lookup table to restore the
# original symbol in every reported table and figure label.
gene_map <- data.frame(safe_name = make.names(rownames(logCPM), unique = TRUE),
                        GeneSymbol = rownames(logCPM),
                        stringsAsFactors = FALSE)
rownames(logCPM) <- gene_map$safe_name


## =============================================================================
## 3. FIGURE 1 — SAMPLE-LEVEL QC (PCA)
## =============================================================================

pca <- prcomp(t(logCPM), scale. = TRUE)
var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)

pca_df <- data.frame(pca$x[, 1:2],
                      Status = meta$Status,
                      Group  = meta$Group,
                      SampleID = meta$SampleID)

fig1 <- ggplot(pca_df, aes(PC1, PC2, color = Status, shape = Group)) +
  geom_point(size = 4, alpha = 0.85) +
  labs(title = "Sample QC: PCA on filtered/normalized expression",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)")) +
  scale_color_manual(values = c(NoRebound = "#3B76B5", Rebound = "#D1495B")) +
  theme_bw(base_size = 13)

ggsave("RF_output/Fig1_PCA_QC.pdf", fig1, width = 6.5, height = 5)
print(fig1)


## =============================================================================
## 4. UNBIASED PERFORMANCE ESTIMATION: NESTED LEAVE-ONE-OUT CROSS-VALIDATION
## =============================================================================
# With n = 20 and ~18,000 candidate genes, any feature selection step MUST be
# repeated *inside* each CV fold using only the training samples, or the
# reported accuracy/AUC will be optimistically biased (selection leakage).
# LOOCV is used (rather than k-fold) because the class sizes are small
# (11 vs 9) and LOOCV maximizes the effective training size per fold.

n <- nrow(meta)
topN_prefilter <- 500     # candidate genes carried into RF at each fold

loocv_pred <- data.frame(SampleID = meta$SampleID,
                          Truth = meta$Status,
                          Prob_Rebound = NA_real_)

for (i in seq_len(n)) {

  train_idx <- setdiff(seq_len(n), i)
  train_expr  <- logCPM[, train_idx]
  test_expr   <- logCPM[, i, drop = FALSE]
  train_status <- meta$Status[train_idx]

  # -- feature pre-selection using ONLY training-fold samples (limma moderated t-test)
  design <- model.matrix(~ train_status)
  fit <- eBayes(lmFit(train_expr, design))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  sel_genes <- rownames(tt)[seq_len(min(topN_prefilter, nrow(tt)))]

  rf_train <- as.data.frame(t(train_expr[sel_genes, , drop = FALSE]))
  rf_train$Status <- train_status
  rf_test  <- as.data.frame(t(test_expr[sel_genes, , drop = FALSE]))

  rf_fold <- randomForest(Status ~ ., data = rf_train, ntree = 1000)
  pred <- predict(rf_fold, newdata = rf_test, type = "prob")

  loocv_pred$Prob_Rebound[i] <- pred[, "Rebound"]
}

loocv_pred$Pred_Class <- factor(ifelse(loocv_pred$Prob_Rebound > 0.5, "Rebound", "NoRebound"),
                                 levels = levels(meta$Status))

conf_mat <- table(Predicted = loocv_pred$Pred_Class, Truth = loocv_pred$Truth)
cat("\nLOOCV confusion matrix:\n"); print(conf_mat)
cat("LOOCV accuracy:", round(mean(loocv_pred$Pred_Class == loocv_pred$Truth), 3), "\n")

roc_obj <- pROC::roc(loocv_pred$Truth, loocv_pred$Prob_Rebound,
                      levels = c("NoRebound", "Rebound"), direction = "<")
cat("LOOCV AUC:", round(pROC::auc(roc_obj), 3), "\n")

## FIGURE 2 — cross-validated ROC curve
pdf("RF_output/Fig2_LOOCV_ROC.pdf", width = 5.5, height = 5.5)
plot(roc_obj, col = "#D1495B", lwd = 3,
     main = paste0("LOOCV ROC — AUC = ", round(pROC::auc(roc_obj), 2)))
abline(a = 1, b = -1, lty = 2, col = "grey60")
dev.off()

write.csv(loocv_pred, "RF_output/LOOCV_Predictions.csv", row.names = FALSE)


## =============================================================================
## 5. BIOMARKER DISCOVERY MODEL (uses all samples — for candidate ranking only;
##    performance is NOT re-estimated here, see Section 4 for that)
## =============================================================================

# 5a. Candidate gene pool: top differentially expressed genes (Rebound vs NoRebound)
design_full <- model.matrix(~ meta$Status)
fit_full <- eBayes(lmFit(logCPM, design_full))
tt_full <- topTable(fit_full, coef = 2, number = Inf, sort.by = "P")
tt_full$Gene <- rownames(tt_full)

n_candidates <- min(2000, nrow(tt_full))
candidate_genes <- rownames(tt_full)[seq_len(n_candidates)]

# 5b. Boruta all-relevant feature selection on the candidate pool
boruta_data <- as.data.frame(t(logCPM[candidate_genes, ]))
boruta_data$Status <- meta$Status

boruta_out   <- Boruta(Status ~ ., data = boruta_data, doTrace = 0,
                        maxRuns = 250, ntree = 1000)
boruta_final <- TentativeRoughFix(boruta_out)

boruta_stats <- attStats(boruta_final)
boruta_stats$Gene <- rownames(boruta_stats)
boruta_stats <- boruta_stats[order(-boruta_stats$meanImp), ]

cat("\nBoruta decisions:\n"); print(table(boruta_stats$decision))

## FIGURE 3 — Boruta importance (all candidate genes, boxplots of shadow vs real)
pdf("RF_output/Fig3_Boruta_Importance.pdf", width = 10, height = 6)
plot(boruta_final, las = 2, cex.axis = 0.5, xlab = "",
     main = "Boruta feature selection: Status ~ candidate genes")
dev.off()

# 5c. Final Random Forest + permutation-based significance on the top-ranked genes
top_for_rf <- boruta_stats$Gene[seq_len(min(100, nrow(boruta_stats)))]

rf_data_top <- as.data.frame(t(logCPM[top_for_rf, ]))
rf_data_top$Status <- meta$Status

rf_final <- randomForest(Status ~ ., data = rf_data_top, ntree = 2000, importance = TRUE)
cat("\nFinal RF OOB confusion matrix:\n"); print(rf_final$confusion)

obs_imp <- importance(rf_final, type = 1, scale = TRUE)[, 1]  # Mean Decrease Accuracy

# Empirical permutation test: shuffle Status 500x, refit RF on the same gene
# set, and ask how often the null importance exceeds the observed importance
n_perm <- 500
imp_null <- matrix(NA_real_, nrow = length(top_for_rf), ncol = n_perm,
                    dimnames = list(top_for_rf, NULL))

for (p in seq_len(n_perm)) {
  perm_data <- rf_data_top
  perm_data$Status <- sample(perm_data$Status)
  rf_perm <- randomForest(Status ~ ., data = perm_data, ntree = 500, importance = TRUE)
  imp_null[, p] <- importance(rf_perm, type = 1, scale = TRUE)[, 1]
}

perm_p <- sapply(names(obs_imp), function(g) mean(imp_null[g, ] >= obs_imp[g]))

## =============================================================================
## 6. ASSEMBLE FINAL BIOMARKER TABLE
## =============================================================================

biomarker_table <- data.frame(
  Gene = names(obs_imp),
  MeanDecreaseAccuracy = as.numeric(obs_imp),
  MeanDecreaseGini = importance(rf_final, type = 2)[names(obs_imp), 1],
  Permutation_p = as.numeric(perm_p)
) %>%
  left_join(boruta_stats %>% select(Gene, Boruta_meanImp = meanImp, Boruta_decision = decision),
            by = "Gene") %>%
  left_join(tt_full %>% select(Gene, logFC, DE_P.Value = P.Value, DE_adj.P.Val = adj.P.Val),
            by = "Gene") %>%
  left_join(gene_map, by = c("Gene" = "safe_name")) %>%
  relocate(GeneSymbol, .after = Gene) %>%
  arrange(Permutation_p, desc(MeanDecreaseAccuracy))

biomarker_table$Permutation_FDR <- p.adjust(biomarker_table$Permutation_p, method = "BH")

# Direction of association:
# The limma coefficient (logFC) comes from the model Rebound vs NoRebound.
#   logFC > 0  means higher expression in Rebound samples.
#   logFC < 0  means higher expression in NoRebound samples.
# These columns and separate CSV files make it explicit which candidate genes
# are associated with rebound and which are associated with no rebound.
biomarker_table <- biomarker_table %>%
  mutate(
    Higher_Expression_In = case_when(
      logFC > 0 ~ "Rebound",
      logFC < 0 ~ "NoRebound",
      TRUE ~ "No clear direction"
    ),
    Associated_With = Higher_Expression_In,
    Association_Interpretation = case_when(
      logFC > 0 ~ "Higher expression in Rebound than NoRebound",
      logFC < 0 ~ "Higher expression in NoRebound than Rebound",
      TRUE ~ "No clear expression difference by direction"
    )
  ) %>%
  relocate(Higher_Expression_In, Associated_With, Association_Interpretation, .after = GeneSymbol)

write.csv(biomarker_table, "RF_output/Biomarker_Table.csv", row.names = FALSE)

rebound_genes <- biomarker_table %>% filter(Associated_With == "Rebound")
no_rebound_genes <- biomarker_table %>% filter(Associated_With == "NoRebound")

write.csv(rebound_genes, "RF_output/Genes_Associated_With_Rebound.csv", row.names = FALSE)
write.csv(no_rebound_genes, "RF_output/Genes_Associated_With_NoRebound.csv", row.names = FALSE)

cat("\nTop candidate biomarkers:\n")
print(head(biomarker_table, 15))

cat("\nGenes associated with Rebound (higher expression in Rebound):\n")
print(rebound_genes %>% select(GeneSymbol, MeanDecreaseAccuracy, logFC, Permutation_p, Permutation_FDR) %>% head(15))

cat("\nGenes associated with NoRebound (higher expression in NoRebound):\n")
print(no_rebound_genes %>% select(GeneSymbol, MeanDecreaseAccuracy, logFC, Permutation_p, Permutation_FDR) %>% head(15))


## =============================================================================
## 7. FIGURE 4 — RANDOM FOREST IMPORTANCE PLOT (TOP GENES, WITH SIGNIFICANCE)
## =============================================================================

top_plot <- biomarker_table %>% slice_head(n = 25) %>%
  mutate(GeneSymbol = factor(GeneSymbol, levels = rev(GeneSymbol)),
         Sig = ifelse(Permutation_FDR < 0.05, "FDR<0.05",
                ifelse(Permutation_p < 0.05, "p<0.05", "n.s.")))

fig4 <- ggplot(top_plot, aes(x = MeanDecreaseAccuracy, y = GeneSymbol, fill = Sig)) +
  geom_col() +
  scale_fill_manual(values = c("FDR<0.05" = "#D1495B", "p<0.05" = "#EDAE49", "n.s." = "grey70")) +
  labs(title = "Top 25 candidate biomarkers — Random Forest importance",
       x = "Mean Decrease in Accuracy (permutation-scaled)", y = NULL,
       fill = "Permutation\ntest") +
  theme_bw(base_size = 12)

ggsave("RF_output/Fig4_RF_Importance_TopGenes.pdf", fig4, width = 7.5, height = 8)
print(fig4)


## =============================================================================
## 8. FIGURE 5 — HEATMAP OF TOP BIOMARKERS
## =============================================================================

top_genes_hm <- head(biomarker_table$Gene, 25)
hm_mat <- logCPM[top_genes_hm, ]
hm_mat_z <- t(scale(t(hm_mat)))   # gene-wise z-score across samples
rownames(hm_mat_z) <- gene_map$GeneSymbol[match(rownames(hm_mat_z), gene_map$safe_name)]

ann_col <- data.frame(Status = meta$Status, Group = meta$Group, SetVL = meta$SetVL,
                       row.names = meta$SampleID)
ann_colors <- list(Status = c(NoRebound = "#3B76B5", Rebound = "#D1495B"))

pdf("RF_output/Fig5_Heatmap_TopBiomarkers.pdf", width = 9, height = 8)
pheatmap(hm_mat_z,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         show_colnames = FALSE,
         color = colorRampPalette(c("#3B76B5", "white", "#D1495B"))(100),
         main = "Top 25 candidate biomarkers (row z-score of log2 CPM)")
dev.off()


## =============================================================================
## 9. FIGURE 6 — EXPRESSION OF TOP 6 BIOMARKERS BY REBOUND STATUS
## =============================================================================

top6 <- head(biomarker_table$Gene, 6)
box_df <- as.data.frame(t(logCPM[top6, ])) %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(meta %>% select(SampleID, Status), by = "SampleID") %>%
  reshape2::melt(id.vars = c("SampleID", "Status"), variable.name = "Gene", value.name = "log2CPM") %>%
  mutate(Gene = gene_map$GeneSymbol[match(Gene, gene_map$safe_name)])

fig6 <- ggplot(box_df, aes(Status, log2CPM, fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.8) +
  facet_wrap(~Gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(NoRebound = "#3B76B5", Rebound = "#D1495B")) +
  labs(title = "Top 6 candidate biomarkers", y = "log2 CPM", x = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "none")

ggsave("RF_output/Fig6_Boxplots_TopGenes.pdf", fig6, width = 9, height = 6)
print(fig6)


## =============================================================================
## 10. SESSION INFO (for reproducibility / methods section)
## =============================================================================
writeLines(capture.output(sessionInfo()), "RF_output/sessionInfo.txt")

cat("\nDone. All figures and tables written to ./RF_output/\n")
################################################################################
