#!/usr/bin/env Rscript
# =============================================================================
# Practical 2 | Exercise 14 - Transmission cluster detection
#
# Run:  Rscript scripts/08_clusters.R sim      # has KNOWN clusters to check against
#       Rscript scripts/08_clusters.R real
#       Rscript scripts/08_clusters.R h1n1
#
# WHAT IS A CLUSTER?
# ------------------
# A group of sequences that are (a) each other's closest relatives, (b) very
# similar, and (c) well supported. This is the R equivalent of TreeCluster.
#
# The definition used here (the standard "max patristic distance" rule):
#     a clade is a cluster if EVERY pair inside it is within THRESHOLD
#     substitutions/site of each other, AND the clade's bootstrap >= MIN_BOOT
#
# THE THRESHOLD IS A CHOICE, NOT A FACT.
#   HIV pol, common convention : 0.015 (1.5%) or 0.045 (4.5%)
#   The number you pick decides how many clusters you "find". Report it.
#
# WHAT A CLUSTER IS *NOT*
#   It is NOT proof that A infected B. It means their viruses share a recent
#   common ancestor. The actual infector may be unsampled. Direction of
#   transmission cannot be read off a tree like this.
# =============================================================================

suppressPackageStartupMessages({
  library(ape); library(dplyr); library(ggplot2); library(readr)
})

args    <- commandArgs(trailingOnly = TRUE)
DATASET <- if (length(args) > 0) args[1] else "sim"

cfg <- list(
  real = list(tree = "results/HIV_alignment.fasta.treefile",
              meta = "data/metadata.csv", group = "country",
              thr = 0.045, tag = "Dataset A - real HIV-1 gag"),
  sim  = list(tree = "results/sim/SIM_alignment.fasta.treefile",
              meta = "data/simulated/SIM_metadata.csv", group = "country",
              thr = 0.02, tag = "Dataset B - SIMULATED"),
  h1n1 = list(tree = "results/h1n1/H1N1_alignment.fasta.treefile",
              meta = "data/h1n1/H1N1_metadata.csv", group = "region",
              thr = 0.001, tag = "Dataset C - 2009 H1N1")
)
if (!DATASET %in% names(cfg)) stop("usage: Rscript scripts/08_clusters.R [real|sim|h1n1]")
CFG <- cfg[[DATASET]]

THRESHOLD <- CFG$thr     # max patristic distance within a cluster (subs/site)
MIN_BOOT  <- 95          # ultrafast bootstrap
MIN_SIZE  <- 2

tree <- read.tree(CFG$tree)
meta <- read.csv(CFG$meta)

# IQ-TREE node labels are "aLRT/UFboot"
ufb <- suppressWarnings(as.numeric(sub(".*/", "", tree$node.label)))

pat <- cophenetic(tree)   # patristic (along-the-tree) distances

cat("=== Exercise 14:", CFG$tag, "===\n")
cat("Threshold      :", THRESHOLD, "subs/site (max within-cluster distance)\n")
cat("Min bootstrap  :", MIN_BOOT, "\n")
cat("Median pairwise distance in tree:", round(median(pat[upper.tri(pat)]), 4), "\n\n")

# --- walk every clade, test the rule ----------------------------------------
found <- list()
for (nd in (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))) {
  tips <- extract.clade(tree, nd)$tip.label
  if (length(tips) < MIN_SIZE) next
  boot <- ufb[nd - Ntip(tree)]
  if (is.na(boot) || boot < MIN_BOOT) next
  if (max(pat[tips, tips]) > THRESHOLD) next
  found[[as.character(nd)]] <- tips
}
# keep only the LARGEST cluster where one is nested inside another
if (length(found) > 1) {
  drop <- rep(FALSE, length(found))
  for (i in seq_along(found)) for (j in seq_along(found))
    if (i != j && all(found[[i]] %in% found[[j]]) &&
        length(found[[i]]) < length(found[[j]])) drop[i] <- TRUE
  found <- found[!drop]
}

cat("CLUSTERS FOUND:", length(found), "\n")
if (length(found) == 0) {
  cat("None at this threshold. That is a RESULT, not a failure - it means no\n")
  cat("group of samples is closely related enough to meet your definition.\n")
}

rows <- list()
for (k in seq_along(found)) {
  tips <- found[[k]]
  nd   <- as.integer(names(found)[k])
  sub  <- pat[tips, tips]
  m    <- meta[match(tips, meta[[1]]), ]
  grp  <- table(m[[CFG$group]])
  cat(sprintf("\nCluster_%s  n=%d  bootstrap=%.0f  max_dist=%.4f\n",
              LETTERS[k], length(tips), ufb[nd - Ntip(tree)], max(sub)))
  cat("  ", paste(names(grp), grp, sep = "=", collapse = ", "), "\n")
  for (t in tips) cat("   -", t, "\n")
  rows[[k]] <- tibble(cluster = paste0("Cluster_", LETTERS[k]),
                      tip = tips, n = length(tips),
                      bootstrap = ufb[nd - Ntip(tree)],
                      max_distance = max(sub),
                      single_group = length(grp) == 1)
}

out <- if (length(rows)) bind_rows(rows) else
  tibble(cluster = character(), tip = character())
write_csv(out, paste0("results/08_clusters_", DATASET, ".csv"))

# --- if we KNOW the truth, mark the answer ----------------------------------
if (DATASET == "sim" && file.exists("data/simulated/SIM_metadata.csv")) {
  truth <- read.csv("data/simulated/SIM_metadata.csv") |>
    filter(cluster != "unclustered")
  cat("\n=== MARKING AGAINST TRUTH (Dataset B only) ===\n")
  cat("True clustered tips :", nrow(truth), "in",
      length(unique(truth$cluster)), "clusters\n")
  if (nrow(out) > 0) {
    tp <- sum(out$tip %in% truth$sample_id)
    cat("Recovered tips      :", nrow(out), "of which", tp, "are truly clustered\n")
    cat("Precision           :", round(tp / nrow(out), 2), "\n")
    cat("Recall              :", round(tp / nrow(truth), 2), "\n")
    cat("\nMissed tips are usually in clusters whose bootstrap fell below",
        MIN_BOOT, "-\ntightening the support cutoff costs you recall. That trade-off is\nthe whole game.\n")
  }
}

# --- how does the threshold change the answer? ------------------------------
# Students should SEE that "how many clusters" depends on a number they chose.
grid <- tibble(threshold = seq(THRESHOLD / 4, THRESHOLD * 2.5, length.out = 12)) |>
  rowwise() |>
  mutate(n_clusters = {
    th <- threshold; cnt <- list()
    for (nd in (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))) {
      tips <- extract.clade(tree, nd)$tip.label
      if (length(tips) < MIN_SIZE) next
      b <- ufb[nd - Ntip(tree)]
      if (is.na(b) || b < MIN_BOOT) next
      if (max(pat[tips, tips]) > th) next
      cnt[[as.character(nd)]] <- tips
    }
    if (length(cnt) > 1) {
      dr <- rep(FALSE, length(cnt))
      for (i in seq_along(cnt)) for (j in seq_along(cnt))
        if (i != j && all(cnt[[i]] %in% cnt[[j]]) &&
            length(cnt[[i]]) < length(cnt[[j]])) dr[i] <- TRUE
      cnt <- cnt[!dr]
    }
    length(cnt)
  }) |> ungroup()

p <- ggplot(grid, aes(threshold, n_clusters)) +
  geom_line(colour = "#2c3e50") + geom_point(size = 2) +
  geom_vline(xintercept = THRESHOLD, linetype = "dashed", colour = "#c0392b") +
  annotate("text", x = THRESHOLD, y = max(grid$n_clusters), hjust = -0.1,
           label = "your threshold", colour = "#c0392b", size = 3) +
  labs(title = paste0("Cluster count depends on the threshold you chose - ", CFG$tag),
       subtitle = paste0("bootstrap cutoff fixed at ", MIN_BOOT,
                         ". There is no 'true' number of clusters independent of this choice."),
       x = "max within-cluster patristic distance (subs/site)",
       y = "number of clusters") +
  theme_bw(base_size = 10) + theme(plot.subtitle = element_text(size = 8))

ggsave(paste0("figures/08_cluster_threshold_", DATASET, ".pdf"), p, width = 8, height = 4.5)
write_csv(grid, paste0("results/08_threshold_sweep_", DATASET, ".csv"))

cat("\nWrote: results/08_clusters_", DATASET, ".csv\n", sep = "")
cat("Wrote: figures/08_cluster_threshold_", DATASET, ".pdf\n", sep = "")

# =============================================================================
# YOUR TURN
#   1. How many clusters did you find? Now look at
#      figures/08_cluster_threshold_*.pdf - how many would you have found with
#      a threshold half as big? Does "3 transmission clusters" mean anything
#      without stating the threshold?
#   2. Are your clusters within ONE country (single_group = TRUE)? Should a
#      real transmission cluster be?
#   3. For Dataset B, compare your clusters to data/simulated/SIM_truth.txt.
#      Which true clusters did you miss, and why?
#   4. Drop MIN_BOOT to 70 and rerun. You find more clusters. Are they real?
# =============================================================================
