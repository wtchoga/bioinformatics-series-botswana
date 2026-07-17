#!/usr/bin/env Rscript
# =============================================================================
# Practical 2 | Exercise 5 - Root-to-tip regression   (TempEst, but in R)
#                Exercise 6 - Interpret it
#
# Run:  Rscript scripts/04_root_to_tip.R real     # Dataset A - real HIV gag
#       Rscript scripts/04_root_to_tip.R sim      # Dataset B - simulated
#       Rscript scripts/04_root_to_tip.R h1n1     # Dataset C - 2009 H1N1
#
# WHAT THIS IS
# ------------
# Plot each tip's genetic distance from the root against its sampling date.
# If the virus evolves like a clock, older samples sit closer to the root and
# newer ones further away, so the points form a rising line.
#
#   slope      = evolutionary rate (substitutions/site/YEAR)
#   x-intercept= tMRCA (the date where divergence = 0)
#   R^2        = how much of the divergence the DATES explain
#
# THE BIG WARNING
# ---------------
# This regression is NOT a statistical test. The points are not independent -
# tips share ancestry, so they share branches. Use it as a DIAGNOSTIC to decide
# whether a clock analysis is worth running, and to spot bad samples. Never
# quote its R^2 as if it were a p-value.
# =============================================================================

suppressPackageStartupMessages({
  library(ape); library(ggplot2); library(dplyr); library(ggrepel)
})

args    <- commandArgs(trailingOnly = TRUE)
DATASET <- if (length(args) > 0) args[1] else "real"

cfg <- list(
  real = list(tree = "results/HIV_alignment.fasta.treefile",
              dates = "data/dates.tsv",
              tag  = "Dataset A - real HIV-1 gag (80 seq, 2005-2025)"),
  sim  = list(tree = "results/sim/SIM_alignment.fasta.treefile",
              dates = "data/simulated/SIM_dates.tsv",
              tag  = "Dataset B - SIMULATED (60 seq, 2005-2025)"),
  h1n1 = list(tree = "results/h1n1/H1N1_alignment.fasta.treefile",
              dates = "data/h1n1/H1N1_dates.tsv",
              tag  = "Dataset C - 2009 pandemic H1N1 (50 genomes, 43 days)")
)
if (!DATASET %in% names(cfg))
  stop("usage: Rscript scripts/04_root_to_tip.R [real|sim|h1n1]")

tree  <- read.tree(cfg[[DATASET]]$tree)
dts   <- read.delim(cfg[[DATASET]]$dates)
dates <- setNames(as.numeric(dts$date), dts$name)

stopifnot(all(tree$tip.label %in% names(dates)))
dates <- dates[tree$tip.label]

# ---------------------------------------------------------------------------
# helper: root-to-tip fit for a given rooted tree
# ---------------------------------------------------------------------------
rtt_fit <- function(rt) {
  d <- node.depth.edgelength(rt)[1:Ntip(rt)]
  # lm() throws the names away, so carry the tip labels explicitly - otherwise
  # your "outliers" come back labelled 1,2,3... and you cannot act on them.
  names(d) <- rt$tip.label
  y   <- dates[rt$tip.label]
  fit <- lm(d ~ y)
  r   <- residuals(fit); names(r) <- rt$tip.label
  list(r2 = summary(fit)$r.squared,
       slope = unname(coef(fit)[2]),
       tmrca = unname(-coef(fit)[1] / coef(fit)[2]),
       corr  = suppressWarnings(cor(d, y)),
       dist  = d, date = y, resid = r, fit = fit)
}

# ---------------------------------------------------------------------------
# 1. BEST-FITTING ROOT  (this is what TempEst's "best fit" button does)
#    Try rooting on every branch, at several positions along it, and keep the
#    root that maximises R^2. It is a search, not magic.
# ---------------------------------------------------------------------------
cat("Searching for the best-fitting root over", nrow(tree$edge), "branches...\n")
best <- list(r2 = -Inf)
for (e in seq_len(nrow(tree$edge))) {
  nd <- tree$edge[e, 2]
  el <- tree$edge.length[e]
  if (is.na(el) || el <= 0) next
  for (f in seq(0.1, 0.9, length.out = 5)) {
    rt <- tryCatch(phytools::reroot(tree, nd, position = f * el),
                   error = function(z) NULL)
    if (is.null(rt)) next
    fitr <- rtt_fit(rt)
    if (is.finite(fitr$r2) && fitr$r2 > best$r2) { best <- fitr; best$tree <- rt }
  }
}

mid  <- rtt_fit(phytools::midpoint.root(tree))

# ---------------------------------------------------------------------------
# 2. REPORT
# ---------------------------------------------------------------------------
band <- function(r2) {
  if (r2 > 0.80) "Excellent - strong molecular clock; proceed confidently"
  else if (r2 > 0.60) "Good - suitable for time-scaled analyses"
  else if (r2 > 0.40) "Moderate - use caution; inspect outliers"
  else if (r2 > 0.20) "Weak - consider removing problematic sequences"
  else "Very weak - clock-based dating is UNRELIABLE"
}

cat("\n==============================================================\n")
cat(cfg[[DATASET]]$tag, "\n")
cat("==============================================================\n")
cat(sprintf("%-22s %10s %10s\n", "", "midpoint", "best-fit"))
cat(sprintf("%-22s %10.3f %10.3f\n", "R^2",          mid$r2,    best$r2))
cat(sprintf("%-22s %10.3f %10.3f\n", "correlation",  mid$corr,  best$corr))
cat(sprintf("%-22s %10.2e %10.2e\n", "slope (subs/site/yr)", mid$slope, best$slope))
cat(sprintf("%-22s %10.1f %10.1f\n", "x-intercept (tMRCA)", mid$tmrca, best$tmrca))
cat(sprintf("%-22s %10.4f %10.4f\n", "residual mean",
            mean(abs(mid$resid)), mean(abs(best$resid))))
cat("\nSampling window   :", sprintf("%.2f - %.2f (%.1f years)",
    min(dates), max(dates), diff(range(dates))), "\n")
cat("INTERPRETATION    :", band(best$r2), "\n")

if (mid$slope < 0)
  cat("\n*** NEGATIVE SLOPE under midpoint rooting: divergence DECREASES with\n",
      "    time. That is biologically impossible under a clock. It means there\n",
      "    is no usable temporal signal, or the dates are wrong.\n")

# --- the trap nobody warns students about ----------------------------------
if (best$r2 - mid$r2 > 0.10)
  cat("\n*** CAUTION: the best-fit root SEARCHES for the root that maximises\n",
      sprintf("    R^2, so its R^2 is optimistically biased (%.3f -> %.3f here).\n",
              mid$r2, best$r2),
      "    You fitted the root TO the dates, then used the dates to judge the\n",
      "    fit. A best-fit R^2 of 0.2-0.3 obtained this way is NOT evidence of\n",
      "    a clock. Trust an outgroup root, or a best-fit root that barely\n",
      "    improves on midpoint.\n")

# ---------------------------------------------------------------------------
# 3. OUTLIERS - large residuals
#    A residual is how far a tip sits from the regression line. A tip far above
#    the line is more diverged than its date predicts; far below = less.
#    Cause is usually: wrong date, bad sequence, recombinant, or contamination.
# ---------------------------------------------------------------------------
res <- tibble(tip = names(best$resid), residual = as.numeric(best$resid),
              distance = best$dist, date = best$date) |>
  mutate(z = as.numeric(scale(residual))) |>
  arrange(desc(abs(z)))

cat("\n=== TOP RESIDUALS (|z| > 2 is worth investigating) ===\n")
print(res |> mutate(across(where(is.numeric), ~round(.x, 4))) |> head(6))

outliers <- res |> filter(abs(z) > 2)
cat("\nFlagged outliers (|z| > 2):", nrow(outliers), "\n")

# ---------------------------------------------------------------------------
# 4. "What if I drop the outliers?" - the honest version of the exercise
# ---------------------------------------------------------------------------
if (nrow(outliers) > 0) {
  keep <- setdiff(tree$tip.label, outliers$tip)
  t2   <- keep.tip(tree, keep)
  b2   <- list(r2 = -Inf)
  for (e in seq_len(nrow(t2$edge))) {
    nd <- t2$edge[e, 2]; el <- t2$edge.length[e]
    if (is.na(el) || el <= 0) next
    rt <- tryCatch(phytools::reroot(t2, nd, position = 0.5 * el),
                   error = function(z) NULL)
    if (is.null(rt)) next
    f <- rtt_fit(rt); if (is.finite(f$r2) && f$r2 > b2$r2) b2 <- f
  }
  cat(sprintf("\nDropping %d outlier(s): R^2 %.3f -> %.3f\n",
              nrow(outliers), best$r2, b2$r2))
  cat("Rule: if removing ONE sequence rescues R^2, that sequence was the\n")
  cat("problem. If R^2 barely moves, the dataset simply has no clock signal\n")
  cat("and deleting samples until it looks good is data dredging - do not.\n")
}

# ---------------------------------------------------------------------------
# 5. PLOT
# ---------------------------------------------------------------------------
df <- tibble(date = best$date, distance = best$dist, tip = names(best$dist)) |>
  left_join(res |> select(tip, z), by = "tip")

p <- ggplot(df, aes(date, distance)) +
  geom_smooth(method = "lm", se = TRUE, colour = "#2c3e50",
              fill = "grey80", linewidth = 0.6, formula = y ~ x) +
  geom_point(aes(colour = abs(z) > 2), size = 2.2, alpha = 0.85) +
  ggrepel::geom_text_repel(data = filter(df, abs(z) > 2),
                           aes(label = tip), size = 2.4, max.overlaps = 12,
                           min.segment.length = 0) +
  scale_colour_manual(values = c(`FALSE` = "#2980b9", `TRUE` = "#c0392b"),
                      labels = c("within expectation", "outlier (|z|>2)"),
                      name = NULL) +
  labs(title = paste0("Root-to-tip regression - ", cfg[[DATASET]]$tag),
       subtitle = sprintf(
         "best-fit root | R2 = %.3f | slope = %.2e subs/site/yr | x-intercept (tMRCA) = %.1f\n%s",
         best$r2, best$slope, best$tmrca, band(best$r2)),
       x = "Sampling date", y = "Root-to-tip distance (subs/site)") +
  theme_bw(base_size = 10) +
  theme(plot.subtitle = element_text(size = 8), legend.position = "top")

ggsave(paste0("figures/04_root_to_tip_", DATASET, ".pdf"), p, width = 8, height = 5.5)
readr::write_csv(res, paste0("results/04_root_to_tip_residuals_", DATASET, ".csv"))
write.tree(best$tree, paste0("results/04_bestroot_", DATASET, ".nwk"))

sink(paste0("results/04_root_to_tip_summary_", DATASET, ".txt"))
cat(cfg[[DATASET]]$tag, "\n")
cat("best-fit root R^2   :", round(best$r2, 4), "\n")
cat("correlation         :", round(best$corr, 4), "\n")
cat("slope (subs/site/yr):", format(best$slope, digits = 3), "\n")
cat("x-intercept (tMRCA) :", round(best$tmrca, 2), "\n")
cat("sampling window     :", sprintf("%.2f - %.2f", min(dates), max(dates)), "\n")
cat("outliers (|z|>2)    :", nrow(outliers), "\n")
cat("interpretation      :", band(best$r2), "\n")
sink()

cat("\nWrote: figures/04_root_to_tip_", DATASET, ".pdf\n", sep = "")
cat("Wrote: results/04_root_to_tip_residuals_", DATASET, ".csv\n", sep = "")
cat("Wrote: results/04_bestroot_", DATASET, ".nwk\n", sep = "")

# =============================================================================
# YOUR TURN
#   1. Run this on all THREE datasets. Write down R^2 for each.
#   2. Dataset A (HIV) and Dataset C (H1N1) both have low R^2 - but for
#      COMPLETELY different reasons. Work out what they are.
#      (Hint: compare the sampling windows, and the divergence on the y-axis.)
#   3. Which dataset would you be willing to date? Defend your answer.
#   4. Compare the "midpoint" and "best-fit" columns. Does the choice of root
#      change your conclusion?
# =============================================================================
