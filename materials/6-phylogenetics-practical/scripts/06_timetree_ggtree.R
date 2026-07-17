#!/usr/bin/env Rscript
# =============================================================================
# Practical 2 | Exercise 9  - Read the TreeTime tree, estimate tMRCA
#                Exercise 10 - Visualise the time tree in R (ggtree)
#                Exercise 12 - Export a publication-quality figure
#                Exercise 13 - Map metadata onto the tree
#
# Run:  Rscript scripts/06_timetree_ggtree.R sim      # start here - it works
#       Rscript scripts/06_timetree_ggtree.R h1n1
#       Rscript scripts/06_timetree_ggtree.R real     # see the warning it prints
# =============================================================================

suppressPackageStartupMessages({
  library(ape); library(treeio); library(ggtree); library(ggplot2)
  library(dplyr); library(lubridate)
})

args    <- commandArgs(trailingOnly = TRUE)
DATASET <- if (length(args) > 0) args[1] else "sim"

cfg <- list(
  real = list(tt = "results/treetime_real/timetree.nexus",
              meta = "data/metadata.csv", colour = "country",
              tag = "Dataset A - real HIV-1 gag"),
  sim  = list(tt = "results/sim/treetime_sim/timetree.nexus",
              meta = "data/simulated/SIM_metadata.csv", colour = "country",
              tag = "Dataset B - SIMULATED epidemic"),
  h1n1 = list(tt = "results/h1n1/treetime_h1n1/timetree.nexus",
              meta = "data/h1n1/H1N1_metadata.csv", colour = "region",
              tag = "Dataset C - 2009 pandemic H1N1")
)
if (!DATASET %in% names(cfg))
  stop("usage: Rscript scripts/06_timetree_ggtree.R [real|sim|h1n1]")
CFG <- cfg[[DATASET]]

# Read R^2 from TreeTime's own output rather than hard-coding it - if you
# rerun TreeTime with different flags, the figure must tell the new truth.
clock_file <- file.path(dirname(CFG$tt), "molecular_clock.txt")
CFG$r2 <- if (file.exists(clock_file)) {
  as.numeric(sub(".*:\\s*", "", grep("r\\^2", readLines(clock_file), value = TRUE)[1]))
} else NA_real_
cat("Root-to-tip R^2 (from TreeTime):", CFG$r2, "\n")

tt   <- read.beast(CFG$tt)
phy  <- as.phylo(tt)
meta <- read.csv(CFG$meta)

# --- 1. Exercise 9: estimate the tMRCA --------------------------------------
# A TreeTime tree is measured in YEARS. The root sits at the tMRCA. To find its
# date: take any tip's sampling date and subtract its distance back to the root.
# (Do NOT trust node names like "NODE_0000000" to be the root - they are not.)
tip_date <- switch(DATASET,
  real = as.numeric(sub(".*_(\\d{4})$", "\\1", phy$tip.label)),
  sim  = as.numeric(sub("^SIM_[A-Z]{3}_[0-9]+_", "", phy$tip.label)),
  h1n1 = as.numeric(sub("^.*[|]", "", phy$tip.label))
)
depth <- node.depth.edgelength(phy)[1:Ntip(phy)]
tmrca <- tip_date - depth

cat("=== Exercise 9: tMRCA ===\n")
cat(sprintf("Estimated tMRCA : %.2f\n", mean(tmrca)))
cat(sprintf("Tree height     : %.2f years\n", max(depth)))
cat(sprintf("Most recent tip : %.2f\n", max(tip_date)))

# A tMRCA without a confidence interval is not a result. TreeTime writes the
# bounds into its OWN dates.tsv (not the nexus), but only when you run it with
# --covariation. The root is the node numbered Ntip+1 in ape's convention -
# do NOT assume "NODE_0000000" is the root, because it usually is not.
root_label <- as_tibble(tt)$label[as_tibble(tt)$node == Ntip(phy) + 1]
dates_file <- file.path(dirname(CFG$tt), "dates.tsv")
ci_txt <- "not available (rerun TreeTime with --covariation)"
if (file.exists(dates_file)) {
  dd <- read.delim(dates_file, comment.char = "#", header = FALSE,
                   col.names = c("node", "date", "numeric_date", "lower", "upper"),
                   fill = TRUE, colClasses = "character")
  # TIP rows carry no bounds, so these columns come back as text. Coerce
  # explicitly rather than trusting read.delim to guess the type.
  for (v in c("numeric_date", "lower", "upper"))
    dd[[v]] <- suppressWarnings(as.numeric(dd[[v]]))
  rr <- dd[dd$node == root_label, ]
  if (nrow(rr) == 1 && !is.na(rr$lower))
    ci_txt <- sprintf("%.2f - %.2f (90%% region)", rr$lower, rr$upper)
  if (nrow(rr) == 1)
    cat(sprintf("Root date (TreeTime): %.2f\n", rr$numeric_date))
}
cat("Root 90% interval  :", ci_txt, "\n")

# --- 2. the honesty check ---------------------------------------------------
if (CFG$r2 < 0.4) {
  cat("\n***********************************************************\n")
  cat("* WARNING: this dataset's root-to-tip R^2 is", CFG$r2, "\n")
  cat("* TreeTime WILL still draw you a beautiful time tree with\n")
  cat("* dates on every node. It will look completely convincing.\n")
  cat("* It is not trustworthy. A time tree is only as good as the\n")
  cat("* temporal signal underneath it, and this one has almost none.\n")
  cat("* We plot it anyway so you can SEE that a pretty figure is\n")
  cat("* not evidence. Do not report these dates.\n")
  cat("***********************************************************\n\n")
}

# --- 3. Exercise 13: attach the metadata ------------------------------------
# %<+% joins a data frame to the tree by matching the FIRST column against tip
# labels. If your colours come out all grey, the names did not match.
stopifnot(all(phy$tip.label %in% meta[[1]]))

# --- 4. Exercise 10 + 12: the figure ----------------------------------------
p <- ggtree(tt, mrsd = NULL, size = 0.45) %<+% meta +
  geom_tippoint(aes(colour = .data[[CFG$colour]]), size = 2.2, alpha = 0.9) +
  scale_colour_brewer(palette = "Dark2",
                      name = tools::toTitleCase(CFG$colour)) +
  labs(title = paste0("Time-scaled phylogeny - ", CFG$tag),
       subtitle = sprintf("TreeTime | estimated tMRCA = %.1f | root-to-tip R2 = %.2f%s",
                          mean(tmrca), CFG$r2,
                          ifelse(CFG$r2 < 0.4, "  <- WEAK SIGNAL: dates unreliable", "")),
       x = "Year", caption = "Tip points coloured by metadata. x-axis is TIME, not substitutions.") +
  theme_tree2() +
  theme(text = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8))

# put a real calendar axis on it: ggtree measures x from the root, so shift it
root_year <- mean(tmrca)
p <- p + scale_x_continuous(
  labels = function(x) round(x + root_year),
  breaks = scales::pretty_breaks(6))

if (Ntip(phy) <= 80)
  p <- p + geom_tiplab(size = 1.5, offset = 0.15)

out <- paste0("figures/06_timetree_", DATASET, ".pdf")
ggsave(out, p, width = 10, height = 8, limitsize = FALSE)
cat("Wrote:", out, "\n")

# --- 5. Exercise 12: the publication version --------------------------------
# Same tree, cleaned up: no tip labels, bigger points, clean theme.
pub <- ggtree(tt, size = 0.5) %<+% meta +
  geom_tippoint(aes(colour = .data[[CFG$colour]]), size = 3, alpha = 0.95) +
  scale_colour_brewer(palette = "Dark2", name = tools::toTitleCase(CFG$colour)) +
  scale_x_continuous(labels = function(x) round(x + root_year),
                     breaks = scales::pretty_breaks(6)) +
  labs(title = CFG$tag, x = "Year") +
  theme_tree2() +
  theme(text = element_text(size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.12, 0.85),
        legend.background = element_rect(fill = "white", colour = "grey80"),
        plot.title = element_text(face = "bold"))

ggsave(paste0("figures/Figure1_", DATASET, ".pdf"), pub, width = 12, height = 8)
cat("Wrote: figures/Figure1_", DATASET, ".pdf  <- Exercise 12 deliverable\n", sep = "")

# --- 6. node dates table ----------------------------------------------------
res <- tibble(tip = phy$tip.label, sampling_date = tip_date,
              root_to_tip_years = depth)
readr::write_csv(res, paste0("results/06_tip_dates_", DATASET, ".csv"))

# =============================================================================
# YOUR TURN
#   1. Exercise 9: what is your tMRCA? What is its confidence interval?
#      Is the interval narrow enough to mean anything?
#   2. For Dataset C (H1N1), compare your tMRCA to the published estimate for
#      the 2009 pandemic (roughly late 2008 / January 2009). How close are you?
#   3. Exercise 13: recolour by a different column (edit CFG$colour above).
#      Does geography explain the tree shape, or not?
#   4. Compare figures/06_timetree_sim.pdf with figures/03_ML_tree_sim.pdf.
#      Same data, same topology - why do they look so different?
# =============================================================================
