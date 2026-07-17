#!/usr/bin/env Rscript
# =============================================================================
# Practical 2 | Exercise 3 - Visualise the ML tree in R, coloured by country
#                Exercise 4 - Rooting (the R equivalent of the FigTree steps)
#
# Run:  Rscript scripts/03_plot_ml_tree.R          # Dataset A (real)
#       Rscript scripts/03_plot_ml_tree.R sim      # Dataset B (simulated)
#
# FigTree does the same job by clicking. Doing it in code means your figure is
# reproducible: you can regenerate it after you fix a sample, without
# remembering which buttons you pressed.
# =============================================================================

suppressPackageStartupMessages({
  library(ape); library(ggtree); library(ggplot2); library(dplyr); library(treeio)
})

args    <- commandArgs(trailingOnly = TRUE)
DATASET <- if (length(args) > 0) args[1] else "real"

cfg <- list(
  real = list(tree = "results/HIV_alignment.fasta.treefile",
              meta = "data/metadata.csv", colour = "country",
              tag  = "Dataset A - real HIV-1 gag, 80 sequences, 2005-2025"),
  sim  = list(tree = "results/sim/SIM_alignment.fasta.treefile",
              meta = "data/simulated/SIM_metadata.csv", colour = "country",
              tag  = "Dataset B - SIMULATED epidemic, 60 sequences"),
  h1n1 = list(tree = "results/h1n1/H1N1_alignment.fasta.treefile",
              meta = "data/h1n1/H1N1_metadata.csv", colour = "region",
              tag  = "Dataset C - 2009 pandemic H1N1, 50 genomes, 43-day window")
)
if (!DATASET %in% names(cfg))
  stop("usage: Rscript scripts/03_plot_ml_tree.R [real|sim|h1n1]")

tree_file <- cfg[[DATASET]]$tree
meta_file <- cfg[[DATASET]]$meta
COLOUR_BY <- cfg[[DATASET]]$colour
tag       <- cfg[[DATASET]]$tag
dir.create("figures", showWarnings = FALSE)

tree <- read.tree(tree_file)
meta <- read.csv(meta_file)

cat("=== TREE SUMMARY ===\n")
cat("Tips        :", Ntip(tree), "\n")
cat("Rooted?     :", is.rooted(tree), "  <- ML trees come out UNROOTED\n")
cat("Total length:", round(sum(tree$edge.length), 4), "subs/site\n")

# --- Q: which tip has the longest branch? -----------------------------------
# A long terminal branch = that sample is unusually different from everything
# else. Could be real divergence, could be a bad sequence. Always check.
term <- which(tree$edge[, 2] <= Ntip(tree))
tb   <- data.frame(
  tip    = tree$tip.label[tree$edge[term, 2]],
  length = tree$edge.length[term]
) |> arrange(desc(length))

cat("\n=== Q: LONGEST TERMINAL BRANCHES (inspect these) ===\n")
print(head(tb, 5))
cat("\nMedian terminal branch:", round(median(tb$length), 5), "subs/site\n")
cat("Longest is", round(max(tb$length) / median(tb$length), 1),
    "x the median.\n")

# --- ROOTING (Exercise 4) ---------------------------------------------------
# An unrooted tree shows relationships but NOT direction of time. You must root
# it before you can say "this lineage came first".
#
# Three ways, in order of preference:
#   1. OUTGROUP  - a sequence you know is outside the group. Best, if you have one.
#   2. MIDPOINT  - put the root at the midpoint of the longest tip-to-tip path.
#                  Assumes a roughly clocklike tree. Cheap and usually sane.
#   3. BEST-FIT  - the root that maximises temporal signal (see script 05).
tree_mid <- phytools::midpoint.root(tree)
cat("\nMidpoint-rooted. Rooted now?:", is.rooted(tree_mid), "\n")
root_out <- file.path(dirname(tree_file),
                      paste0(DATASET, "_ML_tree_midpoint_rooted.nwk"))
write.tree(tree_mid, root_out)
cat("Wrote:", root_out, "\n")

# --- parse the IQ-TREE support labels ---------------------------------------
# IQ-TREE writes node labels as "aLRT/UFboot", e.g. "98.2/100".
supp <- tree_mid$node.label
alrt <- suppressWarnings(as.numeric(sub("/.*", "", supp)))
ufb  <- suppressWarnings(as.numeric(sub(".*/", "", supp)))
cat("\n=== BRANCH SUPPORT ===\n")
cat("Nodes with UFboot >= 95 (strong)   :", sum(ufb >= 95, na.rm = TRUE),
    "/", length(na.omit(ufb)), "\n")
cat("Nodes with UFboot <  70 (ignore)   :", sum(ufb < 70, na.rm = TRUE), "\n")
cat("Nodes strong on BOTH (UFb>=95 & aLRT>=80):",
    sum(ufb >= 95 & alrt >= 80, na.rm = TRUE), "\n")

# --- plot, coloured by country ----------------------------------------------
# Put the support values in a tidy data frame keyed by node number. Do NOT try
# to parse "98/100" inside aes() - ggtree cannot parse that expression.
node_support <- data.frame(
  node   = (Ntip(tree_mid) + 1):(Ntip(tree_mid) + Nnode(tree_mid)),
  ufboot = ufb,
  alrt   = alrt
)

p <- ggtree(tree_mid, size = 0.4) %<+% meta +
  geom_tippoint(aes(colour = .data[[COLOUR_BY]]), size = 2.2, alpha = 0.9) +
  geom_tiplab(size = 1.5, offset = 0.002, align = FALSE) +
  scale_colour_brewer(palette = "Dark2", name = tools::toTitleCase(COLOUR_BY)) +
  geom_treescale(width = 0.01, fontsize = 2.5, offset = 1) +
  labs(title = paste0("Maximum Likelihood tree - ", tag),
       subtitle = "Midpoint-rooted. Black dots = ultrafast bootstrap >= 95.") +
  theme_tree() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 11))

# Mark the strongly-supported nodes. We add them as an explicit layer using the
# coordinates ggtree already computed (p$data), because parsing "98/100" inside
# aes() does not work.
nodes_df <- subset(p$data, !isTip)
nodes_df$ufboot <- ufb[nodes_df$node - Ntip(tree_mid)]
strong <- subset(nodes_df, !is.na(ufboot) & ufboot >= 95)
p <- p + geom_point(data = strong, aes(x = x, y = y), inherit.aes = FALSE,
                    colour = "black", size = 0.9, alpha = 0.7)
cat("Nodes marked on figure (UFboot >= 95):", nrow(strong), "\n")

out <- file.path("figures", paste0("03_ML_tree_", DATASET, ".pdf"))
ggsave(out, p, width = 9, height = 11, limitsize = FALSE)
cat("\nWrote:", out, "\n")

write.csv(tb, file.path("results", paste0("03_terminal_branches_", DATASET, ".csv")),
          row.names = FALSE)

# =============================================================================
# YOUR TURN
#   1. Do sequences from the same country group together? Should they?
#   2. Find the longest branch. Cross-check it against results/01_alignment_qc.csv
#      - is it a badly sequenced sample, or genuinely divergent virus?
#   3. Re-plot with `geom_tiplab(size=2)` and read the tip names. Any surprises?
#   4. Midpoint rooting assumes a clock. Script 05 tests that assumption.
#      Do NOT trust a midpoint root before you have looked at the temporal signal.
# =============================================================================
