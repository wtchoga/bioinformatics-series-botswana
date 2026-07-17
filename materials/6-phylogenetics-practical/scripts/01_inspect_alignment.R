#!/usr/bin/env Rscript
# =============================================================================
# Practical 2 | Exercise 1 - Inspect the alignment
#
# Before you build a tree, look at your data. Almost every strange tree is a
# data problem wearing a costume.
#
# Questions this script answers:
#   How many sequences?  Average length?  Missing data?
#   Which sample is shortest?  Which contains many Ns?  Any duplicates?
#
# Run:  Rscript scripts/01_inspect_alignment.R sim      # Dataset B - simulated
#       Rscript scripts/01_inspect_alignment.R h1n1     # Dataset C - 2009 H1N1
# =============================================================================

suppressPackageStartupMessages({
  library(Biostrings); library(ggplot2); library(dplyr)
})

args    <- commandArgs(trailingOnly = TRUE)
DATASET <- if (length(args) > 0) args[1] else "sim"

aln_file <- switch(DATASET,
  sim  = "data/simulated/SIM_alignment.fasta",
  h1n1 = "data/h1n1/H1N1_alignment.fasta",
  real = "data/HIV_alignment.fasta",     # optional extra dataset
  stop("usage: Rscript scripts/01_inspect_alignment.R [sim|h1n1]"))

if (!file.exists(aln_file))
  stop("Alignment not found: ", aln_file, "\n  Build it first - see run_all.sh step 0.")

cat("Dataset:", DATASET, "->", aln_file, "\n")
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# --- 1. read the alignment --------------------------------------------------
seqs <- readDNAStringSet(aln_file)

cat("\n=== BASIC COUNTS ===\n")
cat("Number of sequences :", length(seqs), "\n")
cat("Alignment width     :", paste(unique(width(seqs)), collapse = ", "), "bp\n")

# If width() returns more than one number the file is NOT aligned. Everything
# downstream assumes column 100 means the same thing in every sequence.
if (length(unique(width(seqs))) > 1)
  stop("Sequences differ in length - this file is not aligned. Align it first.")

# --- 2. per-sequence composition -------------------------------------------
# An aligned sequence is padded with '-'. The UNGAPPED length is how much real
# sequence you actually have for that sample.
mat <- alphabetFrequency(seqs, baseOnly = TRUE)   # A C G T other

qc <- tibble(
  sample_id = names(seqs),
  gaps = as.integer(letterFrequency(seqs, "-")),
  Ns   = as.integer(letterFrequency(seqs, "N")),
  A = mat[, "A"], C = mat[, "C"], G = mat[, "G"], T = mat[, "T"]
) |>
  mutate(
    aligned_width = unique(width(seqs)),
    called_bases  = A + C + G + T,
    ungapped_len  = called_bases + Ns,
    pct_missing   = round(100 * (gaps + Ns) / aligned_width, 2),
    pct_N         = round(100 * Ns / aligned_width, 2),
    gc_content    = round(100 * (G + C) / pmax(called_bases, 1), 2)
  ) |>
  arrange(desc(pct_missing))

cat("Average ungapped length :", round(mean(qc$ungapped_len), 1), "bp\n")
cat("Median  ungapped length :", median(qc$ungapped_len), "bp\n")
cat("Mean missing data       :", round(mean(qc$pct_missing), 2), "%\n")
cat("Mean GC content         :", round(mean(qc$gc_content), 2), "%\n")

# --- 3. answer the questions -----------------------------------------------
cat("\n=== Q: Which sample is SHORTEST? ===\n")
print(qc |> arrange(ungapped_len) |>
        select(sample_id, ungapped_len, pct_missing) |> head(5))

cat("\n=== Q: Which sample contains many Ns? ===\n")
print(qc |> arrange(desc(Ns)) |> select(sample_id, Ns, pct_N) |> head(5))

# --- 4. flag samples for QC -------------------------------------------------
# RULE OF THUMB (see notes/INTERPRETATION_GUIDE.md):
#   >5%  missing -> inspect
#   >10% missing -> usually exclude from a clock analysis
N_THRESHOLD <- 5
flagged <- qc |> filter(pct_missing > N_THRESHOLD)
cat("\n=== FLAGGED (>", N_THRESHOLD, "% missing) ===\n", sep = "")
if (nrow(flagged) == 0) cat("None.\n") else
  print(flagged |> select(sample_id, gaps, Ns, pct_missing))

# --- 5. duplicate sequences -------------------------------------------------
# Identical sequences are NOT evidence of transmission - they are often the
# same patient sequenced twice. Index against names(seqs), NOT qc$sample_id:
# qc has been re-sorted above, so its row order no longer matches the file.
dup_groups <- split(names(seqs), as.character(seqs))
dup_groups <- dup_groups[lengths(dup_groups) > 1]

cat("\n=== EXACT DUPLICATE SEQUENCES ===\n")
if (length(dup_groups) == 0) cat("None.\n") else {
  for (g in dup_groups) cat("  ", paste(g, collapse = "  ==  "), "\n")
  cat("(", length(dup_groups), "duplicate group(s). Keep one per patient.)\n")
}

# --- 6. how much signal is even in here? ------------------------------------
# A column that is identical in every sequence tells you nothing about the
# tree. If almost every column is constant, expect a poorly resolved tree.
cm <- consensusMatrix(seqs, baseOnly = TRUE)
variable <- sum(apply(cm[1:4, , drop = FALSE], 2, function(x) sum(x > 0) > 1))
cat("\n=== SEQUENCE VARIATION ===\n")
cat("Variable sites :", variable, "/", unique(width(seqs)),
    sprintf("(%.1f%%)\n", 100 * variable / unique(width(seqs))))
cat("Constant sites :", unique(width(seqs)) - variable, "\n")
if (variable / unique(width(seqs)) < 0.05)
  cat("-> Very little variation. Expect low bootstrap support: there is\n",
      "   simply not much evidence in the data to resolve the branching.\n")

# --- 7. figure --------------------------------------------------------------
p <- ggplot(qc, aes(x = reorder(sample_id, pct_missing), y = pct_missing)) +
  geom_col(aes(fill = pct_missing > N_THRESHOLD), width = 0.8) +
  scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#c0392b"),
                    labels = c("pass", paste0(">", N_THRESHOLD, "% missing")),
                    name = NULL) +
  coord_flip() +
  labs(title = "Missing data per sample (gaps + Ns)",
       subtitle = paste0(DATASET, ": ", length(seqs), " sequences, ",
                         unique(width(seqs)), " bp alignment"),
       x = NULL, y = "% missing") +
  theme_minimal(base_size = 7) + theme(legend.position = "top")

ggsave(paste0("figures/01_missing_data_", DATASET, ".pdf"), p,
       width = 7, height = 10)
readr::write_csv(qc, paste0("results/01_alignment_qc_", DATASET, ".csv"))

cat("\nWrote: results/01_alignment_qc_", DATASET, ".csv\n", sep = "")
cat("Wrote: figures/01_missing_data_", DATASET, ".pdf\n", sep = "")

# =============================================================================
# YOUR TURN
#   1. Run this on BOTH datasets. Which has more missing data? More variation?
#   2. Dataset C (H1N1) is 13109 bp but only ~2% of sites vary. Predict what
#      that does to bootstrap support in Exercise 2, then check whether you
#      were right.
#   3. Change N_THRESHOLD to 2. Does your list of problem samples change?
#   4. Find any duplicate pair. Which one would you keep, and why?
# =============================================================================
