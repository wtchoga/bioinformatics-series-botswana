#!/usr/bin/env Rscript
# =============================================================================
# Practical 2 | Dataset B - SIMULATE a clocklike dataset with KNOWN TRUTH
#
# WHY THIS EXISTS
# ---------------
# The real HIV gag alignment (Dataset A) has almost no temporal signal
# (R^2 ~ 0.06). That is a genuine property of unlinked, cross-sectional
# sequences from a mature epidemic: two random patients differ by ~14% of
# sites, while 20 years of clock evolution adds only ~2%. The date signal is
# drowned by ancient divergence. Subsetting to one subtype does not fix it.
#
# You cannot learn to READ a molecular clock on data that has no clock.
# So Part 2 of the workshop uses a SIMULATED epidemic where we chose the
# truth ourselves, and can therefore mark the answers.
#
#   *** THIS IS SIMULATED DATA. It is not from any patient.        ***
#   *** Never report these numbers as real epidemiological results. ***
#
# HOW
#   1. Simulate a birth-death epidemic over 1990-2025. Lineages that go
#      extinct = patients sampled (and removed) at that moment, which is what
#      gives us tips at MANY different dates (heterochronous sampling).
#   2. Subsample 60 tips spread across 2005-2025.
#   3. Scale time -> substitutions with a chosen clock rate, then evolve
#      sequences along the tree under GTR.
#   4. Write the truth to SIM_truth.txt so students can be marked.
#
# Run:  Rscript scripts/00_simulate_dataset.R
# =============================================================================

suppressPackageStartupMessages({
  library(ape); library(phangorn)
})

OUT <- "data/simulated"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# --- ground truth we are choosing -------------------------------------------
TRUE_RATE  <- 2.2e-3     # substitutions/site/year - typical published HIV-1 rate
ORIGIN     <- 1990.0     # start of the simulated epidemic
SAMPLE_END <- 2025.0
N_TIP      <- 60
SEQ_LEN    <- 1500

# --- 1. simulate a birth-death epidemic; keep a seed with a wide date span ---
# rlineage() retains extinct lineages, so tips sit at many different times.
found <- FALSE
for (sd in 1:400) {
  set.seed(sd)
  full <- rlineage(birth = 0.35, death = 0.28, Tmax = SAMPLE_END - ORIGIN)
  dep  <- node.depth.edgelength(full)[1:Ntip(full)]
  dates_all <- ORIGIN + dep
  inwin <- which(dates_all >= 2005 & dates_all <= SAMPLE_END)
  if (length(inwin) >= N_TIP + 15) { found <- TRUE; break }
}
if (!found) stop("No seed produced enough tips in the 2005-2025 window.")
cat(sprintf("seed=%d | simulated tips=%d | sampled in 2005-2025: %d\n",
            sd, Ntip(full), length(inwin)))

# --- 2. stratified subsample: spread sampling across the whole window --------
set.seed(20260716)
bins <- cut(dates_all[inwin], breaks = seq(2005, SAMPLE_END, by = 2))
pick <- unlist(lapply(split(inwin, bins), function(ix) {
  if (length(ix) == 0) return(integer(0))
  sample(ix, min(length(ix), ceiling(N_TIP / 10)))
}))
if (length(pick) > N_TIP) pick <- sample(pick, N_TIP)
if (length(pick) < N_TIP) {
  extra <- setdiff(inwin, pick)
  pick  <- c(pick, sample(extra, N_TIP - length(pick)))
}

tr       <- keep.tip(full, full$tip.label[pick])
tip_date <- setNames(ORIGIN + node.depth.edgelength(full)[pick],
                     full$tip.label[pick])[tr$tip.label]

# --- 3. verify the clock is exact on the TRUE tree ---------------------------
d_true <- node.depth.edgelength(tr)[1:Ntip(tr)]
TRUE_TMRCA <- unique(round(tip_date - d_true, 6))
stopifnot(length(TRUE_TMRCA) == 1)          # constant => perfect clock
fit <- lm(d_true ~ tip_date)
cat(sprintf("TRUE tree root-to-tip R^2 = %.4f (must be 1.0000 by construction)\n",
            summary(fit)$r.squared))
cat(sprintf("TRUE tMRCA                = %.2f\n", TRUE_TMRCA))
cat(sprintf("date span achieved        = %.2f - %.2f\n",
            min(tip_date), max(tip_date)))

# --- 4. find the true transmission clusters ---------------------------------
sub_tree <- tr; sub_tree$edge.length <- tr$edge.length * TRUE_RATE
pat <- cophenetic(sub_tree)

clusters <- list()
for (nd in (Ntip(tr) + 1):(Ntip(tr) + Nnode(tr))) {
  tips <- extract.clade(tr, nd)$tip.label
  if (length(tips) < 2 || length(tips) > 8) next   # match script 08
  if (max(pat[tips, tips]) <= 0.02) clusters[[as.character(nd)]] <- tips
}
keep <- rep(TRUE, length(clusters))          # drop nested clusters
for (i in seq_along(clusters)) for (j in seq_along(clusters))
  if (i != j && all(clusters[[i]] %in% clusters[[j]]) &&
      length(clusters[[i]]) < length(clusters[[j]])) keep[i] <- FALSE
clusters <- clusters[keep]
cat(sprintf("TRUE clusters (<=0.02 subs/site, 2-8 tips): %d\n", length(clusters)))

# --- 5. PLANT a metadata error ----------------------------------------------
# One sample gets a badly wrong RECORDED date (a 15-year transcription error).
# The SEQUENCE is simulated from its TRUE date, so it carries the divergence of
# an old sample while claiming to be recent. That is exactly what a real
# metadata error looks like, and root-to-tip regression should catch it.
#
# We plant it HERE, in the simulated data, and not in the real HIV data,
# because you can only SEE a date outlier if the dataset has temporal signal
# in the first place. Dataset A has none, so an outlier planted there would be
# undetectable - which is itself the lesson of Exercise 6.
rec_date <- tip_date
bad      <- which.min(abs(tip_date - 2008))      # an early sample
BAD_TRUE <- tip_date[bad]
rec_date[bad] <- BAD_TRUE + 15                   # "2008" typed as "2023"
cat(sprintf("PLANTED date error: true %.2f -> recorded %.2f\n",
            BAD_TRUE, rec_date[bad]))

# --- 6. informative tip labels ----------------------------------------------
# Labels carry the RECORDED date - that is all an analyst would ever see.
COUNTRY <- c("BWA", "ZAF", "ZMB")
cc <- sample(COUNTRY, Ntip(tr), replace = TRUE, prob = c(0.5, 0.3, 0.2))
for (k in seq_along(clusters)) {            # clusters are local: one country
  cc[match(clusters[[k]], tr$tip.label)] <- COUNTRY[(k - 1) %% 3 + 1]
}
new_lab  <- sprintf("SIM_%s_%03d_%.2f", cc, seq_len(Ntip(tr)), rec_date)
old_lab  <- tr$tip.label
clusters <- lapply(clusters, function(t) new_lab[match(t, old_lab)])
tr$tip.label    <- new_lab
names(tip_date) <- new_lab      # TRUE dates (simulation truth)
names(rec_date) <- new_lab      # RECORDED dates (what the analyst sees)
BAD_TIP <- new_lab[bad]

# --- 7. evolve sequences along the substitution tree -------------------------
sim_tree <- tr
sim_tree$edge.length <- tr$edge.length * TRUE_RATE     # years -> subs/site

set.seed(20260716)
aln <- simSeq(sim_tree, l = SEQ_LEN, type = "DNA",
              Q  = c(1.6, 4.5, 0.9, 1.1, 5.2, 1.0),    # GTR: transitions > transversions
              bf = c(0.36, 0.18, 0.24, 0.22))          # HIV-like base composition

# --- 8. write outputs -------------------------------------------------------
write.dna(as.DNAbin(aln), file.path(OUT, "SIM_alignment.fasta"),
          format = "fasta", nbcol = -1, colsep = "")
write.tree(tr, file.path(OUT, "SIM_true_timetree.nwk"))

write.table(data.frame(name = tr$tip.label,
                       date = sprintf("%.2f", rec_date[tr$tip.label])),
            file.path(OUT, "SIM_dates.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

md <- data.frame(
  sample_id    = tr$tip.label,
  country      = c(BWA = "Botswana", ZAF = "South Africa", ZMB = "Zambia")[cc],
  country_code = cc,
  date         = round(as.numeric(rec_date[tr$tip.label]), 2),
  year         = floor(rec_date[tr$tip.label]),
  cluster      = "unclustered",
  row.names    = NULL
)
for (k in seq_along(clusters))
  md$cluster[md$sample_id %in% clusters[[k]]] <- paste0("Cluster_", LETTERS[k])
write.csv(md, file.path(OUT, "SIM_metadata.csv"), row.names = FALSE)

# --- 9. the answer key ------------------------------------------------------
tx <- c(
  "GROUND TRUTH - SIMULATED dataset (Dataset B)",
  "============================================",
  "*** THIS IS SIMULATED DATA. Not from any patient. Do not report as real. ***",
  "Regenerate with: Rscript scripts/00_simulate_dataset.R",
  "",
  sprintf("TRUE clock rate    : %.2e substitutions/site/year", TRUE_RATE),
  sprintf("TRUE tMRCA (root)  : %.2f", TRUE_TMRCA),
  sprintf("Epidemic origin    : %.1f", ORIGIN),
  sprintf("Sampling window    : %.2f - %.2f", min(tip_date), max(tip_date)),
  sprintf("Tips               : %d", Ntip(tr)),
  sprintf("Alignment length   : %d bp", SEQ_LEN),
  "Substitution model : GTR, no rate heterogeneity, strict clock",
  "",
  "PLANTED METADATA ERROR (Exercise 6 - outlier detection):",
  sprintf("  Tip          : %s", BAD_TIP),
  sprintf("  RECORDED date: %.2f  (what students are given)", rec_date[bad]),
  sprintf("  TRUE date    : %.2f  (used to simulate the sequence)", BAD_TRUE),
  sprintf("  Error        : %+.0f years", rec_date[bad] - BAD_TRUE),
  "  This tip should show a large NEGATIVE root-to-tip residual: it claims to",
  "  be recent but carries the divergence of an old sample. Removing it should",
  "  visibly improve R^2. It is the ONLY planted date error in Dataset B.",
  "",
  sprintf("TRUE transmission clusters : %d", length(clusters)), "")
for (k in seq_along(clusters))
  tx <- c(tx, sprintf("  Cluster_%s (n=%d):", LETTERS[k], length(clusters[[k]])),
          paste0("    ", clusters[[k]]))
tx <- c(tx, "",
  "EXPECTED STUDENT RESULTS (vary slightly - estimation noise is normal):",
  "  Root-to-tip R^2 : ~0.7 - 0.95   (strong temporal signal)",
  "  TreeTime rate   : should bracket 2.2e-3 (roughly 1.7e-3 - 2.7e-3)",
  sprintf("  TreeTime tMRCA  : should bracket %.1f (roughly +/- 5 years)", TRUE_TMRCA),
  "",
  "Marking notes:",
  "  - A rate 10x off usually means the student used Dataset A's dates.tsv.",
  "  - tMRCA is estimated by EXTRAPOLATION back beyond the oldest sample,",
  "    so it always has wider error than the rate. Reward a sensible CI,",
  "    not a bullseye.")
writeLines(tx, file.path(OUT, "SIM_truth.txt"))

cat("\nWrote:", file.path(OUT, "SIM_alignment.fasta"), "\n")
cat("Wrote:", file.path(OUT, "SIM_dates.tsv"), "\n")
cat("Wrote:", file.path(OUT, "SIM_metadata.csv"), "\n")
cat("Wrote:", file.path(OUT, "SIM_true_timetree.nwk"), "\n")
cat("Wrote:", file.path(OUT, "SIM_truth.txt"), " <- ANSWER KEY (instructor only)\n")
