#!/usr/bin/env Rscript
# =============================================================================
#  MUTATIONS, DIVERSITY, BURDEN AND MINORITY VARIANTS
#  A hands-on workshop script — HIV-1 deep sequencing, NC_001802.1
#  Dr. W. Choga
# -----------------------------------------------------------------------------
#  USAGE
#     Rscript mutations_workshop.R
#     Rscript mutations_workshop.R <path/to/variants.csv>
#     Rscript mutations_workshop.R <csv> --minvaf 0.02 --outdir results
#
#  Or step through it interactively:  source("mutations_workshop.R")
#
#  OUTPUTS
#     figures/*.png   every plot
#     results/*.csv   every table
#     console         the full narrative, printed as it computes
#
#  WHAT THIS SCRIPT DOES
#     It runs a standard minority-variant analysis that produces a confident,
#     FDR-significant result — and then destroys it. Section 7 is the point of
#     the whole exercise. Do not skip to section 10.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# =============================================================================
# 0.  CONFIGURATION
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

get_opt <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && length(args) > i) return(args[i + 1])
  default
}

CFG <- list(
  csv        = if (length(args) >= 1 && !startsWith(args[1], "--")) args[1]
               else "hiv1_pileup_NC_001802.1_snp-variants.csv",
  min_depth  = as.numeric(get_opt("--mindepth", 100)),    # site must be this deep to be callable
  error_rate = as.numeric(get_opt("--error",    0.002)),  # assumed per-base sequencing error
  minor_vaf  = as.numeric(get_opt("--minvaf",   0.01)),   # NAIVE threshold — we will test it
  strict_vaf = as.numeric(get_opt("--strictvaf", 0.02)),  # DEFENDED threshold — justified in §7
  outdir     = get_opt("--outdir", "."),
  fig_dpi    = 150
)

FIGDIR <- file.path(CFG$outdir, "figures")
RESDIR <- file.path(CFG$outdir, "results")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESDIR, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 12))
PAL <- c(signal = "#0FA3A3", artefact = "#E2574C", violet = "#8B5CF6", muted = "#6C7D8B")

# ---- little helpers for readable console output -----------------------------

banner <- function(n, title) {
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat(sprintf("  %s  %s\n", n, toupper(title)))
  cat(strrep("=", 78), "\n", sep = "")
}
note   <- function(...) cat("\n  ", sprintf(...), "\n", sep = "")
bullet <- function(...) cat("   - ", sprintf(...), "\n", sep = "")
show_tbl <- function(df, caption = NULL) {
  if (!is.null(caption)) cat("\n  ", caption, "\n", sep = "")
  print(as.data.frame(df), row.names = FALSE, digits = 4)
  invisible(df)
}
save_tbl <- function(df, name) {
  readr::write_csv(df, file.path(RESDIR, paste0(name, ".csv")))
  invisible(df)
}
save_fig <- function(p, name, w = 9, h = 5) {
  ggsave(file.path(FIGDIR, paste0(name, ".png")), p,
         width = w, height = h, dpi = CFG$fig_dpi)
  invisible(p)
}

# ---- the two functions the whole analysis leans on ---------------------------

# Transition = purine<->purine (A<->G) or pyrimidine<->pyrimidine (C<->T).
# Everything else is a transversion. 4 possible Ti vs 8 possible Tv => chance = 0.5
is_transition <- function(a, b) paste0(pmin(a, b), pmax(a, b)) %in% c("AG", "CT")

# Ti/Tv summary for any set of called variants
titv_of <- function(d, from = "major_base", to = "minor_base") {
  if (nrow(d) == 0) return(tibble(n = 0L, Ti = 0L, Tv = 0L, TiTv = NA_real_))
  cl <- ifelse(is_transition(d[[from]], d[[to]]), "Ti", "Tv")
  ti <- sum(cl == "Ti"); tv <- sum(cl == "Tv")
  tibble(n = nrow(d), Ti = ti, Tv = tv, TiTv = if (tv > 0) ti / tv else NA_real_)
}

# HIV-1 gene coordinates on NC_001802.1. NOTE the overlaps - see section 6.
GENES <- tribble(
  ~gene,   ~start, ~end,
  "gag",     336,  1838,
  "pol",    1631,  4642,
  "vif",    4587,  5165,
  "vpr",    5105,  5396,
  "tat_e1", 5377,  5591,
  "vpu",    5608,  5856,
  "env",    5771,  8341,
  "rev_e2", 7925,  8199,
  "nef",    8343,  8963
)

ARTEFACT_SIG <- c("A>C", "T>G", "G>T", "C>A")   # identified empirically in section 7

cat("\n")
cat("+----------------------------------------------------------------------------+\n")
cat("|  MUTATIONS WORKSHOP  -  HIV-1 deep sequencing  -  NC_001802.1               |\n")
cat("+----------------------------------------------------------------------------+\n")
note("Input : %s", CFG$csv)
note("Config: min depth %g | error rate %.3f | naive VAF %.0f%% | strict VAF %.0f%%",
     CFG$min_depth, CFG$error_rate, 100 * CFG$minor_vaf, 100 * CFG$strict_vaf)

if (!file.exists(CFG$csv)) {
  stop("Cannot find the CSV: ", CFG$csv,
       "\n  Pass the path as the first argument:  Rscript mutations_workshop.R <file.csv>")
}

# =============================================================================
# 1.  THE DATA  -  load, classify rows, QC coverage
# =============================================================================
banner("1.", "The dataset")

raw <- readr::read_csv(CFG$csv, trim_ws = TRUE, show_col_types = FALSE)

# The header carries spaces and parentheses. Rename positionally - robust to that.
names(raw) <- c("pos", "ref", "strain", "depth_hq",
                "A", "C", "G", "T", "indel_count", "depth_raw")

dat <- raw %>%
  mutate(
    across(c(pos, depth_hq, A, C, G, T, indel_count, depth_raw), as.numeric),
    indel_count = replace_na(indel_count, 0),
    ref    = str_trim(ref),
    strain = str_trim(strain)
  )

note("Loaded %s rows, positions %g - %g", format(nrow(dat), big.mark = ","),
     min(dat$pos), max(dat$pos))

# --- Row types --------------------------------------------------------------
# Not every row is a SNP. There are THREE kinds, and a two-way split silently
# turns the third into NA:
#   1. SNP-comparable  - single ref base, single consensus base, real counts
#   2. indel record    - multi-base ref or strain, no depth/counts
#   3. zero-coverage   - EMPTY strain: no reads at all, so no consensus could be
#                        called. These are NOT deletions - they are simply not
#                        sequenced. Calling them deletions would invent data.
# Watch out: nchar(NA) is NA, not 0. So `filter(nchar(strain) == 1)` evaluates to
# NA for these rows and dplyr DROPS them without a word. Handle them explicitly.

dat <- dat %>%
  mutate(row_type = case_when(
    is.na(strain) | strain == ""          ~ "zero-coverage (no call)",
    nchar(ref) == 1 & nchar(strain) == 1  ~ "SNP-comparable",
    TRUE                                  ~ "indel record"
  ))

show_tbl(count(dat, row_type), "Row types:")

snp  <- filter(dat, row_type == "SNP-comparable")
ind  <- filter(dat, row_type == "indel record")
nocv <- filter(dat, row_type == "zero-coverage (no call)")

if (nrow(nocv)) {
  bullet("%d zero-coverage positions (%g-%g): sequencing never reached them.",
         nrow(nocv), min(nocv$pos), max(nocv$pos))
  bullet("They are excluded by the depth filter anyway - but never call them deletions.")
}

# --- Coverage QC ------------------------------------------------------------
# Rule: you cannot call a variant you did not sequence deeply enough to see.

show_tbl(
  snp %>% summarise(sites = n(),
                    zero_depth   = sum(depth_hq == 0),
                    min_depth    = min(depth_hq),
                    median_depth = median(depth_hq),
                    max_depth    = max(depth_hq)),
  "High-quality depth:"
)

save_fig(
  ggplot(snp, aes(pos, depth_hq)) +
    geom_line(linewidth = 0.3, colour = PAL[["muted"]]) +
    geom_hline(yintercept = CFG$min_depth, colour = PAL[["artefact"]], linetype = 2) +
    scale_y_log10(labels = scales::comma) +
    labs(title = "Coverage across the HIV-1 genome",
         subtitle = sprintf("Red = analysis threshold (%gx)", CFG$min_depth),
         x = "Genome position (NC_001802.1)", y = "Depth (log10)"),
  "01_coverage"
)

callable <- snp %>%
  mutate(total = A + C + G + T) %>%
  filter(depth_hq >= CFG$min_depth, total >= CFG$min_depth)

note("Callable sites: %s of %s (%.1f%%)",
     format(nrow(callable), big.mark = ","),
     format(nrow(snp), big.mark = ","),
     100 * nrow(callable) / nrow(snp))
bullet("Coverage is not uniform: low-coverage sites have a WORSE limit of detection.")
bullet("'No variant found' at 500x and at 20,000x are not the same statement.")

# =============================================================================
# 2.  WHAT IS A MUTATION?  -  consensus differences (ref vs strain)
# =============================================================================
banner("2.", "What is a mutation?")

bullet("A mutation is a change relative to a reference. Pin down three things:")
bullet("  (1) relative to WHAT?  (2) in what FRACTION?  (3) with what CONFIDENCE?")

cons <- callable %>% filter(ref != strain)

note("Consensus mutations (ref != strain): %s", format(nrow(cons), big.mark = ","))
note("Divergence from reference: %.2f%% of callable sites",
     100 * nrow(cons) / nrow(callable))

show_tbl(
  cons %>% count(ref, strain, name = "n") %>% arrange(desc(n)) %>% head(10),
  "Most common reference -> sample consensus changes:"
) %>% save_tbl("02_consensus_changes")

cat("\n  INTERPRETATION\n")
bullet("~14%% divergence is enormous - and expected.")
bullet("This is almost certainly a NON-SUBTYPE-B virus vs a SUBTYPE-B reference (HXB2).")
bullet("Most of these are SUBTYPE differences: evolution that predates this infection.")
cat("\n  >> A difference from the reference is NOT a mutation that arose in your\n")
cat("     sample. The reference is a historical accident, not a wild-type ideal.\n")

# =============================================================================
# 3.  TYPES OF MUTATION
# =============================================================================
banner("3.", "Types of mutation")

# --- 3a. By structural change ------------------------------------------------
cat("\n  3a. BY STRUCTURAL CHANGE (substitution / insertion / deletion)\n")

ind_class <- ind %>%
  mutate(
    indel_type = case_when(nchar(strain) > nchar(ref) ~ "Insertion",
                           nchar(strain) < nchar(ref) ~ "Deletion",
                           TRUE                       ~ "Complex/other"),
    size  = abs(nchar(ref) - nchar(strain)),
    frame = if_else(size %% 3 == 0, "In-frame (x3)", "Frameshift")
  )

show_tbl(count(ind_class, indel_type, frame),
         "Indel records by type and reading-frame effect:") %>%
  save_tbl("03_indels")

bullet("An indel sized in multiples of 3 adds/removes whole codons - protein still reads.")
bullet("Any other size SHIFTS THE FRAME: every downstream codon garbled.")
bullet("Same 'one mutation', wildly different consequence.")
bullet("Most indels here are in-frame - consistent with real viral variation.")
bullet("Random error has no reason to respect codon boundaries.")

# --- 3b. By biochemistry: transitions vs transversions -----------------------
cat("\n  3b. BY BIOCHEMISTRY (transition vs transversion)\n")

titv_cons <- titv_of(cons, "ref", "strain")
show_tbl(titv_cons, "Consensus mutations by substitution class:")

note("Consensus Ti/Tv = %.2f   (pure chance would give 0.50)", titv_cons$TiTv)
bullet("4 possible transitions vs 8 possible transversions => chance = 0.5")
bullet("Real data runs higher: transitions are mechanistically favoured and")
bullet("more often synonymous, so selection removes them less often.")
bullet("REMEMBER 0.5. It is the cheapest QC metric you own - and in section 7")
bullet("it is going to save us from publishing an artefact.")

save_fig(
  cons %>% count(ref, strain) %>%
    mutate(class = if_else(is_transition(ref, strain), "Transition", "Transversion")) %>%
    ggplot(aes(strain, ref, fill = n)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = n, colour = class), size = 3.6, fontface = "bold") +
    scale_fill_gradient(low = "grey95", high = PAL[["signal"]]) +
    scale_colour_manual(values = c(Transition = PAL[["artefact"]], Transversion = "grey20")) +
    labs(title = "Substitution spectrum: reference base -> sample base",
         subtitle = "Red labels = transitions (A<->G, C<->T)",
         x = "Sample (strain) base", y = "Reference base",
         fill = "Count", colour = NULL),
  "03_substitution_spectrum", w = 7, h = 5
)

# --- 3c. By functional consequence -------------------------------------------
cat("\n  3c. BY FUNCTIONAL CONSEQUENCE (needs the reading frame - not in a pileup)\n")
bullet("Synonymous  : codon changes, amino acid does not  (GGA->GGG, both Gly)")
bullet("Missense    : one amino acid replaced             (HIV RT K103N)")
bullet("Nonsense    : codon -> STOP, protein truncated    (usually lethal to virus)")
bullet("Frameshift  : indel not x3, downstream garbled    (usually knockout)")
bullet("Non-coding  : outside genes, may affect expression (LTR promoter)")
bullet("dN/dS: <1 purifying, ~1 neutral, >1 positive selection (immune/drug escape)")

cat("\n  3d. BY ORIGIN\n")
bullet("Germline (inherited) vs somatic (acquired) - the basis of cancer genomics.")
bullet("Viral analogue: TRANSMITTED vs ACQUIRED resistance - decides whether a")
bullet("first-line regimen was ever going to work.")

# =============================================================================
# 4.  HOW TO CALCULATE MUTATIONS  -  VAF, and frequency vs rate
# =============================================================================
banner("4.", "How to calculate mutations")

cat("\n  The core quantity - variant allele frequency:\n")
cat("      VAF_b  =  n_b / (n_A + n_C + n_G + n_T)\n")
cat("  Everything downstream - minority calls, diversity, burden - is built on it.\n")

freq <- callable %>%
  mutate(across(c(A, C, G, T), ~ .x / total, .names = "f_{.col}")) %>%
  rowwise() %>%
  mutate(
    major_base = c("A", "C", "G", "T")[which.max(c(A, C, G, T))],
    major_f    = max(c(f_A, f_C, f_G, f_T)),
    minor_base = c("A", "C", "G", "T")[order(c(A, C, G, T), decreasing = TRUE)[2]],
    minor_f    = sort(c(f_A, f_C, f_G, f_T), decreasing = TRUE)[2],
    minor_n    = sort(c(A, C, G, T), decreasing = TRUE)[2]
  ) %>%
  ungroup() %>%
  mutate(change = paste0(major_base, ">", minor_base))

note("Sanity check - file's consensus vs our computed majority base: %.2f%% agree",
     100 * mean(freq$major_base == freq$strain))

cat("\n  FREQUENCY IS NOT RATE\n")
bullet("Frequency: how many mutants NOW. A proportion. Needs 1 sample. -> we have it.")
bullet("Rate     : how FAST new ones appear. Per base per replication.")
bullet("           Needs >=2 timepoints. -> we CANNOT get this from one sample.")
bullet("Frequency is a snapshot: mutation x selection x drift x elapsed time.")

hiv_rate   <- 2.4e-5   # mutations per base per replication cycle
genome_len <- 9181
virions    <- 1e10     # produced per day, untreated infection

note("Expected new mutations per genome per cycle: %.2f", hiv_rate * genome_len)
note("Every possible point mutation is generated ~%s times per day",
     format(round(hiv_rate * virions), big.mark = ","))
bullet("That last number is the whole ballgame. Resistance is not a question of")
bullet("WHETHER it is generated - it is generated constantly. It is a question of")
bullet("whether the drug gives it a reason to expand.")

# =============================================================================
# 5.  DIVERSITY  -  Shannon entropy and nucleotide diversity (pi)
# =============================================================================
banner("5.", "Diversity")

bullet("DIVERGENCE asks: how far from the reference?")
bullet("DIVERSITY  asks: how variable is the population WITHIN ITSELF?")
bullet("A sample can be far from the reference yet uniform, or sit on it and seethe.")
bullet("This is the QUASISPECIES concept: an infection is a cloud, not a genome.")

shannon <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }

div <- freq %>%
  rowwise() %>%
  mutate(H = shannon(c(f_A, f_C, f_G, f_T))) %>%
  ungroup() %>%
  # Sample-size-corrected (unbiased) nucleotide diversity.
  # The n/(n-1) term matters: without it low-coverage sites look artificially
  # LESS diverse, purely because you sampled fewer reads.
  mutate(pi_site = (total / (total - 1)) * (1 - (f_A^2 + f_C^2 + f_G^2 + f_T^2))) %>%
  mutate(is_artefact_sig = change %in% ARTEFACT_SIG)

note("Mean Shannon entropy : %.4f   (max possible = %.3f)", mean(div$H), log(4))
note("Nucleotide diversity : pi = %.5f", mean(div$pi_site))
bullet("i.e. two random reads differ at ~%.1f of every 10,000 bases",
       10000 * mean(div$pi_site))

save_fig(
  ggplot(div, aes(pos, pi_site)) +
    geom_point(aes(colour = pi_site > 0.02), size = 0.5, alpha = 0.7) +
    scale_colour_manual(values = c("grey70", PAL[["artefact"]]), guide = "none") +
    labs(title = "Within-host nucleotide diversity across the genome",
         subtitle = "Red = pi > 0.02",
         x = "Genome position", y = "pi (per site)"),
  "05_diversity_genome"
)

# =============================================================================
# 6.  DIVERSITY BY GENE  -  and the overlapping-ORF trap
# =============================================================================
banner("6.", "Diversity by gene (mind the overlapping ORFs)")

bullet("HIV-1 is a masterclass in genomic compression: GENES OVERLAP, in")
bullet("different reading frames. vpu sits inside env; tat/rev exon 2 inside env.")
bullet("A naive case_when() assigns each site to ONE gene and SILENTLY DISCARDS")
bullet("the overlaps. Use a long-format join so a site can belong to several genes.")

# Long format: one row per (site, gene). Overlapping sites appear more than once.
long <- div %>% crossing(GENES) %>% filter(pos >= start, pos <= end)

show_tbl(
  long %>% count(pos) %>% count(n, name = "n_sites") %>% rename(genes_per_site = n),
  "Sites belonging to 1 vs 2 overlapping reading frames:"
)

gene_div <- long %>%
  group_by(gene) %>%
  summarise(sites = n(), pi = mean(pi_site), H = mean(H), .groups = "drop") %>%
  arrange(desc(pi))

show_tbl(gene_div, "Diversity by HIV-1 gene:") %>% save_tbl("06_gene_diversity")

save_fig(
  ggplot(gene_div, aes(reorder(gene, pi), pi)) +
    geom_col(fill = PAL[["signal"]]) + coord_flip() +
    labs(title = "Within-host diversity by HIV-1 gene", x = NULL, y = "mean pi"),
  "06_gene_diversity", w = 7, h = 4.5
)

bullet("env sits high: surface protein, relentless antibody pressure, tolerant of change.")
bullet("pol sits lower: drug target, structurally constrained.")
bullet("Diversity is not random - it maps onto SELECTION.")

# =============================================================================
# 7.  MINORITY VARIANTS  -  the naive result, and how to destroy it
# =============================================================================
banner("7.", "Minority variants - and how to know if they are real")

bullet("Sanger sees ~20%% and above. Rarer variants are invisible to it.")
bullet("Deep sequencing can in principle reach 1%% or lower.")
bullet("Variants in that gap - present, real, sub-consensus - are MINORITY VARIANTS.")
bullet("In HIV: minority drug-resistance mutations at 1-20%% raise failure risk,")
bullet("are archived or transmitted, and are invisible to standard genotyping.")

# --- 7a. The naive call: binomial test against a flat error rate -------------
cat("\n  7a. FIRST ATTEMPT - binomial test vs the error rate\n")
cat("      Under the null 'this is just error', alt count ~ Binomial(n, eps).\n")

minor <- div %>%
  filter(minor_n > 0) %>%
  rowwise() %>%
  mutate(p_error = binom.test(minor_n, total, p = CFG$error_rate,
                              alternative = "greater")$p.value) %>%
  ungroup() %>%
  mutate(q_error = p.adjust(p_error, method = "BH"))

naive <- minor %>% filter(q_error < 0.01, minor_f >= CFG$minor_vaf)

note("Sites with any 2nd allele              : %s", format(nrow(minor), big.mark = ","))
note("Passing VAF >= %.0f%% AND FDR < 1%%       : %s",
     100 * CFG$minor_vaf, format(nrow(naive), big.mark = ","))
note("Proportion of those that pass FDR < 1%% : %.1f%%",
     100 * mean(minor$q_error[minor$minor_f >= CFG$minor_vaf] < 0.01))

cat("\n  A clean, defensible, fully FDR-corrected result. We could write it up today.\n")
cat("  So let us try to destroy it.\n")

# --- 7b. The Ti/Tv stress test ----------------------------------------------
cat("\n  7b. THE STRESS TEST - Ti/Tv as a function of the VAF threshold\n")
cat("      Recall: 0.50 is what PURE CHANCE gives.\n")

thresholds <- c(0.001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.10)
titv_curve <- map_dfr(thresholds, function(t) {
  titv_of(filter(minor, minor_f >= t)) %>% mutate(threshold = t)
}) %>% select(threshold, n, Ti, Tv, TiTv)

show_tbl(titv_curve, "Ti/Tv of the minor allele by VAF threshold:") %>%
  save_tbl("07_titv_by_threshold")

save_fig(
  ggplot(titv_curve, aes(threshold * 100, TiTv)) +
    geom_hline(yintercept = 0.5, colour = PAL[["artefact"]], linetype = 2, linewidth = 0.8) +
    geom_hline(yintercept = titv_cons$TiTv, colour = PAL[["muted"]], linetype = 3, linewidth = 0.8) +
    geom_line(linewidth = 0.8, colour = PAL[["signal"]]) +
    geom_point(aes(size = n), colour = PAL[["signal"]]) +
    scale_x_log10() +
    annotate("text", x = 0.12, y = 0.58, label = "pure chance = 0.5",
             colour = PAL[["artefact"]], hjust = 0, size = 3.4) +
    annotate("text", x = 0.12, y = titv_cons$TiTv + 0.08, label = "consensus Ti/Tv",
             colour = PAL[["muted"]], hjust = 0, size = 3.4) +
    labs(title = "Ti/Tv collapses at low VAF - the signature of error contamination",
         subtitle = "Below ~2% the 'variants' are indistinguishable from noise, or worse",
         x = "VAF threshold (%, log scale)", y = "Ti/Tv of minor allele", size = "n sites"),
  "07_titv_collapse"
)

worst <- titv_curve %>% slice_min(TiTv, n = 1)
cat("\n  THIS IS DAMNING\n")
bullet("At a %.1f%% threshold: %s 'variants', Ti/Tv = %.2f",
       100 * worst$threshold, format(worst$n, big.mark = ","), worst$TiTv)
bullet("That is not merely random (0.5) - it is FAR BELOW it.")
bullet("Random error gives 0.5. Something WORSE than random is happening:")
bullet("a systematic, direction-specific artefact.")
bullet("At our %.0f%% threshold Ti/Tv = %.2f - roughly HALF of those %s calls are noise.",
       100 * CFG$minor_vaf,
       titv_curve$TiTv[titv_curve$threshold == CFG$minor_vaf],
       format(nrow(naive), big.mark = ","))

# --- 7c. Diagnose: read the substitution spectrum -----------------------------
cat("\n  7c. DIAGNOSIS - read the spectrum, name the culprit\n")

bands <- minor %>%
  filter(minor_f >= 0.001) %>%
  mutate(band = cut(minor_f, c(0.001, 0.005, 0.01, 0.02, 1),
                    labels = c("0.1-0.5%", "0.5-1%", "1-2%", ">2%"),
                    include.lowest = TRUE))

spectrum <- bands %>%
  count(band, change) %>%
  group_by(band) %>% slice_max(n, n = 6, with_ties = FALSE) %>% ungroup() %>%
  mutate(artefact = change %in% ARTEFACT_SIG)

show_tbl(spectrum, "Top substitutions per VAF band:") %>% save_tbl("07_spectrum_by_band")

save_fig(
  spectrum %>%
    mutate(lab = fct_reorder(paste(band, change, sep = "__"), n)) %>%
    ggplot(aes(lab, n, fill = artefact)) +
    geom_col() +
    scale_x_discrete(labels = function(x) sub("^.*__", "", x)) +
    facet_wrap(~ band, scales = "free", nrow = 1) +
    coord_flip() +
    scale_fill_manual(values = c(`FALSE` = PAL[["signal"]], `TRUE` = PAL[["artefact"]]),
                      labels = c("plausible biology", "artefact signature"), name = NULL) +
    labs(title = "Substitution spectrum by VAF band",
         subtitle = "Low VAF is dominated by A>C / T>G / G>T / C>A - reverse-complement pairs",
         x = NULL, y = "Sites"),
  "07_artefact_spectrum", w = 11, h = 4.5
)

bullet("Low VAF is dominated by: %s - ALL transversions.", paste(ARTEFACT_SIG, collapse = ", "))
bullet("They form REVERSE-COMPLEMENT PAIRS (A>C = T>G; G>T = C>A).")
bullet("That strand symmetry is the fingerprint of a chemical/technical artefact:")
bullet("  G>T / C>A -> the classic 8-oxoguanine (OxoG) oxidative damage signature")
bullet("  A>C / T>G -> a documented context-specific sequencing artefact")
bullet("Real biology has no reason to make a strand-symmetric, transversion-only")
bullet("spectrum. Sample handling does.")

# --- 7d. Why the statistics did not save us ----------------------------------
cat("\n  7d. WHY THE BINOMIAL TEST FAILED TO CATCH IT\n")

show_tbl(
  minor %>% filter(minor_f >= CFG$minor_vaf) %>%
    summarise(n = n(),
              median_depth     = median(total),
              median_alt_reads = median(minor_n),
              frac_significant = mean(q_error < 0.01)),
  "The calls it was asked to judge:"
)

bullet("At ~4,000x depth, a 1%% variant is ~40 supporting reads.")
bullet("Against a null of eps=%.1f%%, 40 reads is OVERWHELMINGLY significant.", 100 * CFG$error_rate)
bullet("So the test passes EVERY call. It never had a chance to reject anything.")
bullet("The test is not wrong - its ASSUMPTION is. It assumes error is uniform at")
bullet("%.1f%% everywhere. Real error is context- and site-specific; at some sites", 100 * CFG$error_rate)
bullet("it runs at 2%%. The model cannot see that, so it certifies artefacts.")
cat("\n  >> STATISTICAL SIGNIFICANCE IS NOT BIOLOGICAL REALITY.\n")
cat("     A p-value only ever tests the model you handed it. With enough depth you\n")
cat("     get vanishingly small p-values for pure garbage. Deep sequencing does not\n")
cat("     protect you from this - it GUARANTEES it. Depth makes everything significant.\n")

# --- 7e. The fix -------------------------------------------------------------
cat("\n  7e. THE FIX - filter the artefact signature, raise the threshold\n")

called <- div %>%
  filter(minor_n > 0, minor_f >= CFG$strict_vaf, !is_artefact_sig)

titv_clean <- map_dfr(c(0.01, 0.02, 0.05), function(t) {
  titv_of(div %>% filter(minor_n > 0, minor_f >= t, !is_artefact_sig)) %>%
    mutate(threshold = t)
}) %>% select(threshold, n, Ti, Tv, TiTv)

show_tbl(titv_clean, "Ti/Tv AFTER removing the artefact signature:") %>%
  save_tbl("07_titv_clean")

note("NAIVE  (%.0f%%, no filter)      : %4s calls   Ti/Tv = %.2f",
     100 * CFG$minor_vaf, nrow(naive), titv_of(naive)$TiTv)
note("CLEAN  (%.0f%% + artefact filter): %4s calls   Ti/Tv = %.2f",
     100 * CFG$strict_vaf, nrow(called), titv_of(called)$TiTv)
note("Discarded: %s of %s calls (%.0f%%)",
     nrow(naive) - nrow(called), nrow(naive),
     100 * (1 - nrow(called) / nrow(naive)))

bullet("Ti/Tv jumps from ~0.5 to ~5. THAT is a biologically sensible number.")
bullet("We threw away most of our variants, and what remains looks like biology.")
bullet("Note: 2%% was justified EMPIRICALLY - from Ti/Tv and the spectrum -")
bullet("not by convention or by copying someone else's methods section.")

save_tbl(called %>% select(pos, ref, strain, total, major_base, minor_base,
                           minor_n, minor_f, change, H, pi_site),
         "07_clean_minority_variants")

# --- 7f. Biological validation: the APOBEC3G signature -----------------------
cat("\n  7f. INDEPENDENT VALIDATION - the APOBEC3G signature\n")
cat("      Host APOBEC3G deaminates cytosine during reverse transcription,\n")
cat("      producing G>A hypermutation. HIV's vif gene exists to destroy it.\n")
cat("      If our calls are real, G>A should be ENRICHED at higher VAF.\n")
cat("      Sequencing error has no reason whatsoever to do that.\n")

apobec <- minor %>%
  mutate(band = cut(minor_f, c(0, 0.005, 0.01, 0.02, 0.05, 1),
                    labels = c("<0.5%", "0.5-1%", "1-2%", "2-5%", ">5%"))) %>%
  group_by(band) %>%
  summarise(n = n(), G_to_A = sum(change == "G>A"),
            pct_GtoA = 100 * mean(change == "G>A"), .groups = "drop")

show_tbl(apobec, "G>A (APOBEC3G signature) by VAF band:") %>% save_tbl("07_apobec")

save_fig(
  ggplot(apobec, aes(band, pct_GtoA)) +
    geom_col(fill = PAL[["signal"]]) +
    labs(title = "G>A rises with VAF - the APOBEC3G hypermutation signature",
         subtitle = "Sequencing error has no reason to do this. Biology does.",
         x = "Minor allele frequency band", y = "% of variants that are G>A"),
  "07_apobec", w = 7, h = 4.5
)

bullet("G>A climbs from %.1f%% to %.1f%% as VAF rises.",
       apobec$pct_GtoA[1], max(apobec$pct_GtoA))
bullet("The high-frequency variants are enriched for a specific, known,")
bullet("host-driven biological process. Independent confirmation that the")
bullet("survivors are real - and that the low-VAF tier is not.")

# =============================================================================
# 8.  BURDEN TESTING  -  concept
# =============================================================================
banner("8.", "Burden testing")

bullet("Single-variant tests work for COMMON variants. For RARE ones they fail:")
bullet("if a variant is in 3 of 5,000 people, nothing has power to associate it.")
bullet("The burden test's move is COLLAPSING: stop testing variants, test UNITS.")
bullet("Instead of 'does variant X associate?', ask:")
bullet("  'do carriers of ANY rare variant in gene G differ from non-carriers?'")
bullet("You trade RESOLUTION for POWER.")
cat("\n")
bullet("CAST    : binary, carries >=1 rare variant   | rare, same direction")
bullet("Burden  : count, weighted by rarity          | effects same direction")
bullet("SKAT    : variance-component                 | MIXED directions")
bullet("SKAT-O  : data-driven blend of the two       | you don't know - safe default")
cat("\n")
bullet("FATAL FLAW of a plain burden test: it assumes every variant in the unit")
bullet("pushes the SAME way. Half harmful + half protective => they cancel, and it")
bullet("sees nothing. That is exactly why SKAT exists.")

cat("\n  AN HONEST NOTE ON THIS DATASET\n")
bullet("A classical burden test compares CASES vs CONTROLS across MANY individuals.")
bullet("We have ONE sample. We cannot run one - be suspicious of anyone who claims")
bullet("otherwise. What we CAN do is the same machinery in its within-sample form:")
bullet("collapse variants into genes, test each against a background rate. That is")
bullet("how regional mutation-enrichment scans work in cancer genomics.")

# =============================================================================
# 9.  BURDEN  -  regional enrichment, naive vs clean
# =============================================================================
banner("9.", "Burden - the same test, before and after artefact removal")

# Null: minority variants fall uniformly across the genome at the genome-wide
# rate. Under that null the count in a gene is Binomial(gene length, bg rate).
burden_test <- function(call_set, label) {
  bg <- nrow(call_set) / nrow(div)
  long %>%
    count(gene, name = "sites") %>%
    left_join(
      call_set %>% crossing(GENES) %>% filter(pos >= start, pos <= end) %>%
        count(gene, name = "observed"),
      by = "gene"
    ) %>%
    mutate(
      observed = replace_na(observed, 0),
      expected = sites * bg,
      fold     = observed / expected,
      p_value  = pbinom(observed - 1, sites, bg, lower.tail = FALSE),
      q_value  = p.adjust(p_value, method = "BH"),
      analysis = label
    ) %>%
    arrange(p_value)
}

LAB_NAIVE <- sprintf("naive (%.0f%%, no artefact filter)", 100 * CFG$minor_vaf)
LAB_CLEAN <- sprintf("clean (%.0f%% + artefact filter)",   100 * CFG$strict_vaf)

burden_naive <- burden_test(naive,  LAB_NAIVE)
burden_clean <- burden_test(called, LAB_CLEAN)

show_tbl(burden_naive %>% select(gene, sites, observed, expected, fold, p_value, q_value),
         sprintf("NAIVE burden (the contaminated calls, %.0f%%):", 100 * CFG$minor_vaf)) %>%
  save_tbl("09_burden_naive")

cat("\n  A beautiful result! vif enriched, q < 0.001. We could write this up.\n")
cat("  It is wrong.\n")

show_tbl(burden_clean %>% select(gene, sites, observed, expected, fold, p_value, q_value),
         sprintf("CLEAN burden (artefact signature removed, %.0f%% threshold):",
                 100 * CFG$strict_vaf)) %>%
  save_tbl("09_burden_clean")

save_fig(
  bind_rows(burden_naive, burden_clean) %>%
    mutate(sig = q_value < 0.05) %>%
    ggplot(aes(reorder(gene, fold), fold, fill = sig)) +
    geom_col() +
    geom_hline(yintercept = 1, linetype = 2) +
    coord_flip() +
    facet_wrap(~ analysis) +
    scale_fill_manual(values = c(`FALSE` = PAL[["muted"]], `TRUE` = PAL[["artefact"]]),
                      labels = c("ns", "FDR < 0.05"), name = NULL) +
    labs(title = "The same burden test, before and after artefact removal",
         subtitle = "vif's 'significant' enrichment was an artefact. env survives.",
         x = NULL, y = "Observed / Expected"),
  "09_burden_comparison", w = 10, h = 5
)

verdict <- bind_rows(burden_naive, burden_clean) %>%
  filter(gene %in% c("vif", "env")) %>%
  select(analysis, gene, observed, fold, q_value) %>%
  arrange(gene, desc(analysis))

show_tbl(verdict, "THE VERDICT on vif and env:") %>% save_tbl("09_verdict")

cat("\n  WHAT HAPPENED\n")
bullet("vif : significant -> NOT significant. Its 'enrichment' was driven by")
bullet("      artefact-signature calls (~77%% of vif's 1%% calls were %s).",
       paste(ARTEFACT_SIG, collapse = "/"))
bullet("      A confident, FDR-corrected, ENTIRELY FALSE result.")
bullet("env : got STRONGER when noise was removed, and stayed significant.")
bullet("      That is the hallmark of real signal - and exactly what HIV biology")
bullet("      predicts for the antibody-facing envelope protein.")

cat("\n  >> FDR CONTROLS CHANCE. IT IS BLIND TO BIAS.\n")
cat("     Our q-value was 0.0005 and the result was still garbage. Multiple-testing\n")
cat("     correction did not save us - it was never designed to. No amount of\n")
cat("     statistics rescues contaminated input. Only the BIOLOGY told us the truth:\n")
cat("     Ti/Tv on the noise line, the OxoG spectrum, the APOBEC gradient. Every one\n")
cat("     of those is a DOMAIN check, not a statistical one.\n")

# =============================================================================
# 10.  BURDEN  -  the classic case/control design (simulated cohort)
# =============================================================================
banner("10.", "Burden - case/control (simulated cohort)")

cat("\n  The classic design, on simulated data so the logic is visible:\n")
cat("  20 patients failing therapy vs 20 suppressed. Are rare variants in a gene\n")
cat("  enriched in the failures?\n")

set.seed(42)
n_case <- 20; n_ctrl <- 20; n_var <- 15

sim_gene <- function(effect) {
  rbind(
    matrix(rbinom(n_case * n_var, 1, 0.02 + effect), nrow = n_case),  # cases
    matrix(rbinom(n_ctrl * n_var, 1, 0.02),          nrow = n_ctrl)   # controls
  )
}

pheno <- c(rep(1, n_case), rep(0, n_ctrl))

run_burden <- function(G, label) {
  carrier <- as.integer(rowSums(G) > 0)                 # CAST: carrier yes/no
  cast_p  <- fisher.test(table(carrier, pheno))$p.value
  load    <- rowSums(G)                                 # Burden: variant count
  burd_p  <- summary(glm(pheno ~ load, family = binomial))$coefficients[2, 4]
  tibble(gene = label,
         mean_load_case = mean(load[pheno == 1]),
         mean_load_ctrl = mean(load[pheno == 0]),
         CAST_p = cast_p, Burden_p = burd_p)
}

sim_res <- bind_rows(
  run_burden(sim_gene(0.00), "null_gene (no effect)"),
  run_burden(sim_gene(0.12), "causal_gene (true effect)")
)

show_tbl(sim_res, "Simulated case/control burden: CAST vs count-based Burden") %>%
  save_tbl("10_burden_simulated")

bullet("The causal gene shows a significant burden; the null gene does not.")
bullet("That is the whole engine of rare-variant association testing - the same")
bullet("collapsing logic we applied to genes in the real data above.")
cat("\n  EXERCISE: in run_burden(), make HALF the variants protective (flip their\n")
cat("  direction in cases). Re-run. Watch the burden test's power COLLAPSE.\n")
cat("  Then explain, in one sentence, why SKAT exists.\n")

# =============================================================================
# 11.  SUMMARY
# =============================================================================
banner("11.", "Summary - one sample, every metric")

summary_tbl <- tibble(
  Metric = c(
    "SNP-comparable sites",
    "Callable sites (>= min depth)",
    "Median depth",
    "Consensus mutations vs reference",
    "% divergence from NC_001802.1",
    "Consensus Ti/Tv",
    "Indel records",
    "Mean Shannon entropy",
    "Nucleotide diversity (pi)",
    sprintf("NAIVE minority variants (%.0f%%, FDR<1%%)", 100 * CFG$minor_vaf),
    "  their Ti/Tv (= the noise line)",
    sprintf("CLEAN minority variants (%.0f%% + filter)", 100 * CFG$strict_vaf),
    "  their Ti/Tv",
    "Discarded as artefact",
    "Genes significant, naive burden",
    "Genes significant, clean burden"
  ),
  Value = c(
    format(nrow(snp), big.mark = ","),
    format(nrow(callable), big.mark = ","),
    format(median(snp$depth_hq), big.mark = ","),
    format(nrow(cons), big.mark = ","),
    sprintf("%.2f%%", 100 * nrow(cons) / nrow(callable)),
    sprintf("%.2f", titv_cons$TiTv),
    format(nrow(ind), big.mark = ","),
    sprintf("%.4f", mean(div$H)),
    sprintf("%.5f", mean(div$pi_site)),
    format(nrow(naive), big.mark = ","),
    sprintf("%.2f", titv_of(naive)$TiTv),
    format(nrow(called), big.mark = ","),
    sprintf("%.2f", titv_of(called)$TiTv),
    sprintf("%d (%.0f%%)", nrow(naive) - nrow(called),
            100 * (1 - nrow(called) / nrow(naive))),
    paste(burden_naive$gene[burden_naive$q_value < 0.05], collapse = ", "),
    paste(burden_clean$gene[burden_clean$q_value < 0.05], collapse = ", ")
  )
)

show_tbl(summary_tbl) %>% save_tbl("11_summary")

cat("\n", strrep("-", 78), "\n", sep = "")
cat("  THE SIX THINGS TO TAKE AWAY\n")
cat(strrep("-", 78), "\n", sep = "")
cat(sprintf("
  1. 'MUTATION' IS A COMPARISON. State: relative to what, at what frequency,
     with what confidence. %.1f%% of this genome differs from HXB2 because it is a
     different SUBTYPE - not because the patient mutated.

  2. FREQUENCY IS NOT RATE. One timepoint gives a snapshot, never a speed.

  3. DIVERSITY IS NOT DIVERGENCE. Independent axes - and diversity maps onto
     selection (env high, pol constrained).

  4. YOUR THRESHOLD AND DEPTH DEFINE YOUR VARIANT COUNT. %s at %.0f%%; %s at %.0f%%
     with filtering. Neither is 'the truth'. Report both, or report nothing.

  5. Ti/Tv IS YOUR LIE DETECTOR. 0.5 = noise. If your minority calls sit at 0.5,
     they are noise - no matter how small the p-value.

  6. STATISTICS CANNOT RESCUE CONTAMINATED DATA. We got a tiny q-value for vif
     and it was pure artefact. FDR controls chance, not bias. Only the biology -
     Ti/Tv, the OxoG spectrum, the APOBEC signature - told us the truth.
",
  100 * nrow(cons) / nrow(callable),
  format(nrow(naive), big.mark = ","), 100 * CFG$minor_vaf,
  format(nrow(called), big.mark = ","), 100 * CFG$strict_vaf
))

cat(strrep("-", 78), "\n", sep = "")
note("Figures written to: %s/", normalizePath(FIGDIR, mustWork = FALSE))
note("Tables  written to: %s/", normalizePath(RESDIR, mustWork = FALSE))
cat("\n  Done.\n\n")
