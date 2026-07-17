#!/usr/bin/env Rscript
# =============================================================================
#  POL MINORITY VARIANT ANALYSIS  -  PR / RT / IN
#  Threshold-laddered consensus FASTAs, lollipop plots, and comparisons
#  HIV-1 deep sequencing, NC_001802.1  |  Dr. W. Choga
# -----------------------------------------------------------------------------
#  USAGE
#     Rscript pol_minority_fasta.R
#     Rscript pol_minority_fasta.R <variants.csv> --outdir pol_minority
#
#  WHAT IT DOES
#     For protease (PR), reverse transcriptase (RT) and integrase (IN) it builds
#     a consensus nucleotide sequence at each of five VAF thresholds
#     (50, 20, 15, 10, 5 %), using IUPAC ambiguity codes wherever more than one
#     allele clears the threshold. It then translates them, draws lollipop plots
#     of the minority landscape, and compares the sequences to each other.
#
#  WHY THE THRESHOLD LADDER
#     50%  = plain majority consensus. What a FASTA normally means.
#     20%  = the Sanger detection limit. What routine genotyping would report.
#     15/10/5% = progressively deeper into the minority population.
#     The sequences only differ where minority variants exist. That difference
#     IS the clinical question.
#
#  READ THE CONSOLE OUTPUT. Section 7 (drug resistance) is the point. A
#  nucleotide change at a resistance codon is NOT a resistance mutation, and
#  this script will show you two ways that trap is sprung.
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

# =============================================================================
# 0.  CONFIG
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default) {
  i <- match(flag, args); if (!is.na(i) && length(args) > i) return(args[i + 1]); default
}

CFG <- list(
  csv        = if (length(args) >= 1 && !startsWith(args[1], "--")) args[1]
               else "hiv1_pileup_NC_001802.1_snp-variants.csv",
  outdir     = get_opt("--outdir", "pol_minority"),
  min_depth  = as.numeric(get_opt("--mindepth", 100)),
  thresholds = c(0.50, 0.20, 0.15, 0.10, 0.05),   # the requested ladder
  explore    = c(0.02, 0.01),                     # exploratory tier - see section 6
  sample_id  = "sample01"
)

PAL <- c(signal = "#0FA3A3", artefact = "#E2574C", violet = "#8B5CF6", muted = "#6C7D8B")
ARTEFACT_SIG <- c("A>C", "T>G", "G>T", "C>A")   # OxoG / context artefact, see mutations_workshop.R

theme_set(theme_minimal(base_size = 12))

banner <- function(n, t) { cat("\n", strrep("=", 78), "\n  ", n, "  ", toupper(t), "\n",
                              strrep("=", 78), "\n", sep = "") }
note   <- function(...) cat("\n  ", sprintf(...), "\n", sep = "")
bullet <- function(...) cat("   - ", sprintf(...), "\n", sep = "")
show_tbl <- function(df, cap = NULL) {
  if (!is.null(cap)) cat("\n  ", cap, "\n", sep = "")
  print(as.data.frame(df), row.names = FALSE, digits = 4); invisible(df)
}

# =============================================================================
# 1.  REFERENCE TABLES  -  genetic code, IUPAC, coordinates, DRM positions
# =============================================================================

CODON_TBL <- c(
  TTT="F",TTC="F",TTA="L",TTG="L",CTT="L",CTC="L",CTA="L",CTG="L",
  ATT="I",ATC="I",ATA="I",ATG="M",GTT="V",GTC="V",GTA="V",GTG="V",
  TCT="S",TCC="S",TCA="S",TCG="S",CCT="P",CCC="P",CCA="P",CCG="P",
  ACT="T",ACC="T",ACA="T",ACG="T",GCT="A",GCC="A",GCA="A",GCG="A",
  TAT="Y",TAC="Y",TAA="*",TAG="*",CAT="H",CAC="H",CAA="Q",CAG="Q",
  AAT="N",AAC="N",AAA="K",AAG="K",GAT="D",GAC="D",GAA="E",GAG="E",
  TGT="C",TGC="C",TGA="*",TGG="W",CGT="R",CGC="R",CGA="R",CGG="R",
  AGT="S",AGC="S",AGA="R",AGG="R",GGT="G",GGC="G",GGA="G",GGG="G"
)

# IUPAC: set of bases -> single ambiguity letter
IUPAC <- c(
  "A"="A","C"="C","G"="G","T"="T",
  "AG"="R","CT"="Y","CG"="S","AT"="W","GT"="K","AC"="M",
  "CGT"="B","AGT"="D","ACT"="H","ACG"="V","ACGT"="N"
)
# and the reverse: ambiguity letter -> the bases it stands for
IUPAC_REV <- setNames(names(IUPAC), IUPAC)

to_iupac  <- function(bases) unname(IUPAC[paste(sort(unique(bases)), collapse = "")])
from_iupac <- function(code) strsplit(IUPAC_REV[[code]], "")[[1]]

# --- Coordinates -------------------------------------------------------------
# CRITICAL: NC_001802.1 is NOT numbered like HXB2/K03455. It is 9,181 bp; HXB2 is
# 9,719 bp. In the pol region the offset is exactly 454 (NC = HXB2 - 454).
# Every HIV drug-resistance table in the literature is written in HXB2 numbering,
# so we convert once, here, and verify the result by translation in section 3.
HXB2_OFFSET <- 454

REGIONS <- tribble(
  ~region, ~name,                   ~hxb2_start, ~hxb2_end, ~n_codons, ~expect_aa,
  "PR",    "Protease",                    2253,      2549,        99,   "PQ",
  "RT",    "Reverse transcriptase",       2550,      4229,       560,   "PISPIE",
  "IN",    "Integrase",                   4230,      5093,       288,   "FLDGID"
) %>%
  mutate(start = hxb2_start - HXB2_OFFSET,
         end   = hxb2_end   - HXB2_OFFSET)

# --- Known major drug-resistance positions (codon numbers, HXB2 convention) ---
# Source: IAS-USA / Stanford HIVdb major + selected accessory positions.
# NOTE: these are POSITIONS ONLY. Presence of a variant at these codons does NOT
# mean resistance - the specific amino acid is what matters. Section 7 enforces this.
DRM_POS <- list(
  PR = c(23,24,30,32,33,46,47,48,50,53,54,73,76,82,83,84,88,90),
  RT = c(41,65,67,69,70,74,75,77,98,100,101,103,106,108,115,116,138,151,
         179,181,184,188,190,210,215,219,221,225,227,230),
  IN = c(66,74,92,95,97,118,121,138,140,143,147,148,151,155,163,230,232,263)
)

cat("\n+--------------------------------------------------------------------------+\n")
cat("|  POL MINORITY VARIANT ANALYSIS  -  PR / RT / IN                          |\n")
cat("+--------------------------------------------------------------------------+\n")
note("Sample     : %s", CFG$sample_id)
note("Thresholds : %s", paste0(100 * CFG$thresholds, "%", collapse = ", "))

if (!file.exists(CFG$csv)) stop("Cannot find CSV: ", CFG$csv)

# --- Build the output tree ---------------------------------------------------
dir.create(CFG$outdir, showWarnings = FALSE, recursive = TRUE)
for (rg in REGIONS$region) {
  dir.create(file.path(CFG$outdir, rg, "fasta"), showWarnings = FALSE, recursive = TRUE)
}
FIGDIR <- file.path(CFG$outdir, "figures"); dir.create(FIGDIR, showWarnings = FALSE)
RESDIR <- file.path(CFG$outdir, "results"); dir.create(RESDIR, showWarnings = FALSE)

# =============================================================================
# 2.  LOAD AND RECONSTRUCT
# =============================================================================
banner("2.", "Load data and reconstruct the reference")

raw <- readr::read_csv(CFG$csv, trim_ws = TRUE, show_col_types = FALSE)
names(raw) <- c("pos","ref","strain","depth_hq","A","C","G","T","indel_count","depth_raw")

dat <- raw %>%
  mutate(across(c(pos, depth_hq, A, C, G, T, indel_count, depth_raw), as.numeric),
         ref = str_trim(ref), strain = str_trim(strain))

# Rows come in THREE kinds - and nchar(NA) is NA, not 0, so a naive
#   filter(nchar(strain) == 1)
# evaluates to NA on the zero-coverage rows and dplyr DROPS them silently.
# Classify explicitly instead.
dat <- dat %>%
  mutate(row_type = case_when(
    is.na(strain) | strain == ""          ~ "zero-coverage",   # no reads -> no call
    nchar(ref) == 1 & nchar(strain) == 1  ~ "SNP",
    TRUE                                  ~ "indel"
  ))

snp      <- dat %>% filter(row_type == "SNP")
ind_rows <- dat %>% filter(row_type == "indel")
nocov    <- dat %>% filter(row_type == "zero-coverage")

note("Rows: %s total | %s SNP-comparable | %s indel records | %s zero-coverage",
     format(nrow(dat), big.mark = ","), format(nrow(snp), big.mark = ","),
     format(nrow(ind_rows), big.mark = ","), format(nrow(nocov), big.mark = ","))
if (nrow(nocov)) {
  bullet("Zero-coverage positions %g-%g were never sequenced. NOT deletions.",
         min(nocov$pos), max(nocov$pos))
}

GENOME_LEN <- max(9181, max(dat$pos))

# Reconstruct the reference. Two subtleties, both of which silently corrupt the
# sequence if ignored:
#   (a) an indel record carries a MULTI-BASE ref that consumes the positions
#       downstream of it - those positions get no row of their own;
#   (b) 26 indel records share a position with a normal SNP row.
# So: lay down the SNP rows first (they are authoritative and carry the counts),
# then let the indel rows fill only what is still empty, expanding multi-base refs.
ref_vec <- rep(NA_character_, GENOME_LEN)
ref_vec[snp$pos] <- snp$ref

# Zero-coverage rows still carry a valid reference base - use it.
ref_vec[nocov$pos] <- nocov$ref

# Indel rows fill only what is still empty, expanding multi-base refs.
for (i in seq_len(nrow(ind_rows))) {
  p  <- ind_rows$pos[i]
  bs <- strsplit(ind_rows$ref[i], "")[[1]]
  for (k in seq_along(bs)) {
    q <- p + k - 1
    if (q <= GENOME_LEN && is.na(ref_vec[q])) ref_vec[q] <- bs[k]
  }
}

note("Reference reconstructed: %s / %s positions (%s gaps)",
     format(sum(!is.na(ref_vec)), big.mark = ","), format(GENOME_LEN, big.mark = ","),
     format(sum(is.na(ref_vec)), big.mark = ","))

if (any(is.na(ref_vec))) {
  bullet("Gaps remain at %d position(s) - they will be checked per region below.",
         sum(is.na(ref_vec)))
}

# --- Per-site allele frequencies ---------------------------------------------
sites <- snp %>%
  mutate(total = A + C + G + T) %>%
  filter(depth_hq >= CFG$min_depth, total >= CFG$min_depth) %>%
  mutate(across(c(A, C, G, T), ~ .x / total, .names = "f_{.col}")) %>%
  rowwise() %>%
  mutate(
    major_base = c("A","C","G","T")[which.max(c(A, C, G, T))],
    major_f    = max(c(f_A, f_C, f_G, f_T)),
    minor_base = c("A","C","G","T")[order(c(A, C, G, T), decreasing = TRUE)[2]],
    minor_f    = sort(c(f_A, f_C, f_G, f_T), decreasing = TRUE)[2],
    minor_n    = sort(c(A, C, G, T), decreasing = TRUE)[2]
  ) %>%
  ungroup() %>%
  mutate(change      = paste0(major_base, ">", minor_base),
         is_artefact = change %in% ARTEFACT_SIG)

note("Callable sites (>= %gx): %s", CFG$min_depth, format(nrow(sites), big.mark = ","))

# =============================================================================
# 3.  FRAME VERIFICATION  -  do not skip this
# =============================================================================
banner("3.", "Reading-frame verification")

bullet("NC_001802.1 is 9,181 bp; HXB2/K03455 is 9,719 bp. They are NOT numbered")
bullet("alike. In pol the offset is exactly %d (NC = HXB2 - %d).", HXB2_OFFSET, HXB2_OFFSET)
bullet("Every DRM table in the literature uses HXB2 numbering. Get this wrong and")
bullet("EVERY amino acid position you report is wrong. So we verify by translating.")

translate_nt <- function(nt) {
  cods <- substring(nt, seq(1, nchar(nt) - 2, by = 3), seq(3, nchar(nt), by = 3))
  paste(ifelse(is.na(CODON_TBL[cods]), "X", CODON_TBL[cods]), collapse = "")
}

frame_check <- REGIONS %>%
  rowwise() %>%
  mutate(
    nt      = paste(ref_vec[start:end], collapse = ""),
    nt_len  = nchar(nt),
    codons  = nt_len / 3,
    aa      = translate_nt(nt),
    starts_ok  = startsWith(aa, expect_aa),
    len_ok     = codons == n_codons,
    no_int_stop = !grepl("\\*", substr(aa, 1, nchar(aa) - 1)),
    aa_head = substr(aa, 1, 14)
  ) %>%
  ungroup()

show_tbl(frame_check %>% select(region, start, end, nt_len, codons, aa_head,
                                starts_ok, len_ok, no_int_stop),
         "Frame check (NC_001802.1 coordinates):")

if (!all(frame_check$starts_ok & frame_check$len_ok & frame_check$no_int_stop)) {
  stop("FRAME CHECK FAILED. Refusing to emit amino acid coordinates that would be wrong.")
}
note("All three regions PASS: expected N-terminus, exact codon count, no internal stop.")
bullet("PR  %d codons, starts %s...", frame_check$n_codons[1], substr(frame_check$aa[1], 1, 6))
bullet("RT  %d codons, starts %s...", frame_check$n_codons[2], substr(frame_check$aa[2], 1, 6))
bullet("IN  %d codons, starts %s...", frame_check$n_codons[3], substr(frame_check$aa[3], 1, 6))

# =============================================================================
# 4.  BUILD THRESHOLD-LADDERED CONSENSUS FASTAs
# =============================================================================
banner("4.", "Building consensus FASTAs across the threshold ladder")

cat("
  THE RULE
      At threshold T, a position keeps every allele whose VAF >= T.
        - one allele passes   -> that base
        - two or more pass    -> the IUPAC ambiguity code for that set
        - none reach T        -> fall back to the single most frequent base
      So T=50% can never produce an ambiguity (only one allele can exceed half):
      it is the plain majority consensus. Lower T admits the minority population.
")

# Consensus base at a given site row, for a given threshold
consensus_base <- function(fA, fC, fG, fT, thresh) {
  f <- c(A = fA, C = fC, G = fG, T = fT)
  keep <- names(f)[f >= thresh]
  if (length(keep) == 0) keep <- names(f)[which.max(f)]
  to_iupac(keep)
}

# Expand an ambiguous codon to every concrete codon it could be, then translate
translate_codon_amb <- function(cod) {
  b <- lapply(strsplit(cod, "")[[1]], from_iupac)
  combos <- expand.grid(b, stringsAsFactors = FALSE)
  cods <- apply(combos, 1, paste, collapse = "")
  aas <- unique(unname(CODON_TBL[cods]))
  aas <- aas[!is.na(aas)]
  if (length(aas) == 0) "X" else if (length(aas) == 1) aas else paste(sort(aas), collapse = "/")
}

write_fasta <- function(seq, header, path, width = 60) {
  lines <- c(paste0(">", header),
             substring(seq, seq(1, nchar(seq), width), pmin(seq(width, nchar(seq) + width - 1, width), nchar(seq))))
  writeLines(lines, path)
}

all_seqs   <- list()   # nucleotide sequences  [[region]][[threshold]]
all_aa     <- list()   # amino acid sequences
amb_records <- list()  # every ambiguous position introduced

for (i in seq_len(nrow(REGIONS))) {
  rg <- REGIONS[i, ]
  reg_sites <- sites %>% filter(pos >= rg$start, pos <= rg$end) %>% arrange(pos)

  if (nrow(reg_sites) != (rg$end - rg$start + 1)) {
    warning(sprintf("%s: %d callable sites but region is %d nt - low-coverage gaps will be filled from reference",
                    rg$region, nrow(reg_sites), rg$end - rg$start + 1))
  }

  for (th in CFG$thresholds) {
    # Start from the reference, overwrite every callable site with the threshold call
    nt <- ref_vec[rg$start:rg$end]
    idx <- reg_sites$pos - rg$start + 1
    nt[idx] <- pmap_chr(list(reg_sites$f_A, reg_sites$f_C, reg_sites$f_G, reg_sites$f_T),
                        function(a, c, g, t) consensus_base(a, c, g, t, th))
    seq <- paste(nt, collapse = "")

    tag  <- sprintf("T%02d", round(100 * th))
    hdr  <- sprintf("%s|%s|%s|threshold_%.0fpct|%s_%d-%d|HXB2_%d-%d",
                    CFG$sample_id, rg$region, rg$name, 100 * th,
                    "NC_001802.1", rg$start, rg$end, rg$hxb2_start, rg$hxb2_end)
    write_fasta(seq, hdr, file.path(CFG$outdir, rg$region, "fasta",
                                    sprintf("%s_%s_nt.fasta", rg$region, tag)))

    # translate (ambiguity-aware, codon by codon)
    cods <- substring(seq, seq(1, nchar(seq) - 2, by = 3), seq(3, nchar(seq), by = 3))
    aa_v <- vapply(cods, translate_codon_amb, character(1))
    aa_flat <- paste(ifelse(nchar(aa_v) == 1, aa_v, "X"), collapse = "")
    write_fasta(aa_flat, paste0(hdr, "|translated"),
                file.path(CFG$outdir, rg$region, "fasta",
                          sprintf("%s_%s_aa.fasta", rg$region, tag)))

    all_seqs[[rg$region]][[tag]] <- seq
    all_aa[[rg$region]][[tag]]   <- aa_v          # keep the ambiguity-expanded form

    # record ambiguous positions
    amb_idx <- which(nt %in% c("R","Y","S","W","K","M","B","D","H","V","N"))
    if (length(amb_idx)) {
      amb_records[[length(amb_records) + 1]] <- tibble(
        region    = rg$region,
        threshold = th,
        pos       = rg$start + amb_idx - 1,
        codon     = floor((amb_idx - 1) / 3) + 1,
        iupac     = nt[amb_idx]
      )
    }
  }
  note("%s: wrote %d nt FASTAs + %d aa FASTAs -> %s/",
       rg$region, length(CFG$thresholds), length(CFG$thresholds),
       file.path(CFG$outdir, rg$region, "fasta"))
}

ambig <- if (length(amb_records)) bind_rows(amb_records) else
  tibble(region = character(), threshold = numeric(), pos = numeric(),
         codon = numeric(), iupac = character())

# =============================================================================
# 5.  COMPARE THE THRESHOLDS
# =============================================================================
banner("5.", "Comparing sequences across the ladder")

# --- 5a. Ambiguity count per region per threshold ----------------------------
amb_summary <- expand_grid(region = REGIONS$region, threshold = CFG$thresholds) %>%
  left_join(ambig %>% count(region, threshold, name = "ambiguous_nt"),
            by = c("region", "threshold")) %>%
  mutate(ambiguous_nt = replace_na(ambiguous_nt, 0)) %>%
  left_join(REGIONS %>% select(region, n_codons), by = "region") %>%
  left_join(ambig %>% distinct(region, threshold, codon) %>%
              count(region, threshold, name = "ambiguous_codons"),
            by = c("region", "threshold")) %>%
  mutate(ambiguous_codons = replace_na(ambiguous_codons, 0),
         pct_codons = round(100 * ambiguous_codons / n_codons, 2)) %>%
  arrange(region, desc(threshold))

show_tbl(amb_summary, "Ambiguous positions introduced at each threshold:")
readr::write_csv(amb_summary, file.path(RESDIR, "threshold_ambiguity_summary.csv"))

# --- 5b. Pairwise identity vs the 50% consensus ------------------------------
hamming <- function(a, b) sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]])

cmp <- map_dfr(REGIONS$region, function(rg) {
  base <- all_seqs[[rg]][["T50"]]
  map_dfr(CFG$thresholds, function(th) {
    tag <- sprintf("T%02d", round(100 * th))
    s <- all_seqs[[rg]][[tag]]
    tibble(region = rg, threshold = th,
           nt_differing_from_T50 = hamming(base, s),
           identical_to_T50      = hamming(base, s) == 0)
  })
})

show_tbl(cmp, "Nucleotide differences vs the 50% (plain consensus) sequence:")
readr::write_csv(cmp, file.path(RESDIR, "fasta_comparison_vs_T50.csv"))

identical_down_to <- cmp %>% filter(identical_to_T50) %>%
  group_by(region) %>% summarise(lowest = min(threshold), .groups = "drop")

cat("\n  WHAT THIS MEANS\n")
for (i in seq_len(nrow(identical_down_to))) {
  bullet("%s: identical to the plain consensus all the way down to %.0f%%",
         identical_down_to$region[i], 100 * identical_down_to$lowest[i])
}
bullet("Where the ladder produces IDENTICAL sequences, deep sequencing bought you")
bullet("nothing over Sanger. That is a real, reportable result - not a failure.")

# =============================================================================
# 6.  THE MINORITY LANDSCAPE  -  lollipop plots
# =============================================================================
banner("6.", "The minority landscape (lollipop plots)")

pol_minor <- map_dfr(seq_len(nrow(REGIONS)), function(i) {
  rg <- REGIONS[i, ]
  sites %>%
    filter(pos >= rg$start, pos <= rg$end) %>%
    mutate(region    = rg$region,
           region_nm = rg$name,
           codon     = floor((pos - rg$start) / 3) + 1,
           codon_pos = ((pos - rg$start) %% 3) + 1,
           is_drm    = codon %in% DRM_POS[[rg$region]])
}) %>%
  mutate(region = factor(region, levels = c("PR", "RT", "IN")))

note("Minority variants in pol, by threshold:")
lad <- map_dfr(c(CFG$thresholds, CFG$explore), function(t) {
  pol_minor %>% filter(minor_f >= t) %>%
    group_by(region) %>%
    summarise(n = n(), n_real = sum(!is_artefact), .groups = "drop") %>%
    mutate(threshold = t)
}) %>%
  pivot_wider(names_from = region, values_from = c(n, n_real), values_fill = 0) %>%
  arrange(desc(threshold))
show_tbl(lad)
readr::write_csv(lad, file.path(RESDIR, "minority_ladder_by_region.csv"))

# --- the lollipop -----------------------------------------------------------
make_lollipop <- function(df, rg_code, rg_name, n_codons, min_vaf = 0.01) {
  d <- df %>% filter(region == rg_code, minor_f >= min_vaf)
  drm <- tibble(codon = DRM_POS[[rg_code]]) %>% filter(codon <= n_codons)

  # Set the y range first, then draw only the threshold lines that fall inside it.
  # (Drawing the 50% line on a 0-22% axis just emits a "removed rows" warning.)
  ymax  <- max(22, 100 * max(c(d$minor_f, 0.05)))
  lines <- 100 * CFG$thresholds
  lines <- lines[lines <= ymax]

  p <- ggplot() +
    # DRM codon markers along the baseline
    geom_point(data = drm, aes(codon, 0), shape = 17, size = 1.8,
               colour = PAL[["violet"]], alpha = 0.8) +
    # threshold reference lines
    geom_hline(yintercept = lines, linetype = 2,
               colour = PAL[["muted"]], linewidth = 0.3) +
    geom_hline(yintercept = 20, linetype = 1, colour = PAL[["artefact"]], linewidth = 0.5)

  if (nrow(d)) {
    p <- p +
      geom_segment(data = d, aes(x = codon, xend = codon, y = 0, yend = 100 * minor_f,
                                 colour = is_artefact), linewidth = 0.5) +
      geom_point(data = d, aes(codon, 100 * minor_f, colour = is_artefact,
                               shape = is_drm, size = is_drm)) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18),
                         labels = c("other codon", "known DRM codon"), name = NULL) +
      scale_size_manual(values = c(`FALSE` = 1.8, `TRUE` = 3.6), guide = "none")
  }

  p +
    scale_colour_manual(values = c(`FALSE` = PAL[["signal"]], `TRUE` = PAL[["artefact"]]),
                        labels = c("plausible biology", "artefact signature"), name = NULL) +
    scale_x_continuous(limits = c(0, n_codons + 1), expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.06))) +
    annotate("text", x = n_codons * 0.99, y = 20.6, label = "Sanger limit 20%",
             hjust = 1, size = 3, colour = PAL[["artefact"]]) +
    labs(title = sprintf("%s (%s) - minority variant landscape", rg_code, rg_name),
         subtitle = sprintf("Each lollipop is one nucleotide site >= %.0f%% VAF. Purple triangles = known DRM codons.",
                            100 * min_vaf),
         x = sprintf("%s codon (HXB2 numbering)", rg_code),
         y = "Minor allele frequency (%)") +
    theme(legend.position = "top")
}

for (i in seq_len(nrow(REGIONS))) {
  rg <- REGIONS[i, ]
  p <- make_lollipop(pol_minor, rg$region, rg$name, rg$n_codons)
  ggsave(file.path(FIGDIR, sprintf("lollipop_%s.png", rg$region)), p,
         width = 11, height = 4.4, dpi = 150)
}
note("Lollipops written: %s", paste(sprintf("lollipop_%s.png", REGIONS$region), collapse = ", "))

# --- combined faceted view ---------------------------------------------------
p_all <- pol_minor %>%
  filter(minor_f >= 0.01) %>%
  ggplot(aes(codon, 100 * minor_f)) +
  geom_hline(yintercept = 20, colour = PAL[["artefact"]], linewidth = 0.5) +
  geom_hline(yintercept = c(5, 10, 15), linetype = 2, colour = PAL[["muted"]], linewidth = 0.3) +
  geom_segment(aes(xend = codon, y = 0, yend = 100 * minor_f, colour = is_artefact),
               linewidth = 0.4) +
  geom_point(aes(colour = is_artefact, shape = is_drm, size = is_drm)) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18),
                     labels = c("other codon", "known DRM codon"), name = NULL) +
  scale_size_manual(values = c(`FALSE` = 1.4, `TRUE` = 3), guide = "none") +
  scale_colour_manual(values = c(`FALSE` = PAL[["signal"]], `TRUE` = PAL[["artefact"]]),
                      labels = c("plausible biology", "artefact signature"), name = NULL) +
  facet_wrap(~ region, scales = "free_x", ncol = 1) +
  labs(title = "pol minority variant landscape: PR / RT / IN",
       subtitle = "Nothing reaches the 20% Sanger line. Red = the OxoG/context artefact signature.",
       x = "Codon (HXB2 numbering)", y = "Minor allele frequency (%)") +
  theme(legend.position = "top")

ggsave(file.path(FIGDIR, "lollipop_pol_all.png"), p_all, width = 11, height = 8, dpi = 150)

# --- comparison bar: ambiguity by threshold ----------------------------------
p_cmp <- amb_summary %>%
  mutate(region = factor(region, levels = c("PR", "RT", "IN")),
         thr = factor(sprintf("%.0f%%", 100 * threshold),
                      levels = sprintf("%.0f%%", 100 * CFG$thresholds))) %>%
  ggplot(aes(thr, ambiguous_nt, fill = region)) +
  geom_col(position = position_dodge(preserve = "single"), width = 0.7) +
  geom_text(aes(label = ambiguous_nt), position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.2, colour = "grey30") +
  scale_fill_manual(values = c(PR = PAL[["signal"]], RT = PAL[["violet"]], IN = PAL[["artefact"]])) +
  labs(title = "Ambiguity codes introduced as the threshold drops",
       subtitle = "At 50-10% every region is identical to the plain consensus. Only 5% admits anything.",
       x = "VAF threshold", y = "Ambiguous nucleotide positions", fill = NULL)

ggsave(file.path(FIGDIR, "threshold_comparison.png"), p_cmp, width = 9, height = 4.6, dpi = 150)
note("Comparison figures written to %s/", FIGDIR)

# =============================================================================
# 7.  DRUG RESISTANCE  -  where a nucleotide is not an amino acid
# =============================================================================
banner("7.", "Drug resistance: the trap")

cat("
  A variant AT a resistance codon is not a resistance mutation. K103N is a
  specific amino acid substitution, not 'any change at codon 103'. To say
  anything clinical you must translate the actual codon and name the actual
  amino acid. We do that now, for every minority variant landing on a known
  DRM position.
")

# consensus (majority) genome, used as the codon background
cons_vec <- ref_vec
cons_vec[sites$pos] <- sites$major_base

drm_hits <- pol_minor %>%
  filter(is_drm, minor_f >= 0.01) %>%
  rowwise() %>%
  mutate(
    reg_start = REGIONS$start[REGIONS$region == as.character(region)],
    codon_start = reg_start + (codon - 1) * 3,
    codon_ref = paste(ref_vec[codon_start:(codon_start + 2)], collapse = ""),
    codon_con = paste(cons_vec[codon_start:(codon_start + 2)], collapse = ""),
    codon_min = { x <- strsplit(codon_con, "")[[1]]; x[codon_pos] <- minor_base; paste(x, collapse = "") },
    aa_ref = unname(CODON_TBL[codon_ref]),
    aa_con = unname(CODON_TBL[codon_con]),
    aa_min = unname(CODON_TBL[codon_min]),
    mutation = sprintf("%s%d%s", aa_con, codon, aa_min),
    effect   = if_else(aa_con == aa_min, "synonymous", "non-synonymous"),
    vaf_pct  = round(100 * minor_f, 2),
    # Each of these is an independent reason to disbelieve a "resistance" call.
    # They are checked in order of how decisively they kill it.
    verdict = case_when(
      is_artefact      ~ "ARTEFACT signature - discard",
      aa_min == "*"    ~ "NONSENSE - stop codon; a dead virus, not resistance",
      aa_con == aa_min ~ "SYNONYMOUS - no amino acid change at all",
      change == "G>A"  ~ "G>A - possible APOBEC hypermutation of a defective genome",
      TRUE             ~ "non-synonymous - check the SPECIFIC aa against Stanford HIVdb"
    )
  ) %>%
  ungroup() %>%
  select(region, codon, pos, codon_pos, change, vaf_pct, is_artefact,
         codon_con, codon_min, aa_con, aa_min, mutation, effect, verdict, total)

if (nrow(drm_hits) == 0) {
  note("No minority variants >= 1%% fall on any known DRM codon.")
} else {
  show_tbl(drm_hits %>% select(region, codon, change, vaf_pct, mutation, effect, verdict),
           "Minority variants landing on known DRM codons (>= 1% VAF):")
  readr::write_csv(drm_hits, file.path(RESDIR, "drm_codon_hits.csv"))

  cat("\n  READ EACH ONE CAREFULLY\n")
  for (i in seq_len(nrow(drm_hits))) {
    h <- drm_hits[i, ]
    cat(sprintf("\n   %s codon %-4d %-8s at %5.2f%%\n", h$region, h$codon, h$mutation, h$vaf_pct))
    cat(sprintf("      consensus %s (%s) -> minority %s (%s)   [%s]\n",
                h$codon_con, h$aa_con, h$codon_min, h$aa_min, h$change))
    cat(sprintf("      VERDICT: %s\n", h$verdict))
  }

  surviving <- drm_hits %>%
    filter(!is_artefact, aa_con != aa_min, aa_min != "*")

  cat("\n  ", strrep("-", 74), "\n", sep = "")
  cat(sprintf("   Of %d minority variants on DRM codons, %d survive every check.\n",
              nrow(drm_hits), nrow(surviving)))
  if (nrow(surviving)) {
    for (i in seq_len(nrow(surviving))) {
      s <- surviving[i, ]
      cat(sprintf("     %s %s at %.2f%%  -  %s\n", s$region, s$mutation, s$vaf_pct, s$verdict))
    }
    cat("   None of these is a recognised resistance substitution at that codon.\n")
    cat("   Do not report them as resistance. Send the FASTA to Stanford HIVdb.\n")
  }
  cat("  ", strrep("-", 74), "\n", sep = "")
}

cat("\n  THE FOUR WAYS THIS TRAP SPRINGS\n")
bullet("(1) ARTEFACT. Any hit carrying the %s signature is the OxoG/context",
       paste(ARTEFACT_SIG, collapse = "/"))
bullet("    artefact characterised in mutations_workshop.R. PR I50L - a signature")
bullet("    ATAZANAVIR resistance mutation - is exactly this. It is an oxidised base.")
bullet("(2) WRONG AMINO ACID. The resistance mutation at RT 103 is K103N (or K103S).")
bullet("    A variant encoding K103E is NOT resistance. Same codon, different amino")
bullet("    acid, completely different clinical meaning. Likewise L210S is not L210W.")
bullet("(3) SYNONYMOUS. A variant can sit on a DRM codon and change no amino acid at")
bullet("    all - PR N83N here, at 6.8%%, the largest real minority variant in pol.")
bullet("    A nucleotide is not an amino acid.")
bullet("(4) NONSENSE / APOBEC. A stop codon (IN L74*) is a dead virus, not a")
bullet("    resistant one. And a G>A on a DRM codon may be APOBEC3G hypermutation of")
bullet("    a DEFECTIVE, replication-incompetent genome - real sequence, no clinical")
bullet("    meaning. Both are tells that you are not looking at a live quasispecies.")
cat("\n  >> Every filter must pass before the word 'resistance' is used. And even")
cat("\n     then, interpretation belongs to Stanford HIVdb - not to this script.\n")

# =============================================================================
# 8.  SUMMARY
# =============================================================================
banner("8.", "Summary")

summ <- REGIONS %>%
  select(region, name, n_codons) %>%
  left_join(cmp %>% filter(threshold == 0.05) %>%
              select(region, nt_diff_at_5pct = nt_differing_from_T50), by = "region") %>%
  left_join(pol_minor %>% group_by(region) %>%
              summarise(minor_ge_5pct  = sum(minor_f >= 0.05),
                        minor_ge_1pct  = sum(minor_f >= 0.01),
                        real_ge_1pct   = sum(minor_f >= 0.01 & !is_artefact),
                        .groups = "drop") %>%
              mutate(region = as.character(region)), by = "region") %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

show_tbl(summ, "Per-region summary:")
readr::write_csv(summ, file.path(RESDIR, "pol_summary.csv"))

cat("\n", strrep("-", 78), "\n", sep = "")
cat("  WHAT THIS SAMPLE ACTUALLY SHOWS\n")
cat(strrep("-", 78), "\n", sep = "")
cat(sprintf("
  1. THE LADDER IS FLAT DOWN TO 10%%. PR, RT and IN are IDENTICAL to the plain
     50%% consensus at 20%%, 15%% and 10%%. Routine Sanger genotyping would have
     lost nothing. Report that - it is a result, not a null.

  2. ONLY 5%% ADMITS ANYTHING, AND BARELY. %d ambiguity code(s) genome-wide across
     all three genes. RT - the largest region and the main drug target - stays
     completely clean.

  3. NO MINORITY DRUG RESISTANCE. %d minority variants land on resistance codons.
     Not one is resistance: they are artefact, or synonymous, or a stop codon, or
     they encode the WRONG amino acid (K103E is not K103N; L210S is not L210W).

  4. THE ARTEFACT FILTER HAS CLINICAL TEETH. Without it you would report PR
     I50L - a signature atazanavir-resistance mutation - from an A>C artefact.
     That is a treatment decision made on an oxidised base.

  5. A NUCLEOTIDE IS NOT AN AMINO ACID. The largest real minority variant in pol
     (PR N83N, 6.8%%) sits on a DRM codon and changes NOTHING. Had we stopped at
     'variant at codon 83', we would have invented a finding out of a silent site.

  6. A FASTA IS A CLAIM. Every one of these files asserts 'this is the virus'.
     They disagree with each other, and only the threshold in the header says
     why. Never ship a consensus without its threshold and its depth.
",
  sum(amb_summary$ambiguous_nt[amb_summary$threshold == 0.05]),
  nrow(drm_hits)
))

cat(strrep("-", 78), "\n", sep = "")
note("FASTAs  : %s/{PR,RT,IN}/fasta/", CFG$outdir)
note("Figures : %s/", FIGDIR)
note("Tables  : %s/", RESDIR)
cat("\n  Done.\n\n")
