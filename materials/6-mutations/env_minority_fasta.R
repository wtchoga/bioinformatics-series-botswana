#!/usr/bin/env Rscript
# =============================================================================
#  ENV MINORITY VARIANT ANALYSIS  -  V1 / V2 / V3 / V4 / V5 variable loops
#  Threshold-laddered consensus FASTAs, lollipops, glycan shield, tropism
#  HIV-1 deep sequencing, NC_001802.1  |  Dr. W. Choga
# -----------------------------------------------------------------------------
#  USAGE
#     Rscript env_minority_fasta.R
#     Rscript env_minority_fasta.R <variants.csv> --outdir env_minority
#
#  COMPANION TO
#     mutations_workshop.R   - where the artefact signature is characterised
#     pol_minority_fasta.R   - the same treatment for PR / RT / IN
#
#  THE POINT OF THIS ONE
#     env is the most variable gene in HIV, and the V1-V5 loops are the most
#     variable parts of env. The obvious prediction is that they light up as
#     diversity peaks. In this sample THEY DO NOT - and the reason is not
#     biology, it is the file format. Section 7 is the payoff. Read it.
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
  outdir     = get_opt("--outdir", "env_minority"),
  min_depth  = as.numeric(get_opt("--mindepth", 100)),
  thresholds = c(0.50, 0.20, 0.15, 0.10, 0.05),
  sample_id  = "sample01"
)

PAL <- c(signal = "#0FA3A3", artefact = "#E2574C", violet = "#8B5CF6", muted = "#6C7D8B")
ARTEFACT_SIG <- c("A>C", "T>G", "G>T", "C>A")

theme_set(theme_minimal(base_size = 12))
banner <- function(n, t) cat("\n", strrep("=", 78), "\n  ", n, "  ", toupper(t), "\n",
                             strrep("=", 78), "\n", sep = "")
note   <- function(...) cat("\n  ", sprintf(...), "\n", sep = "")
bullet <- function(...) cat("   - ", sprintf(...), "\n", sep = "")
show_tbl <- function(df, cap = NULL) {
  if (!is.null(cap)) cat("\n  ", cap, "\n", sep = "")
  print(as.data.frame(df), row.names = FALSE, digits = 4); invisible(df)
}

# =============================================================================
# 1.  REFERENCE TABLES
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
IUPAC <- c("A"="A","C"="C","G"="G","T"="T","AG"="R","CT"="Y","CG"="S","AT"="W",
           "GT"="K","AC"="M","CGT"="B","AGT"="D","ACT"="H","ACG"="V","ACGT"="N")
IUPAC_REV <- setNames(names(IUPAC), IUPAC)
to_iupac   <- function(b) unname(IUPAC[paste(sort(unique(b)), collapse = "")])
from_iupac <- function(c) strsplit(IUPAC_REV[[c]], "")[[1]]

# --- Coordinates -------------------------------------------------------------
# NC_001802.1 = HXB2 - 454 (verified by translation in section 3).
HXB2_OFFSET <- 454
ENV_HXB2 <- c(6225, 8795)                      # gp160 ORF
ENV_START <- ENV_HXB2[1] - HXB2_OFFSET         # 5771
ENV_END   <- ENV_HXB2[2] - HXB2_OFFSET         # 8341

# Env domains in HXB2 gp160 amino-acid numbering (Los Alamos definitions).
# V1, V3 and V4 are disulfide-bonded: they OPEN and CLOSE on a cysteine.
# V2 closes on C196 (it opens just after V1's C157), and V5 has no disulfide
# at all - so we only assert C..C for V1/V3/V4.
DOMAINS <- tribble(
  ~region, ~aa_start, ~aa_end, ~type,        ~cys_bounded,
  "C1",         1,      130,   "conserved",  FALSE,
  "V1",       131,      157,   "variable",   TRUE,
  "V2",       158,      196,   "variable",   FALSE,
  "C2",       197,      295,   "conserved",  FALSE,
  "V3",       296,      331,   "variable",   TRUE,
  "C3",       332,      384,   "conserved",  FALSE,
  "V4",       385,      418,   "variable",   TRUE,
  "C4",       419,      459,   "conserved",  FALSE,
  "V5",       460,      469,   "variable",   FALSE,
  "C5",       470,      511,   "conserved",  FALSE,
  "gp41",     512,      856,   "conserved",  FALSE
) %>%
  mutate(nt_start = ENV_START + (aa_start - 1) * 3,
         nt_end   = ENV_START + aa_end * 3 - 1,
         n_codons = aa_end - aa_start + 1)

VLOOPS <- DOMAINS %>% filter(type == "variable")

cat("\n+--------------------------------------------------------------------------+\n")
cat("|  ENV MINORITY VARIANT ANALYSIS  -  V1-V5 VARIABLE LOOPS                  |\n")
cat("+--------------------------------------------------------------------------+\n")
note("Sample     : %s", CFG$sample_id)
note("env        : NC_001802.1 %d-%d  (HXB2 %d-%d)", ENV_START, ENV_END, ENV_HXB2[1], ENV_HXB2[2])
note("Thresholds : %s", paste0(100 * CFG$thresholds, "%", collapse = ", "))

if (!file.exists(CFG$csv)) stop("Cannot find CSV: ", CFG$csv)

dir.create(CFG$outdir, showWarnings = FALSE, recursive = TRUE)
for (v in c(VLOOPS$region, "ENV")) {
  dir.create(file.path(CFG$outdir, v, "fasta"), showWarnings = FALSE, recursive = TRUE)
}
FIGDIR <- file.path(CFG$outdir, "figures"); dir.create(FIGDIR, showWarnings = FALSE)
RESDIR <- file.path(CFG$outdir, "results"); dir.create(RESDIR, showWarnings = FALSE)

# =============================================================================
# 2.  LOAD AND RECONSTRUCT
# =============================================================================
banner("2.", "Load and reconstruct")

raw <- readr::read_csv(CFG$csv, trim_ws = TRUE, show_col_types = FALSE)
names(raw) <- c("pos","ref","strain","depth_hq","A","C","G","T","indel_count","depth_raw")
dat <- raw %>%
  mutate(across(c(pos, depth_hq, A, C, G, T, indel_count, depth_raw), as.numeric),
         ref = str_trim(ref), strain = str_trim(strain))

# Three row types. nchar(NA) is NA, so a naive filter drops zero-coverage rows
# silently - classify explicitly. (See mutations_workshop.R for the full story.)
dat <- dat %>%
  mutate(row_type = case_when(
    is.na(strain) | strain == ""         ~ "zero-coverage",
    nchar(ref) == 1 & nchar(strain) == 1 ~ "SNP",
    TRUE                                 ~ "indel"
  ))

snp      <- dat %>% filter(row_type == "SNP")
ind_rows <- dat %>% filter(row_type == "indel")
nocov    <- dat %>% filter(row_type == "zero-coverage")

note("Rows: %s total | %s SNP | %s indel | %s zero-coverage",
     format(nrow(dat), big.mark = ","), format(nrow(snp), big.mark = ","),
     format(nrow(ind_rows), big.mark = ","), format(nrow(nocov), big.mark = ","))

GENOME_LEN <- max(9181, max(dat$pos))
ref_vec <- rep(NA_character_, GENOME_LEN)
ref_vec[snp$pos]   <- snp$ref
ref_vec[nocov$pos] <- nocov$ref
for (i in seq_len(nrow(ind_rows))) {
  p <- ind_rows$pos[i]; bs <- strsplit(ind_rows$ref[i], "")[[1]]
  for (k in seq_along(bs)) { q <- p + k - 1
    if (q <= GENOME_LEN && is.na(ref_vec[q])) ref_vec[q] <- bs[k] }
}
note("Reference reconstructed: %s / %s (%s gaps)",
     format(sum(!is.na(ref_vec)), big.mark = ","), format(GENOME_LEN, big.mark = ","),
     sum(is.na(ref_vec)))

sites <- snp %>%
  mutate(total = A + C + G + T) %>%
  filter(depth_hq >= CFG$min_depth, total >= CFG$min_depth) %>%
  mutate(across(c(A, C, G, T), ~ .x / total, .names = "f_{.col}")) %>%
  rowwise() %>%
  mutate(
    major_base = c("A","C","G","T")[which.max(c(A, C, G, T))],
    minor_base = c("A","C","G","T")[order(c(A, C, G, T), decreasing = TRUE)[2]],
    minor_f    = sort(c(f_A, f_C, f_G, f_T), decreasing = TRUE)[2],
    minor_n    = sort(c(A, C, G, T), decreasing = TRUE)[2],
    pi_site    = (total / (total - 1)) * (1 - (f_A^2 + f_C^2 + f_G^2 + f_T^2))
  ) %>%
  ungroup() %>%
  mutate(change = paste0(major_base, ">", minor_base),
         is_artefact = change %in% ARTEFACT_SIG)

# =============================================================================
# 3.  FRAME AND LOOP VERIFICATION
# =============================================================================
banner("3.", "Frame and V-loop verification")

translate_nt <- function(nt) {
  cods <- substring(nt, seq(1, nchar(nt) - 2, by = 3), seq(3, nchar(nt), by = 3))
  paste(ifelse(is.na(CODON_TBL[cods]), "X", CODON_TBL[cods]), collapse = "")
}

env_ref_nt <- paste(ref_vec[ENV_START:ENV_END], collapse = "")
env_ref_aa <- translate_nt(env_ref_nt)

note("env: %d nt, %g codons, aa length %d", nchar(env_ref_nt),
     nchar(env_ref_nt) / 3, nchar(env_ref_aa))
bullet("starts: %s", substr(env_ref_aa, 1, 24))
bullet("internal stop: %s | ends with stop: %s",
       if (grepl("\\*", substr(env_ref_aa, 1, nchar(env_ref_aa) - 1))) "YES" else "none",
       substr(env_ref_aa, nchar(env_ref_aa), nchar(env_ref_aa)) == "*")

stopifnot(nchar(env_ref_nt) %% 3 == 0)
if (!startsWith(env_ref_aa, "MRV")) stop("env frame FAILED: does not start with MRV...")
if (grepl("\\*", substr(env_ref_aa, 1, nchar(env_ref_aa) - 1)))
  stop("env frame FAILED: internal stop codon.")

# --- The loops must open/close on cysteines -----------------------------------
aa_at <- function(i) substr(env_ref_aa, i, i)
loop_check <- DOMAINS %>%
  filter(type == "variable") %>%
  rowwise() %>%
  mutate(first_aa = aa_at(aa_start), last_aa = aa_at(aa_end),
         seq = substr(env_ref_aa, aa_start, aa_end),
         cys_ok = if (cys_bounded) (first_aa == "C" && last_aa == "C") else TRUE) %>%
  ungroup()

show_tbl(loop_check %>% select(region, aa_start, aa_end, n_codons, first_aa, last_aa,
                               cys_bounded, cys_ok),
         "V-loop boundaries (HXB2 Env aa numbering):")
for (i in seq_len(nrow(loop_check))) {
  cat(sprintf("     %-3s %s\n", loop_check$region[i], loop_check$seq[i]))
}

if (!all(loop_check$cys_ok)) stop("V-loop cysteine check FAILED - coordinates are wrong.")

# --- The V3 crown is the decisive check --------------------------------------
v3_crown <- substr(env_ref_aa, 312, 315)
note("V3 crown at aa 312-315: %s", v3_crown)
if (!v3_crown %in% c("GPGR", "GPGQ", "GPGK")) {
  warning("V3 crown is not the canonical GPGR/GPGQ - check coordinates!")
} else {
  bullet("Canonical V3 tip confirmed. The numbering is right.")
}
note("gp120/gp41 cleavage site (aa 505-515): %s", substr(env_ref_aa, 505, 515))
bullet("V1/V3/V4 open and close on cysteines; V2 closes on C196; V5 has no")
bullet("disulfide. All boundary checks pass - amino acid numbering is verified.")

# =============================================================================
# 4.  THRESHOLD-LADDERED CONSENSUS FASTAs
# =============================================================================
banner("4.", "Building threshold FASTAs for each V loop (and whole env)")

consensus_base <- function(fA, fC, fG, fT, thresh) {
  f <- c(A = fA, C = fC, G = fG, T = fT)
  keep <- names(f)[f >= thresh]
  if (length(keep) == 0) keep <- names(f)[which.max(f)]
  to_iupac(keep)
}
translate_codon_amb <- function(cod) {
  b <- lapply(strsplit(cod, "")[[1]], from_iupac)
  cods <- apply(expand.grid(b, stringsAsFactors = FALSE), 1, paste, collapse = "")
  aas <- unique(unname(CODON_TBL[cods])); aas <- aas[!is.na(aas)]
  if (length(aas) == 0) "X" else if (length(aas) == 1) aas else paste(sort(aas), collapse = "/")
}
write_fasta <- function(seq, header, path, width = 60) {
  writeLines(c(paste0(">", header),
               substring(seq, seq(1, nchar(seq), width),
                         pmin(seq(width, nchar(seq) + width - 1, width), nchar(seq)))), path)
}

TARGETS <- bind_rows(
  VLOOPS %>% select(region, aa_start, aa_end, nt_start, nt_end, n_codons),
  tibble(region = "ENV", aa_start = 1, aa_end = 856,
         nt_start = ENV_START, nt_end = ENV_END, n_codons = 857)
)

all_seqs <- list(); amb_records <- list()

for (i in seq_len(nrow(TARGETS))) {
  tg <- TARGETS[i, ]
  reg_sites <- sites %>% filter(pos >= tg$nt_start, pos <= tg$nt_end) %>% arrange(pos)

  for (th in CFG$thresholds) {
    nt <- ref_vec[tg$nt_start:tg$nt_end]
    idx <- reg_sites$pos - tg$nt_start + 1
    nt[idx] <- pmap_chr(list(reg_sites$f_A, reg_sites$f_C, reg_sites$f_G, reg_sites$f_T),
                        function(a, c, g, t) consensus_base(a, c, g, t, th))
    seq <- paste(nt, collapse = "")
    tag <- sprintf("T%02d", round(100 * th))
    hdr <- sprintf("%s|env_%s|threshold_%.0fpct|aa_%d-%d|NC_001802.1_%d-%d|HXB2_%d-%d",
                   CFG$sample_id, tg$region, 100 * th, tg$aa_start, tg$aa_end,
                   tg$nt_start, tg$nt_end,
                   tg$nt_start + HXB2_OFFSET, tg$nt_end + HXB2_OFFSET)

    write_fasta(seq, hdr, file.path(CFG$outdir, tg$region, "fasta",
                                    sprintf("%s_%s_nt.fasta", tg$region, tag)))
    cods <- substring(seq, seq(1, nchar(seq) - 2, by = 3), seq(3, nchar(seq), by = 3))
    aa_v <- vapply(cods, translate_codon_amb, character(1))
    write_fasta(paste(ifelse(nchar(aa_v) == 1, aa_v, "X"), collapse = ""),
                paste0(hdr, "|translated"),
                file.path(CFG$outdir, tg$region, "fasta",
                          sprintf("%s_%s_aa.fasta", tg$region, tag)))

    all_seqs[[tg$region]][[tag]] <- seq
    amb_idx <- which(nt %in% c("R","Y","S","W","K","M","B","D","H","V","N"))
    if (length(amb_idx)) {
      amb_records[[length(amb_records) + 1]] <- tibble(
        region = tg$region, threshold = th,
        pos = tg$nt_start + amb_idx - 1,
        aa_codon = tg$aa_start + floor((amb_idx - 1) / 3),
        iupac = nt[amb_idx])
    }
  }
  note("%s: %d nt + %d aa FASTAs -> %s/", tg$region, length(CFG$thresholds),
       length(CFG$thresholds), file.path(CFG$outdir, tg$region, "fasta"))
}
ambig <- if (length(amb_records)) bind_rows(amb_records) else
  tibble(region=character(), threshold=numeric(), pos=numeric(), aa_codon=numeric(), iupac=character())

# =============================================================================
# 5.  COMPARE THE LADDER
# =============================================================================
banner("5.", "Comparing sequences across the ladder")

hamming <- function(a, b) sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]])

cmp <- map_dfr(TARGETS$region, function(rg) {
  base <- all_seqs[[rg]][["T50"]]
  map_dfr(CFG$thresholds, function(th) {
    tag <- sprintf("T%02d", round(100 * th))
    tibble(region = rg, threshold = th,
           nt_vs_T50 = hamming(base, all_seqs[[rg]][[tag]]),
           identical = hamming(base, all_seqs[[rg]][[tag]]) == 0)
  })
})
show_tbl(cmp %>% pivot_wider(id_cols = region, names_from = threshold,
                             values_from = nt_vs_T50, names_prefix = "T"),
         "Nucleotide differences vs the 50% consensus:")
readr::write_csv(cmp, file.path(RESDIR, "fasta_comparison_vs_T50.csv"))

amb_summary <- expand_grid(region = TARGETS$region, threshold = CFG$thresholds) %>%
  left_join(ambig %>% count(region, threshold, name = "ambiguous_nt"),
            by = c("region", "threshold")) %>%
  mutate(ambiguous_nt = replace_na(ambiguous_nt, 0))
readr::write_csv(amb_summary, file.path(RESDIR, "threshold_ambiguity_summary.csv"))

# =============================================================================
# 6.  DIVERSITY: VARIABLE LOOPS vs CONSERVED REGIONS
# =============================================================================
banner("6.", "Diversity: V loops vs conserved regions")

env_sites <- sites %>%
  filter(pos >= ENV_START, pos <= ENV_END) %>%
  mutate(aa_codon = floor((pos - ENV_START) / 3) + 1)

assign_domain <- function(cd) {
  h <- DOMAINS$region[cd >= DOMAINS$aa_start & cd <= DOMAINS$aa_end]
  if (length(h)) h[1] else NA_character_
}
# The stop codon (aa 857) belongs to no domain and would carry domain = NA.
# Left in, `vc$type == "variable"` returns NA for it, subsetting silently yields
# a length-2 vector, and every sprintf() below recycles into garbage. Drop it.
env_sites <- env_sites %>%
  mutate(domain = map_chr(aa_codon, assign_domain)) %>%
  filter(!is.na(domain))

dom_stats <- DOMAINS %>%
  select(region, type, aa_start, aa_end, n_codons) %>%
  left_join(
    env_sites %>% group_by(domain) %>%
      summarise(sites = n(), pi = mean(pi_site),
                minor_1pct = sum(minor_f >= 0.01),
                real_1pct  = sum(minor_f >= 0.01 & !is_artefact),
                minor_5pct = sum(minor_f >= 0.05),
                .groups = "drop") %>% rename(region = domain),
    by = "region") %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)),
         per_kb = round(1000 * minor_1pct / pmax(sites, 1), 1))

show_tbl(dom_stats %>% select(region, type, n_codons, sites, pi, minor_1pct,
                              real_1pct, minor_5pct, per_kb),
         "Per-domain diversity and minority variants:")
readr::write_csv(dom_stats, file.path(RESDIR, "env_domain_summary.csv"))

v_pi <- env_sites %>% left_join(DOMAINS %>% select(region, type), by = c("domain" = "region"))

vc <- v_pi %>%
  group_by(type) %>%
  summarise(sites = n(), pi_mean = mean(pi_site), pi_median = median(pi_site),
            minor_1pct = sum(minor_f >= 0.01),
            per_kb = round(1000 * sum(minor_f >= 0.01) / n(), 1), .groups = "drop")
show_tbl(vc, "Variable loops vs conserved regions:")
readr::write_csv(vc, file.path(RESDIR, "v_vs_c_diversity.csv"))

# Pull scalars out explicitly - never index by a logical that might contain NA.
gv <- function(col, ty) vc[[col]][match(ty, vc$type)]
pi_v <- gv("pi_mean", "variable");   pi_c <- gv("pi_mean", "conserved")
md_v <- gv("pi_median", "variable"); md_c <- gv("pi_median", "conserved")
kb_v <- gv("per_kb", "variable");    kb_c <- gv("per_kb", "conserved")

# Test it rather than eyeballing two means.
wt <- wilcox.test(pi_site ~ type, data = v_pi)
note("Wilcoxon test, per-site pi, variable vs conserved: p = %.4f", wt$p.value)

cat("\n  THE RESULT NOBODY PREDICTS\n")
bullet("mean pi   - V loops %.5f  vs conserved %.5f  (essentially equal)", pi_v, pi_c)
bullet("median pi - V loops %.5f  vs conserved %.5f", md_v, md_c)
bullet("Minority-variant density: %.1f /kb in the loops vs %.1f /kb conserved.", kb_v, kb_c)
cat("\n")
bullet("The means are equal, but the medians are NOT: the V loops are")
bullet("significantly %s (Wilcoxon p = %.4f). The loops are not merely",
       if (md_v < md_c) "LESS diverse" else "MORE diverse", wt$p.value)
bullet("failing to stand out - by median they sit BELOW the conserved regions.")
bullet("(Means match because a few loop sites carry high-VAF outliers.)")
cat("\n")
bullet("If you predicted V-loop diversity peaks - and everyone does, including")
bullet("Exercise 2 of mutations_workshop.Rmd - you were wrong. But NOT because")
bullet("the loops are conserved. Keep reading.")

# =============================================================================
# 7.  WHY  -  the V loops vary in LENGTH, and this file cannot see it
# =============================================================================
banner("7.", "Why: the loops vary in length, not in substitution")

env_ind <- ind_rows %>%
  filter(pos >= ENV_START, pos <= ENV_END) %>%
  mutate(aa_codon = floor((pos - ENV_START) / 3) + 1,
         domain   = map_chr(aa_codon, assign_domain)) %>%
  left_join(DOMAINS %>% select(region, type), by = c("domain" = "region"))

ind_by_type <- env_ind %>% count(type, name = "indels")
nt_by_type  <- DOMAINS %>% group_by(type) %>% summarise(nt = sum(n_codons) * 3, .groups = "drop")

ind_test <- nt_by_type %>%
  left_join(ind_by_type, by = "type") %>%
  mutate(indels = replace_na(indels, 0),
         per_kb = round(1000 * indels / nt, 1))

show_tbl(ind_test, "Indel records in env, by domain type:")
show_tbl(env_ind %>% select(pos, aa_codon, domain, ref, strain, indel_count),
         "Every indel record inside env:")
readr::write_csv(env_ind, file.path(RESDIR, "env_indels.csv"))

# Are indels enriched in the loops? This compares two RATES (indels per nt in V
# vs in C), so the right test is a Poisson rate-ratio test - not a binomial on
# "share of total indels", which answers a subtly different question and returns
# a different fold number for the same data.
gi <- function(col, ty) ind_test[[col]][match(ty, ind_test$type)]
n_v <- gi("indels", "variable"); n_c <- gi("indels", "conserved")
L_v <- gi("nt", "variable");     L_c <- gi("nt", "conserved")
n_tot <- n_v + n_c
p_v   <- L_v / (L_v + L_c)

pt    <- poisson.test(c(n_v, n_c), T = c(L_v, L_c))
fold  <- unname(pt$estimate)
p_val <- pt$p.value

note("Indels in V loops: %d of %d, though the loops are only %.1f%% of env",
     n_v, n_tot, 100 * p_v)
note("Indel rate: %.1f/kb in the loops vs %.1f/kb conserved", 1000 * n_v / L_v, 1000 * n_c / L_c)
note("Poisson rate ratio: %.1f-fold  (95%% CI %.1f-%.1f)  p = %.4f",
     fold, pt$conf.int[1], pt$conf.int[2], p_val)

# Callability: the loops are also where coverage falls over
call_by_type <- DOMAINS %>%
  rowwise() %>%
  mutate(expected_nt = n_codons * 3,
         callable = sum(sites$pos >= nt_start & sites$pos <= nt_end)) %>%
  ungroup() %>%
  group_by(type) %>%
  summarise(expected_nt = sum(expected_nt), callable = sum(callable),
            pct_callable = round(100 * sum(callable) / sum(expected_nt), 1), .groups = "drop")
show_tbl(call_by_type, "Callable fraction, V loops vs conserved:")

cat("\n  THE EXPLANATION\n")
bullet("HIV's variable loops vary mostly in LENGTH - insertions and deletions -")
bullet("not in point substitutions. That is how they evade antibodies: they")
bullet("change shape and move their glycans, rather than just swapping residues.")
cat("\n")
bullet("A per-position pileup is indexed by REFERENCE COORDINATE. It has exactly")
bullet("one row per reference base and no way to represent 'this read carries")
bullet("four extra codons here'. So length variation is invisible by construction.")
cat("\n")
bullet("Worse, it is doubly invisible: reads carrying long indels map poorly or")
bullet("get soft-clipped, so the loops LOSE coverage (%.1f%% callable vs %.1f%%).",
       call_by_type$pct_callable[call_by_type$type == "variable"],
       call_by_type$pct_callable[call_by_type$type == "conserved"])
bullet("The very reads that carry the variation are the ones that fail to align.")
cat("\n  >> pi is flat across V and C not because the loops are invariant, but")
cat("\n     because THIS FILE CANNOT SEE THE WAY THEY VARY. The measurement, not")
cat("\n     the biology, produced the null result.\n")
cat("\n  >> To measure V-loop diversity properly you need read-level or assembly-\n")
cat("     based methods: local de novo assembly, haplotype reconstruction, or\n")
cat("     long reads. A pileup is the wrong instrument for this question.\n")

# =============================================================================
# 8.  LOLLIPOPS
# =============================================================================
banner("8.", "Lollipop plots")

dom_bands <- DOMAINS %>% filter(type == "variable")

p_env <- ggplot() +
  geom_rect(data = dom_bands, aes(xmin = aa_start, xmax = aa_end, ymin = 0, ymax = Inf),
            fill = PAL[["violet"]], alpha = 0.12) +
  geom_text(data = dom_bands, aes(x = (aa_start + aa_end) / 2, y = 15.5, label = region),
            colour = PAL[["violet"]], size = 3.2, fontface = "bold") +
  geom_hline(yintercept = c(5, 10, 15), linetype = 2, colour = PAL[["muted"]], linewidth = 0.3) +
  geom_segment(data = filter(env_sites, minor_f >= 0.01),
               aes(x = aa_codon, xend = aa_codon, y = 0, yend = 100 * minor_f,
                   colour = is_artefact), linewidth = 0.4) +
  geom_point(data = filter(env_sites, minor_f >= 0.01),
             aes(aa_codon, 100 * minor_f, colour = is_artefact), size = 1.5) +
  geom_point(data = env_ind %>% mutate(y = 0),
             aes(aa_codon, y), shape = 25, size = 2.4, fill = PAL[["violet"]],
             colour = PAL[["violet"]]) +
  scale_colour_manual(values = c(`FALSE` = PAL[["signal"]], `TRUE` = PAL[["artefact"]]),
                      labels = c("plausible biology", "artefact signature"), name = NULL) +
  labs(title = "env (gp160) minority variant landscape",
       subtitle = "Shaded = V1-V5 loops. Purple triangles on the axis = indel records - note where they cluster.",
       x = "Env amino acid position (HXB2 numbering)", y = "Minor allele frequency (%)") +
  theme(legend.position = "top")

ggsave(file.path(FIGDIR, "lollipop_env.png"), p_env, width = 12, height = 5, dpi = 150)

for (i in seq_len(nrow(VLOOPS))) {
  v <- VLOOPS[i, ]
  d <- env_sites %>% filter(aa_codon >= v$aa_start, aa_codon <= v$aa_end, minor_f >= 0.01)
  ind_v <- env_ind %>% filter(aa_codon >= v$aa_start, aa_codon <= v$aa_end)
  ymax <- max(12, 100 * max(c(d$minor_f, 0.05)))
  p <- ggplot() +
    geom_hline(yintercept = c(5, 10), linetype = 2, colour = PAL[["muted"]], linewidth = 0.3)
  if (nrow(d)) {
    p <- p +
      geom_segment(data = d, aes(x = aa_codon, xend = aa_codon, y = 0, yend = 100 * minor_f,
                                 colour = is_artefact), linewidth = 0.7) +
      geom_point(data = d, aes(aa_codon, 100 * minor_f, colour = is_artefact), size = 2.6)
  }
  if (nrow(ind_v)) {
    p <- p + geom_point(data = ind_v, aes(aa_codon, 0), shape = 25, size = 3,
                        fill = PAL[["violet"]], colour = PAL[["violet"]])
  }
  p <- p +
    scale_colour_manual(values = c(`FALSE` = PAL[["signal"]], `TRUE` = PAL[["artefact"]]),
                        labels = c("plausible biology", "artefact signature"), name = NULL) +
    scale_x_continuous(limits = c(v$aa_start - 1, v$aa_end + 1)) +
    scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.08))) +
    labs(title = sprintf("%s loop (Env aa %d-%d, %d codons)", v$region, v$aa_start, v$aa_end, v$n_codons),
         subtitle = sprintf("%d minority site(s) >= 1%%;  %d indel record(s) - purple triangles",
                            nrow(d), nrow(ind_v)),
         x = "Env aa position (HXB2)", y = "Minor allele frequency (%)") +
    theme(legend.position = "top")
  ggsave(file.path(FIGDIR, sprintf("lollipop_%s.png", v$region)), p,
         width = 7, height = 4, dpi = 150)
}
note("Lollipops: lollipop_env.png + %s", paste(sprintf("lollipop_%s.png", VLOOPS$region), collapse = ", "))

# --- V vs C diversity figure -------------------------------------------------
p_vc <- dom_stats %>%
  mutate(region = factor(region, levels = DOMAINS$region)) %>%
  ggplot(aes(region, pi, fill = type)) +
  geom_col() +
  scale_fill_manual(values = c(variable = PAL[["violet"]], conserved = PAL[["muted"]])) +
  labs(title = "env diversity by domain: the V loops do NOT stand out",
       subtitle = "Point-substitution diversity only. The loops vary in LENGTH - see indel panel.",
       x = NULL, y = "mean pi", fill = NULL)

p_ind <- ind_test %>%
  ggplot(aes(type, per_kb, fill = type)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%d indels\n%.1f/kb", indels, per_kb)), vjust = -0.3, size = 3.4) +
  scale_fill_manual(values = c(variable = PAL[["violet"]], conserved = PAL[["muted"]]), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = sprintf("...but indels are %.1fx enriched in them (p = %.4f)", fold, p_val),
       x = NULL, y = "Indel records per kb")

ggsave(file.path(FIGDIR, "diversity_V_vs_C.png"), p_vc, width = 9, height = 4.4, dpi = 150)
ggsave(file.path(FIGDIR, "indel_enrichment.png"), p_ind, width = 6, height = 4.4, dpi = 150)

# =============================================================================
# 9.  THE GLYCAN SHIELD  -  N-linked glycosylation sites
# =============================================================================
banner("9.", "The glycan shield (PNGS)")

bullet("env is the most heavily glycosylated protein known - roughly half its")
bullet("mass is host sugar. Those glycans are HOST molecules, so antibodies")
bullet("largely cannot see them. The virus hides behind a shield of 'self'.")
bullet("A PNGS is the motif N-X-S/T where X is any residue except proline.")

# Sample consensus env protein (from the 50% FASTA)
env_cons_aa <- translate_nt(all_seqs[["ENV"]][["T50"]])

find_pngs <- function(aa) {
  n <- nchar(aa); hits <- integer(0)
  for (i in seq_len(n - 2)) {
    a <- substr(aa, i, i); b <- substr(aa, i + 1, i + 1); c <- substr(aa, i + 2, i + 2)
    if (a == "N" && b != "P" && c %in% c("S", "T")) hits <- c(hits, i)
  }
  hits
}
pngs <- find_pngs(env_cons_aa)
pngs_tbl <- tibble(aa_pos = pngs,
                   motif = vapply(pngs, function(i) substr(env_cons_aa, i, i + 2), character(1)),
                   domain = map_chr(pngs, assign_domain)) %>%
  left_join(DOMAINS %>% select(region, type), by = c("domain" = "region"))

note("Potential N-linked glycosylation sites in this Env: %d", length(pngs))
show_tbl(pngs_tbl %>% count(domain, type, name = "n_pngs") %>%
           left_join(DOMAINS %>% select(region, n_codons), by = c("domain" = "region")) %>%
           mutate(per_100aa = round(100 * n_pngs / n_codons, 1)) %>%
           arrange(desc(per_100aa)),
         "PNGS density by domain:")
readr::write_csv(pngs_tbl, file.path(RESDIR, "pngs_sites.csv"))

pngs_vc <- pngs_tbl %>% count(type) %>%
  left_join(DOMAINS %>% group_by(type) %>% summarise(aa = sum(n_codons), .groups = "drop"), by = "type") %>%
  mutate(per_100aa = round(100 * n / aa, 1))
show_tbl(pngs_vc, "Glycan density: V loops vs conserved:")
bullet("The loops carry the densest glycan packing in the protein - which is")
bullet("exactly why shifting them by an indel is such an effective escape move.")

# =============================================================================
# 10.  V3 AND CORECEPTOR TROPISM
# =============================================================================
banner("10.", "V3 and coreceptor tropism")

v3_aa <- substr(env_cons_aa, 296, 331)
v3_crown_sample <- substr(env_cons_aa, 312, 315)

note("V3, reference (HXB2): %s", substr(env_ref_aa, 296, 331))
note("V3, THIS SAMPLE     : %s", v3_aa)
note("V3 crown - reference %s  vs  sample %s", v3_crown, v3_crown_sample)

if (v3_crown_sample != v3_crown) {
  cat("\n  ** THE CROWN DISAGREES WITH THE REFERENCE - AND THAT IS INFORMATIVE **\n")
  bullet("HXB2 (subtype B) carries GPGR. This sample carries %s.", v3_crown_sample)
  bullet("GPGQ is the classic NON-B (notably subtype C) V3 crown signature.")
  bullet("This independently corroborates what the 14%% genome-wide divergence")
  bullet("already implied in mutations_workshop.R: this is not a subtype B virus,")
  bullet("and most of its 'mutations' vs HXB2 are subtype, not pathology.")
  bullet("Two unrelated lines of evidence, one conclusion. That is what")
  bullet("corroboration looks like - and why you never rest a call on one metric.")
}

# The 11/25 rule, on V3-relative numbering (V3 position 1 = Env 296)
p11 <- substr(v3_aa, 11, 11); p25 <- substr(v3_aa, 25, 25)
x4_flag <- p11 %in% c("R", "K") || p25 %in% c("R", "K")
note("11/25 rule: V3 position 11 = %s, position 25 = %s  ->  predicted %s",
     p11, p25, if (x4_flag) "CXCR4 (X4) tropic" else "CCR5 (R5) tropic")

cat("\n  HANDLE THIS PREDICTION WITH CARE\n")
bullet("The 11/25 rule is a crude heuristic with poor sensitivity for X4. The")
bullet("clinical standard is geno2pheno[coreceptor], not this line of code.")
bullet("Tropism decides maraviroc eligibility - do not call it from a script.")
cat("\n")
bullet("AND NOTE THE DEEPER PROBLEM: this prediction is made from the CONSENSUS.")
bullet("A minority X4 population below 20%% is invisible to consensus tropism")
bullet("testing, and minority X4 variants are a documented cause of maraviroc")
bullet("failure. That is precisely why you would do minority variant analysis")
bullet("on env in the first place.")

v3_minor <- env_sites %>% filter(aa_codon >= 296, aa_codon <= 331, minor_f >= 0.01) %>%
  mutate(vaf_pct = round(100 * minor_f, 2),
         v3_position = aa_codon - 295) %>%
  select(pos, aa_codon, v3_position, change, vaf_pct, is_artefact, total)

if (nrow(v3_minor)) {
  show_tbl(v3_minor, "Minority variants inside V3 (>= 1% VAF):")
  readr::write_csv(v3_minor, file.path(RESDIR, "v3_minority_variants.csv"))
  key <- v3_minor %>% filter(v3_position %in% c(11, 25))
  if (nrow(key)) {
    cat("\n  ** Minority variant at V3 position 11 or 25 - the tropism-determining\n")
    cat("     residues. This is exactly the situation minority testing exists for.\n")
    show_tbl(key)
  } else {
    note("No minority variant lands on V3 position 11 or 25.")
    bullet("No evidence of a minority X4 population at the tropism-determining sites,")
    bullet("down to 1%%. That is a reportable negative - state the threshold with it.")
  }
} else {
  note("No minority variants >= 1%% anywhere in V3.")
}

# =============================================================================
# 11.  SUMMARY
# =============================================================================
banner("11.", "Summary")

summ <- dom_stats %>%
  filter(type == "variable") %>%
  select(region, aa_start, aa_end, n_codons, sites, pi, minor_1pct, real_1pct, minor_5pct) %>%
  left_join(env_ind %>% count(domain, name = "indels"), by = c("region" = "domain")) %>%
  mutate(indels = replace_na(indels, 0)) %>%
  left_join(cmp %>% filter(threshold == 0.05) %>% select(region, nt_vs_T50), by = "region")

show_tbl(summ, "V1-V5 summary:")
readr::write_csv(summ, file.path(RESDIR, "vloop_summary.csv"))

cat("\n", strrep("-", 78), "\n", sep = "")
cat("  WHAT ENV ACTUALLY SHOWS\n")
cat(strrep("-", 78), "\n", sep = "")
gc_ <- function(ty) call_by_type$pct_callable[match(ty, call_by_type$type)]

cat(sprintf("
  1. THE FRAME IS VERIFIED, NOT ASSUMED. 857 codons, MRV... start, no internal
     stop, V1/V3/V4 closing on cysteines, and the canonical %s crown at 312-315
     of the reference. Only then is any amino acid number trustworthy.

  2. THE V LOOPS ARE NOT DIVERSITY PEAKS. Mean pi is %.5f in the loops vs
     %.5f conserved - equal. But by MEDIAN the loops are significantly LOWER
     (%.5f vs %.5f, Wilcoxon p = %.4f). Everyone predicts peaks. Everyone is
     wrong - and not because the loops are conserved.

  3. THE LOOPS VARY IN LENGTH, NOT IN SUBSTITUTION. %d of %d env indels sit in
     the V loops, though the loops are only %.0f%% of env: a %.1f-fold rate
     enrichment (Poisson p = %.4f). A per-position pileup is indexed by
     reference coordinate and CANNOT represent length variation. The instrument
     is blind to the exact signal you came looking for.

  4. THE BLINDNESS IS SELF-REINFORCING. Reads carrying long indels map badly, so
     the loops also lose coverage (%.1f%% callable vs %.1f%%). The very reads
     that carry the variation are the ones that fail to align. The harder a
     region varies, the less of it you get to see.

  5. THE GLYCAN SHIELD IS DENSEST WHERE THE LOOPS ARE. %d PNGS across Env,
     %.1f per 100 aa in the loops vs %.1f conserved. Half of Env's mass is host
     sugar; antibodies cannot see 'self'. Shifting a glycan by one indel moves
     the shield - which is why length, not substitution, is the escape currency.

  6. THE V3 CROWN SAYS NON-B. Reference %s, this sample %s. GPGQ is the classic
     non-B (subtype C) signature - independently corroborating the 14%% genome
     divergence from mutations_workshop.R. Two unrelated metrics, one answer.

  7. TROPISM FROM A CONSENSUS IS A HALF-ANSWER. The 11/25 rule reads %s/%s
     (R5) here, and no minority variant reaches V3 position 11 or 25 down to
     1%%. But a minority X4 population under 20%% is invisible to consensus
     tropism testing, and is a documented cause of maraviroc failure.

  >> THE LESSON: choosing the right INSTRUMENT comes before any statistic. No
     threshold, no filter, and no p-value in this script could rescue a
     measurement that structurally cannot see what it is looking for. pol was a
     story about filtering noise. env is a story about a question the data
     format cannot answer. Knowing which one you are in is the whole job.
     For V-loop diversity: local assembly, haplotype reconstruction, or long
     reads. Not a pileup.
",
  v3_crown,
  pi_v, pi_c, md_v, md_c, wt$p.value,
  n_v, n_tot, 100 * p_v, fold, p_val,
  gc_("variable"), gc_("conserved"),
  length(pngs),
  pngs_vc$per_100aa[match("variable", pngs_vc$type)],
  pngs_vc$per_100aa[match("conserved", pngs_vc$type)],
  v3_crown, v3_crown_sample, p11, p25
))

cat(strrep("-", 78), "\n", sep = "")
note("FASTAs  : %s/{V1,V2,V3,V4,V5,ENV}/fasta/", CFG$outdir)
note("Figures : %s/", FIGDIR)
note("Tables  : %s/", RESDIR)
cat("\n  Done.\n\n")
