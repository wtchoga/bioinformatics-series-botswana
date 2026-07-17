#!/usr/bin/env Rscript
# =============================================================================
# Practical 2 | Exercise 11 - The molecular clock: read the rate, judge it
#
# Run:  Rscript scripts/07_clock_summary.R      # reads all three datasets
#
# Pulls the rate and R^2 out of every TreeTime run, puts them side by side, and
# compares each against PUBLISHED rates for that organism. The point of the
# exercise is not to produce a number - it is to decide whether the number is
# believable.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
})

runs <- list(
  list(key = "sim",  dir = "results/sim/treetime_sim",
       organism = "SIMULATED",    tag = "B: simulated, 60 seq, 20 yr"),
  list(key = "h1n1", dir = "results/h1n1/treetime_h1n1",
       organism = "Influenza A H1N1pdm09", tag = "C: 2009 H1N1, 50 genomes, 43 d"),
  # Dataset A (real HIV gag) is OPTIONAL - see notes/WORKSHOP_MANUAL.md.
  # If you build it, it reappears here automatically.
  list(key = "real", dir = "results/treetime_real",
       organism = "HIV-1 (gag)",  tag = "A: real HIV-1 gag, 80 seq, 20 yr")
)
runs <- Filter(function(r) dir.exists(r$dir), runs)   # skip what was not run

# --- published rate ranges (substitutions/site/year) -------------------------
# These are the numbers you sanity-check YOUR estimate against. An estimate
# outside the published range is not a discovery - it is a bug, until proven
# otherwise.
published <- tribble(
  ~organism,                ~low,   ~high,  ~typical, ~source_note,
  "HIV-1 (gag)",            1e-3,   3e-3,   2.1e-3,   "typical HIV-1 within-subtype estimates",
  "Influenza A H1N1pdm09",  3e-3,   5e-3,   4.0e-3,   "2009 pandemic H1N1 whole genome",
  "SARS-CoV-2",             5e-4,   1e-3,   8e-4,     "early-pandemic genomes",
  "HBV",                    1e-5,   1e-4,   5e-5,     "long-term HBV estimates",
  "SIMULATED",              1.7e-3, 2.7e-3, 2.2e-3,   "TRUE simulated rate = 2.2e-3"
)

parse_clock <- function(dir) {
  f <- file.path(dir, "molecular_clock.txt")
  if (!file.exists(f)) return(NULL)
  tx <- readLines(f)
  rate_line <- grep("--rate:", tx, value = TRUE)[1]
  r2_line   <- grep("r\\^2", tx, value = TRUE)[1]
  # "--rate:\t2.116e-03 +/- 1.26e-04 (one std-dev)"
  nums <- regmatches(rate_line, gregexpr("[0-9.]+e?-?[0-9]*", rate_line))[[1]]
  rate <- as.numeric(nums[1])
  sd   <- if (length(nums) > 1) as.numeric(nums[2]) else NA_real_
  r2   <- as.numeric(sub(".*:\\s*", "", r2_line))
  list(rate = rate, sd = sd, r2 = r2)
}

get_tmrca <- function(dir) {
  f <- file.path(dir, "dates.tsv")
  if (!file.exists(f)) return(list(tmrca = NA, lo = NA, hi = NA))
  dd <- read.delim(f, comment.char = "#", header = FALSE, fill = TRUE,
                   colClasses = "character",
                   col.names = c("node", "date", "num", "lo", "hi"))
  for (v in c("num", "lo", "hi")) dd[[v]] <- suppressWarnings(as.numeric(dd[[v]]))
  dd <- dd[!is.na(dd$lo), ]
  root <- dd[which.min(dd$num), ]     # the root is the OLDEST dated node
  list(tmrca = root$num, lo = root$lo, hi = root$hi)
}

rows <- list()
for (r in runs) {
  cl <- parse_clock(r$dir)
  if (is.null(cl)) next
  tm <- get_tmrca(r$dir)
  pub <- published |> filter(organism == r$organism)
  verdict <- if (is.na(cl$rate)) "no estimate"
    else if (cl$rate < pub$low)  "BELOW published range"
    else if (cl$rate > pub$high) "ABOVE published range"
    else "consistent with published range"
  rows[[r$key]] <- tibble(
    dataset = r$tag, organism = r$organism,
    rate = cl$rate, rate_sd = cl$sd, r2 = cl$r2,
    tmrca = tm$tmrca, tmrca_lo = tm$lo, tmrca_hi = tm$hi,
    tmrca_ci_width = tm$hi - tm$lo,
    pub_low = pub$low, pub_high = pub$high, verdict = verdict
  )
}
tab <- bind_rows(rows)

cat("\n=========================== EXERCISE 11 ===========================\n")
for (i in seq_len(nrow(tab))) {
  x <- tab[i, ]
  cat("\n", x$dataset, "\n", sep = "")
  cat(sprintf("  rate           : %.3e", x$rate))
  if (!is.na(x$rate_sd)) cat(sprintf(" +/- %.2e", x$rate_sd))
  cat(sprintf("   [published %.1e - %.1e]\n", x$pub_low, x$pub_high))
  cat(sprintf("  verdict        : %s\n", x$verdict))
  cat(sprintf("  root-to-tip R2 : %.2f\n", x$r2))
  cat(sprintf("  tMRCA          : %.2f  (90%% region %.2f - %.2f, width %.1f yr)\n",
              x$tmrca, x$tmrca_lo, x$tmrca_hi, x$tmrca_ci_width))
  cat(sprintf("  precision      : %s\n",
      ifelse(x$tmrca_ci_width < 5, "good - a usable interval",
      ifelse(x$tmrca_ci_width < 15, "poor - wide interval, treat as indicative",
             "USELESS - the interval is so wide it excludes nothing"))))
}

cat("\n---------------------------------------------------------------\n")
cat("LOW R2 IS A SYMPTOM, NOT A DIAGNOSIS. ASK *WHY* IT IS LOW.\n\n")
cat("  Dataset C (H1N1) has R2 ~ 0.2 - 'weak' by the usual lookup table - and\n")
cat("  yet it recovers a textbook pandemic rate (~4e-3) AND a tMRCA of ~2009.0\n")
cat("  to within a few weeks, matching the published estimate for the 2009\n")
cat("  pandemic. Its R2 is low only because 43 days of sampling gives the\n")
cat("  x-axis almost no spread to lever against. The clock is real, and this\n")
cat("  analysis is USABLE.\n\n")
cat("  Compare that with Dataset B, where R2 ~ 0.8 because we simulated 20\n")
cat("  years of sampling. Same organism speed, four times the R2 - the\n")
cat("  difference is the SAMPLING WINDOW, not the biology.\n\n")
cat("  So: a low R2 caused by a short window can still give a good rate.\n")
cat("  A low R2 caused by ancient, non-clocklike divergence (as in the\n")
cat("  optional real HIV dataset, where the tMRCA lands centuries before HIV\n")
cat("  existed) cannot. Read rate, R2, sampling window and tMRCA interval\n")
cat("  TOGETHER - never the R2 alone.\n")
cat("---------------------------------------------------------------\n")

write_csv(tab, "results/07_clock_summary.csv")

# --- figure: estimate vs published range ------------------------------------
plt <- tab |> filter(!is.na(rate))
p <- ggplot(plt, aes(y = reorder(dataset, rate))) +
  geom_linerange(aes(xmin = pub_low, xmax = pub_high), linewidth = 6,
                 colour = "grey85") +
  geom_point(aes(x = rate, colour = verdict), size = 3.5) +
  geom_errorbarh(aes(xmin = pmax(rate - rate_sd, 1e-6), xmax = rate + rate_sd),
                 height = 0.15, na.rm = TRUE) +
  scale_x_log10() +
  scale_colour_manual(values = c("consistent with published range" = "#27ae60",
                                 "BELOW published range" = "#c0392b",
                                 "ABOVE published range" = "#c0392b"),
                      name = NULL) +
  labs(title = "Estimated clock rate vs published range",
       subtitle = "Grey bar = published range for that organism. Point = your TreeTime estimate (+/- 1 sd).",
       x = "substitutions / site / year (log scale)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 8))

ggsave("figures/07_clock_rates.pdf", p, width = 9, height = 4)

# --- figure: tMRCA intervals ------------------------------------------------
p2 <- ggplot(plt, aes(y = reorder(dataset, -tmrca_ci_width))) +
  geom_errorbarh(aes(xmin = tmrca_lo, xmax = tmrca_hi), height = 0.2,
                 colour = "#2c3e50") +
  geom_point(aes(x = tmrca), size = 3, colour = "#c0392b") +
  labs(title = "Estimated tMRCA and its 90% interval",
       subtitle = "A tMRCA is only as good as its interval. Wide interval = no answer.",
       x = "Year", y = NULL) +
  theme_bw(base_size = 10) +
  theme(plot.subtitle = element_text(size = 8))

ggsave("figures/07_tmrca_intervals.pdf", p2, width = 9, height = 4)

cat("\nWrote: results/07_clock_summary.csv\n")
cat("Wrote: figures/07_clock_rates.pdf\n")
cat("Wrote: figures/07_tmrca_intervals.pdf\n")

# =============================================================================
# YOUR TURN
#   1. Fill in the table in your report: rate, R2, tMRCA, CI, for each dataset.
#   2. Dataset B's TRUE rate is 2.2e-3 (see data/simulated/SIM_truth.txt).
#      Did your estimate's error bar actually contain the truth?
#   3. Dataset A's tMRCA is centuries before HIV-1 is known to have emerged in
#      humans (~1920s for group M). What does that tell you about the estimate?
#   4. Which of these three would you put in a paper? Which would you not?
# =============================================================================
