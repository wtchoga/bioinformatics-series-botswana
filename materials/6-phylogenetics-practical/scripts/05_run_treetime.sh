#!/usr/bin/env bash
# =============================================================================
# Practical 2 | Exercise 7 - Prepare TreeTime
#                Exercise 8 - Time-scaled phylogeny
#
# Run:  bash scripts/05_run_treetime.sh real     # Dataset A - real HIV gag
#       bash scripts/05_run_treetime.sh sim      # Dataset B - simulated
#       bash scripts/05_run_treetime.sh h1n1     # Dataset C - 2009 H1N1
#
# TreeTime takes THREE things:
#     1. a tree       (the ML tree from IQ-TREE)
#     2. an alignment (the same one you built the tree from)
#     3. dates        (a TSV: name <TAB> date)
# and rescales the tree from "substitutions" into "years".
#
# dates.tsv MUST look like this, and the names MUST match the tree tips exactly:
#     name<TAB>date
#     BWA_WS001_2015<TAB>2015
#     BWA_WS002_2016<TAB>2016
#
# Decimal dates (2021.23) are better than years (2021) when you have them.
# 2021.23 means ~23% of the way through 2021, i.e. late March.
# =============================================================================
set -euo pipefail

DATASET="${1:-real}"

case "$DATASET" in
  real) TREE="results/HIV_alignment.fasta.treefile"
        ALN="results/HIV_alignment.fasta"
        DATES="data/dates.tsv"; OUT="results/treetime_real" ;;
  sim)  TREE="results/sim/SIM_alignment.fasta.treefile"
        ALN="results/sim/SIM_alignment.fasta"
        DATES="data/simulated/SIM_dates.tsv"; OUT="results/sim/treetime_sim" ;;
  h1n1) TREE="results/h1n1/H1N1_alignment.fasta.treefile"
        ALN="results/h1n1/H1N1_alignment.fasta"
        DATES="data/h1n1/H1N1_dates.tsv"; OUT="results/h1n1/treetime_h1n1" ;;
  *)    echo "usage: bash scripts/05_run_treetime.sh [real|sim|h1n1]"; exit 1 ;;
esac

command -v treetime >/dev/null || { echo "ERROR: treetime not found."; exit 1; }

# Start from a clean output directory. If you re-run TreeTime into a directory
# that already has results, it can pick up the old state and report numbers
# from your PREVIOUS run - we hit exactly that while building this workshop
# (a rate silently changed from 3.4e-3 to 1.1e-3 with no error message).
# If your numbers move between runs for no reason, this is why.
rm -rf "$OUT"

# --- sanity check BEFORE running: do the names match? -----------------------
# The single most common TreeTime failure is a name mismatch between the tree
# and dates.tsv. Check it yourself - TreeTime will silently drop what it
# cannot match, and then quietly give you a wrong answer.
python3 - "$TREE" "$DATES" <<'EOF'
import re, sys
tree = open(sys.argv[1]).read()
tips = set(re.findall(r'[(,]\s*([^(),:]+)\s*:', tree))
dates = {l.split('\t')[0] for l in open(sys.argv[2]).read().splitlines()[1:] if l.strip()}
missing = tips - dates
print(f"tips in tree      : {len(tips)}")
print(f"names in dates.tsv: {len(dates)}")
print(f"tips WITHOUT a date: {len(missing)}")
if missing:
    print("  e.g.", list(missing)[:3])
    sys.exit("FATAL: fix the name mismatch before running TreeTime.")
print("name check passed.\n")
EOF

# -----------------------------------------------------------------------------
# FLAGS
#   --reroot best     find the root that best fits the dates. Convenient, but
#                     remember it OPTIMISES the fit - see script 04's warning.
#   --coalescent const  use a constant-size coalescent prior on node times.
#                     Regularises the time tree. 'skyline' allows changing Ne.
#   --confidence      report confidence intervals on node dates. ALWAYS ask
#                     for these - a tMRCA without a CI is not a result.
#   --clock-filter 3  flag tips >3 interquartile ranges off the clock and drop
#                     them. Set to 0 to keep everything (do this first, so you
#                     SEE your outliers before deciding). Dropped tips are
#                     listed in outliers.tsv.
#   --rng-seed N      TreeTime is STOCHASTIC (it resolves polytomies randomly
#                     and does marginal inference). Without a fixed seed your
#                     rate and tMRCA shift slightly on every run. We pin it so
#                     everyone's numbers match the notes. In a real study, run
#                     it several times WITHOUT a seed and check your
#                     conclusions are stable across runs - if they are not,
#                     that instability IS your result.
#   --covariation     THE IMPORTANT ONE. Tips share ancestry, so their
#                     root-to-tip distances are NOT independent. A naive
#                     regression ignores that and gives an over-confident,
#                     biased rate. --covariation accounts for it, and in
#                     exchange gives you a rate WITH a standard deviation and
#                     date confidence bounds in dates.tsv. Without it you get
#                     a point estimate and no honest error bar.
# -----------------------------------------------------------------------------
treetime \
  --tree "$TREE" \
  --aln "$ALN" \
  --dates "$DATES" \
  --outdir "$OUT" \
  --reroot best \
  --rng-seed 20260716 \
  --coalescent const \
  --confidence \
  --covariation \
  --clock-filter 3

echo
echo "=============================================================="
echo "OUTPUTS in $OUT/"
echo "  timetree.nexus            time-scaled tree -> open in FigTree"
echo "  timetree.pdf              quick look"
echo "  molecular_clock.txt       the RATE and R^2   <- Exercise 11"
echo "  root_to_tip_regression.pdf  TreeTime's own version of Exercise 5"
echo "  ancestral_sequences.fasta reconstructed ancestors"
echo "  dates.tsv                 inferred date + 90% bounds for EVERY node"
echo "  outliers.tsv              tips dropped by --clock-filter"
echo "  branch_mutations.txt      which mutations happened on which branch"
echo "=============================================================="
echo
echo "--- molecular_clock.txt ---"
cat "$OUT/molecular_clock.txt" 2>/dev/null | head -6
echo
echo "QUESTIONS"
echo "  1. What rate did TreeTime estimate? Compare it to the published range"
echo "     for this organism (see notes/INTERPRETATION_GUIDE.md)."
echo "  2. Find the root's date in $OUT/dates.tsv - that is your tMRCA."
echo "  3. Did --clock-filter drop any tips? Which, and do you agree?"
