#!/usr/bin/env bash
# =============================================================================
# Practical 2 | Exercise 2 - Build the Maximum Likelihood tree with IQ-TREE2
#
# Run from the Phylogenetics_Workshop/ folder:
#     bash scripts/02_run_iqtree.sh real     # Dataset A - real HIV gag
#     bash scripts/02_run_iqtree.sh sim      # Dataset B - simulated
#     bash scripts/02_run_iqtree.sh h1n1     # Dataset C - 2009 H1N1
#
# Runtime: ~1 minute for 80 sequences on a laptop.
# =============================================================================
set -euo pipefail

DATASET="${1:-real}"

case "$DATASET" in
  real) ALN="data/HIV_alignment.fasta";           OUTDIR="results" ;;
  sim)  ALN="data/simulated/SIM_alignment.fasta"; OUTDIR="results/sim" ;;
  h1n1) ALN="data/h1n1/H1N1_alignment.fasta";     OUTDIR="results/h1n1" ;;
  *)    echo "usage: bash scripts/02_run_iqtree.sh [real|sim|h1n1]"; exit 1 ;;
esac

command -v iqtree2 >/dev/null || { echo "ERROR: iqtree2 not found on PATH."; exit 1; }
mkdir -p "$OUTDIR"
cp "$ALN" "$OUTDIR"/
BASE="$OUTDIR/$(basename "$ALN")"

# -----------------------------------------------------------------------------
# WHAT EACH FLAG DOES  (this is the bit students should actually understand)
#
#   -s        the input ALIGNMENT. Must be aligned - IQ-TREE will not align it.
#   -m MFP    ModelFinder Plus. Tests substitution models and picks the best by
#             BIC, then builds the tree with it. Do not just assume GTR.
#   -B 1000   Ultrafast bootstrap, 1000 replicates. Branch support.
#             Read these as: >=95 = strong support.
#   -alrt 1000  SH-aLRT branch test, 1000 replicates. A SECOND, different
#             support measure. Read as: >=80 = strong support.
#             Trust a clade most when BOTH are high (>=95 UFboot AND >=80 aLRT).
#   -T AUTO   use a sensible number of CPU threads.
#   -redo     overwrite previous results.
#   -seed N   IQ-TREE starts its tree search from RANDOM starting trees. Without
#             a fixed seed you get a slightly different tree every run, and
#             every number downstream (R^2, rate, tMRCA) shifts a little. That
#             is normal stochasticity, not a bug - but for a workshop we pin it
#             so everyone's numbers match the notes. In a real study, run it
#             several times WITHOUT a seed and check your conclusions are stable.
# -----------------------------------------------------------------------------
iqtree2 \
  -s "$BASE" \
  -m MFP \
  -B 1000 \
  -alrt 1000 \
  -T AUTO \
  -seed 20260716 \
  -redo

echo
echo "=============================================================="
echo "KEY OUTPUT FILES ( in $OUTDIR/ )"
echo "  *.treefile  the ML tree, Newick. Node labels are 'aLRT/UFboot'."
echo "  *.iqtree    human-readable report: chosen model, likelihood, stats."
echo "  *.log       full run log."
echo "  *.contree   consensus tree."
echo "=============================================================="
echo
echo "Best-fit model chosen by ModelFinder:"
grep "Best-fit model" "$BASE.iqtree" || true
echo
echo "QUESTIONS (Exercise 2)"
echo "  1. Which model was selected, and what do its terms mean?"
echo "     e.g. GTR+F+R6 = General Time Reversible, empirical base"
echo "     Frequencies, 6-category FreeRate rate heterogeneity."
echo "  2. Open the .treefile. How many clusters can you see?"
echo "  3. Which tip has the longest branch? Is it real, or bad sequence?"
