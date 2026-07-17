#!/usr/bin/env bash
# =============================================================================
# Practical 2 - run the ENTIRE workshop end to end.
#
#   bash run_all.sh            # everything (~10 min on a laptop)
#   bash run_all.sh sim        # just one dataset
#
# Instructors: run this once before the session to generate every result and
# figure, so you can show students the finished output and so the /results
# folder is populated if their own runs fail.
#
# Students: you normally run the scripts ONE AT A TIME and read the output.
# This script is the "show me everything" button, not the exercise.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")"

# Two datasets, deliberately chosen to teach OPPOSITE lessons:
#   sim  - simulated, strong clock, KNOWN truth -> "what good looks like"
#   h1n1 - real 2009 pandemic flu, 43-day window -> "low R^2 that is still fine"
[[ $# -gt 0 ]] && DATASETS=("$@") || DATASETS=(sim h1n1)

FAILURES=0
hr() { printf '\n\033[1m=== %s ===\033[0m\n' "$1"; }
fail() { echo "!!! FAILED: $1"; FAILURES=$((FAILURES+1)); }
need() { command -v "$1" >/dev/null || { echo "MISSING TOOL: $1"; exit 1; }; }

need Rscript; need iqtree2; need treetime; need python3

hr "0. Build the datasets"
Rscript  scripts/00_simulate_dataset.R             # Dataset B - simulated
python3  scripts/00b_prepare_h1n1.py               # Dataset C - 2009 H1N1

for D in "${DATASETS[@]}"; do
  hr "1. Inspect the alignment - $D (Exercise 1)"
  Rscript scripts/01_inspect_alignment.R "$D" || fail "inspect ($D)"
done

for D in "${DATASETS[@]}"; do
  hr "2. IQ-TREE ML tree - $D (Exercise 2)"
  bash scripts/02_run_iqtree.sh "$D" || fail "IQ-TREE ($D)"

  hr "3. Plot ML tree + rooting - $D (Exercises 3-4)"
  Rscript scripts/03_plot_ml_tree.R "$D" || fail "plot ($D)"

  hr "4. Root-to-tip regression - $D (Exercises 5-6)"
  Rscript scripts/04_root_to_tip.R "$D" || fail "root-to-tip ($D)"

  hr "5. TreeTime - $D (Exercises 7-8)"
  bash scripts/05_run_treetime.sh "$D" || fail "TreeTime ($D)"

  hr "6. Time tree + metadata + publication figure - $D (Ex 9,10,12,13)"
  Rscript scripts/06_timetree_ggtree.R "$D" || fail "timetree ($D)"

  hr "8. Clusters - $D (Exercise 14)"
  Rscript scripts/08_clusters.R "$D" || fail "clusters ($D)"
done

hr "7. Clock rate comparison (Exercise 11)"
Rscript scripts/07_clock_summary.R || fail "clock summary"

hr "DONE"
if [[ $FAILURES -gt 0 ]]; then
  echo "*** $FAILURES step(s) FAILED - scroll up. Results are incomplete. ***"
  echo "*** Common cause: no disk space. Check with: df -h . ***"
fi
echo "Results : $(ls results/*.csv results/*.txt 2>/dev/null | wc -l | tr -d ' ') tables"
echo "Figures : $(ls figures/*.pdf 2>/dev/null | wc -l | tr -d ' ') PDFs"
echo
echo "Start reading here:"
echo "  notes/WORKSHOP_MANUAL.md       the practical itself"
echo "  notes/INTERPRETATION_GUIDE.md  how to read the numbers"
echo "  results/07_clock_summary.csv   the headline table"
exit $(( FAILURES > 0 ))
