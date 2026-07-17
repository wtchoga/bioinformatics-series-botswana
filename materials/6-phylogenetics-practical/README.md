# Practical 2 — Molecular Evolution and Phylogenetics

A beginner's workshop taking you from a sequence alignment to a time-scaled
phylogeny — and, more importantly, to knowing whether to believe it.

**Start here:** [`notes/WORKSHOP_MANUAL.md`](notes/WORKSHOP_MANUAL.md)
**Keep open:** [`notes/INTERPRETATION_GUIDE.md`](notes/INTERPRETATION_GUIDE.md)
**Instructors:** [`notes/INSTRUCTOR_KEY.md`](notes/INSTRUCTOR_KEY.md) *(answers — don't hand out early)*

---

## Quick start

```bash
cd Phylogenetics_Workshop
bash run_all.sh          # regenerates every result and figure (~10 min)
```

Then work through the exercises one at a time — the printed output *is* the
lesson.

**Needs:** `Rscript`, `iqtree2`, `treetime`, `python3`
**R packages:** Biostrings, ape, phangorn, phytools, treeio, ggtree, ggplot2,
dplyr, tidyr, readr, ggrepel, scales, lubridate

---

## Two datasets that disagree

| | **B — simulated** | **C — 2009 pandemic H1N1** |
|---|---|---|
| 60 seq × 1,500 bp | 50 genomes × 13,109 bp |
| Sampled over **20 years** | Sampled over **43 days** |
| Truth **known** (rate 2.2e-3, tMRCA 1992.7) | Compare against published estimates |
| R² ≈ **0.80** → "trust it" | R² ≈ **0.20** → "weak signal" |
| Recovers the truth | **Recovers the truth anyway** |

Dataset C scores "very weak" on every R² lookup table, and then produces a
textbook pandemic H1N1 rate and a tMRCA matching the published 2009 origin to
within weeks. Its R² is low only because 43 days of sampling gives the
regression nothing to lever against.

**That contradiction is the workshop.** Low R² is a symptom, not a diagnosis —
and if you read the table without asking *why*, you throw away the good analysis.

An optional third dataset (real HIV-1 *gag*) is the negative control: it dates
its own ancestor to **1793**, centuries before HIV existed in humans, without a
single warning from the software. Build it with
`python3 scripts/00a_prepare_hiv.py`.

---

## Scripts

| Script | Exercise | Does |
|---|---|---|
| `00_simulate_dataset.R` | — | Builds Dataset B + the answer key |
| `00b_prepare_h1n1.py` | — | Builds Dataset C (source is NEXUS, not FASTA) |
| `00a_prepare_hiv.py` | — | Builds optional Dataset A |
| `01_inspect_alignment.R` | 1 | QC: length, Ns, duplicates, variable sites |
| `02_run_iqtree.sh` | 2 | ML tree + bootstrap + SH-aLRT |
| `03_plot_ml_tree.R` | 3–4 | Plot, colour by metadata, root |
| `04_root_to_tip.R` | 5–6 | Root-to-tip regression (TempEst in R) |
| `05_run_treetime.sh` | 7–8 | Time-scaled tree |
| `06_timetree_ggtree.R` | 9,10,12,13 | tMRCA, ggtree, publication figure |
| `07_clock_summary.R` | 11 | Rates vs published ranges |
| `08_clusters.R` | 14 | Cluster detection + threshold sweep |

All scripts take a dataset argument: `sim`, `h1n1`, or `real`.

```bash
Rscript scripts/04_root_to_tip.R h1n1
```

Run everything from this folder, not from inside `scripts/`.

---

## Layout

```
data/simulated/   Dataset B (+ SIM_truth.txt = answers)
data/h1n1/        Dataset C
data/README_DATA.md   <- provenance: what is real, what is simulated. Read it.
scripts/          00-08
results/          tables + trees (generated)
figures/          PDFs (generated)
notes/            manual, interpretation guide, instructor key
report/           report template
```

---

## The one thing to take away

> The software will **always** give you a tree. It will be beautiful, fully
> annotated, and completely confident. It will hand you a tMRCA of 1793 for a
> virus that emerged in the 20th century without complaint.
>
> **Your job is not to produce the tree. It is to decide whether to believe it.**
