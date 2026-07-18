# Practical 2 — Molecular Evolution and Phylogenetics
## From Sequence Alignment to Time-Scaled Phylogenetic Analysis

**Duration:** 4–6 hours · **Level:** beginner · **You need:** R, IQ-TREE2, TreeTime, (optionally FigTree)

---

## Learning objectives

By the end you should be able to:

1. Perform Maximum Likelihood phylogenetic analysis
2. Assess branch support
3. Root trees correctly
4. **Evaluate temporal signal — and know when it is lying to you**
5. Estimate a molecular clock
6. Infer time-scaled phylogenies
7. Export publication-quality trees

Objective 4 is the one that matters. Anyone can make a tree. This practical is
about deciding whether to believe it.

---

## The datasets, and why there are two

Most workshops hand you one dataset that works, you follow the steps, everything
succeeds, and you learn nothing about judgement. We use **two real analyses that
disagree with each other**.

| | **Dataset B — simulated** | **Dataset C — 2009 pandemic H1N1** |
|---|---|---|
| What | A simulated HIV-like epidemic | **Real** influenza A(H1N1)pdm09 whole genomes |
| Size | 60 sequences × 1,500 bp | 50 genomes × 13,109 bp |
| Sampling | 2005–2025 (**20 years**) | Apr–May 2009 (**43 days**) |
| Truth | **Known** — we chose it | Unknown, but published estimates exist |
| Teaches | What good looks like | Why "weak signal" ≠ "bad analysis" |

**Why a simulated dataset?** Because you cannot learn to read a molecular clock
on data that has no clock, and you cannot mark your own answer without knowing
the truth. For Dataset B we *chose* the rate (2.2 × 10⁻³) and the tMRCA (1992.7),
so at the end you can check whether your estimate actually contains the right
answer. Everything is in `data/simulated/SIM_truth.txt` — **don't peek until
Exercise 11.**

> ⚠️ **Dataset B is simulated. It is not from any patient. Never report its
> numbers as real epidemiological results.**

**Why real H1N1?** Because the 2009 pandemic is one of the best-characterised
outbreaks in history, so you can check your answer against the literature. It is
also a genuinely awkward dataset — 6 weeks of sampling — and it will score
"weak temporal signal" while being perfectly usable. Learning to tell that apart
is the point.

**Optional Dataset A — real HIV-1 gag** (80 sequences, 6 African countries,
2005–2025). It is the *negative control*: unlinked cross-sectional HIV has
almost no temporal signal (R² ≈ 0.06) and TreeTime dates its tMRCA to **1793**,
centuries before HIV existed in humans. If you want it, run
`python3 scripts/00a_prepare_hiv.py` (needs `GAG_sdb_n417.fasta`), then pass
`real` to any script. Highly recommended if you have time — nothing teaches
scepticism faster.

---

## Folder structure

```
Phylogenetics_Workshop/
├── data/
│   ├── simulated/        Dataset B  (+ SIM_truth.txt = the answer key)
│   └── h1n1/             Dataset C
├── scripts/              00–08, run in order
├── results/              tables + tree files (generated)
├── figures/              PDFs (generated)
├── notes/
│   ├── WORKSHOP_MANUAL.md      <- you are here
│   ├── INTERPRETATION_GUIDE.md <- how to read the numbers. Keep it open.
│   └── INSTRUCTOR_KEY.md       <- instructors only
└── report/               your write-up
```

**Run everything from the `Phylogenetics_Workshop/` folder**, not from inside
`scripts/`. Paths are relative to the top.

Fastest way to see it all work: `bash run_all.sh` (~10 min). But do the
exercises one at a time — the output is the lesson.

---

## Exercise 1 — Inspect the alignment

```bash
Rscript scripts/01_inspect_alignment.R sim
Rscript scripts/01_inspect_alignment.R h1n1
```

Before you build a tree, look at your data. Almost every strange tree is a data
problem wearing a costume.

**Answer these:**
- How many sequences? Average length? How much missing data?
- Which sample is shortest? Which has the most Ns?
- Any exact duplicates?
- **What fraction of sites actually vary?**

**Look carefully at that last one for H1N1.** It is 13,109 bp — sounds
enormous — but **98.1% of sites are constant** and only **47 sites are
parsimony-informative**. Predict now what that will do to your bootstrap values
in Exercise 2, then check whether you were right.

> **Concept — the alignment is an assumption.** Column 100 must mean "the same
> position" in every sequence. If `width()` returns more than one number, your
> file isn't aligned and everything downstream is nonsense.

---

## Exercise 2 — Build the ML tree (IQ-TREE2)

```bash
bash scripts/02_run_iqtree.sh sim
bash scripts/02_run_iqtree.sh h1n1
```

```
iqtree2 -s <alignment> -m MFP -B 1000 -alrt 1000
```

- `-m MFP` — ModelFinder tests substitution models and picks the best by BIC.
  **Don't just assume GTR.**
- `-B 1000` — ultrafast bootstrap. Read as: **≥95 = strong**.
- `-alrt 1000` — SH-aLRT, a *different* support test. Read as: **≥80 = strong**.
  Trust a clade most when **both** are high.

**Answer these:**
- Which model was chosen? What do its terms mean?
  (`TIM3+F` = transition model 3 + empirical base Frequencies. `+I` = a
  proportion of Invariant sites. `+R6`/`+G4` = rate variation across sites.)
- How many clusters can you see? Which tip has the longest branch?
- **How many nodes reach UFboot ≥ 95 in each dataset?**

You should find ~43/58 for the simulated data but only **~5/48** for H1N1.
Same software, same settings. The difference is *information in the alignment* —
those 47 informative sites. **You cannot bootstrap your way to evidence that
isn't there.**

---

## Exercise 3 — Visualise the ML tree in R

```bash
Rscript scripts/03_plot_ml_tree.R sim
Rscript scripts/03_plot_ml_tree.R h1n1
```

Colour tips by country (B) or region (C). Output: `figures/03_ML_tree_*.pdf`.

**Answer these:**
- Do sequences from the same place group together? Should they?
- Find the longest terminal branch. Cross-check it against
  `results/01_alignment_qc_*.csv` — bad sequence, or genuinely divergent virus?

---

## Exercise 4 — Rooting

An unrooted tree shows relationships but **not the direction of time**. You
cannot say "this came first" until you root it.

Three options, best first:
1. **Outgroup** — a sequence you know sits outside the group. Best, if you have one.
2. **Midpoint** — root at the middle of the longest tip-to-tip path. Assumes a
   roughly clocklike tree.
3. **Best-fit** — the root that maximises temporal signal (Exercise 5). **Read
   the warning there before you trust it.**

Script 03 writes a midpoint-rooted tree. In **FigTree** you can do the same by
clicking: open the `.treefile`, *Midpoint root*, turn on node labels to show
support, collapse weak clades, export PDF.

> **Concept.** Rooting is a *hypothesis*, not an observation. Midpoint rooting
> assumes a clock — so don't use it to "prove" there is one. That's circular.

---

## Exercise 5 — Root-to-tip regression (temporal signal)

```bash
Rscript scripts/04_root_to_tip.R sim
Rscript scripts/04_root_to_tip.R h1n1
```

This is TempEst, reimplemented in R so everyone can run it. (If you have
TempEst, load the `.treefile`, set dates from tip labels, click *Best-fitting
root*, and compare.)

Plot each tip's distance from the root against its sampling date:
- **slope** = evolutionary rate (subs/site/year)
- **x-intercept** = tMRCA
- **R²** = how much of the divergence the dates explain

**Record for both datasets:** R², correlation, slope, x-intercept, residual
mean, best root.

> ⚠️ **This regression is not a statistical test.** Tips share ancestry, so the
> points are not independent. Use it as a *diagnostic*, never as a p-value.

---

## Exercise 6 — Interpret the regression

Fill this in:

| Statistic | Dataset B (sim) | Dataset C (H1N1) |
|---|---|---|
| Slope | | |
| R² | | |
| Residual mean | | |
| Best root | | |

**Now the real questions:**

1. Dataset C's R² is ~0.2. The standard table calls that "very weak — dating
   unreliable". **Is that verdict correct here?** Look at the sampling window
   before you answer.
2. **What happens if R² = 0.03?** (Dataset A scores ~0.06 if you built it.)
   Is that the *same* problem as Dataset C's 0.19, or a different one?
3. Compare the `midpoint` and `best-fit` columns the script prints. How much of
   your R² was bought by *searching* for a root?
4. **Outliers.** Which tips have |z| > 2? Dataset B has exactly one dramatic
   outlier (z ≈ −5.4). What does the script report when you drop it?

> **This is the heart of the practical.** See §1b and §1c of
> `INTERPRETATION_GUIDE.md`. Two datasets both score "weak R²" for completely
> different reasons; one is fine and one is junk. If you only read the R² table,
> you get it exactly backwards.

---

## Exercise 7 — Prepare TreeTime

TreeTime needs three things: the **tree**, the **alignment** it came from, and
the **dates**.

`dates.tsv` — names must match the tree tips *exactly*:

```
name	date
SIM_BWA_001_2006.42	2006.42
SIM_BWA_002_2005.64	2005.64
```

Decimal dates beat year-only dates: `2021.23` ≈ late March 2021. Year-only
dates say "January 1st", which quietly injects error up to a year.

> **The #1 TreeTime failure is a name mismatch.** It silently drops tips it
> cannot match, then confidently gives you a wrong answer. Script 05 checks this
> for you *before* running and refuses to continue if names don't match. Do this
> check yourself, always.

---

## Exercise 8 — Time-scaled phylogeny

```bash
bash scripts/05_run_treetime.sh sim
bash scripts/05_run_treetime.sh h1n1
```

Key flags (read the comments in the script):
- `--reroot best` — see the Exercise 5 warning.
- `--coalescent const` — a prior on node times.
- `--confidence` + `--covariation` — **always**. Tips are not independent;
  `--covariation` accounts for that and gives you a rate **with an error bar**
  and dates **with intervals**. Without it you get false confidence.
- `--clock-filter 3` — flags tips far off the clock (see `outliers.tsv`).

**Outputs:** `timetree.nexus` · `molecular_clock.txt` · `dates.tsv` ·
`ancestral_sequences.fasta` · `branch_mutations.txt`

> **Reproducibility trap we hit ourselves:** re-running TreeTime into an
> existing output directory can silently reuse old state — we watched a rate
> change from 3.4e-3 to 1.1e-3 with no error message. Script 05 wipes the
> directory first. If your numbers move for no reason, that's why.

---

## Exercise 9 — Read the time tree, estimate tMRCA

```bash
Rscript scripts/06_timetree_ggtree.R sim
Rscript scripts/06_timetree_ggtree.R h1n1
```

Open `timetree.nexus` in FigTree too.

**Answer these:**
- What is your tMRCA, **and its interval**?
- For H1N1, compare to the published 2009 pandemic origin (≈ late 2008 /
  early 2009). How close are you?
- Most recent transmission? Longest lineage?

> **Find the root properly.** `NODE_0000000` is usually **not** the root. Take
> a tip's date minus its root-to-tip distance, or read the oldest dated node in
> TreeTime's `dates.tsv`.

---

## Exercise 10 — Visualise the time tree (ggtree)

Script 06 does this. The x-axis is now **years**, not substitutions.

Compare `figures/06_timetree_sim.pdf` with `figures/03_ML_tree_sim.pdf`. Same
data, same topology — why do they look so different?

**Note what script 06 prints for a weak dataset.** TreeTime will hand you a
beautiful, fully-dated, utterly convincing tree for data with no temporal signal
whatsoever. **A pretty figure is not evidence.**

---

## Exercise 11 — The molecular clock

```bash
Rscript scripts/07_clock_summary.R
```

Compares every dataset's rate against **published ranges** for that organism,
and plots the tMRCA intervals.

**Answer these:**
1. What rate did you get for each dataset? Is it in the published range?
2. **Now open `data/simulated/SIM_truth.txt`.** The true rate is 2.2 × 10⁻³ and
   the true tMRCA is 1992.7. **Did your error bars actually contain the truth?**
3. Why is the H1N1 rate well estimated while its R² is "weak"?
4. Which of these would you put in a paper? Which would you not?

---

## Exercise 12 — Publication figure

Script 06 writes `figures/Figure1_*.pdf`. Rebuild it yourself with `ggsave`:

```r
ggsave("Figure1.pdf", width = 12, height = 8)
```

A good figure states: what the tree is, what the axis means, what the colours
mean, and how well supported it is. If a reader can't tell whether to trust it,
it's decoration.

---

## Exercise 13 — Metadata mapping

`%<+%` joins a data frame to a tree by matching the first column to tip labels.

```r
p <- ggtree(tree) %<+% metadata + geom_tippoint(aes(colour = country))
```

Recolour by different columns (`country`, `region`, `cluster`). Does geography
explain the tree shape, or not?

> **If your colours come out all grey, the names didn't match.** That's the same
> bug as the TreeTime one, wearing a different hat.

---

## Exercise 14 — Cluster detection

```bash
Rscript scripts/08_clusters.R sim
Rscript scripts/08_clusters.R h1n1
```

A cluster = a clade where every pair is within a **distance threshold** and the
clade is **well supported**.

**Answer these:**
1. How many clusters did you find?
2. Open `figures/08_cluster_threshold_*.pdf`. **How many would you have found
   with a threshold half as big?** Does "3 transmission clusters" mean anything
   without stating the threshold?
3. For Dataset B, the script marks you against the known truth. What are your
   precision and recall? Which true clusters did you miss, and why?
4. Drop `MIN_BOOT` to 70 and rerun. You find more clusters — are they real?

> **A cluster is not proof that A infected B.** The real infector may be
> unsampled, and direction of transmission cannot be read off this tree. These
> are real people; inferred linkage has been misused in court. Say "consistent
> with", not "proves".

---

## Exercise 15 — Final interpretation

Answer for **each** dataset, with `INTERPRETATION_GUIDE.md` open:

1. Is there temporal signal? **And if it's weak — why?**
2. Estimated substitution rate? Consistent with published values?
3. Estimated tMRCA? **Is its interval narrow enough to mean anything?**
4. Which sequences cluster? At what threshold and support?
5. Evidence of transmission — how strong, and what would change your mind?
6. Outliers — which, and what's the likely cause?
7. Which samples should be removed, and **why**? What happens to your
   conclusions if you keep them?

---

## Deliverables

1. ML tree (`.treefile`) for both datasets
2. Rooted tree, FigTree or ggtree PDF
3. Root-to-tip regression plot + interpretation (TempEst screenshot if used)
4. TreeTime time-scaled tree (`.nexus`)
5. Publication-quality tree (ggtree)
6. Metadata-annotated phylogeny
7. **A 2–3 page report** covering evolutionary relationships, temporal signal,
   and molecular clock estimates — including, explicitly:
   - whether your Dataset B estimates **contained the known truth**
   - why Dataset C's R² is low but its rate is still good
   - the thresholds you chose, and what changes if you choose differently

Template: `report/report_template.md`

---

## Extensions

- **Model selection** — compare the models ModelFinder chose for B vs C. Why different?
- **Support** — where do SH-aLRT and UFboot disagree, and what does that mean?
- **Tree comparison** — ML tree vs time-scaled tree: what changed, what didn't?
- **Ancestral reconstruction** — interpret `ancestral_sequences.fasta` and
  `branch_mutations.txt`.
- **Phylogeography** — use metadata to infer geographic spread. What confounds this?
- **BEAST comparison** — contrast TreeTime with a Bayesian analysis: assumptions,
  computational cost, when each is appropriate.
- **The negative control** — build Dataset A (`scripts/00a_prepare_hiv.py`) and
  explain why real cross-sectional HIV dates its own ancestor to 1793.
- **Sampling window** — resimulate Dataset B with `SAMPLE_END <- 2010` (5-year
  window) in `scripts/00_simulate_dataset.R`. Same clock, same organism. Watch
  R² collapse. **That is Exercise 6's lesson, on demand.**

---

## Common problems

| Symptom | Cause |
|---|---|
| "0 sequences" from a `.fasta` | It may not be FASTA. Dataset C is NEXUS with old-Mac CR line endings — see `scripts/00b_prepare_h1n1.py`. |
| TreeTime drops tips silently | Name mismatch between tree and `dates.tsv`. |
| Numbers change between identical runs | Stale TreeTime output directory. Delete it. |
| All tips grey in ggtree | Metadata names don't match tip labels. |
| Negative slope in root-to-tip | No temporal signal, or wrong dates. Stop; do not date. |
| `Error in width(seqs)` / multiple widths | Your file isn't aligned. |
| Very low bootstrap everywhere | Not enough variable sites. Not fixable by rerunning. |

---

## The one thing to take away

> The software will **always** give you a tree. It will be beautiful, fully
> annotated, and completely confident. It will give you a tMRCA of 1793 for a
> virus that emerged in the 20th century without a single warning.
>
> **Your job is not to produce the tree. It is to decide whether to believe it.**
