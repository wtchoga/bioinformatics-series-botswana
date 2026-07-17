# Final Interpretation Guide
### How to read the numbers — Practical 2

Students can nearly always *generate* a phylogeny. What they struggle with is
knowing whether the numbers mean anything. This guide gives you the ranges, and
then tells you when the ranges lie.

> **Read this first.** Every threshold below is a **rule of thumb, not a law**.
> The single most important skill in this practical is asking *why* a number
> came out the way it did — not looking it up in a table. Section 1b exists
> because this workshop's own data breaks the table on page one.

---

## 1. Is there temporal signal?

**Where:** root-to-tip regression — TempEst, or `scripts/04_root_to_tip.R`.

Look at: R², the slope, and the scatter plot. **Always look at the plot.**

### R² — the standard table

| R² | Interpretation | Usual recommendation |
|---|---|---|
| > 0.80 | Excellent temporal signal | Strong clock; proceed confidently |
| 0.60 – 0.80 | Good | Suitable for time-scaled analyses |
| 0.40 – 0.60 | Moderate | Use caution; inspect outliers |
| 0.20 – 0.40 | Weak | Investigate before trusting |
| < 0.20 | Very weak | Clock dating usually unreliable |

### Correlation

| r | Interpretation |
|---|---|
| > 0.80 | Excellent |
| 0.60 – 0.80 | Good |
| 0.40 – 0.60 | Moderate |
| < 0.40 | Poor |

### The slope is a sanity check you can't skip

The slope **is** the evolutionary rate (substitutions/site/year).

- **Negative slope** → divergence *decreases* with time. Biologically
  impossible under a clock. It means no temporal signal, or wrong dates.
  **Do not date this tree.** No amount of downstream cleverness fixes it.
- Slope wildly outside the published range for your organism (Section 2) →
  treat the whole analysis as broken until you find out why.

---

## 1b. ⚠️ Why the R² table will mislead you

**Low R² is a symptom, not a diagnosis. You must ask what caused it.**

There are two completely different diseases with the same symptom:

**Cause A — short sampling window.** Your samples span weeks, not years, so the
x-axis has almost no spread. R² collapses even though the clock is perfectly
real and the rate is perfectly estimable.

**Cause B — deep, non-clocklike divergence.** Your sequences differ mostly
because of ancient splits, not because of the time between sampling. The dates
explain nothing. The clock is meaningless here.

**This workshop demonstrates both.** From the actual runs:

| Dataset | Window | R² | Rate estimated | tMRCA (90% region) | Verdict |
|---|---|---|---|---|---|
| **B — simulated** | 20 years | **~0.80** | ~2.11e-3 ± 0.13e-3 (**true = 2.2e-3**) | ~1994 (±~2 yr); **true = 1992.7** | Trust it |
| **C — 2009 H1N1** | **43 days** | **~0.20** | ~3.5–3.9e-3 ± ~1.4e-3 | ~2009.1 (2008.9–2009.2) | **Trust it anyway** |

*(Your exact numbers will differ slightly — see "Why your numbers wobble" below.)*

Dataset C scores 0.19 — "very weak, unreliable" by the table. And yet it
recovers a **textbook pandemic H1N1 rate (~4e-3)** and a tMRCA of **~2009.1**,
which matches the published estimate for the 2009 pandemic to within weeks.
Its R² is low *only* because 43 days of sampling gives the regression nothing
to lever against.

If you had obeyed the table, you would have thrown away a perfectly good
analysis.

**The diagnostic question is not "is R² big?" It is:**

> Is my sampling window long enough, relative to how fast this organism
> evolves, for divergence to have visibly accumulated across it?

**Quick check:** multiply `rate × window`. If that product is not much bigger
than your sequencing/phylogenetic noise, R² will be low no matter how good your
clock is.
- H1N1: 4e-3 × 0.12 yr ≈ 5e-4 subs/site ≈ **6 substitutions** across 13 kb.
  Real, but tiny. Low R² is *expected* and *not* a problem.
- Simulated: 2.2e-3 × 20 yr ≈ 4.4e-2 subs/site ≈ **66 substitutions** across
  1.5 kb. Huge. R² is high.

**And the reverse trap:** an R² that looks fine can still be junk if it was
obtained by *searching* for the root that maximises it (below).

---

## 1c. ⚠️ The best-fitting root inflates R²

TempEst's "best fit" and TreeTime's `--reroot best` **search over roots to
maximise R²**. You then judge the fit using the very dates you fitted to. The
resulting R² is optimistically biased.

Measured in this workshop on the optional real HIV dataset:

```
midpoint root : R² = 0.001
best-fit root : R² = 0.214      <- gained entirely by searching
```

Nothing about the data improved. Only the root moved.

**Rules:**
- Prefer an **outgroup root** when you have a credible outgroup.
- If best-fit R² *massively* beats midpoint R², be suspicious, not pleased.
- A best-fit R² of 0.2–0.3 is **not** evidence of a clock.
- Always report which rooting you used. "R² = 0.21" without that is meaningless.

---

## 1d. Why your numbers wobble between identical runs

Re-run the whole pipeline and your numbers will move slightly. This is **not** a
mistake, and it is worth understanding:

- **IQ-TREE** starts its search from random trees, and `-T AUTO` varies the
  thread count (which changes the order of floating-point summation).
- **TreeTime** resolves polytomies randomly and does marginal inference.

The workshop scripts pin both seeds, so variation is small — but not zero.
Observed across repeated runs here:

| | Dataset B (simulated) | Dataset C (H1N1) |
|---|---|---|
| Rate | 2.108 – 2.116e-3 (**stable**) | 3.5 – 3.9e-3 (**±5%**) |
| R² | 0.78 – 0.82 | 0.19 – 0.22 |
| tMRCA | 1993.5 – 1994.1 | 2009.05 – 2009.15 |

**Read the pattern.** Dataset B barely moves; Dataset C wobbles by 5%. That
difference is informative: C has only 47 parsimony-informative sites, so the
likelihood surface is nearly flat and many trees fit almost equally well. The
instability *is* a measurement of how little information the data contains.

**In a real study:** run the analysis several times **without** a fixed seed. If
your conclusions change between runs, that instability is your result, and you
must report it. A single seeded run can hide this completely.

---

## 2. Estimated substitution rate

From TreeTime (`molecular_clock.txt`) or BEAST. Units: **substitutions/site/year**.

Always run TreeTime with `--covariation` — tips share ancestry, so they are not
independent observations. Without it you get an over-confident point estimate
and **no honest error bar**.

### HIV-1

| Rate | Interpretation |
|---|---|
| 5×10⁻⁴ – 1×10⁻³ | Slightly slow |
| **1×10⁻³ – 3×10⁻³** | **Typical HIV-1 evolution** |
| > 3×10⁻³ | Rapid — investigate data quality |

Typical published value ≈ **2.1 × 10⁻³** (gene- and region-dependent; *env* is
faster than *gag*/*pol*).

### Influenza A (H1N1pdm09)

| Rate | Interpretation |
|---|---|
| 3×10⁻³ – 5×10⁻³ | Typical whole-genome pandemic H1N1 |
| > 6×10⁻³ | Suspicious — check alignment and dates |

*Workshop result: ~3.5–3.9 × 10⁻³ ± ~1.4 × 10⁻³ — squarely in range.*

### SARS-CoV-2

| Rate | Interpretation |
|---|---|
| 5×10⁻⁴ – 1×10⁻³ | Typical (≈8×10⁻⁴ early pandemic) |
| > 1×10⁻³ | Faster than expected |

### HBV

| Rate | Interpretation |
|---|---|
| 1×10⁻⁵ – 1×10⁻⁴ | Typical HBV evolution |
| > 5×10⁻⁴ | Usually a calibration or data problem |

> HBV is ~100× slower than HIV. A 20-year sampling window that gives HIV
> excellent signal gives HBV almost none. **The window you need depends on the
> organism's rate.**

**How to report it:**

> "TreeTime estimated a rate of 3.7 × 10⁻³ (± 1.4 × 10⁻³) substitutions/site/year,
> consistent with published estimates for pandemic H1N1."

A rate inside the published range does **not** prove your clock is good — see
Section 1b. It is a necessary check, not a sufficient one.

---

## 3. Estimated tMRCA

TreeTime reports a date for every node; the root's date is your tMRCA.

**Find the root properly.** Do **not** assume `NODE_0000000` is the root — in
our runs it usually wasn't. Take a tip's date minus its root-to-tip distance in
years, or read the oldest dated node in TreeTime's `dates.tsv`.

### A tMRCA without a confidence interval is not a result

| CI width | Interpretation |
|---|---|
| < 5 years | Good precision — usable |
| 5 – 15 years | Poor — indicative only |
| > 15 years | Useless — the interval excludes nothing |

Worked examples from this workshop:

| Dataset | tMRCA | 90% region | Width | Verdict |
|---|---|---|---|---|
| B (simulated) | ~1994 | ~1991.8 – 1996.0 | ~4 yr | Good — **and it contains the true 1992.7** |
| C (H1N1) | ~2009.1 | ~2008.9 – 2009.2 | ~0.3 yr | Good — matches published |
| A (real HIV, optional) | **1793** | **1679 – 1875** | **196 yr** | **Junk** |

That HIV row is the lesson. HIV-1 group M did not exist in humans in 1793 — the
estimate is centuries out, and its own CI is two centuries wide. The software
printed it without complaint.

**tMRCA is always extrapolated back beyond your oldest sample, so it carries
more error than the rate.** Expect it, report the interval, and sanity-check
against what is biologically possible.

---

## 4. Which sequences cluster?

A **transmission cluster** = a group of sequences that are each other's closest
relatives, very similar, and well supported.

Look for: short branches • high support • same location • similar dates.

### Bootstrap

| UFboot | Interpretation |
|---|---|
| > 95 | Very strong support |
| 80 – 95 | Strong |
| 70 – 80 | Moderate |
| 50 – 70 | Weak |
| < 50 | Ignore |

IQ-TREE gives **two** measures — read them together:
- **UFboot ≥ 95** *and* **SH-aLRT ≥ 80** → trust the clade.
- One high, one low → treat as unresolved.

> **Low support is often the data's fault, not yours.** Dataset C has 13,109 bp
> but only **47 parsimony-informative sites** (98.1% of sites are constant) —
> so only **5 of 48** nodes reach UFboot ≥ 95. You cannot bootstrap your way to
> evidence that isn't in the alignment.

### The threshold defines the answer

A cluster is normally defined by a maximum within-cluster patristic distance.
Common HIV conventions: **1.5%** (0.015) or **4.5%** (0.045) — *different papers
use different numbers*.

**"We found 3 transmission clusters" is meaningless without the threshold and
the support cutoff.** See `figures/08_cluster_threshold_*.pdf`: the same tree
yields a different cluster count at every threshold. Always report:

> "Clusters were defined as clades with ≥95% ultrafast bootstrap support and a
> maximum within-cluster patristic distance of 0.02 substitutions/site."

### Precision vs recall — you cannot have both

Workshop result on Dataset B, where the true clusters are known:

```
True clusters : 15
Found         : 12   (precision 0.90, recall 0.64)
```

At UFboot ≥ 95 we got 90% precision but missed a third of the real clusters.
Loosen to ≥70 and you find more clusters — some of them false. **That trade-off
is the whole game.** Choose it deliberately, and say which way you chose.

---

## 5. Is there evidence of transmission?

**Strong evidence:** bootstrap > 95 • very short branches • samples months
apart • same district • epidemiological link independent of the sequences.

**Weak evidence:** long branches • low support • different continents • large
sampling interval.

### What a cluster does *not* tell you

- **Not** "A infected B." The real infector may be unsampled.
- **Direction of transmission cannot be read off this kind of tree.**
- Clustering shows shared recent ancestry — that is all.

Phylogenetics can **exclude** a suspected link far more safely than it can
**confirm** one. Say "consistent with" and "cannot exclude", not "proves".

> **Ethics.** These are people. A cluster is not evidence of wrongdoing, and
> inferred linkage has been misused in criminal prosecution. Handle
> accordingly.

---

## 6. Outliers

TempEst/root-to-tip flags sequences far from the regression line (large
residual). `scripts/04_root_to_tip.R` reports a **z-score** per tip.

| |z| | Action |
|---|---|
| < 2 | Expected scatter |
| 2 – 3 | Inspect |
| > 3 | Strong outlier — investigate |

**Possible causes:** wrong sampling date • contamination • recombinant • poor
sequence quality • mislabelled sample.

Worked example — the planted date error in Dataset B:

```
SIM_ZMB_014_2022.97   residual -0.0381   z = -5.39     <- next worst is 1.78
   truth: really sampled 2007.97, recorded as 2022.97 (+15 yr error)
```

z = −5.39 with the next worst at 1.78 is unmistakable. It carries the
divergence of a 2008 sample while claiming to be from 2023.

> **The catch:** you can only detect a date outlier if the dataset *has*
> temporal signal. The same error planted in the real HIV dataset scored just
> z = −1.88 and was **not flagged**. In a dataset with no clock, outlier
> detection quietly stops working.

---

## 7. Which samples should be removed?

Remove for a **reason you can state**, not because removal improves your number.

**1. Poor sequence quality.** Many Ns or large gaps.
- \> 5% missing → inspect. > 10% → usually exclude from clock analyses.

**2. Incorrect dates.** e.g. a sequence collected in 2024 that clusters with
1995 viruses → suspect a metadata error. **Fix the metadata if you can; only
delete as a last resort.**

**3. Recombinants.** Detect with RDP, jpHMM, or SimPlot. A recombinant has two
different histories and belongs to neither tree. Usually analysed separately
rather than deleted — depends on your study question.

**4. Long branches.** May indicate contamination, sequencing artefact, or
genuine accelerated evolution. **Investigate before excluding** — a long branch
can be your most interesting sample.

**5. Duplicates.** Identical sequences from the same patient: keep one. But
identical sequences from *different* patients are a real finding, not a
QC problem — check before deleting.

### The "removing one sequence rescues R²" test

If dropping a single sequence moves R² from 0.35 → 0.82, that sequence was the
problem. Workshop example (Dataset B):

```
R² 0.806 -> 0.904 after dropping SIM_ZMB_014 (the planted date error)
```

**But apply the test honestly.** If R² barely moves when you drop your
outliers, the dataset has no clock signal, and deleting samples until the plot
looks nice is **data dredging**. In the real HIV dataset, dropping the flagged
outlier moved R² by 0.001 — the problem was never one bad sequence.

**Rule: state how many sequences you removed and why. Report the analysis both
with and without them.**

---

## 8. Final interpretation table

| Question | A good answer looks like |
|---|---|
| Temporal signal | R² > 0.6 with positive slope — **or** a defensible explanation of why R² is low (short window) backed by a plausible rate |
| Substitution rate | Consistent with published rates for the organism, **quoted with its error** |
| tMRCA | Biologically plausible, with a **narrow enough CI to mean something** |
| Clusters | UFboot ≥ 95 (or ≥80) and short branches — **with the threshold stated** |
| Transmission | Genetically close *and* epidemiologically consistent; phrased as "consistent with", never "proves" |
| Outliers | Large root-to-tip residuals, long branches, or inconsistent metadata — **investigated, not just deleted** |
| Samples removed | Poor quality, wrong dates, contamination, duplicates — each with a stated reason, and the analysis reported both ways |

---

## 9. Example student conclusion

> The simulated dataset (B) exhibited strong temporal signal (root-to-tip
> R² ≈ 0.80, best-fit root). TreeTime estimated an evolutionary rate of
> 2.11 × 10⁻³ ± 0.13 × 10⁻³ substitutions/site/year and a tMRCA of ≈1994
> (90% region ≈1991.8–1996.0). Both intervals contain the true simulated values
> (2.2 × 10⁻³ and 1992.7), so the method is behaving correctly on data where
> the answer is known.
>
> One sequence (SIM_ZMB_014) was a marked root-to-tip outlier (z = −5.4);
> removing it improved R² from 0.81 to 0.90. Its recorded date was inconsistent
> with its divergence, indicating a metadata error rather than unusual biology.
>
> The 2009 H1N1 dataset (C) returned a much lower R² (0.19), but this reflects
> its 43-day sampling window rather than an absent clock: the estimated rate
> (≈3.7 × 10⁻³ ± 1.4 × 10⁻³) is consistent with published pandemic H1N1
> estimates, and the tMRCA (≈2009.1, 90% region ≈2008.9–2009.2) matches the
> published origin of the 2009 pandemic. Twelve clusters were identified in
> Dataset B at ≥95% bootstrap and a 0.02 substitutions/site threshold
> (precision 0.90, recall 0.64 against the known truth), demonstrating the
> trade-off between cluster support and sensitivity.
>
> These interpretations integrate tree topology, branch support, temporal
> signal and organism-specific expectations, rather than relying on any single
> statistic.

---

## 10. The one-page version

1. **Plot it.** Never interpret R² without looking at the scatter.
2. **Ask why**, not just how big. Low R² from a short window ≠ low R² from
   ancient divergence.
3. **Check the slope's sign** and magnitude against published rates.
4. **Report intervals.** A tMRCA without a CI is not a result.
5. **State your thresholds** — root, cluster distance, bootstrap. They *are*
   the answer.
6. **Sanity-check against biology.** A 1793 tMRCA for HIV is a bug, not a discovery.
7. **Justify every deletion**, and report the analysis with and without.
8. Software will always give you a beautiful, confident, fully-annotated tree.
   **That is not evidence.**
