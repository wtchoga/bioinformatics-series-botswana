# Practical 2 — Report Template

**Name:** ________________  **Date:** ________________

2–3 pages. Marks are for **reasoning**, not for producing numbers. A correct
number with a mechanical interpretation scores less than a well-argued reading
of an awkward result.

Keep `notes/INTERPRETATION_GUIDE.md` open while you write.

---

## 1. Data and methods (½ page)

State, for each dataset:
- what it is, how many sequences, alignment length, **sampling window**
- substitution model chosen by ModelFinder, and what its terms mean
- how you rooted the tree, and why
- software and versions

> Someone should be able to reproduce you from this section alone.

---

## 2. Results table

| | Dataset B (simulated) | Dataset C (H1N1) |
|---|---|---|
| Sequences / length | | |
| **Sampling window** | | |
| Variable sites | | |
| Model (BIC) | | |
| Nodes UFboot ≥ 95 | | |
| Root-to-tip R² (midpoint) | | |
| Root-to-tip R² (best-fit) | | |
| Slope (subs/site/yr) | | |
| TreeTime rate ± sd | | |
| **tMRCA** | | |
| **tMRCA 90% region** | | |
| Outliers (\|z\| > 2) | | |
| Clusters found (threshold?) | | |

---

## 3. Temporal signal (½–1 page) — **the important section**

For **each** dataset:

1. Is there temporal signal? Quote R², correlation, slope.
2. **If R² is low — why?** Short sampling window, or non-clocklike divergence?
   These are different diagnoses with different prognoses. Justify your answer
   using the sampling window and the rate, not just the R² value.
3. How much of your R² came from *searching* for the best-fit root? Compare the
   midpoint and best-fit values. What does that comparison tell you?
4. Would you date this tree? Defend it.

> Dataset C scores "very weak" on the standard R² table. Decide whether the
> table is right, and argue your case. This is the question the practical is
> actually about.

---

## 4. Molecular clock and tMRCA (½ page)

1. Your rate estimate ± error, per dataset. Consistent with published values?
2. Your tMRCA **with its interval**. Is the interval narrow enough to mean
   anything?
3. **Dataset B: open `data/simulated/SIM_truth.txt`.** True rate 2.2 × 10⁻³,
   true tMRCA 1992.70. **Did your intervals actually contain the truth?** If
   they did, say so plainly. If not, work out why.
4. **Dataset C:** compare your tMRCA to the published origin of the 2009
   pandemic. How close are you, and is that surprising given your R²?

---

## 5. Clusters and transmission (½ page)

1. How many clusters, at what distance threshold and what support cutoff?
2. What happens to that count if you halve the threshold? (See
   `figures/08_cluster_threshold_*.pdf`.) Is "N clusters" a fact about the
   virus or about you?
3. Dataset B: your precision and recall against the known truth. Which real
   clusters did you miss, and why?
4. What would you actually claim epidemiologically — and what would you refuse
   to claim?

> Phrase conclusions as "consistent with" / "cannot exclude". A cluster is not
> proof that A infected B: the true infector may be unsampled, and direction of
> transmission cannot be read off this tree.

---

## 6. Outliers and quality control (½ page)

1. Which sequences are outliers? Give tip names and z-scores.
2. What is the likely cause of each — wrong date, poor quality, recombinant,
   contamination, genuine divergence?
3. What happened to R² when you removed them?
4. **Which would you remove, and why?** Report your conclusions **with and
   without** them.

> "R² improved" is not a reason to delete a sample. Deleting sequences until the
> plot looks good is data dredging. State a cause, or keep the sample.

---

## 7. Conclusion (¼ page)

Pull it together: relationships, temporal signal, clock, and what you would and
would not report. Integrate topology, support, temporal signal and biological
plausibility — don't lean on one statistic.

---

## Checklist before you submit

- [ ] Every R² is accompanied by the plot and the rooting method
- [ ] Every tMRCA is accompanied by an interval
- [ ] Every cluster count is accompanied by its threshold and support cutoff
- [ ] Every removed sample has a stated reason (not "it improved R²")
- [ ] Dataset B results are explicitly compared against the known truth
- [ ] Dataset B is never described as real epidemiological data
- [ ] No claim that the tree "proves" a transmission event
- [ ] You explained **why** Dataset C's R² is low, not just that it is

---

## Deliverables to attach

1. ML tree (`.treefile`), both datasets
2. Rooted tree PDF (FigTree or ggtree)
3. Root-to-tip regression plot (+ TempEst screenshot if used)
4. TreeTime time tree (`.nexus`)
5. Publication-quality tree (`Figure1_*.pdf`)
6. Metadata-annotated phylogeny
