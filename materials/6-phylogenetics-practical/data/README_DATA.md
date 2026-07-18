# Data provenance — Practical 2

Read this before you quote any number out of this folder.

---

## Dataset B — `simulated/` — SIMULATED

> ⚠️ **This data is simulated. It is not from any patient, anywhere.
> Never report these numbers as real epidemiological results.**

| File | What it is |
|---|---|
| `SIM_alignment.fasta` | 60 sequences × 1,500 bp, simulated under GTR |
| `SIM_dates.tsv` | sampling dates, decimal (**contains one planted error**) |
| `SIM_metadata.csv` | country, date, true cluster membership |
| `SIM_true_timetree.nwk` | the true time-scaled tree |
| `SIM_truth.txt` | **ANSWER KEY** — true rate, tMRCA, clusters, planted error |

**Built by:** `scripts/00_simulate_dataset.R` (seeded — identical every run).

**How:** a birth–death epidemic is simulated over 1990–2025; extinct lineages
give tips at many different dates (heterochronous sampling). 60 tips are
subsampled across 2005–2025. Branch lengths in years are multiplied by a chosen
clock rate to get substitutions, then sequences are evolved along that tree.

**Chosen truth:** rate **2.2 × 10⁻³** subs/site/year · tMRCA **1992.70** ·
strict clock · GTR, no rate heterogeneity.

Countries (Botswana / South Africa / Zambia) are **labels for colouring only** —
no geography was simulated. Don't read anything into them.

**Why simulated data is in a workshop about real analysis:** you cannot learn to
read a molecular clock on data that has no clock, and you cannot check your
answer without knowing the truth. Dataset B is the only place here where "is my
estimate correct?" has a definite answer. It is deliberately kinder than
reality — that's what Dataset C is for.

---

## Dataset C — `h1n1/` — REAL

Influenza A(H1N1)pdm09 — the 2009 pandemic.

| File | What it is |
|---|---|
| `H1N1_alignment.fasta` | 50 whole genomes × 13,109 bp |
| `H1N1_dates.tsv` | decimal sampling dates |
| `H1N1_metadata.csv` | strain, location, region, date |

**Built by:** `scripts/00b_prepare_h1n1.py`
**Source:** `New Folder With Items 3/H1N1pdm_2009.fasta`

**Note the source file is misnamed.** Despite the `.fasta` extension it is a
**NEXUS** file with classic-Mac **CR (`\r`) line endings** — which is why `grep`
and most tools report *0 sequences*. The prep script normalises the line
endings and parses the NEXUS matrix. It verifies the parse against the header's
own `ntax=50 nchar=13109` and fails loudly on mismatch.

**Sampling window: 2009.285 – 2009.403 — about 43 days.** This matters more
than anything else about this dataset. Keep it in mind for Exercise 6.

**Dates** come from the original taxon names
(`A/Beijing/01/2009_135_2009.37` → strain, day-of-year, decimal date) and are
carried into the tip labels so TempEst/TreeTime can read them directly.

**`region`** is a coarse continent grouping derived from the strain's location
field, for colouring only.

**Caveat:** influenza is segmented, and this is an 8-segment concatenation.
Reassortment means segments can have different evolutionary histories. Across
this 6-week window the pdm09 lineage is effectively clonal, so a single tree is
defensible — but state the assumption if you write this up.

---

## Dataset A — real HIV-1 gag — OPTIONAL, not built by default

**Build with:** `python3 scripts/00a_prepare_hiv.py`
(requires `GAG_sdb_n417.fasta` in `~/Downloads`)

Produces `HIV_alignment.fasta` (80 seq × 1,503 bp), `reference.fasta` (HXB2,
K03455), `dates.tsv`, `metadata.csv`.

**Real** HIV-1 *gag* sequences from 6 African countries (Botswana, Kenya,
Rwanda, Uganda, South Africa, Zambia), sampled 2005–2025, subsampled from 417
sequences and stratified by country and year.

### What is real and what is not

| Column | Status |
|---|---|
| `sample_id` | **Re-issued** — `<ISO3>_WS###_<year>`. No original patient or study identifier is carried through. |
| `country` | **REAL** |
| `year` | **REAL**, but **year-precision only** — the source headers carry no month or day |
| `gene`, `length_bp` | **REAL** |
| `province` | ⚠️ **SIMULATED** — plausible admin-1 units, randomly assigned |
| `risk_group` | ⚠️ **SIMULATED** — randomly assigned |
| `subtype` | ⚠️ **SIMULATED** — drawn from regionally plausible frequencies, **not** determined from the sequence |

**The simulated columns exist so Exercise 13 has something to colour by. They
are fabricated. Do not analyse them as data, and do not let them leave this
workshop attached to these sequences.** If you need real subtypes, call them
with COMET/REGA/jpHMM; if you need real risk-group data, get it from the study.

**Year-precision dates** mean up to a year of error per tip — worth remembering
when you look at its temporal signal.

### Why it's optional, and why it's worth building

Dataset A is the **negative control**. Unlinked cross-sectional HIV from a
mature epidemic has almost no temporal signal: two random subtype-C patients
differ by ~0.145 subs/site, while 20 years of clock evolution adds only ~0.02.
The dates explain essentially nothing (R² ≈ 0.06).

TreeTime dates its tMRCA to **1793** (90% region 1679–1875) — centuries before
HIV-1 group M entered humans — and reports it without a single warning.

Nothing in this workshop teaches scepticism faster.

**It contains three deliberately planted teaching artefacts.** They are listed
in `notes/INSTRUCTOR_KEY.md`. Instructors: don't hand that out early.

---

## Ethics

Dataset A derives from real patient-derived viral sequences. Sample identifiers
have been re-issued so that no original patient or study identifier appears in
workshop materials. The fabricated `province` / `risk_group` / `subtype` columns
must never be presented as real attributes of these samples.

More generally: a phylogenetic cluster is **not** evidence that one person
infected another, and inferred linkage has been misused in criminal
prosecutions. Teach that alongside the methods.
