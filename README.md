# Bioinformatics Series in Botswana

Workshop website for **Bioinformatics, AI and Data Science Training**, hosted in
Gaborone by the Botswana Harvard Partnership (BHP) and the Harvard T.H. Chan
School of Public Health (HSPH).

## Instructions on contributing to the website

1. clone the repo locally
2. make changes - you can preview changes locally by clicking "Render Website" in the Build tab in RStudio, or by running `quarto preview` in the terminal.
3. commit changes and push to the `bhl-customisation` branch
4. run `quarto publish gh-pages` to render the website and update the github.io content

### Requirements

- [Quarto](https://quarto.org/docs/get-started/) 1.4 or later
- R 4.3.1 or later

```r
install.packages(c("tidyverse", "here", "knitr", "rmarkdown", "quarto"))
```

Session materials need more — see [installation-instructions.qmd](installation-instructions.qmd)
for the full list, including the Bioconductor packages (`DESeq2`, `edgeR`,
`limma`, `mixOmics`, `ggtree`) and the standalone tools (AliView, IQ-TREE 3,
FigTree 1.4.4, BEAST 1.10.5).

## Repository layout

```
├── _quarto.yml                       site config, navbar and sidebar
├── index.qmd                         landing page: overview, schedule, instructors
├── about.qmd                         teaching team, partners, contact
├── past_workshops.qmd                Sept 2025 Virus Evolution & Genomics Workshop
├── installation-instructions.qmd
├── exercises.qmd                     Day 1 + Day 2 exercises
└── materials/
    ├── 1-welcome/                    welcome session
    ├── 2-ml-session/                 Random Forest — HIV/SIV rebound
    ├── 3-bulk-rnaseq/                DESeq2 + GSEA
    ├── 1-intro-bioinfo-ai/           what is bioinformatics / AI / data science
    ├── 2-ai-application/             ML foundations → LLMs (5 modules)
    ├── 3-bioinformatics-application/ QC → phylogenetics (4 modules)
    ├── 4-hands-on/                   mtDNA case study (synthetic data)
    ├── 5-survey/                     feedback
    ├── 6-mutations/                  HIV-1 mutations & minority variants
    └── 6-phylogenetics-practical/    alignment → time-scaled tree
```

## A note on participant data

Participant-level data is **not** in this repository and must not be added. The
mtDNA case study now generates a **synthetic** dataset inline — same structure
and same sex/HIV confounding as the real Botswana cohort — so it renders anywhere
with no private data. The real `datasets/mtdna_haplogroups.csv` is not
distributed and is kept out of the repo (the gitignore entry remains as a
guard). It contains HIV status for a real cohort; obtain it from the workshop
organizers if you need to work with the actual data.

The Welcome slides are served view-only via Google Slides rather than hosted
here for download.

## Acknowledgements

The September 2025 Virus Evolution and Genomics Workshop was supported by
SANTHE and the Botswana Harvard Partnership.
