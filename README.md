# sleepMenopause

Analysis code for the study **"Multimodal mapping associates age-related
androgen decline with perimenopausal sleep-onset disruption."**

The project combines population and clinical cohort analyses, causal inference
and single-cell transcriptomics to dissect perimenopausal sleep disruption and
the role of androgens in sleep onset. This repository contains the analysis
scripts only — **no raw or participant-level data are included** (see
[Data availability](#data-availability)).

---

## Repository structure

```
sleepMenopause/
├── README.md                  # this file (index)
├── LICENSE                    # academic / non-commercial use license
├── CITATION.cff               # how to cite this work
├── .gitignore                 # excludes data, outputs and system files
└── scripts/
    ├── 01_population_age_sleep_comparison.py
    ├── 02_cohort_preprocessing_and_regression.Rmd
    └── 03_hormone_cell_atlas_overview.R
```

## Scripts index

| # | Script | Language | What it does | Manuscript |
|---|--------|----------|--------------|------------|
| 01 | `scripts/01_population_age_sleep_comparison.py` | Python | Age-stratified comparison of sleep outcomes in women (<40 vs ≥40 y); Mann–Whitney U, Welch's t-tests and ordinal models — the population baseline of age-related sleep degradation. | Figure 1 |
| 02 | `scripts/02_cohort_preprocessing_and_regression.Rmd` | R (R Markdown) | Cleans the perimenopausal sleep-survey responses, derives the analysis dataset, and models associations between hormonal use, clinical triggers and sleep outcomes via ordinal logistic regression (`MASS::polr`) and Spearman correlations. | Figures 2–4 |
| 03 | `scripts/03_hormone_cell_atlas_overview.R` | R | Summarises the single-cell Hormone Cell Atlas subset (ovary, adrenal, breast): donor counts, cell-type counts and tissue composition for the steroidogenic-enzyme and androgen-receptor analyses. | Figure 6 / overview |

> The Mendelian randomisation analysis reported in the manuscript (Figure 5)
> used external GWAS summary statistics and is not included here.

## Analysis workflow

The scripts follow the manuscript narrative:

1. **Population baseline (01)** — establish age-related sleep degradation in
   women over 40.
2. **Perimenopausal cohort (02)** — deep-phenotype an independent cohort and map
   clinical triggers and hormonal use onto specific sleep domains.
3. **Mechanism (03)** — characterise age-related steroidogenic decline and
   androgen-receptor expression in the single-cell Hormone Cell Atlas.

## Data availability

Raw and participant-level data are **not** distributed with this repository for
privacy and ethical reasons. Each script documents the input files it expects
in its header. Survey data were collected under the study's ethical approval;
the single-cell analyses draw on the Hormone Cell Atlas resource (Fei et al.,
2026). Requests for derived or summary data should be directed to the
corresponding authors.

## Requirements

**R** (≥ 4.4 recommended)

```r
install.packages(c("dplyr", "tidyr", "caret", "ggplot2", "stringr",
                   "ggbiplot", "MASS", "broom", "data.table", "scales",
                   "cowplot"))
```

**Python** (≥ 3.9)

```bash
pip install pandas numpy scipy statsmodels
```

## Reproduction

1. Obtain the input data (see *Data availability*) and place them where each
   script expects them, updating any absolute file paths.
2. Run `scripts/01_population_age_sleep_comparison.py`.
3. Knit `scripts/02_cohort_preprocessing_and_regression.Rmd` (e.g. in RStudio).
4. Run `scripts/03_hormone_cell_atlas_overview.R`.

## Citation

If you use this code, please cite the study (see `CITATION.cff`).

## Authors

Suleyman Noordeen, Lijiang Fei, Mariam Bihnam, Isabel Huang-Doran, Dan Reisel,
Kari Nightingale, Louise Newson and Yizhou Yu.

Correspondence: Yizhou Yu and Louise Newson.

## License

Free for academic, research and other non-commercial use with attribution.
**Commercial use requires a separate license** from the authors — see
[`LICENSE`](LICENSE).
