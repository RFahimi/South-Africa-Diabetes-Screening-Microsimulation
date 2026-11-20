# South Africa diabetes screening: individual-level microsimulation

R analysis code for a **20-year, individual-level health-economic microsimulation** comparing alternative strategies for type 2 diabetes detection in South Africa.

## Research context

This code accompanies a manuscript titled *Cost-effectiveness of fast-and-frugal tree-based and machine learning-based strategies for type 2 diabetes screening in South Africa: an individual-level microsimulation*, led by Reza Fahimi at Warwick Medical School, University of Warwick. **Manuscript status: accepted for publication in Frontiers in Clinical Diabetes and Healthcare (November 2026).** Final publication details pending.

The script implements a four-state (no diabetes; undiagnosed diabetes; diagnosed diabetes; death), annual-cycle stochastic microsimulation. It compares a no-active-screening baseline with machine-learning targeting and three fast-and-frugal-tree (FFT) strategies. It estimates per-person costs and quality-adjusted life years (QALYs), incremental cost-effectiveness ratios (ICERs), the cost-effectiveness frontier, deterministic sensitivity analyses, probabilistic sensitivity analyses (1,000 parameter draws), cost-effectiveness acceptability analyses, and expected value of perfect information (EVPI). It also implements common random numbers across strategies to reduce Monte Carlo noise.

## What is included

- `Objective5_full_analysis_v2.R`: annotated R analysis script. Its calculations are unchanged from the uploaded script; only two header comments were updated for clarity.
- `.gitignore`: rules intended to prevent accidental addition of datasets, generated outputs and local R files.
- `README.md`: this file.

**No participant-level data, prediction datasets, model outputs or figures are included.** The public code can be inspected, but it cannot be run end-to-end without separately authorised access to the required input files.

## Data requirements

Place the following three files in the working directory **only in a private, access-controlled environment**, not in this public repository:

| Required input | Purpose |
| --- | --- |
| `df_t2dm.csv` | Prepared South African DHS analytical cohort, including expected column names |
| `ml_predictions_sa_fixed.csv` | Precomputed ML predictions, joined by individual `id` |
| `fft_predictions_sa.csv` | Precomputed FFT rule flags, joined by individual `id` |

The code assigns cohort IDs after filtering for South Africa; the two prediction files must use matching row/ID definitions. Exact input schemas, survey access rights and prediction-generation procedures must be confirmed against the research team's authorised analysis environment. **Do not substitute arbitrary public datasets and present their outputs as study results.**

Data are not redistributed here because access, use and sharing are subject to the relevant original data-provider and institutional terms. Please obtain data through the appropriate authorised routes.

## Software and running the analysis

The source script documents **R >= 4.0** and uses the `tidyverse` and `scales` packages. The following is an illustration of how the author ran the analysis; it has **not been executed or independently validated in preparing this repository package**.

```r
install.packages(c("tidyverse", "scales"))
```

Run from the working directory that contains the three authorised input CSVs:

```bash
Rscript Objective5_full_analysis_v2.R
```

The script creates local output files under `output/` and produces `Figure 1_corrected` graphics. These outputs are excluded by `.gitignore` and should not be committed without separate permission and review. The manuscript's Figures 2 and 3 are produced by a separate Python file mentioned in the source comments; that Python file is **not** included in this release.

## Scientific and reproducibility notes

- Several transition, screening, cost and health-utility inputs are **calibrated or assumed rather than individually estimated from observed longitudinal transitions**. The source code labels these modelling assumptions; the resulting uncertainty distributions should not be interpreted as fully empirical parameter uncertainty.
- The ML-based targeting scores are described in the analysis comments as **in-sample predictions**. That limits the interpretation of prospective deployment performance.
- The EVPI calculation uses an illustrative eligible-population size of **35 million** with an explicit source-code note that a cited estimate is required. Do not treat this number as a verified population estimate without checking its source.
- The program contains internal checks and expected-result assertions. These do not replace external methodological validation, independent code review or access to the original study data.
- `Figure 1_corrected.pdf` may require a working Cairo graphics device in the local R environment.

## Availability, citation and reuse

This repository provides code-only methodological transparency. It does **not** claim full public reproducibility without the three authorised input files. For scholarly citation, refer to the corresponding manuscript when its bibliographic details are final. Until the paper is published, cite this repository with its URL and access date if appropriate.

No open-source licence has been assigned here. Public visibility alone should not be interpreted as permission to redistribute, modify or commercially reuse the code. Reuse permissions and authorship/IP arrangements should be clarified with the research team and institution before public release.

## Contact

Reza Fahimi  
Warwick Medical School, University of Warwick  
GitHub: https://github.com/RFahimi
