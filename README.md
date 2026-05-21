# Reproducibility Analysis and Audit of the vdLW Implementation

This repository contains the executable notebooks, figures, and supporting material accompanying the rebuttal paper:

> *Rebuttal to “Comment on ‘A theoretical upper limit for offshore wind energy extraction’ by Simão Ferreira et al. (2026)”: reproducibility analysis, audit of the vdLW implementation, and clarification of finite wind-farm correction methodologies* 

by 

Jens Nørkær Sørensen, Gunner Chr. Larsen, and Carlos Simão Ferreira

first released on 2026 May 16th 

at https://wes.copernicus.org/preprints/wes-2026-59/  CC1: 'Comment on wes-2026-59', Jens Nørkær Sørensen, 16 May 2026

The purpose of this repository is to provide a transparent and executable reproducibility framework for the analyses presented in the rebuttal paper. The notebooks reproduce the principal figures and comparisons discussed in the manuscript and document the progressive correction of implementation inconsistencies identified during the audit of the vdLW methodology.

The pdf of the rebuttal document is also in the repository.

---

# Repository Contents

## Main Notebook

- `Generation_figures_1_2_3_5_and_6.ipynb`

This notebook reproduces the analyses associated with:

- Figures 1–3:
  Type 1 analytical implementation errors and reproducibility analysis

- Figures 5–6:
  Type 2 geometry-processing and finite wind-farm methodology analyses

The notebook combines two originally independent reproducibility workflows into a single executable document for archival and publication purposes.

---

# Structure of the Notebook

The notebook is organized into two main parts corresponding to the structure of the rebuttal paper.

## Part I — Type 1 Errors

This section reproduces the analysis presented in Section 3.1 of the rebuttal paper and investigates the impact of:

- removal of the hub-height wind-speed correction,
- incorrect cut-in and cut-out assumptions,
- use of fixed latitude values,
- and incorrect implementation of the analytical capacity-factor equation.

The notebook progressively corrects these implementation inconsistencies and demonstrates restoration of agreement with the published SLS results.

This section corresponds primarily to Figures 1–3 of the rebuttal paper.

---

## Part II — Type 2 Errors

This section reproduces the analyses associated with geometry-processing inconsistencies and finite wind-farm preprocessing effects discussed in Section 3.2 of the rebuttal paper.

The analyses include examples associated with:

- neighboring wind-farm handling,
- artificial connection of disconnected wind-farm regions,
- suppression of clean inflow,
- edge-detection inconsistencies,
- and finite wind-farm correction methodologies.

This section corresponds primarily to Figures 5–6 of the rebuttal paper.

---

# Note on Structure and Repetition

The two notebook sections were originally developed as separate standalone reproducibility workflows. They were later merged into a single notebook to simplify publication and public access through NBViewer.

As a consequence, some repetition of:

- imports,
- helper functions,
- explanatory text,
- variable definitions,
- and plotting configuration

has intentionally been preserved in order to maintain the integrity and readability of the original standalone workflows.

---

# Reproducibility

The notebook is intended to function as an executable scientific audit and reproducibility document.

For best reproducibility:

1. Restart the notebook kernel
2. Run all cells sequentially
3. Verify that all figures reproduce correctly

---

# Rendered Version

A rendered version of the notebook can be viewed through NBViewer:

[(NBViewer link here)](https://nbviewer.org/github/csimaoferreira/offshore-wind-limit-reproducibility-audit/blob/main/Generation_figures_1_2_3_5_and_6.ipynb)

---

# Related Publication

The repository accompanies the rebuttal manuscript submitted to *Wind Energy Science Discussions*.

---

# License

(Add license information here)
