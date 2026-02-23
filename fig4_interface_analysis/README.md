# Figure 4 — Interface Construction and Comparative Analysis

This directory contains the construction and analysis of tumor–immune interfaces used in Figure 4 and related supplementary panels.

---

## Step 0 — Interface Construction

`step0_construct_interfaces.ipynb`

Interfaces were defined from spatial cell coordinates and lineage annotations. This step establishes the spatial framework for all downstream comparisons.

---

## Step 1 — MSI vs MSS Interface Comparison

`step1_MSI_vs_MSS_interfaces.ipynb`

This notebook compares interface properties between MSI and MSS samples, corresponding to Figure 4 panels.

---

## Step 2 — Hub+ vs Hub− Interface Comparison

`step2_hubPos_vs_hubNeg.ipynb`

This notebook analyzes differences between hub-positive and hub-negative interfaces.

---

## Supplementary and Extended Analyses

- `supplementary_figure_4.ipynb`
- `Annotated_Step2_analysis_compare_hubPos_vs_hubNeg.ipynb`
- `Gene_patterning_around_MSI_interfaces.ipynb`
- `Spatial_gene_expression_hub_vs_nonhub_vs_MSS.ipynb`

These notebooks generate supplementary panels and extended gene-patterning analyses related to Figure 4.

---

### Execution Notes

Interface analyses require:
- Harmonized cell embeddings (Figure 1)
- Region definitions (Figure 3)