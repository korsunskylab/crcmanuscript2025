# Figure 3 — Tessera Tiling and Tissue Region Identification

This directory contains the workflow used to partition tissue into spatially coherent regions and identify region-level transcriptional programs.

---

## Tessera Tiling

Tissue was divided into spatial tiles using the Tessera framework. These tiles were clustered to identify spatially distinct tissue regions.

Utility functions supporting tiling and preprocessing are located in the `utils/` subdirectory.

---

## Region Characterization

Cluster-level differential expression analyses were performed to define region-specific gene expression signatures. These results form the basis of Figure 3 panels.

---

### Execution Notes

Region identification depends on harmonized cell embeddings from Figure 1.



---
# Manuscript Figures

| Figure | Description |
|---|---|
| `figure_3.qmd` | Figure 3 |
| `figure_s4.qmd` | Figure S4 |