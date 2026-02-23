# Figure 1 — Segmentation and Cell Typing Pipeline

This directory contains the complete workflow used to generate Figure 1, starting from raw MERFISH output and ending with harmonized, labeled single-cell profiles and marker gene identification.

The pipeline proceeds through segmentation, hierarchical cell typing, harmonization, and statistical modeling. Each stage is organized below.

---

## Step 1 — Baysor Segmentation (`step1_baysor_segmentation/`)

Spatial transcript counts were segmented into single cells using Baysor. The notebooks in this folder reconstruct the segmentation workflow used to define cell boundaries and generate cell-level expression matrices.

---

## Step 2 — Cell Typing (`step2_cell_typing/`)

### Step 2A — Coarse Typing (`step2a_coarse_typing/`)

Broad lineage identities were assigned independently per sample.  
`coarse_typing_example.ipynb` demonstrates the full workflow for a representative sample. The same procedure was applied to all samples.

### Step 2B — scRNA-seq Reference Preparation (`step2b_scRNA_references/`)

Single-cell RNA-seq reference datasets were processed to define lineage-specific transcriptional profiles used for label transfer and fine typing.

### Step 2C — Fine Typing (`step2c_fine_typing/`)

Within each lineage, weighted KNN-based approaches were used to refine cell-type annotations.

### Step 2D — Harmonization (`step2d_harmonize_samples/`)

Cell embeddings were harmonized across samples to enable joint analysis and cross-sample comparisons.

### Step 2E — GLMM Modeling and Marker Identification (`step2e_glmms_and_markers/`)

Generalized linear mixed models (GLMMs) were used to identify cluster-specific marker genes within samples, followed by meta-analysis across samples to define robust transcriptional programs.

---

### Execution Order

1. Segmentation  
2. Coarse typing  
3. Reference preparation  
4. Fine typing  
5. Harmonization  
6. GLMM marker analysis