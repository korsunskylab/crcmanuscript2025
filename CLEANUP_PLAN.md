# Repository Cleanup Plan
*Audit date: 2026-02-23 — reflects current state of `main` after 4 reorganization commits + unstaged working-tree changes*

---

## 1. Coarse Typing Notebooks

### Per-sample notebooks
The original per-sample notebooks (one per sample ID) have been removed from the repo during reorganization. The full list that existed was:

| Notebook | Cells | Notes |
|---|---|---|
| C107.ipynb | 101 | Had extra cells vs standard template |
| C110.ipynb | 95 | Standard |
| C123.ipynb | 95 | Standard |
| C163.ipynb | 95 | Standard |
| C164.ipynb | 95 | Standard |
| C167.ipynb | 95 | Standard |
| G4209.ipynb | 96 | Extra cell vs standard |
| G4423.ipynb | 96 | Extra cell; used `/n/data1/` paths (earlier run) |
| G4554.ipynb | 95 | Standard |
| G4595.ipynb | 95 | Standard |
| G4630.ipynb | 95 | Standard |
| G4659-CP-MET_Beta8.ipynb | 94 | Met sample; slightly different path scheme |
| G4659-CP-MET_VMSC04701.ipynb | 94 | Met sample; slightly different path scheme |
| G4659_Beta8.ipynb | 95 | Standard |
| G4669_reg0.ipynb | 95 | Standard |
| G4669_reg1.ipynb | 95 | Standard |
| G4671.ipynb | 96 | Extra cell vs standard |
| G4712_Beta10.ipynb | 95 | Standard |
| G4738_Beta10_06.ipynb | 95 | Standard |
| G4738_Beta10_08.ipynb | 95 | Standard |
| Jax001.ipynb | 95 | Standard |
| compare_knn_and_geneformer.ipynb | 147 | Comparison/diagnostic; not a per-sample typing notebook |
| Untitled.ipynb | 16 | Scratch; deleted |

### Representative example (current)
**`fig1_segmentation_and_cell_typing/step2_cell_typing/step2a_coarse_typing/coarse_typing_example.ipynb`**

- 82 cells, R kernel
- Has a clean configuration block (Cell 1) with `SAMPLE_ID`, `RAW_DATA_PATH`, and `REFERENCE_PATH` as placeholder variables
- All subsequent data loading uses those variables (e.g., `readr::read_rds(RAW_DATA_PATH)`)
- No hardcoded absolute paths remain
- **Status: DONE — this is the cleaned representative notebook**

The simplest and most self-contained template for reference would have been any of the 95-cell standard notebooks (e.g., `G4554.ipynb`); `G4554` was structurally identical to the majority and used the consistent `/n/scratch3/` path scheme. The current `coarse_typing_example.ipynb` was derived from `G4209` (one of the 96-cell notebooks) and has been further simplified.

---

## 2. Duplicate Files

### Previously duplicated — both now resolved

**`collect_harmonized_embeddings_20241209.ipynb`**
- Old copy 1: `Cell typing/Step 2 - Fine_typing/Fine typing pipeline scripts/collect_harmonized_embeddings_20241209.ipynb`
- Old copy 2: `Cell typing/Step 2 - Fine_typing/Fine_typing_with_weighted_KNN/collect_harmonized_embeddings_20241209.ipynb`
- **Current state:** Single copy at `fig1_segmentation_and_cell_typing/step2_cell_typing/step2e_glmms_and_markers/collect_harmonized_embeddings_20241209.ipynb` — still has active hardcoded paths (see Section 6).

**`process_for_niches.ipynb`**
- Old copy 1: `Cell typing/Step 2 - Fine_typing/Fine typing pipeline scripts/process_for_niches.ipynb`
- Old copy 2: `Cell typing/Step 2 - Fine_typing/Fine_typing_with_weighted_KNN/process_for_niches.ipynb`
- **Current state:** Deleted (unstaged deletion in working tree, not yet committed). **Action needed: stage and commit this deletion.**

---

## 3. Files to Delete

### Already deleted (committed)
- `Imputed counts/` (entire folder) — removed in reorganization commits
- `Interface analysis/supplementary_figure_4_DEFUNCT.ipynb` — removed
- `Cell typing/Supplementary Figure 1 - not used/` (entire folder) — removed
- `Cell typing/Step 1 - Coarse_typing/Untitled.ipynb` — removed
- `Cell typing/Step 2 - Fine_typing/Fine_typing_with_weighted_KNN/Strom/Untitled.ipynb` — removed
- `Cell typing/Step 3 - Harmonize all MERFISH samples/harmonize_merfish-Copy1.ipynb` — removed
- Duplicate `Fine_typing_with_weighted_KNN/collect_harmonized_embeddings_20241209.ipynb` — removed
- Per-sample coarse typing notebooks (all 21) — removed

### Deleted in working tree but not yet committed
These show as unstaged deletions in `git status` — need to be staged and committed:

| File | Rationale |
|---|---|
| `fig1_.../step2a_coarse_typing/compare_knn_and_geneformer.ipynb` | Diagnostic comparison notebook, not part of pipeline |
| `fig1_.../step2b_scRNA_references/Bcells_process_reference_for_fine_typing.ipynb` | Cell-type-specific reference prep; redundant with representative |
| `fig1_.../step2b_scRNA_references/Mast_process_reference_for_fine_typing.ipynb` | Same |
| `fig1_.../step2b_scRNA_references/Myeloid_process_reference_for_fine_typing.ipynb` | Same |
| `fig1_.../step2b_scRNA_references/Myeloid_subcluster_reference_clusters.ipynb` | Same |
| `fig1_.../step2b_scRNA_references/Plasma_process_reference_for_fine_typing.ipynb` | Same |
| `fig1_.../step2b_scRNA_references/Stromal_process_reference_for_fine_typing.ipynb` | Same |
| `fig1_.../step2c_fine_typing/B_fine_typing_results.ipynb` | Cell-type-specific fine typing; 4 of 5 removed, only TNKILC remains |
| `fig1_.../step2c_fine_typing/Epi_fine_typing_results.ipynb` | Same |
| `fig1_.../step2c_fine_typing/Epi_subcluster_reference_clusters.ipynb` | Same |
| `fig1_.../step2c_fine_typing/Myeloid_fine_typing_results.ipynb` | Same |
| `fig1_.../step2c_fine_typing/Plasma_fine_typing_results.ipynb` | Same |
| `fig1_.../step2c_fine_typing/fine_typing_Stromal_cells_in_G4209.ipynb` | Sample-specific sub-analysis |
| `fig1_.../step2c_fine_typing/stromal_fine_typing_results.ipynb` | Cell-type-specific fine typing |
| `fig1_.../step2c_fine_typing/fine_typing_B.sh` | HPC submission script |
| `fig1_.../step2c_fine_typing/fine_typing_Epi.sh` | Same |
| `fig1_.../step2c_fine_typing/fine_typing_Myeloid.sh` | Same |
| `fig1_.../step2c_fine_typing/fine_typing_Plasma.sh` | Same |
| `fig1_.../step2c_fine_typing/fine_typing_strom.sh` | Same |
| `fig1_.../step2c_fine_typing/fine_typing_TNKILC.sh` | Same — kept initially, deleted during cleanup |
| `fig1_.../step2e_glmms_and_markers/glmm_coarse_20250703.ipynb` | Superseded by `glmm_coarse.ipynb` |
| `fig1_.../step2e_glmms_and_markers/process_for_niches.ipynb` | Deleted (was duplicate) |
| `fig1_segmentation_and_cell_typing/step1_baysor_segmentation/Plotting convenience - Micron coordinates of FOVs/README.md` | Utility note folder removed |

### Still present — should be deleted

| File | Rationale |
|---|---|
| `./.DS_Store` | macOS metadata |
| `./fig1_segmentation_and_cell_typing/.DS_Store` | macOS metadata |
| `./fig1_segmentation_and_cell_typing/step2_cell_typing/.DS_Store` | macOS metadata |
| `./fig1_segmentation_and_cell_typing/step2_cell_typing/step2a_coarse_typing/.DS_Store` | macOS metadata |
| `./fig1_segmentation_and_cell_typing/step2_cell_typing/step2d_harmonize_samples/.DS_Store` | macOS metadata |

**Also confirm `.DS_Store` is in `.gitignore` to prevent recurrence.**

---

## 4. Data Files

**No data files are present in the repository.** The search for `.h5ad`, `.h5`, `.csv`, `.parquet`, `.rds`, `.loom`, `.pkl`, `.npy`, `.npz`, `.zarr`, `.mtx`, and `.tsv` returned no results (excluding `.git` and `.claude`). The repo contains only code.

---

## 5. TODO / Collaboration Comments in Notebooks

Two matches found. Neither is a blocking TODO — both are cosmetic.

### `fig1_segmentation_and_cell_typing/step2_cell_typing/step2b_scRNA_references/TNKILCs_process_reference_for_fine_typing.ipynb`
- **Cell 88** (code): Contains the string "Pelka" as a data column label:
  ```r
  cbind(table(scRNA_TNKILC@meta.data$cleaned_fine_types) %>% as.matrix(),
        table(scRNA_TNKILC@meta.data$merged_fine_types) %>% as.matrix()) %>%
    as.data.frame() %>%
    rename("Cleaned Fine Types" = V1, "Pelka Labels" = V2)
  ```
  This is a meaningful column name referring to the source annotation (not a TODO). Consider renaming `"Pelka Labels"` to something generic like `"Reference Labels"` for public release.

### `fig2_pathology_annotations/figure2.ipynb`
- **Cell 6** (markdown): Section header reads:
  ```
  ## load pathology regions as transferred by pelka lab
  ```
  Internal lab provenance note. Consider revising to `## load pathology regions` for public release.

---

## 6. Hardcoded Absolute Paths

### Files with ACTIVE (non-commented) hardcoded paths — require fixing

**`fig1_segmentation_and_cell_typing/step1_baysor_segmentation/Step 1 Baysor segmentation results/baysor_segmentation_script.R`**
| Line | Path |
|---|---|
| 7 | `/n/data1/bwh/medicine/korsunsky/lab/data/BroadCancer/Nghia_Crescendo_cell_segmentations_C110/...` (3 paths in `dirs` vector) |
| 13 | `/home/mup728/.julia/bin/baysor` (inside function body) |
| 51 | `/home/mup728/.julia/bin/baysor` (top-level assignment) |

**`fig1_segmentation_and_cell_typing/step2_cell_typing/step2c_fine_typing/fine_typing_TNKILC.sh`**
| Line | Path |
|---|---|
| 12 | `source /home/mup728/jupytervenv/bin/activate` |

**`fig1_segmentation_and_cell_typing/step2_cell_typing/step2c_fine_typing/TNKILC_fine_typing_results.ipynb`**
| Cell | Path |
|---|---|
| 3 | `list.files('/n/scratch/users/m/mup728//mup728/Pelka_Baysor_segmentation/coarse_typing_with_weighted_knn/...')` |
| 28 | `readr::read_rds(glue('/n/scratch/users/m/mup728//mup728/Pelka_Baysor_segmentation/coarse_typing_with_weighted_knn/...'))` |

**`fig1_segmentation_and_cell_typing/step2_cell_typing/step2d_harmonize_samples/harmonize_merfish.ipynb`**
| Cell | Path |
|---|---|
| 2 | `/n/data1/bwh/medicine/korsunsky/lab/mup728/CRC_MERFISH/figure2/pathology_regions_postQC.csv` |
| 4 | `/n/data1/bwh/medicine/korsunsky/lab/mup728/CRC_MERFISH_niches/labeled_seurat_objects/renamed_cell_states/` |
| 9 | `/n/data1/bwh/medicine/korsunsky/lab/mup728/CRC_MERFISH_niches/labeled_seurat_objects/renamed_cell_states/renamed_merfish_cell_types.csv` |

**`fig1_segmentation_and_cell_typing/step2_cell_typing/step2e_glmms_and_markers/collect_harmonized_embeddings_20241209.ipynb`**
| Cell | Path |
|---|---|
| 6 | `/n/data1/bwh/medicine/korsunsky/lab/mup728/CRC_MERFISH/Pelka_Baysor_segmentation/coarse_typing_with_weighted_knn/Coarse_typing_with_weighted_knn/MSI/` |
| 9 | `/n/data1/bwh/medicine/korsunsky/lab/mup728/CRC_MERFISH/Pelka_Baysor_segmentation/data_and_ingest/` |

**`fig1_segmentation_and_cell_typing/step2_cell_typing/step2e_glmms_and_markers/meta_analysis_of_fine_type_markers.ipynb`**
| Cell | Path |
|---|---|
| 67 | `/n/data1/bwh/medicine/korsunsky/lab/mup728/CRC_MERFISH_niches/labeled_seurat_objects/annotated_merged_merfish.rds` |

**`fig4_interface_analysis/Gene_patterning_around_MSI_interfaces.ipynb`**
| Cell | Path |
|---|---|
| 5 | `/n/data1/bwh/medicine/korsunsky/lab/mup728/CRC_MERFISH/tissue_regions_with_tessera_20241105/agg_metadata_2025-07-22.rds` |
| 5 | `/n/data1/bwh/medicine/korsunsky/lab/mup728/CRC_MERFISH/tissue_regions_with_tessera_20241105/tile_metadata_2025-07-22.rds` |
| 12 | `/n/data1/bwh/medicine/korsunsky/lab/mup728/CRC_MERFISH/CRC_Figure_1/harmonized_merfish_20241105.rds` |

### Files with only COMMENTED-OUT hardcoded paths — cosmetically unclean

These contain old paths in commented lines; they don't affect execution but are untidy.

| File | Cells with commented paths |
|---|---|
| `fig1_.../step2b_scRNA_references/ingest_pelka_crc_V2.ipynb` | cell 3: `#setwd('/n/data1/...')` |
| `fig1_.../step2e_glmms_and_markers/meta_analysis_of_fine_type_markers.ipynb` | cells 28, 39, 49, 56: commented `list.files('/n/data1/...')` |
| `fig2_pathology_annotations/figure2.ipynb` | cells 168, 195, 221: commented `# metadata = readr::read_rds('/n/data1/...')` |

---

## 7. README Files

| File | Status | Notes |
|---|---|---|
| `README.md` (root) | Present — sparse | 5 lines: title + one sentence describing intent. No figure index, no citation/DOI, no data access note. Intentionally omits software requirements (repo documents code, not a fully reproducible pipeline). **Needs expansion before public release.** |
| `fig1_segmentation_and_cell_typing/README.md` | ✅ Complete | Covers all 5 pipeline steps, execution order, and per-step descriptions. |
| `fig2_pathology_annotations/README.md` | Empty (1 line) | Intentional — collaborators will fill in. No action needed now. |
| `fig3_tessera_tissue_regions/README.md` | ✅ Complete | Describes tiling, region characterization, and dependency on Fig 1. |
| `fig4_interface_analysis/README.md` | ✅ Complete | Lists all notebooks (including supplementary), execution dependencies on Fig 1 + 3. |
| `fig5_colocalization/README.md` | ✅ Complete | Lists both notebooks, dependencies on Fig 3 + 4. |

**Action needed:** Expand root `README.md` with at minimum: paper title/citation placeholder, brief description of repository scope, and a figure-by-figure directory index. Software/environment requirements to be deferred.

---

## Summary of Remaining Actions

| Priority | Action | Status |
|---|---|---|
| **Must do** | Delete 5 `.DS_Store` files; `.DS_Store` already in `.gitignore` | ✅ Done |
| **Must do** | Stage and commit all unstaged deletions (working tree ahead of index) | ⏳ Pending commit |
| **Must do** | Fix active hardcoded paths in 7 files (config block + placeholder variables) | ✅ Done |
| **Should do** | Delete commented-out absolute paths from 3 files | ✅ Done |
| **Should do** | Rename `"Pelka Labels"` → `"Reference Labels"` in `TNKILCs_process_reference_for_fine_typing.ipynb` | ✅ Done |
| **Should do** | Revise `## load pathology regions as transferred by pelka lab` in `figure2.ipynb` | ✅ Done |
| **Should do** | Delete `fine_typing_TNKILC.sh` (kept initially; removed for consistency with other .sh deletions) | ✅ Done |
| **Deferred** | `fig2/README.md` — collaborators will fill in | — |
| **Deferred** | Root `README.md` software requirements — not needed given repo scope | — |
| **Pending** | Expand root `README.md` with citation placeholder and figure index | ⏳ Not yet done |
