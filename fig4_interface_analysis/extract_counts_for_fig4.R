#!/usr/bin/env Rscript
# =============================================================================
# Build fig4_counts_slim.rds from counts_complete.rds
#
# counts_complete.rds is a DENSE data.table, 479 genes (rows) x 6.79M cells
# (columns), ~14.5 GB. as.matrix() on it allocates a further ~26 GB plus a
# sparse copy, which is what kills a 60 GB session. Each column is one cell, a
# numeric vector of length 479, so one pass over the columns yields both the
# per-cell total (over ALL 479 genes) and the gene rows the analysis needs.
#
#   sbatch --mem=40G --time=2:00:00 --job-name=fig4_counts \
#     --output=fig4_counts_%j.log \
#     --wrap "source \$(conda info --base)/etc/profile.d/conda.sh && \
#             conda activate r-analysis-env && \
#             /usr/bin/time -v Rscript extract_counts_for_fig4.R"
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(readxl)
})

COUNTS_RDS      <- "../Labeled MERFISH data/counts/counts_complete.rds"
GENE_PANEL_XLSX <- "../supplementary_tables/Table S2.xlsx"
OUT_RDS         <- "../Labeled MERFISH data/counts/fig4_counts_slim.rds"

## ---- genes -----------------------------------------------------------------
GENES_AXIS   <- c("CXCL9", "CXCL10", "CXCL11")                       # y axis
GENES_B      <- c("ITGAE", "CCR7", "CXCL13")                         # panel B
GENES_C      <- c("ISG15", "STAT1", "CD274", "GBP1", "TAP1", "IDO1") # panel C
GENES_C_CTRL <- c("CXCL9", "CXCL10")                                 # panel C controls
GENES_MARKER <- c("EPCAM", "KRT8", "CDH1", "COL1A1", "DCN", "PTPRC",
                  "CD3E", "CD3D", "LYZ", "CD68", "C1QA", "MS4A1", "PECAM1")

GENES_REQUIRED <- unique(c(GENES_AXIS, GENES_B, GENES_C, GENES_C_CTRL))
GENES_WANTED   <- unique(c(GENES_REQUIRED, GENES_MARKER))

t0 <- Sys.time()
message("reading ", COUNTS_RDS)
dt <- readRDS(COUNTS_RDS)
message("read in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
message("class: ", paste(class(dt), collapse = ", "))
message("dim: ", nrow(dt), " genes x ", ncol(dt), " cells")

cell_ids <- names(dt)
message("example barcodes: ", paste(head(cell_ids, 3), collapse = ", "))
if (anyDuplicated(cell_ids)) stop("duplicate column names in the counts object")

## ---- gene names come from row ORDER, via the panel -------------------------
gene_panel <- readxl::read_excel(GENE_PANEL_XLSX, col_names = TRUE)[[1]]
gene_panel <- gene_panel[!is.na(gene_panel) & gene_panel != ""]
message("gene panel entries: ", length(gene_panel))

if (length(gene_panel) != nrow(dt)) {
    stop("panel has ", length(gene_panel), " genes but the matrix has ",
         nrow(dt), " rows — row order cannot be assumed.")
}

## every REQUIRED gene must exist; markers are best-effort
missing_req <- setdiff(GENES_REQUIRED, gene_panel)
if (length(missing_req)) {
    stop("required genes absent from the 479-gene panel: ",
         paste(missing_req, collapse = ", "),
         "\nThe analysis cannot proceed without these.")
}
missing_mk <- setdiff(GENES_MARKER, gene_panel)
if (length(missing_mk)) {
    message("markers not on the panel (skipped): ",
            paste(missing_mk, collapse = ", "))
}

need     <- intersect(GENES_WANTED, gene_panel)
need_idx <- match(need, gene_panel)
message("extracting ", length(need), " gene rows")

## ---- single pass over columns ---------------------------------------------
message("one pass over ", ncol(dt), " columns")
t1 <- Sys.time()
out <- vapply(dt,
              function(col) c(sum(col), col[need_idx]),
              numeric(length(need_idx) + 1L))
message("pass done in ", round(difftime(Sys.time(), t1, units = "mins"), 1), " min")

total_all <- out[1, ]
names(total_all) <- cell_ids

counts_small <- out[-1, , drop = FALSE]
rownames(counts_small) <- need
colnames(counts_small) <- cell_ids

rm(out, dt); gc()
counts_small <- Matrix::Matrix(counts_small, sparse = TRUE)

message("median depth per cell: ", median(total_all))
message("depth range: ", min(total_all), " - ", max(total_all))
message("nonzero fraction per gene:")
print(round(Matrix::rowMeans(counts_small > 0), 3))

saveRDS(list(
    counts      = counts_small,
    total       = total_all,
    cells       = cell_ids,
    genes       = need,
    n_genes_all = length(gene_panel),
    panel_src   = normalizePath(GENE_PANEL_XLSX),
    source      = normalizePath(COUNTS_RDS),
    created     = Sys.time()
), OUT_RDS, compress = FALSE)

message("done. size on disk: ", round(file.info(OUT_RDS)$size / 1e6, 1), " MB")
message("total elapsed: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
