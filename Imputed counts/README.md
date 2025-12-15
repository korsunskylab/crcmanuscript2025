# Imputed counts for CRC MERFISH samples:

## Description & Method

This folder contains imputed data for CRC MERFISH samples. Briefly, we inferred gene expression values for all MERFISH cells by Harmonization with the Pelka single-cell atlas and pooling counts from the 10 nearest scRNA neighbors. This is exactly the same procedure as described in the lung immunity hubs paper (https://www.nature.com/articles/s41590-024-01792-2). 

## Folder structure:

- One folder per sample
- - Subfolders for each lineage
- - - Imputed counts in Matrix Market format. 
- - - For lineages with > 100000 cells, the data is split into smaller chunks. Each chunk is written to its own Matrix Market files, labeled like so: sampleID_lineageName_counts_inferred_chunk_IDNumber.mtx, plus additional files for row and column names.
- metadata
- - One file per sample, with a table of per-cell metadata including cell state labels and ids.

## File format:

The data is presented in Matrix Market format, which can be handled by specialized libraries in both R and Python. Here is an example in R, using the Matrix library: https://stat.ethz.ch/R-manual/R-patched/RHOME/library/Matrix/html/externalFormats.html . I believe scipy in Python can also handle these files, as described here: https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.mmread.html#mmread 

**Matrices for each sample should be concatenated to recreate the imputed count matrix for the whole sample.**

```
require(Matrix)
counts = readMM('G4209_Epi_counts_inferred_chunk_14.mtx')
dim(counts) # 1605 127411
rownames(counts) = readLines('G4209_Epi_counts_inferred_chunk_14_colnames.txt')
colnames(counts) = readLines('G4209_Epi_counts_inferred_chunk_14_rownames.txt)
counts[1:10, 1:10]
```
