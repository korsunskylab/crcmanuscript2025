require(baysorrr)
require(tidyverse)
require(glue)

# ---- CONFIGURATION ----
# Set these paths before running
baysor_binpath = "path/to/baysor"          # path to Baysor binary (e.g. ~/.julia/bin/baysor)
dirs = c(
    "path/to/sample1/crescendo_segmentation/",  # directory containing cellpose_tx_annotation.csv
    "path/to/sample2/crescendo_segmentation/",
    "path/to/sample3/crescendo_segmentation/"
)
# -----------------------

for (dir in dirs){
    tx = readr::read_delim(glue::glue(dir, '/cellpose_tx_annotation.csv')) %>% rename('x' = global_x, 'y' = global_y, 'z' = global_z, 'cell' = cellID) %>% select(x, y, gene, cell,)
    head(tx)

    output_dir = dir
    n_clusters = 10
    prior_segmentation_confidence = 0.7
    max_tx_per_tile = 5e6
    scale = 5 
    no_ncv_estimation = FALSE
    min_molecules_per_cell = 10 ## doesn't seem to work? 
    remove_temp_files = TRUE ## only set to FALSE for debugging purposes 
    max_attempts = 2 ## retry Baysor up to max_attempts-1 times. max_attempts=1 means no retries.

    stopifnot(all(colnames(tx) == c('x', 'y', 'gene', 'cell')))
    stopifnot(file.exists(baysor_binpath))

    ## split transcripts into tiles
    #environment(baysorrr:::split_tx_files) <- environment()
    ntiles <- baysorrr:::split_tx_files(output_dir, max_tx_per_tile) 

    ncv_str = ''
    #if (no_ncv_estimation) ncv_str = '--no-ncv-estimation'

    ## run baysor in each tile
    cmds = purrr::map_chr(glue::glue('{output_dir}/g{1:ntiles}/'), function(outdir) {
        as.character(glue::glue('{baysor_binpath} run --x-column x --y-column y --gene-column gene --scale={scale} --save-polygons=GeoJSON --min-molecules-per-cell={min_molecules_per_cell} {ncv_str} --prior-segmentation-confidence={prior_segmentation_confidence} --n-clusters={n_clusters}  -o {outdir} {outdir}/tx_baysor.csv :cell'))
    }) 

    print(cmds)
    
    fileConn<-file(glue::glue(dir, 'baysor_cmds.txt'))
    writeLines(cmds, fileConn)
    close(fileConn)

}


 # finish and clean up 


output_dir = dir
n_clusters = 10
prior_segmentation_confidence = 0.7
max_tx_per_tile = 5e6
scale = 5 
no_ncv_estimation = FALSE
min_molecules_per_cell = 10 ## doesn't seem to work? 
remove_temp_files = FALSE ## only set to FALSE for debugging purposes 
max_attempts = 2 ## retry Baysor up to max_attempts-1 times. max_attempts=1 means no retries.

for (dir in dirs){
    output_dir = dir
    ## stitch and summarize results 
    tx <- baysorrr:::baysor.collect_tx(output_dir)
    counts <- baysorrr:::tx_to_counts(tx$gene, tx$cell, remove_bg = TRUE)
    #environment(baysor.collect_cells) <- environment()
    cells = baysorrr:::baysor.collect_cells(output_dir, no_ncv_estimation)  
    
    ## QC: remove low count cells 
    ##     why doesn't Baysor do this internally? 
    i_keep = which(cells$n_transcripts >= min_molecules_per_cell)
    tx$cell[!tx$cell %in% i_keep] = 0
    cells = cells[i_keep, ]
    counts = counts[, i_keep]    

    ## cache results
    #environment(baysor.finish) <- environment()
    baysorrr:::baysor.finish(remove_temp_files) 
}