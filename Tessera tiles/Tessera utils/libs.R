suppressPackageStartupMessages({
    
    # library(TissueSegmentation) ## useful st_ functions
    
    ## This is my custom library for various helper functions for spatial pre-processing 
    # devtools::install_github('korsunskylab/spatula', dependencies = FALSE)    
    # devtools::install_github('korsunskylab/spatula')
    library(harmony)
    library(uwot)
    library(singlecellmethods)
    library(sfarrow)
    library(mclust)
    library(geojsonsf)
    library(Rcpp)
    library(sf)
    library(spatula)
    library(purrr)
    # library(dendextend)
    library(igraph)
    library(stars)
    library(furrr)
    library(circlize)
    
    library(ComplexHeatmap)
    library(ggrepel)
    library(future)
    library(scales)
    library(glue)
    library(data.table)
    library(spatstat)
    library(tidyr)
    library(dplyr)
    library(data.table)
    library(presto)
    library(ggplot2)
    library(ggthemes)
    library(patchwork)
    library(Matrix)
    library(stringr)
    library(ggridges)    
    
    # library(ggraph)
    # library(tidygraph)
    
})

fig.size <- function(h, w) {
    options(repr.plot.height = h, repr.plot.width = w)
}
