renumber_clusters <- function(clusters) {
    factor(as.integer(factor(clusters, names(sort(table(clusters), TRUE)))))        
}


## This code provides Leiden with random restarts 
## Leiden iterations work by feeding in the last solution into the next one
do_leiden_one <- function(adj, resolution, n_starts, n_iterations = 2, verbose = FALSE) {
    graph = igraph::from_adjacency()$fun(adj, mode = 'undirected')     
    graph = igraph::from_adjacency()$fun(adj, mode = 'undirected', weighted = TRUE)     
    res_best = NULL
    res_new = NULL
    if (n_starts == 1) {
        plan(sequential)
    } else {
        plan(multicore)
    }
    for (iter in seq_len(n_iterations)) {
        if (verbose) message(iter)
        solns = future_map(seq_len(n_starts), function(i) {
            set.seed(i)
            soln = igraph::cluster_leiden(
                graph = graph, 
                objective_function = 'modularity', 
                n_iterations = 1, ## default = 2
                resolution_parameter = resolution,
                initial_membership = res_best
            )
            return(list(membership = soln$membership, quality = soln$quality))
            # factor(res$membership, sort(unique(res$membership))) 
        }, .options = furrr_options(seed = 2L)) 
        scores = map_dbl(solns, 'quality')
        if (verbose) message(glue('best mod={max(scores)}'))
        res_new = solns[[which.max(scores)]]$membership
        if (!is.null(res_best) & all(res_new == res_best)) {
            res_best = res_new
            break 
        }
        res_best = res_new
    }
    
    ## order result by cluster size and return factor 
    # res_best <- factor(as.integer(factor(res_best, names(sort(table(res_best), TRUE)))))
    res_best <- renumber_clusters(res_best)
    return(res_best) 
}


delete_clusters <- function(clusters, graph, clusters_remove=NULL, clusters_block=c(''), min_size=NULL, renumber=FALSE) {
    if (is.null(clusters_remove)) {
        if (is.null(min_size)) {
            stop('Must specify either clusters to remove or minimum cluster size')
        }
        clusters_remove = names(which(table(clusters) < min_size))
    }
    
    i_reassign = which(clusters %in% clusters_remove)
    clusters_block = union(clusters_block, clusters_remove)
    for (i in i_reassign) {
        idx = seq(graph@p[i] + 1, graph@p[i+1]) 
        i_nn = graph@i[idx] + 1
        p = graph@x[idx]
        p[which(clusters[i_nn] %in% clusters_block)] = 0
        if (max(p) == 0) {
        # if (max(graph[i, i_keep]) == 0) {
            ## Case 0: set disconnected points to NA 
            clusters[i] = NA
        } else {
            ## Case I: point is connected, reassign to most connected neighbor 
            # clusters[i] = clusters[i_keep[which.max(graph[i, i_keep])]]
            clusters[i] = split(p, droplevels(clusters[i_nn])) %>% map_dbl(sum) %>% which.max() %>% names()            
        }
    }
    clusters <- droplevels(clusters)
    if (renumber) clusters <- renumber_clusters(clusters)
    return(clusters)
}


merge_clusters <- function(clusters, clusters_merge, renumber=FALSE) {
    clusters <- as.character(clusters)
    clusters[which(clusters %in% clusters_merge)] <- clusters_merge[1] 
    if (renumber) {
        clusters <- factor(as.integer(factor(clusters, names(sort(table(clusters), TRUE)))))        
        # clusters <- factor(as.integer(factor(clusters, names(sort(table(clusters), TRUE)))))        
    } else {
        clusters <- factor(clusters)
    }
    return(clusters)
}


split_clusters <- function(clusters, cluster_split, graph, resolution, n_starts, n_iterations = 2, verbose = FALSE, min_size = 0) {
    if (length(cluster_split) > 1) stop('can only split one cluster at a time') 
    clusters <- as.character(clusters)
    i_split = which(clusters %in% cluster_split)
    clusters_new = do_leiden_one(graph[i_split, i_split], resolution, 1, 1, FALSE)
    clusters[i_split] <- paste0(cluster_split, '_', clusters_new)
    clusters <- factor(clusters)
    
    ## order clusters better
    cluster_levels = levels(clusters)
    o = order(as.integer(gsub('^(.*)_.*$', '\\1', cluster_levels)))
    cluster_levels = cluster_levels[o]
    clusters = factor(clusters, cluster_levels)

    ## merge back small clusters 
    if (min_size > 0) 
        clusters <- delete_clusters(clusters, graph, min_size = min_size, renumber = FALSE)
    
    return(clusters)
}


get_markers <- function(counts, meta_data, cluster_colname) {
    meta_data$cluster <- meta_data[[cluster_colname]]
    pb <- presto::collapse_counts(
        counts, meta_data, 
        c('library', 'cluster'), 
        min_cells_per_group = 10, 
        get_norm = TRUE
    )  
    # pb$meta_data$cluster <- paste0('C', pb$meta_data$cluster)

    system.time({
        suppressWarnings({
            presto_res <- presto.presto(
                y ~ 1 + (1|cluster) + (1|cluster:library) + (1|library) + offset(logUMI), 
                pb$meta_data, 
                pb$counts_mat,
                size_varname = 'logUMI', 
                effects_cov = c('cluster'),
                ncore = 20, 
                min_sigma = .05, 
                family = 'poisson',
                nsim = 1000
            ) 
        })
    })

    contrasts_mat <- make_contrast.presto(presto_res, 'cluster')
    effects_marginal <- contrasts.presto(presto_res, contrasts_mat, one_tailed = TRUE) %>% 
        dplyr::mutate(cluster = contrast) %>% 
        dplyr::mutate(
            ## convert stats to log2 for interpretability 
            logFC = sign(beta) * log2(exp(abs(beta))),
            SD = log2(exp(sigma)),
            zscore = logFC / SD
        ) %>% 
        dplyr::select(cluster, feature, logFC, SD, zscore, pvalue) %>% 
        arrange(pvalue)
    effects_marginal$fdr <- p.adjust(effects_marginal$pvalue, method = 'BH')
    return(effects_marginal)
}


tt <- function(effects_marginal, n=10, fdr_max = .2, logFC_min = .2) {
    cluster_max = max(as.integer(gsub('.*?(\\d+)', '\\1', unique(effects_marginal$cluster))))    
    effects_marginal %>% 
        subset(fdr < fdr_max & logFC > logFC_min) %>% 
        # dplyr::mutate(cluster = factor(cluster, paste0('C', 0:cluster_max))) %>% 
        split(.$cluster) %>% 
        map(function(.SD) {
            .SD %>% 
                # dplyr::arrange(-zscore) %>% 
                dplyr::arrange(-logFC) %>% 
                head(n) %>% 
                dplyr::select(feature, cluster) %>% 
                tibble::rowid_to_column('rank')
        }) %>% 
        bind_rows() %>% 
        tidyr::spread(cluster, feature) %>% 
        dplyr::select(-rank) %>% 
        t() %>% 
        identity()
}


plot_markers <- function(genes, clusters, cluster) {
    U2$embedding %>% 
        data.table() %>% 
        cbind(CLUSTER = clusters) %>% 
        cbind(as.matrix(t(normalizeData(counts_all, 100, 'log')[genes, ]))) %>% 
        subset(CLUSTER %in% cluster) %>%  ## CXCL cluster
        tidyr::gather(key, val, all_of(genes)) %>% 
        dplyr::arrange(val) %>% 
        ggplot(aes(V1, V2, color = val)) + 
            geom_point(shape = '.', alpha = .6) + 
            scale_color_gradient2(mid = 'lightgrey') + 
            facet_wrap(~key) + 
            NULL
}


## A: graph 
# ## X: embeddings
# prune_graph_global <- function(A, X) {
#     dists = sqrt(rowSums((X[A@i + 1, ] - X[rep(1:nrow(A), diff(A@p)), ]) ^ 2))    
#     idx_remove = which(dists > quantile(dists, .68) * 3) ## 3 SDs 
#     # idx_remove = which(dists > quantile(dists, .68) * 2) ## 2 SDs
        
#     ## prune edges
#     A@x[idx_remove] <- 0
#     A <- Matrix::drop0(A)           

#     ## remove unpaired edges
#     at_mask = t(A)
#     at_mask@x = rep(1, length(at_mask@x))
#     A <- A * at_mask
    
#     return(A)
# }


