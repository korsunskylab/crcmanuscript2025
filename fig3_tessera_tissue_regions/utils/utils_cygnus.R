get_hubs = function(obj) {
    require(lwgeom)
    stopifnot('cluster' %in% colnames(obj$aggs$meta_data))
    obj$aggs$meta_data$hub_id = NA
    obj$aggs$edges$cluster_from = obj$aggs$meta_data$cluster[obj$aggs$edges$from]
    obj$aggs$edges$cluster_to = obj$aggs$meta_data$cluster[obj$aggs$edges$to]
    
    for (.cluster in unique(obj$aggs$meta_data$cluster)) {
        E = obj$aggs$edges[cluster_from == .cluster & cluster_to == .cluster]
        N = nrow(obj$aggs$meta_data)
        adj = sparseMatrix(i =E$from, j = E$to, x = 1, dims = c(N, N))
        adj = igraph::from_adjacency()$fun(adj, 'max') ## mode='max' makes it undirected 
        i = which(obj$aggs$meta_data$cluster == .cluster)
        obj$aggs$meta_data$hub_id[i] = paste0(.cluster, '_', igraph::components(adj)$membership[i])
    }
    obj$aggs$meta_data$hub_id = as.integer(factor(obj$aggs$meta_data$hub_id))

    ## INITIALIZE HUBS OBJECT 
    agg_to_hub = sparseMatrix(i=obj$aggs$meta_data$id, j=obj$aggs$meta_data$hub_id, x=1)
    hub_shapes = obj$aggs$meta_data %>% 
        split(.$hub_id) %>% 
        imap(function(.SD, .hub_id) {
            st_sf(id = .hub_id, shape = st_union(c(.SD$shape)))
        }) %>% 
        bind_rows()
    stopifnot(all(hub_shapes$id == seq_len(nrow(hub_shapes))))
    obj$hubs = list(
        meta_data = hub_shapes, 
        counts = obj$aggs$counts %*% agg_to_hub
    )
    obj$aggs$edges$hub_from = obj$aggs$meta_data$hub_id[obj$aggs$edges$from]
    obj$aggs$edges$hub_to = obj$aggs$meta_data$hub_id[obj$aggs$edges$to]
    obj$hubs$edges = unique(obj$aggs$edges[cluster_from != cluster_to][, .(from=hub_from, to=hub_to)])
    obj$hubs$edges = data.table(
        from = pmin(obj$hubs$edges$from, obj$hubs$edges$to), 
        to = pmax(obj$hubs$edges$from, obj$hubs$edges$to)
    )
    obj$hubs$meta_data$area = st_area(obj$hubs$meta_data$shape)
    obj$hubs$meta_data$perimeter = st_perimeter(st_boundary(obj$hubs$meta_data$shape))
    obj$hubs$meta_data$naggs = colSums(agg_to_hub)
    obj$hubs$meta_data$npts = as.numeric(matrix(obj$aggs$meta_data$npts, nrow = 1) %*% agg_to_hub)
    hub_degrees = table(c(obj$hubs$edges$from, obj$hubs$edges$to))
    obj$hubs$meta_data$degree = 0
    obj$hubs$meta_data$degree[as.integer(names(hub_degrees))] = as.numeric(hub_degrees)
    obj$hubs$meta_data = cbind(obj$hubs$meta_data, st_coordinates(st_centroid(obj$hubs$meta_data$shape)))
    D = Matrix::sparse.model.matrix(~0+cluster, obj$aggs$meta_data)
    D = t(D) %*% agg_to_hub
    obj$hubs$meta_data$cluster = factor(D@i+1)
    obj$hubs$edges$cluster_from = obj$hubs$meta_data$cluster[obj$hubs$edges$from]
    obj$hubs$edges$cluster_to = obj$hubs$meta_data$cluster[obj$hubs$edges$to]
    return(obj)
}
    