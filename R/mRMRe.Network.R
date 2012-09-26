## Definition

setClass("mRMRe.Network", representation(topologies = "list", mi_matrix = "matrix", causality_matrix = "list",
                feature_names = "character", target_indices = "integer"))

## Wrappers

## FIXME: Add wrappers for network

## initialize

setMethod("initialize", signature("mRMRe.Network"), function(.Object, data, prior_weight, target_indices, levels,
                layers, ..., mi_threshold = -Inf, causality_threshold = Inf)
{
    if (missing(layers))
        layers <- 1L
    
    .Object@mi_matrix <- matrix(nrow = featureCount(data), ncol = featureCount(data))
    .Object@feature_names <- featureNames(data)
    .Object@target_indices <- as.integer(target_indices)
    .Object@topologies <- list()
    
    lapply(seq(layers), function(layer)
    {
        filter <- new("mRMRe.Filter", data = data, prior_weight = prior_weight, target_indices = target_indices,
                levels = levels, ...)
        solutions <- solutions(filter, mi_threshold = mi_threshold, causality_threshold = causality_threshold)
        lapply(names(solutions), function(i) .Object@topologies[[i]] <<- solutions[[i]])
        
        # mi matrix
                        
        # causality
                        
        new_target_indices <- unique(unlist(solutions))
        target_indices <<- new_target_indices[!as.character(new_target_indices) %in% names(.Object@topologies)]
    })

    ## Perform last-layer linking
    
    filter <- new("mRMRe.Filter", data = data, prior_weight = prior_weight, target_indices = target_indices,
            levels = levels, ...)
    solutions <- solutions(filter, mi_threshold = mi_threshold, causality_threshold = causality_threshold)
    
    lapply(target_indices, function(target_index)
    {
        solution <- solutions[[as.character(target_index)]]
        new_solutions <- apply(solution, c(1, 2), function(feature_index)
                    ifelse(as.character(feature_index) %in% names(.Object@topologies), feature_index, 0))
        
        if (sum(new_solutions == 0) > 0)
            .Object@topologies[[as.character(target_index)]] <<- new_solutions
    })
    
    return(.Object)
})

## show

setMethod("show", signature("mRMRe.Network"), function(object)
{
    str(object)
})

## featureNames

setMethod("featureNames", signature("mRMRe.Network"), function(object)
{
    return(object@feature_names)
})

## mim

setMethod("mim", signature("mRMRe.Network"), function(object)
{

})

## causality

setMethod("causality", signature("mRMRe.Network"), function(object)
{

})

## adjacencyMatrix

setMethod("adjacencyMatrix", signature("mRMRe.Network"), function(object)
{
    
})

## visualize

setMethod("visualize", signature("mRMRe.Network"), function(object)
{
    ## FIXME : Cannot find a way to display vertex names...
    
    #adjacency <- adjacencyMatrix(object)
    #graph <- graph.adjacency(adjacency, mode = "undirected", add.rownames = TRUE)
    #V(graph)$name <- object@feature_names
    
    #return(plot.igraph(graph))
})