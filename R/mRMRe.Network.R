setClass("mRMRe.Network", representation(topologies = "list", target_indices = "integer", levels = "integer"))

setMethod("initialize", "mRMRe.Network", function(.Object, data, prior_weight, target_indices, levels, layers, ...)
{
    if (missing(layers))
        layers <- 1L
    
    topologies <- list()
    length(topologies) <- getFeatureCount(data)
    
    lapply(seq(layers), function(layer)
    {
        target_indices <<- unlist(lapply(target_indices, function(target_index)
        {
            filter <- new("mRMRe.Filter", data = data, prior_weight = prior_weight, target_index = target_index,
                    levels = levels, ...)

            topologies[[target_index]] <<- getSolutions(filter)

            return(as.vector(getSolutions(filter)))
        }))

        target_indices <<- intersect(target_indices, which(sapply(topologies, is.null)))
    })

    .Object@topologies <- topologies

    return(.Object)
})

setMethod("getAdjacencyMatrix", "mRMRe.Network", function(.Object)
{
    matrix <- sapply(seq(.Object@topologies), function(i) sapply(seq(.Object@topologies), function(j)
    {
        if (i %in% .Object@topologies[[j]])
            return(1L)
        else
            return(0L)
    }))

    return(matrix)
})

setMethod("visualize", "mRMRe.Network", function(.Object)
{
    network <- getAdjacencyMatrix(.Object)
    
    return(plot.igraph(graph.adjacency(network)))
})