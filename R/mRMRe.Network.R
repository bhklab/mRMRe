## Definition

setClass("mRMRe.Network", representation(topologies = "list", feature_names = "character", target_indices = "integer",
                levels = "integer"))

## Wrappers

## FIXME: Add wrappers for network

## initialize

setMethod("initialize", signature("mRMRe.Network"), function(.Object, data, prior_weight, target_indices, levels,
                layers)
{
    if (missing(layers))
        layers <- 1L
    
    topologies <- list() #array(dim = c(featureCount(data), prod(levels), 3))
    length(topologies) <- featureCount(data)
    
    # z = 1 is for solutions
    # z = 2 is for MIs
    # z = 3 is for causality
    
    lapply(seq(layers), function(layer)
    {
        target_indices <<- unlist(lapply(target_indices, function(target_index)
        {
            filter <- new("mRMRe.Filter", data = data, prior_weight = prior_weight, target_index = target_index,
                    levels = levels)

            solutions <- solutions(filter)
            
            topologies[[target_index]] <<- solutions

            return(as.vector(solutions))
        }))

        target_indices <<- intersect(target_indices, which(sapply(topologies, is.null)))
    })

    .Object@topologies <- topologies
    .Object@feature_names <- featureNames(data)

    return(.Object)
})

## adjacencyMatrix

setMethod("adjacencyMatrix", signature("mRMRe.Network"), function(.Object)
{
    matrix <- sapply(seq(.Object@topologies), function(i) sapply(seq(.Object@topologies), function(j)
    {
        if (i %in% .Object@topologies[[j]])
            return(1L)
        else
            return(0L)
    }))

    rownames(matrix) <- .Object@feature_names
    colnames(matrix) <- .Object@feature_names

    return(t(matrix))
})

## visualize

setMethod("visualize", signature("mRMRe.Network"), function(.Object)
{
    network <- adjacencyMatrix(.Object)
    
    return(plot.igraph(graph.adjacency(network)))
})