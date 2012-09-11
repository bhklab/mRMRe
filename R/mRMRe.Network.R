## Definition

setClass("mRMRe.Network", representation(topologies = "list", mi_matrix = "matrix", causality_cube = "array",
                feature_names = "character"))

## Wrappers

## FIXME: Add wrappers for network

## initialize

setMethod("initialize", signature("mRMRe.Network"), function(.Object, data, prior_weight, target_indices, levels,
                layers)
{
    if (missing(layers))
        layers <- 1L
    
    .Object@topologies <- list()
    .Object@mi_matrix <- matrix(ncol = featureCount(data), ncol = featureCount(data))
    .Object@causality_cube <- array(dim = rep(featureCount(data), 3))
    .Object@feature_names <- featureNames(data)
    
    length(.Object@topologies) <- featureCount(data)
    
    lapply(seq(layers), function(layer)
    {
        target_indices <<- unlist(lapply(target_indices, function(target_index)
        {
            filter <- new("mRMRe.Filter", data = data, prior_weight = prior_weight, target_index = target_index,
                    levels = levels)

            solutions <- shrink(filter)
            
            .Object@topologies[[target_index]] <<- solutions

            return(as.vector(solutions))
        }))

        target_indices <<- intersect(target_indices, which(sapply(.Object@topologies, is.null)))
    })

    return(.Object)
})

## featureNames

setMethod("featureNames", signature("mRMRe.Network"), function(.Object)
{
    return(.Object@feature_names)
})

## mim

setMethod("mim", signature("mRMRe.Network"), function(.Object)
{
    matrix <- .Object@mi_matrix
    rownames(matrix) <- featureNames(.Object)
    colnames(matrix) <- featureNames(.Object)
    
    return(matrix)
})

## causality

setMethod("causality", signature("mRMRe.Network"), function(.Object)
{
    cube <- .Object@causality_cube
    dimnames(cube) <- as.list(rep(featureNames(.Object), 3))
    
    return(cube)
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