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
    .Object@mi_matrix <- matrix(nrow = featureCount(data), ncol = featureCount(data))
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

## show

setMethod("show", signature("mRMRe.Network"), function(object)
{
    ## FIXME : Implement show method for this S4 class
    
    stop("No show method!")
})

## featureNames

setMethod("featureNames", signature("mRMRe.Network"), function(object)
{
    return(object@feature_names)
})

## mim

setMethod("mim", signature("mRMRe.Network"), function(object)
{
    matrix <- object@mi_matrix
    rownames(matrix) <- featureNames(object)
    colnames(matrix) <- featureNames(object)
    
    return(matrix)
})

## causality

setMethod("causality", signature("mRMRe.Network"), function(object)
{
    cube <- object@causality_cube
    dimnames(cube) <- as.list(rep(featureNames(object), 3))
    
    return(cube)
})

## adjacencyMatrix

setMethod("adjacencyMatrix", signature("mRMRe.Network"), function(object)
{
    matrix <- sapply(seq(object@topologies), function(i) sapply(seq(object@topologies), function(j)
    {
        if (i %in% object@topologies[[j]])
            return(1L)
        else
            return(0L)
    }))

    rownames(matrix) <- object@feature_names
    colnames(matrix) <- object@feature_names

    return(t(matrix))
})

## visualize

setMethod("visualize", signature("mRMRe.Network"), function(object)
{
    network <- adjacencyMatrix(object)
    
    return(plot.igraph(graph.adjacency(network)))
})