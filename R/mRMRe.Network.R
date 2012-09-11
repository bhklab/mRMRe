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
    .Object@mi_matrix <- matrix(NA, nrow = featureCount(data), ncol = featureCount(data))
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
            mi_matrix <- mim(filter, method = "cor")
            
            .Object@topologies[[target_index]] <<- solutions

            ## FIXME : Is there a more elegant/efficient way to combine MI matrices?

            lapply(seq(featureCount(data)^2), function(i)
            {
                if (is.na(.Object@mi_matrix[[i]]))
                    .Object@mi_matrix[[i]] <<- mi_matrix[[i]]
            })
            
    
            ## FIXME : Don't know yet if it is necessary to ensure symmetry in all
            ## three dimensions of the cube (-> symmetricCube)
    
            .Object@causality_cube[target_index, , ] <<- causality(filter)
            
            return(unlist(solutions))
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
    dimnames(cube)[[1]] <- featureNames(object)
    dimnames(cube)[[2]] <- featureNames(object)
    dimnames(cube)[[3]] <- featureNames(object)
    
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