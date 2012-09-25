## Definition

setClass("mRMRe.Network", representation(topologies = "list", mi_matrix = "matrix", causality_cube = "array",
                feature_names = "character"))

## Wrappers

## FIXME: Add wrappers for network

## initialize

setMethod("initialize", signature("mRMRe.Network"), function(.Object, data, prior_weight, target_indices, levels,
                layers, ..., mi_threshold = -Inf, causality_threshold = -Inf)
{
    if (missing(layers))
        layers <- 1L
    
    .Object@topologies <- list()
    .Object@mi_matrix <- matrix(NA, nrow = featureCount(data), ncol = featureCount(data))
    #.Object@causality_cube <- array(dim = rep(featureCount(data), 3))
    .Object@feature_names <- featureNames(data)
    
    length(.Object@topologies) <- featureCount(data)
    
    lapply(seq(layers), function(layer)
    {
        target_indices <<- unlist(lapply(target_indices, function(target_index)
        {
            filter <- new("mRMRe.Filter", data = data, prior_weight = prior_weight, target_index = target_index,
                    levels = levels, ...)

            ## FIXME : Code might not handle too brutal of a cutoff (mi_threshold and causality_threshold...)
            
            solutions <- shrink(filter, mi_threshold = mi_threshold, causality_threshold = causality_threshold)
            #mi_matrix <- mim(filter, method = "cor")
            
            .Object@topologies[[target_index]] <<- solutions

            ## FIXME : Is there a more elegant/efficient way to combine MI matrices?

            #screen <- sapply(seq(featureCount(data)^2), function(i) is.na(.Object@mi_matrix[[i]]))
            #.Object@mi_matrix[screen] <<- mi_matrix[screen]
            
            ## FIXME : Don't know yet if it is necessary to ensure symmetry in all
            ## three dimensions of the cube (-> symmetricCube)
    
            #.Object@causality_cube[target_index, , ] <<- causality(filter)
            
            return(unlist(solutions))
        }))
    
        target_indices <<- intersect(target_indices, which(sapply(.Object@topologies, is.null)))
    })

    ## Perform last-layer linking
    
    
    ## FIXME: Validity checks and efficiency assessment needed here
    
    lapply(target_indices, function(target_index)
    {
        filter <- new("mRMRe.Filter", data = data, prior_weight = prior_weight, target_index = target_index,
                levels = levels, ...)
        
        solutions <- shrink(filter, mi_threshold = mi_threshold, causality_threshold = causality_threshold)
        
        new_solutions <- lapply(solutions, function(solution) as.list(unlist(lapply(solution, function(feature)
        {
            if (!is.null(.Object@topologies[[feature]]))
                return(feature)
            else
                return(NULL)
        }))))   
            
        .Object@topologies[[target_index]] <<- new_solutions
    })

    return(.Object)
})

## show

setMethod("show", signature("mRMRe.Network"), function(object)
{
    ## FIXME: More detailed show method?
    
    print(str(object))
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
                            ifelse(i %in% unlist(object@topologies[[j]]), 1L, 0L)))

    rownames(matrix) <- object@feature_names
    colnames(matrix) <- object@feature_names

    return(t(matrix))
})

## visualize

setMethod("visualize", signature("mRMRe.Network"), function(object)
{
    ## FIXME : Cannot find a way to display vertex names...
    
    adjacency <- adjacencyMatrix(object)
    graph <- graph.adjacency(adjacency, mode = "undirected", add.rownames = TRUE)
    V(graph)$name <- object@feature_names
    
    return(plot.igraph(graph))
})