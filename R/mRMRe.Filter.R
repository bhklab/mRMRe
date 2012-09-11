## Definition

setClass("mRMRe.Filter", representation(solutions = "matrix", mi_matrix = "matrix", feature_names = "character",
                target_index = "integer", levels = "integer", causality_matrix = "matrix"))

## Wrappers

`mRMR.ensemble` <- function(data, prior_weight, target_index, solution_count, feature_count, uses_ranks, outX,
        bootstrap_count)
{
    return(new("mRMRe.Filter", data = data, prior_weight = prior_weight, target_index = target_index,
                    levels = c(solution_count, rep(1, feature_count - 1)), uses_ranks = uses_ranks, outX = outX,
                    bootstrap_count = bootstrap_count))
}

`mRMR.classic` <- function(data, prior_weight, target_index, feature_count, uses_ranks, outX, bootstrap_count)
{
    return(mRMR.ensemble(data = data, prior_weight = prior_weight, target_index = target_index, solution_count = 1,
                    feature_count = feature_count, uses_ranks = uses_ranks, outX = outX,
                    bootstrap_count = bootstrap_count))
}

## initialize

setMethod("initialize", signature("mRMRe.Filter"),
        function(.Object, data, prior_weight, target_index, levels, uses_ranks = TRUE, outX = TRUE,
                bootstrap_count = 0)
{
    if (class(data) != "mRMRe.Data")
        stop("data must be of type mRMRe.Data")
    
    ## Prior processing
    
    if (length(priors(data)) != 0)
    {
        if (missing(prior_weight))
            stop("prior weight must be provided if there are priors")
        else if  (prior_weight < 0 || prior_weight > 1)
            stop("prior weight must be a value ranging from 0 to 1")
    }
    else
        prior_weight <- 0
    
    ## Target processing
    
    if (target_index < 1 || target_index > featureCount(data))
        stop("target_index must be a value ranging from 1 to the amount of features in data")
            
    ## Level processing
    
    ## FIXME : Not sure algorithm for combinatorial prediction is ok
    if (missing(levels))
        stop("levels must be provided")
    else if (prod(levels) - 1 > gamma(featureCount(data)) / gamma(featureCount(data) - length(levels)))
        stop("user cannot request for more solutions than is possible given the data set")
    
    .Object@target_index <- as.integer(c(target_index))
    .Object@levels <- as.integer(c(levels))
    
    target_index <- as.integer(expandFeatureIndices(data, target_index)) - 1
    
    wrap <- function(i) t(matrix(i[length(i):1], nrow = length(levels), ncol = length(i) / length(levels)))
    
    ## Filter
    
    filter <- .Call(mRMRe:::.C_export_filter, .Object@levels, as.vector(data@data), as.vector(data@priors),
            as.numeric(prior_weight), data@strata, data@weights, data@feature_types, nrow(data@data), ncol(data@data),
            as.integer(length(unique(data@strata))), target_index, uses_ranks, outX, bootstrap_count)
    filter$solutions <- wrap(filter$solutions)
    filter$mi_matrix <- matrix(filter$mi_matrix, ncol = sqrt(length(filter$mi_matrix)),
            nrow = sqrt(length(filter$mi_matrix)))
    
    .Object@solutions <- matrix(compressFeatureIndices(data, as.vector(filter$solutions) + 1),
            nrow = nrow(filter$solutions), ncol = ncol(filter$solutions))
    .Object@mi_matrix <- compressFeatureMatrix(data, filter$mi_matrix)
    .Object@feature_names <- featureNames(data)

    return(.Object)
})

## featureNames

setMethod("featureNames", signature("mRMRe.Filter"), function(.Object)
{
    return(.Object@feature_names)
})

## shrink

setMethod("shrink", signature("mRMRe.Filter"), function(.Object, mi_threshold, causality_threshold)
{
    solutions <- .Object@solutions
    
    if (!missing(mi_threshold))
    {
        ## FIXME: Not sure which direction priors are in, so you may have to inverse target_index and J
        
        solutions <- apply(solutions, 1, function(solution)
        {
            screen <- sapply(solution, function(feature) mi_threshold >= -.5 * log(1 -
                                        (.Object@mi_matrix[.Object@target_index, feature])))
            
            return(as.list(solution[screen]))
        })
    }
    
    if (!missing(causality_threshold))
    {
        causality_matrix <- causalityMatrix(.Object)
        
        solutions <- apply(solutions, 1, function(solution)
        {
            screen <- sapply(solution, function(feature) causality_threshold >=
                                max(causality_matrix[feature, solution]))
            
            return(as.list(solution[screen]))
        })
    }
    
    return(solutions)
})

## solutions

setMethod("solutions", signature("mRMRe.Filter"), function(.Object)
{
    return(.Object@solutions)
})

## mim

setMethod("mim", signature("mRMRe.Filter"), function(.Object)
{
    return(.Object@mi_matrix)
})

## causalityMatrix

setMethod("causalityMatrix", signature("mRMRe.Filter"), function(.Object)
{
    if (length(.Object@causality_matrix) == 0)
    {
        target_index <- .Object@target_index
        matrix <- matrix(NA, ncol = ncol(.Object@mi_matrix), nrow = ncol(.Object@mi_matrix))
        
        apply(.Object@solutions, 1, function(row)
        {
            pairs <- combn(row, 2)
            
            apply(pairs, 2, function(pair)
            {
                i <- pair[[1]]
                j <- pair[[2]]
                
                if (is.na(matrix[i, j]))
                {
                    coefficient <- -1/2 * log(((1 - .Object@mi_matrix[i, j]^2) * (1 - .Object@mi_matrix[i, target_index]^2)
                                        * (1 - .Object@mi_matrix[j, target_index]^2)) / (1 + 2 * .Object@mi_matrix[i, j] *
                                        .Object@mi_matrix[i, target_index] * .Object@mi_matrix[j, target_index] -
                                        .Object@mi_matrix[i, j]^2 - .Object@mi_matrix[i, target_index]^2 -
                                        .Object@mi_matrix[j, target_index]^2))
                    
                    matrix[i, j] <<- coefficient
                    matrix[j, i] <<- coefficient
                }
            })
        })

        .Object@causality_matrix <- matrix
    }

    return(.Object@causality_matrix)
})

## targetIndex

setMethod("targetIndex", signature("mRMRe.Filter"), function(.Object)
{
    return(.Object@target_index)
})
