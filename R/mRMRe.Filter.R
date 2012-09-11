## Definition

setClass("mRMRe.Filter", representation(solutions = "matrix", mi_matrix = "matrix", feature_names = "character",
                target_index = "integer", levels = "integer", causality_matrix = "matrix"))

## Wrappers

`mRMR.ensemble` <- function(solution_count, feature_count, ...)
{
    return(new("mRMRe.Filter", levels = c(solution_count, rep(1, feature_count - 1)), ...))
}

`mRMR.classic` <- function(feature_count, ...)
{
    return(new("mRMRe.Filter", levels = rep(1, feature_count), ...))
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

## show

setMethod("show", signature("mRMRe.Filter"), function(object)
{
    ## FIXME : Implement show method for this S4 class
    
    stop("No show method!")
})

## featureNames

setMethod("featureNames", signature("mRMRe.Filter"), function(object)
{
    return(object@feature_names)
})

## shrink

setMethod("shrink", signature("mRMRe.Filter"), function(object, mi_threshold, causality_threshold)
{
    solutions <- object@solutions
    
    if (!missing(mi_threshold))  
    {
        ## FIXME: Not sure which direction priors are in, so you may have to inverse target_index and J
        
        solutions <- apply(solutions, 1, function(solution)
        {
            screen <- sapply(solution, function(feature) mi_threshold >= -.5 * log(1 -
                                        (mim(object, method = "cor")[object@target_index, feature])))
            
            return(as.list(solution[screen]))
        })
    }
                                                                   
    if (!missing(causality_threshold))
    {
        solutions <- apply(solutions, 1, function(solution)
        {
            screen <- sapply(solution, function(feature) causality_threshold >=
                                max(causality(object)[feature, solution]))
            
            return(as.list(solution[screen]))
        })
    }
    
    if (!missing(mi_threshold) && !missing(causality_threshold))
        solutions <- apply(solutions, 1, as.list)
    
    return(solutions)
})

## solutions

setMethod("solutions", signature("mRMRe.Filter"), function(object)
{
    return(object@solutions)
})

## mim

setMethod("mim", signature("mRMRe.Filter"), function(object)
{
    return(object@mi_matrix)
})

## causality

setMethod("causality", signature("mRMRe.Filter"), function(object)
{
    if (length(object@causality_matrix) == 0)
    {
        target_index <- object@target_index
        matrix <- matrix(NA, ncol = ncol(object@mi_matrix), nrow = ncol(object@mi_matrix))
        
        apply(object@solutions, 1, function(row)
        {
            pairs <- combn(row, 2)
            
            apply(pairs, 2, function(pair)
            {
                i <- pair[[1]]
                j <- pair[[2]]
                
                if (is.na(matrix[i, j]))
                {
                    coefficient <- -1/2 * log(((1 - object@mi_matrix[i, j]^2) * (1 - object@mi_matrix[i, target_index]^2)
                                        * (1 - object@mi_matrix[j, target_index]^2)) / (1 + 2 * object@mi_matrix[i, j] *
                                        object@mi_matrix[i, target_index] * object@mi_matrix[j, target_index] -
                                        object@mi_matrix[i, j]^2 - object@mi_matrix[i, target_index]^2 -
                                        object@mi_matrix[j, target_index]^2))
                    
                    matrix[i, j] <<- coefficient
                    matrix[j, i] <<- coefficient
                }
            })
        })

        object@causality_matrix <- matrix
    }

    return(object@causality_matrix)
})

## target

setMethod("target", signature("mRMRe.Filter"), function(object)
{
    return(object@target_index)
})
