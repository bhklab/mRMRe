## Definition

setClass("mRMRe.Filter", representation(solutions = "list", mi_matrix = "matrix", feature_names = "character",
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
        function(.Object, data, prior_weight, target_index, levels, continuous_estimator = "pearson", outX = TRUE,
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
    
    if (missing(levels))
        stop("levels must be provided")
    else if ((prod(levels) - 1) > choose(featureCount(data) - 1, length(levels)))
        stop("user cannot request for more solutions than is possible given the data set")
    
    .Object@target_index <- as.integer(c(target_index))
    .Object@levels <- as.integer(c(levels))
    
    target_index <- as.integer(expandFeatureIndices(data, target_index)) - 1
    
    ## Filter

    mi_matrix <- as.numeric(matrix(NA, ncol = ncol(data@data), nrow = ncol(data@data)))
    
    solutions <- .Call(mRMRe:::.C_export_filter, as.integer(.Object@levels), as.numeric(data@data),
            as.numeric(data@priors), as.numeric(prior_weight), as.integer(data@strata), as.numeric(data@weights),
            as.integer(data@feature_types), as.integer(nrow(data@data)), as.integer(ncol(data@data)),
            as.integer(length(unique(data@strata))), as.integer(target_index),
            as.integer(mRMRe:::.map.continuous.estimator(continuous_estimator)), as.integer(outX),
            as.integer(bootstrap_count), mi_matrix)
    
    solutions <- matrix(compressFeatureIndices(data, as.integer(rev(solutions)) + 1), nrow = length(levels),
            ncol = length(solutions) / length(levels))
    solutions <- split(solutions, 1:ncol(solutions))
    names(solutions) <- NULL

    .Object@solutions <- solutions
    .Object@mi_matrix <- compressFeatureMatrix(data, matrix(mi_matrix, ncol = ncol(data@data), nrow = ncol(data@data)))
    .Object@feature_names <- featureNames(data)

    return(.Object)
})

## show

setMethod("show", signature("mRMRe.Filter"), function(object)
{
    str(object)
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
    
    if (!missing(mi_threshold) && mi_threshold != -Inf)  
    {
        solutions <- lapply(solutions, function(solution)
        {
            if (length(solution) == 0)
                return(solution)
                
            screen <- sapply(solution, function(feature) mi_threshold <= 
                                -.5 * log(1 - object@mi_matrix[feature, object@target_index]))
            
            return(solution[screen])
        })
    }
                                                                   
    if (!missing(causality_threshold) && causality_threshold != -Inf)
    {
        causality_matrix <- causality(object)
        
        solutions <- lapply(solutions, function(solution)
        {
            if (length(solution) == 0)
                return(solution)
            
            screen <- sapply(solution, function(feature) causality_threshold <=
                                max(causality_matrix[feature, solution], na.rm = TRUE))
            
            return(solution[screen])
        })
    }
    
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

## FIXME : cache matrix in Filter object
setMethod("causality", signature("mRMRe.Filter"), function(object)
{
    if (length(object@causality_matrix) == 0)
    {
        target_index <- object@target_index
        object@causality_matrix <- matrix(NA, ncol = ncol(object@mi_matrix), nrow = ncol(object@mi_matrix))
        
        lapply(object@solutions, function(row)
        {
            pairs <- combn(row, 2)
            
            apply(pairs, 2, function(pair)
            {
                j <- pair[[1]]
                i <- pair[[2]]
                
                if (is.na(object@causality_matrix[i, j]))
                {   
                    cor_ij <- max(object@mi_matrix[i, j], object@mi_matrix[j, i])
                    
                    coefficient <- -1/2 * log(((1 - cor_ij^2) * (1 - object@mi_matrix[i, target_index]^2)
                                        * (1 - object@mi_matrix[j, target_index]^2)) / (1 + 2 * cor_ij *
                                        object@mi_matrix[i, target_index] * object@mi_matrix[j, target_index] -
                                        cor_ij^2 - object@mi_matrix[i, target_index]^2 -
                                        object@mi_matrix[j, target_index]^2))
                    
                    object@causality_matrix[i, j] <<- coefficient
                    object@causality_matrix[j, i] <<- coefficient
                }
            })
        })
    }

    return(object@causality_matrix)
})

## target

setMethod("target", signature("mRMRe.Filter"), function(object)
{
    return(object@target_index)
})
