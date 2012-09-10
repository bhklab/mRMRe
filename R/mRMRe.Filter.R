## Definition

setClass("mRMRe.Filter", representation(solutions = "matrix", mi_matrix = "matrix", feature_names = "character",
                target_index = "integer", levels = "integer"))

## initialize

setMethod("initialize", signature("mRMRe.Filter"),
        function(.Object, data, prior_weight, target_index, levels, uses_ranks = TRUE, outX = TRUE,
                bootstrap_count = 0)
{
    if (class(data) != "mRMRe.Data")
        stop("data must be of type mRMRe.Data")
    
    ## Prior processing
    
    if (length(getPriors(data)) != 0)
    {
        if (missing(prior_weight))
            stop("prior weight must be provided if there are priors")
        else if  (prior_weight < 0 || prior_weight > 1)
            stop("prior weight must be a value ranging from 0 to 1")
    }
    else
        prior_weight <- 0
    
    ## Target processing
    
    if (target_index < 1 || target_index > getFeatureCount(data))
        stop("target_index must be a value ranging from 1 to the amount of features in data")
            
    ## Level processing
    
    if (missing(levels))
        stop("levels must be provided")
    else if (prod(levels) - 1 > gamma(getFeatureCount(data)) / gamma(getFeatureCount(data) - length(levels)))
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
    .Object@feature_names <- getFeatureNames(data)

    return(.Object)
})

## getFeatureNames

setMethod("getFeatureNames", signature("mRMRe.Filter"), function(.Object)
{
    return(.Object@feature_names)
})

## getSolutions

setMethod("getSolutions", signature("mRMRe.Filter"), function(.Object)
{
    return(.Object@solutions)
})

## getCausalityMatrix

setMethod("getCausalityMatrix", signature("mRMRe.Filter"), function(.Object)
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
    
    return(matrix)
})

## getMutualInformationMatrix

setMethod("getMutualInformationMatrix", signature("mRMRe.Filter"), function(.Object)
{
    return(.Object@mi_matrix)
})

## getTargetIndex

setMethod("getTargetIndex", signature("mRMRe.Filter"), function(.Object)
{
    return(.Object@target_index)
})

## getLevels

setMethod("getLevels", signature("mRMRe.Filter"), function(.Object)
{
    return(.Object@levels)
})
