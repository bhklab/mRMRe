## Definition

setClass("mRMRe.Filter", representation(filters = "array", mi_matrix = "matrix", causality_matrix = "matrix",
                feature_names = "character", target_indices = "integer", levels = "integer"))

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
        function(.Object, data, prior_weight, target_indices, levels, continuous_estimator = "pearson", outX = TRUE,
                bootstrap_count = 0)
{
    if (class(data) != "mRMRe.Data")
        stop("data must be of type mRMRe.Data")
    
    ## Prior Processing
    
    if (length(priors(data)) != 0)
    {
        if (missing(prior_weight))
            stop("prior weight must be provided if there are priors")
        else if  (prior_weight < 0 || prior_weight > 1)
            stop("prior weight must be a value ranging from 0 to 1")
    }
    else
        prior_weight <- 0
    
    ## Target Processing
    
    if (sum(sapply(target_indices, function(index) index < 1 || index > featureCount(data))) > 1)
        stop("target_indices must only contain values ranging from 1 to the amount of features in data")
            
    ## Level Processing
    
    if (missing(levels))
        stop("levels must be provided")
    else if ((prod(levels) - 1) > choose(featureCount(data) - 1, length(levels)))
        stop("user cannot request for more solutions than is possible given the data set")
    
    .Object@target_indices <- as.integer(c(target_indices))
    .Object@levels <- as.integer(c(levels))
    
    target_indices <- as.integer(expandFeatureIndices(data, target_indices)) - 1
    
    ## Filter; Mutual Information and Causality Matrix

    mi_matrix <- as.numeric(matrix(NA, ncol = ncol(data@data), nrow = ncol(data@data)))
    filters <- vector(mode = "integer", length = length(target_indices) * prod(levels) * length(levels))
    
    .Call(mRMRe:::.C_export_filters, as.integer(.Object@levels), as.numeric(data@data),
            as.numeric(data@priors), as.numeric(prior_weight), as.integer(data@strata), as.numeric(data@weights),
            as.integer(data@feature_types), as.integer(nrow(data@data)), as.integer(ncol(data@data)),
            as.integer(length(unique(data@strata))), as.integer(target_indices),
            as.integer(mRMRe:::.map.continuous.estimator(continuous_estimator)), as.integer(outX),
            as.integer(bootstrap_count), mi_matrix, filters)
    
    .Object@feature_names <- featureNames(data)
    
    .Object@filters <- array(compressFeatureIndices(data, filters + 1), dim = c(length(levels), prod(levels),
                    length(target_indices)))
    .Object@filters <- apply(.Object@filters, c(2, 3), rev) # C code has feature vectors compiled from the bottom-up
    dimnames(.Object@filters) <- list(NULL, NULL, .Object@feature_names[.Object@target_indices])
    
    .Object@mi_matrix <- compressFeatureMatrix(data, matrix(mi_matrix, ncol = ncol(data@data), nrow = ncol(data@data)))
    
    .Object@causality_matrix <- matrix(ncol = length(target_indices), nrow = ncol(.Object@mi_matrix),
            dimnames = list(.Object@feature_names, .Object@feature_names[.Object@target_indices]))
    
    lapply(seq(target_indices), function(target_index_index)
    {
        target_index <- .Object@target_indices[[target_index_index]]
        
        apply(.Object@filters[, , target_index_index, drop = FALSE], 2, function(solution)
        {
            apply(combn(solution, 2), 2, function(pair)
            {
                i <- pair[[1]]
                j <- pair[[2]]
                
                cor_ij <- .Object@mi_matrix[i, j]
                
                if (abs(cor_ij) < abs(.Object@mi_matrix[j, i]))
                    cor_ij <- .Object@mi_matrix[j, i]
                
                coefficient <- -.5 * log(((1 - cor_ij^2) * (1 - .Object@mi_matrix[i, target_index]^2)
                                    * (1 - .Object@mi_matrix[j, target_index]^2)) / (1 + 2 * cor_ij *
                                    .Object@mi_matrix[i, target_index] * .Object@mi_matrix[j, target_index] -
                                    cor_ij^2 - .Object@mi_matrix[i, target_index]^2 -
                                    .Object@mi_matrix[j, target_index]^2))
                
                .Object@causality_matrix[i, target_index_index] <<- min(.Object@causality_matrix[i,
                                target_index_index], coefficient, na.rm = TRUE)
                
                .Object@causality_matrix[j, target_index_index] <<- min(.Object@causality_matrix[j, 
                                target_index_index], coefficient, na.rm = TRUE)
            })
        })
    })

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

setMethod("solutions", signature("mRMRe.Filter"), function(object, mi_threshold = -Inf, causality_threshold = Inf)
{
    # filters[, solution, target] is a vector of selected features
    # in a solution for a target. Features denoted by a missing value
    # have been purged by shrinkage
    
    lapply(seq(object@target_indices), function(target_index_index)
    {
        target_index <- object@target_indices[[target_index_index]]
        
        lapply(seq(prod(object@levels)), function(solution_index)
        {
            lapply(seq(object@levels), function(feature_index_index)
            {
                feature_index <- object@filters[feature_index_index, solution_index, target_index_index]
                
                if (mi_threshold > -.5 * log(1 - object@mi_matrix[feature_index, target_index]) ||
                        causality_threshold < object@causality_matrix[feature_index, target_index])
                    object@filters[feature_index_index, solution_index, target_index_index] <<- NA
            })
        })
    })
    
    return(object@filters)
})

## mim

setMethod("mim", signature("mRMRe.Filter"), function(object)
{
    # mi_matrix[i, j] contains the biased correlation between
    # features i and j (i -> j directionality)
    
    return(object@mi_matrix)
})

## causality

setMethod("causality", signature("mRMRe.Filter"), function(object)
{
    # causality_matrix[feature, target] contains the causality coefficient
    # between feature and target (feature -> target directionality)
    
    return(object@causality_matrix)
})
    
## target

setMethod("target", signature("mRMRe.Filter"), function(object)
{
    return(object@target_indices)
})
