setClass("mRMRe.Filter", representation(solutions = "matrix", mi_matrix = "matrix", prior_weight = "numeric",
                target_index = "integer", levels = "integer", uses_ranks = "logical", outX = "logical",
                bootstrap_count = "integer"))

setMethod("initialize", "mRMRe.Filter",
        function(.Object, data, prior_weight, target_index, levels, uses_ranks, outX, bootstrap_count)
{
    if (class(data) != "mRMRe.Data")
        stop("data must be of type mRMRe.Data")
    
    if (length(getPriors(data)) != 0)
    {
        if (missing(prior_weight))
            stop("prior weight must be provided if there are priors")
        else if  (prior_weight < 0 || prior_weight > 1)
            stop("prior weight must be a value ranging from 0 to 1")
    }
    else
        prior_weight <- 0
    
    if (target_index < 1 || target_index > getFeatureCount(data))
        stop("target_index must be a value ranging from 1 to the amount of features in data")
            
    if (missing(levels))
        stop("levels must be provided")
    else if (prod(levels) > gamma(getFeatureCount(data) + 1) - 1)
        stop("user cannot request for more solutions than is possible given the data set")
    
    if (missing(uses_ranks))
        uses_ranks <- TRUE
    
    if (missing(outX))
        outX <- TRUE

    if (missing(bootstrap_count))
        bootstrap_count <- 0
    
    .Object@prior_weight <- as.numeric(prior_weight) ## No getter method
    .Object@target_index <- as.integer(c(target_index))
    .Object@levels <- as.integer(c(levels))
    .Object@uses_ranks <- as.logical(uses_ranks) ## No getter method
    .Object@outX <- as.logical(outX) ## No getter method
    .Object@bootstrap_count <- as.integer(bootstrap_count) ## No getter method
    
    target_index <- expandFeatureIndices(data, target_index)
    
    wrap <- function(i) t(matrix(i[length(i):1], nrow=length(levels), ncol=length(i)/length(levels)))
    
    filter <- .Call(mRMRe:::.C_export_filter, .Object@levels, as.vector(data@data), as.vector(data@priors),
            .Object@prior_weight, data@strata, data@weights, data@feature_types, nrow(data@data),
            ncol(data@data), as.integer(length(unique(data@strata))), .Object@target_index - 1, .Object@uses_ranks,
            .Object@outX, .Object@bootstrap_count)
    filter$solutions <- wrap(filter$solutions)
    filter$mi_matrix <- matrix(filter$mi_matrix, ncol=sqrt(length(filter$mi_matrix)),
            nrow=sqrt(length(filter$mi_matrix)))
    
    .Object@solutions <- matrix(compressFeatureIndices(data, as.vector(filter$solutions) + 1),
            nrow=nrow(filter$solutions), ncol=ncol(filter$solutions))
    .Object@mi_matrix <- compressFeatureMatrix(data, filter$mi_matrix)

    return(.Object)
})

setMethod("getSolutions", "mRMRe.Filter", function(.Object)
{
    return(.Object@solutions)
})

setMethod("getMutualInformationMatrix", "mRMRe.Filter", function(.Object)
{
    return(.Object@mi_matrix)
})

setMethod("getTargetIndex", "mRMRe.Filter", function(.Object)
{
    return(.Object@target_index)
})

setMethod("getLevels", "mRMRe.Filter", function(.Object)
{
    return(.Object@levels)
})
