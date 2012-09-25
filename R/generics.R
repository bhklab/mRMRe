setGeneric("featureData", function(object) standardGeneric("featureData"))

setGeneric("subsetData", function(object, ...) standardGeneric("subsetData"))

setGeneric("sampleCount", function(object) standardGeneric("sampleCount"))

setGeneric("featureCount", function(object) standardGeneric("featureCount"))

setGeneric("featureNames", function(object) standardGeneric("featureNames"))

setGeneric("sampleStrata", function(object) standardGeneric("sampleStrata"))

setGeneric("sampleStrata<-", function(object, value) standardGeneric("sampleStrata<-"))

setGeneric("sampleWeights", function(object) standardGeneric("sampleWeights"))

setGeneric("sampleWeights<-", function(object, value) standardGeneric("sampleWeights<-"))

setGeneric("priors", function(object) standardGeneric("priors"))

setGeneric("priors<-", function(object, value) standardGeneric("priors<-"))

setGeneric("mim", function(object, method = c("MI", "cor"), ...)
{
    method <- match.arg(method)
    matrix <- standardGeneric("mim")
    
    if (method == "MI")
        matrix <- -.5 * log(1 - (matrix^2))
    
    return(matrix)
})

setGeneric("expandFeatureMatrix", function(object, ...) standardGeneric("expandFeatureMatrix"))

setGeneric("compressFeatureMatrix", function(object, ...) standardGeneric("compressFeatureMatrix"))

setGeneric("expandFeatureIndices", function(object, ...) standardGeneric("expandFeatureIndices"))

setGeneric("compressFeatureIndices", function(object, ...) standardGeneric("compressFeatureIndices"))

setGeneric("shrink", function(object, ...) standardGeneric("shrink"))

setGeneric("solutions", function(object) standardGeneric("solutions"))

setGeneric("causality", function(object) standardGeneric("causality"))

setGeneric("target", function(object) standardGeneric("target"))

setGeneric("adjacencyMatrix", function(object) standardGeneric("adjacencyMatrix"))

setGeneric("visualize", function(object) standardGeneric("visualize"))

`.map.continuous.estimator` <- function(continuous_estimator)
{
    value <- switch(continuous_estimator, "pearson" = 0L, "spearman" = 1L, "kendall" = 2L, "frequency" = 3L, -1L)
    
    if (value < 0L || value > 4L || !is.character(continuous_estimator))
        stop("estimator must be of the following: pearson, spearman, kendall, frequency")
    
    return(value)
}

`correlate` <- function(X, Y, method = "pearson", strata, weights, outX = TRUE, bootstrap_count = 0)
{
    if (method == "pearson" || method == "spearman" || method == "kendall")
    {
        X <- as.numeric(X)
        Y <- as.numeric(Y)
    }
    else if (method == "cramersv")
    {
        X <- as.factor(X)
        Y <- as.factor(Y)
    }
    else if (method != "cindex")
        stop("estimator must be of the following: pearson, spearman, kendall, frequency, cramersv, cindex")
    
    if (!missing(strata) && !missing(weights))
        data <- mRMR.data(data = data.frame(X, Y), strata = strata, weights = weights)
    else if (!missing(strata))
        data <- mRMR.data(data = data.frame(X, Y), strata = strata)
    else if (!missing(weights))
        data <- mRMR.data(data = data.frame(X, Y), weights = weights)
    else
        data <- mRMR.data(data = data.frame(X, Y))
    
    if (method == "cindex")
    {
        empty <- vector(mode = "numeric", length = 0)
        
        if (length(data@feature_types) == 2)
            input <- list(data@data[, 1], data@data[, 2], empty, empty)
        else if (length(data@feature_types) == 3)
        {
            if (data@feature_types[[1]] == 2)
                input <- list(data@data[, 1], data@data[, 3], data@data[, 2], empty)
            else if (data@feature_types[[2]] == 2)
                input <- list(data@data[, 2], data@data[, 1], data@data[, 3], empty)
        }
        else if (length(data@feature_types) == 4)
            input <- list(data@data[, 1], data@data[, 3], data@data[, 2], data@data[, 4])
        
        out <- vector(mode = "numeric", length = 5)
        
        .Call(mRMRe:::.C_export_concordance_index, as.numeric(input[[1]]), as.numeric(input[[2]]),
                as.numeric(input[[3]]), as.numeric(input[[4]]), as.integer(data@strata), as.numeric(data@weights),
                as.integer(length(unique(data@strata))), outX, out)
        
        names(out) <- c("statistic", "concordant_weight", "discordant_weight", "uninformative_weight",
                "relevant_weight")
        
        return(out)
    }
    else if (method == "cramersv")
        return(list(statistic = mim(data, method = "cor", outX = outX, bootstrap_count = bootstrap_count)[1, 2]))
    else
        return(list(statistic = mim(data, method = "cor", continuous_estimator = method, outX = outX,
                            bootstrap_count = bootstrap_count)[1, 2]))
}

`get.thread.count` <- function()
{
    thread_count <- vector(mode = "integer", length = 1)
    
    .Call(mRMRe:::.C_get_thread_count, thread_count)
    
    return(thread_count)
}

`set.thread.count` <- function(thread_count)
{
    thread_count <- as.integer(thread_count)
    
    .Call(mRMRe:::.C_set_thread_count, thread_count)
    
    return(thread_count)
}