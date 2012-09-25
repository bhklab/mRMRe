setGeneric("featureData", function(object) standardGeneric("featureData"))

setGeneric("subsetData", function(object, row_indices, column_indices) standardGeneric("subsetData"))

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

setGeneric("expandFeatureMatrix", function(object, matrix) standardGeneric("expandFeatureMatrix"))

setGeneric("compressFeatureMatrix", function(object, matrix) standardGeneric("compressFeatureMatrix"))

setGeneric("expandFeatureIndices", function(object, indices) standardGeneric("expandFeatureIndices"))

setGeneric("compressFeatureIndices", function(object, indices) standardGeneric("compressFeatureIndices"))

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
        stop("please provide one of the following continuous estimators: pearson, spearman, kendall, frequency")
    
    return(value)
}

`correlate` <- function(X, Y, method = "pearson", strata, weights, outX, bootstrap_count)
{
    data <- mRMR.data(data = data.frame(X, Y), strata = strata, weights = weights)
    
    if (method == "cindex")
    {
        empty <- vector(mode = "numeric", length = 0)
        
        a <- empty
        b <- empty
        c <- empty
        d <- empty
        
        if (length(data@feature_types) == 2)
        {
            a <- data@data[, 1]
            b <- data@data[, 2]
        }
        else if (length(data@feature_types) == 3)
        {
            if (data@feature_types[[1]] == 2)
            {
                a <- data@data[, 1]
                b <- data@data[, 3]
                c <- data@data[, 2]
            }
            else if (data@feature_types[[2]] == 2)
            {
                a <- data@data[, 2]
                b <- data@data[, 1]
                c <- data@data[, 3]
            }
        }
        else if (length(data@feature_types) == 4)
        {
            a <- data@data[, 1]
            b <- data@data[, 3]
            c <- data@data[, 2]
            d <- data@data[, 4]
        }
        
        out <- vector(mode = "numeric", length = 5)
        
        .Call(mRMRe:::.C_export_concordance_index, as.numeric(a), as.numeric(b), as.numeric(c), as.numeric(d),
                as.integer(data@strata), as.numeric(data@weights), as.integer(length(unique(data@strata))),
                outX, out)
        
        names(out) <- c("statistic", "concordant_weight", "discordant_weight", "uninformative_weight",
                "relevant_weight")
        
        return(out)
    }
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