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

`set.thread.count` <- function(thread_count)
{
    .Call(mRMRe:::.C_set_thread_count, as.integer(thread_count))
}