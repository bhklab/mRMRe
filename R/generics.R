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

## FIXME : find a new name for this function 
setGeneric("freqMim", function(object) standardGeneric("freqMim"))

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

`correlate` <- function()
{
    ## FIXME : Recode correlate function ...
    
    #        double statistic; out[0]
    #        double concordant_weight; out[1]
    #        double discordant_weight; out[2]
    #        double uninformative_weight; out[3]
    #        double relevant_weight; out[4]
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