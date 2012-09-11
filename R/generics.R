setGeneric("featureData", function(.Object) standardGeneric("featureData"))

setGeneric("subsetData", function(.Object, row_indices, column_indices) standardGeneric("subsetData"))

setGeneric("sampleCount", function(.Object) standardGeneric("sampleCount"))

setGeneric("featureCount", function(.Object) standardGeneric("featureCount"))

setGeneric("featureNames", function(.Object) standardGeneric("featureNames"))

setGeneric("sampleStrata", function(.Object) standardGeneric("sampleStrata"))

setGeneric("sampleStrata<-", function(.Object, value) standardGeneric("sampleStrata<-"))

setGeneric("sampleWeights", function(.Object) standardGeneric("sampleWeights"))

setGeneric("sampleWeights<-", function(.Object, value) standardGeneric("sampleWeights<-"))

setGeneric("priors", function(.Object) standardGeneric("priors"))

setGeneric("mim", function(.Object, method = c("MI", "cor"), ...)
{
    method <- match.arg(method)
    matrix <- standardGeneric("mim")
    
    if (method == "MI")
        matrix <- -.5 * log(1 - (matrix^2))
    
    return(matrix)
})

setGeneric("expandFeatureMatrix", function(.Object, matrix) standardGeneric("expandFeatureMatrix"))

setGeneric("compressFeatureMatrix", function(.Object, matrix) standardGeneric("compressFeatureMatrix"))

setGeneric("expandFeatureIndices", function(.Object, indices) standardGeneric("expandFeatureIndices"))

setGeneric("compressFeatureIndices", function(.Object, indices) standardGeneric("compressFeatureIndices"))

setGeneric("shrink", function(.Object, mi_threshold, causality_threshold) standardGeneric("shrink"))

setGeneric("solutions", function(.Object) standardGeneric("solutions"))

setGeneric("causality", function(.Object) standardGeneric("causality"))

setGeneric("target", function(.Object) standardGeneric("target"))

setGeneric("adjacencyMatrix", function(.Object) standardGeneric("adjacencyMatrix"))

setGeneric("visualize", function(.Object) standardGeneric("visualize"))
