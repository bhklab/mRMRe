setGeneric("featureData", function(.Object) standardGeneric("featureData"))

setGeneric("sampleCount", function(.Object) standardGeneric("sampleCount"))

setGeneric("featureCount", function(.Object) standardGeneric("featureCount"))

setGeneric("featureNames", function(.Object) standardGeneric("featureNames"))

setGeneric("priors", function(.Object) standardGeneric("priors"))

setGeneric("mim", function(.Object, method = c("cor", "MI"), ...)
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

setGeneric("solutions", function(.Object) standardGeneric("solutions"))

setGeneric("causalityMatrix", function(.Object) standardGeneric("causalityMatrix"))

setGeneric("targetIndex", function(.Object) standardGeneric("targetIndex"))

setGeneric("adjacencyMatrix", function(.Object) standardGeneric("adjacencyMatrix"))

setGeneric("visualize", function(.Object) standardGeneric("visualize"))
