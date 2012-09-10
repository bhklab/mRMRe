setGeneric("getData", function(.Object) standardGeneric("getData"))

setGeneric("getSampleCount", function(.Object) standardGeneric("getSampleCount"))

setGeneric("getFeatureCount", function(.Object) standardGeneric("getFeatureCount"))

setGeneric("getFeatureNames", function(.Object) standardGeneric("getFeatureNames"))

setGeneric("getPriors", function(.Object) standardGeneric("getPriors"))

setGeneric("getMutualInformationMatrix", function(.Object, method = c("cor", "MI"), ...)
{
    method <- match.arg(method)
    matrix <- standardGeneric("getMutualInformationMatrix")
    
    if (method == "MI")
        matrix <- -.5 * log(1 - (matrix^2))
    
    return(matrix)
})

setGeneric("expandFeatureMatrix", function(.Object, matrix) standardGeneric("expandFeatureMatrix"))

setGeneric("compressFeatureMatrix", function(.Object, matrix) standardGeneric("compressFeatureMatrix"))

setGeneric("expandFeatureIndices", function(.Object, indices) standardGeneric("expandFeatureIndices"))

setGeneric("compressFeatureIndices", function(.Object, indices) standardGeneric("compressFeatureIndices"))

setGeneric("getSolutions", function(.Object) standardGeneric("getSolutions"))

setGeneric("getTargetIndex", function(.Object) standardGeneric("getTargetIndex"))

setGeneric("getLevels", function(.Object) standardGeneric("getLevels"))

setGeneric("getAdjacencyMatrix", function(.Object) standardGeneric("getAdjacencyMatrix"))

setGeneric("visualize", function(.Object) standardGeneric("visualize"))
