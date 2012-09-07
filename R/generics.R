## mRMRe.Data.R

# initialize

setGeneric("getData", function(.Object) standardGeneric("getData"))

setGeneric("getSampleCount", function(.Object) standardGeneric("getSampleCount"))

setGeneric("getFeatureCount", function(.Object) standardGeneric("getFeatureCount"))

setGeneric("getPriors", function(.Object) standardGeneric("getPriors"))

setGeneric("getMutualInformationMatrix", function(.Object, ...) standardGeneric("getMutualInformationMatrix"))

setGeneric("expandFeatureMatrix", function(.Object, matrix) standardGeneric("expandFeatureMatrix"))

setGeneric("compressFeatureMatrix", function(.Object, matrix) standardGeneric("compressFeatureMatrix"))

setGeneric("expandFeatureIndices", function(.Object, indices) standardGeneric("expandFeatureIndices"))

setGeneric("compressFeatureIndices", function(.Object, indices) standardGeneric("compressFeatureIndices"))

## mRMRe.Filter.R

# initialize

setGeneric("getSolutions", function(.Object) standardGeneric("getSolutions"))

setGeneric("getMutualInformationMatrix", function(.Object) standardGeneric("getMutualInformationMatrix"))

setGeneric("getTargetIndex", function(.Object) standardGeneric("getTargetIndex"))

setGeneric("getLevels", function(.Object) standardGeneric("getLevels"))

## mRMRe.Network.R

# initialize

setGeneric("getAdjacencyMatrix", function(.Object) standardGeneric("getAdjacencyMatrix"))

setGeneric("visualize", function(.Object) standardGeneric("visualize"))