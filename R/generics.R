setGeneric("getData", function(.Object) standardGeneric("getData"))

setGeneric("getSampleCount", function(.Object) standardGeneric("getSampleCount"))

setGeneric("getFeatureCount", function(.Object) standardGeneric("getFeatureCount"))

setGeneric("getPriors", function(.Object) standardGeneric("getPriors"))

setGeneric("expandFeatureMatrix", function(.Object, matrix) standardGeneric("expandFeatureMatrix"))

setGeneric("compressFeatureMatrix", function(.Object, matrix) standardGeneric("compressFeatureMatrix"))

setGeneric("expandFeatureIndices", function(.Object, indices) standardGeneric("expandFeatureIndices"))

setGeneric("compressFeatureIndices", function(.Object, indices) standardGeneric("compressFeatureIndices"))

setGeneric("getSolutions", function(.Object) standardGeneric("getSolutions"))

setGeneric("getMutualInformationMatrix", function(.Object) standardGeneric("getMutualInformationMatrix"))

setGeneric("getTargetIndex", function(.Object) standardGeneric("getTargetIndex"))

setGeneric("getLevels", function(.Object) standardGeneric("getLevels"))
