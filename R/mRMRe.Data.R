## Definition

setClass("mRMRe.Data", representation(feature_names = "character", feature_types = "numeric", data = "matrix",
                strata = "numeric", weights = "numeric", priors = "matrix", mi_matrix = "matrix"))

## Wrapper

`mRMR.data` <- function(...)
{
    return(new("mRMRe.Data", ...))
}

## initialize

setMethod("initialize", signature("mRMRe.Data"), function(.Object, data, strata, weights, priors)
{
    ## Data processing
    
    if (!is.data.frame(data))
        stop("data must be of type data frame")
    
    feature_types <- sapply(data, function(feature) paste(class(feature), collapse = "_"))
    
    if (any(!is.element(feature_types, c("numeric", "ordered_factor", "Surv"))))
        stop("data columns must be either of numeric, ordered factor or Surv type")
    
    .Object@feature_names <- colnames(data)
    .Object@feature_types <- unlist(lapply(feature_types, switch, "Surv" = c(2, 3), "ordered_factor" = 1, 0))
    names(.Object@feature_types) <- NULL
    
    .Object@data <- do.call(cbind, lapply(seq(feature_types), function(i) switch(feature_types[[i]],
                                "Surv" = cbind(event = data[, i][, "status"], time = data[, i][, "time"]),
                                "ordered_factor" = as.numeric(as.integer(data[, i]) - 1),
                                as.numeric(data[, i]))))
    
    rownames(.Object@data) <- rownames(data)
    colnames(.Object@data)[!.Object@feature_types %in% c(2, 3)] <- colnames(data)[feature_types != "Surv"]
    colnames(.Object@data)[.Object@feature_types %in% c(2, 3)] <- paste(rep(colnames(data)[feature_types == "Surv"],
                    each = 2), rep(c("event", "time"), sum(feature_types == "Surv")), sep = "_")
    
    ## Sample stratum processing
    
    if (missing(strata)) 
        .Object@strata <- rep.int(0, nrow(data))
    else
        strata(.Object) <- strata
    
    ## Sample weight processing
    
    if (missing(weights)) 
        .Object@weights <- rep(1, nrow(data))
    else
        weights(.Object) <- weights

    ## Prior feature matrix processing
    
    if (!missing(priors))
    {
        if (ncol(priors) != ncol(data) || nrow(priors) != ncol(data))
            stop("priors matrix must be a symmetric matrix containing as many features as data")
        else
            .Object@priors <- expandFeatureMatrix(.Object, priors)
    }
    
    return(.Object)
})

## show

setMethod("show", signature("mRMRe.Data"), function(object)
{
    ## FIXME : Implement show method for this S4 class
    
    stop("No show method!")
})

## featureData

setMethod("featureData", signature("mRMRe.Data"), function(object)
{
    data <- lapply(seq(object@feature_types), function(i) switch(as.character(object@feature_types[[i]]),
                        "3" = Surv(time = object@data[, i], event = object@data[, i - 1]),
                        "2" = NULL,
                        "1" = object@data[, i] + 1,
                        "0" = object@data[, i],
                        NULL))
    data <- data.frame(data[!sapply(data, is.null)])
    colnames(data) <- object@feature_names
    
    return(data)
})

## subsetData

setMethod("subsetData", signature("mRMRe.Data"), function(object, row_indices, column_indices)
{
    ## FIXME : Inefficient, as everything is compressed and then re-expanded
    
    data <- featureData(object)[row_indices, column_indices, drop=FALSE]
    strata <- sampleStrata(object)[row_indices]
    weights <- sampleWeights(object)[row_indices]
    priors <- priors(object)[row_indices, column_indices, drop=FALSE]
    
    return(new("mRMRe.Data", data = data, strata = strata, weights = weight, priors = priors))
})

## sampleCount

setMethod("sampleCount", signature("mRMRe.Data"), function(object)
{
    return(nrow(object@data))
})

## featureCount

setMethod("featureCount", signature("mRMRe.Data"), function(object)
{
    return(length(object@feature_names))
})

## featureNames

setMethod("featureNames", signature("mRMRe.Data"), function(object)
{
    return(object@feature_names)
})

## sampleStrata

setMethod("sampleStrata", signature("mRMRe.Data"), function(object)
{
    strata <- object@strata
    names(strata) <- rownames(object@data)
    
    return(strata)
})

## sampleStrata<-

setReplaceMethod("sampleStrata", signature("mRMRe.Data"), function(object, value)
{
    if (length(value) != nrow(object@data))
        stop("data and strata must contain the same number of samples")
    else if (!is.factor(value))
        stop("strata must be provided as factors")
    else
        object@strata <- as.integer(value) - 1
})

## sampleWeights

setMethod("sampleWeights", signature("mRMRe.Data"), function(object)
{
    weights <- object@weights
    names(weights) <- rownames(object@data)
    
    return(weights)
})

## sampleWeights<-

setReplaceMethod("sampleWeights", signature("mRMRe.Data"), function(object, value)
{
    if (length(value) != nrow(object@data))
        stop("data and weight must contain the same number of samples")
    else
        object@weights <- as.numeric(value)
})

## priors

setMethod("priors", signature("mRMRe.Data"), function(object)
{
    if (length(object@priors) == 0)
        return(object@priors)
    else
        return(compressFeatureMatrix(object, object@priors))
})

## mim

setMethod("mim", signature("mRMRe.Data"),
        function(object, prior_weight = 0, uses_ranks = TRUE, outX = TRUE, bootstrap_count = 0)
{
    if (length(object@mi_matrix) == 0)
    {
        if (length(object@priors) != 0)
        {
            if (missing(prior_weight))
                stop("prior weight must be provided if there are priors")
            else if  (prior_weight < 0 || prior_weight > 1)
                stop("prior weight must be a value ranging from 0 to 1")
        }
        else
            prior_weight <- 0
        
        mi_matrix <- .Call(mRMRe:::.C_export_mim, as.vector(object@data), as.vector(object@priors),
                as.numeric(prior_weight), object@strata, object@weights, object@feature_types, nrow(object@data),
                ncol(object@data), as.integer(length(unique(object@strata))), as.integer(uses_ranks), as.integer(outX),
                as.integer(bootstrap_count))
        mi_matrix <- matrix(mi_matrix, nrow = ncol(object@data), ncol = ncol(object@data))
        
        object@mi_matrix <- compressFeatureMatrix(object, mi_matrix)
    }
    
    return(object@mi_matrix)
})

## expandFeatureMatrix

setMethod("expandFeatureMatrix", signature("mRMRe.Data"), function(object, matrix)
{
    adaptor <- which(object@feature_types != 3)
    matrix <- do.call(cbind, lapply(seq(adaptor), function(i)
    {
        column <- do.call(rbind, lapply(seq(adaptor), function(j)
        {
            item <- matrix[j, i]
            
            if (object@feature_types[[adaptor[[j]]]] == 2)
                return(rbind(item, item, deparse.level = 0))
            else
                return(item)
        }))
        
        if (object@feature_types[[adaptor[[i]]]] == 2)
            return(cbind(column, column, deparse.level = 0))
        else
            return(column)
    }))

    return(matrix)
})

## compressFeatureMatrix

setMethod("compressFeatureMatrix", signature("mRMRe.Data"), function(object, matrix)
{
    adaptor <- which(object@feature_types != 3)
    matrix <- matrix[adaptor, adaptor]
    colnames(matrix) <- object@feature_names
    rownames(matrix) <- object@feature_names
    
    return(matrix)
})

## expandFeatureIndices

setMethod("expandFeatureIndices", signature("mRMRe.Data"), function(object, indices)
{
    adaptor <- which(object@feature_types == 3)
    if(length(adaptor) != 0)
    {
        indices <- sapply(indices, function(i) i + sum(sapply(seq(adaptor), function(j) i >= (adaptor[[j]] - j + 1))))
    }

    return(indices)
})

## compressFeatureIndices

setMethod("compressFeatureIndices", signature("mRMRe.Data"), function(object, indices)
{
    adaptor <- which(object@feature_types == 3)
    if(length(adaptor) != 0)
    {
        indices <- sapply(indices, function(i) i - sum(i >= adaptor))
    }
    
    return(indices)
})