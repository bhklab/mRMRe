setClass("mRMRe.Data", representation(feature_names = "character", feature_types = "numeric", data = "matrix",
                strata = "numeric", weights = "numeric", priors = "matrix"))

setMethod("initialize", "mRMRe.Data", function(.Object, data, strata, weights, priors)
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
    else if (is.factor(strata))
        .Object@strata <- as.integer(strata) - 1
    else if (length(strata) != nrow(data))
        stop("data and strata must contain the same number of samples")
    else
        stop("strata must be provided as factors")
    
    ## Sample weight processing
    
    if (missing(weights)) 
        .Object@weights <- rep(1, nrow(data))
    else if (length(weights) == nrow(data))
        .Object@weights <- as.numeric(weights)
    else
        stop("data and weight must contain the same number of samples")

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

setMethod("getData", "mRMRe.Data", function(.Object)
{
    data <- lapply(seq(.Object@feature_types), function(i) switch(as.character(.Object@feature_types[[i]]),
                        "3" = Surv(time = .Object@data[, i], event = .Object@data[, i - 1]),
                        "2" = NULL,
                        "1" = .Object@data[, i] + 1,
                        "0" = .Object@data[, i],
                        NULL))
    data <- data.frame(data[!sapply(data, is.null)])
    colnames(data) <- .Object@feature_names
    
    return(data)
})

setMethod("getSampleCount", "mRMRe.Data", function(.Object)
{
    return(nrow(.Object@data))
})

setMethod("getFeatureCount", "mRMRe.Data", function(.Object)
{
    return(length(.Object@feature_names))
})

setMethod("getPriors", "mRMRe.Data", function(.Object)
{
    if (length(.Object@priors) == 0)
        return(.Object@priors)
    else
        return(compressFeatureMatrix(.Object, .Object@priors))
})

setMethod("expandFeatureMatrix", "mRMRe.Data", function(.Object, matrix)
{
    adaptor <- which(.Object@feature_types != 3)
    matrix <- do.call(cbind, lapply(seq(adaptor), function(i)
    {
        column <- do.call(rbind, lapply(seq(adaptor), function(j)
        {
            item <- matrix[j, i]
            
            if (.Object@feature_types[[adaptor[[j]]]] == 2)
                return(rbind(item, item, deparse.level = 0))
            else
                return(item)
        }))
        
        if (.Object@feature_types[[adaptor[[i]]]] == 2)
            return(cbind(column, column, deparse.level = 0))
        else
            return(column)
    }))

    return(matrix)
})

setMethod("compressFeatureMatrix", "mRMRe.Data", function(.Object, matrix)
{
    adaptor <- which(.Object@feature_types != 3)
    matrix <- matrix[adaptor, adaptor]
    colnames(matrix) <- .Object@feature_names
    rownames(matrix) <- .Object@feature_names
    
    return(matrix)
})

setMethod("expandFeatureIndices", "mRMRe.Data", function(.Object, indices)
{
    adaptor <- which(.Object@feature_types == 3)
    indices <- sapply(indices, function(i) i + sum(sapply(seq(adaptor), function(j) i >= (adaptor[[j]] - j + 1))))

    return(indices)
})

setMethod("compressFeatureIndices", "mRMRe.Data", function(.Object, indices)
{
    adaptor <- which(.Object@feature_types == 3)
    indices <- sapply(indices, function(i) i - sum(i >= adaptor))
    
    return(indices)
})