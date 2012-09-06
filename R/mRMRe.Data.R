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
            .Object@priors <- do.call(cbind, lapply(seq(feature_types), function(i)
            {
                column <- do.call(rbind, lapply(seq(feature_types), function(j)
                {
                    item <- priors[j, i]

                    if (feature_types[[j]] == "Surv")
                        return(rbind(item, item, deparse.level = 0))
                    else
                        return(item)
                }))
                
                if (feature_types[[i]] == "Surv")
                    return(cbind(column, column, deparse.level = 0))
                else
                    return(column)
            }))
    }
    
    return(.Object)
})

setMethod("getSampleCount", "mRMRe.Data", function(.Object) nrow(.Object@data))

setMethod("getFeatureCount", "mRMRe.Data", function(.Object) length(.Object@feature_names))

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

setMethod("getPriors", "mRMRe.Data", function(.Object)
{
    indices <- which(.Object@feature_types != 3)
    priors <- .Object@priors[indices, indices]
    colnames(priors) <- .Object@feature_names
    rownames(priors) <- .Object@feature_names
    
    return(priors)
})