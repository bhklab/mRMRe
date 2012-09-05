setClass("mRMRe.Data", representation(feature_names = "character", feature_types = "numeric", data = "matrix",
                strata = "numeric", weights = "numeric", priors = "matrix"))

setMethod("initialize", "mRMRe.Data", function(.Object, data, strata, weights, priors)
{
    if (!is.data.frame(data))
        stop("data must be of type data frame")
    
    feature_types <- sapply(data, function(feature) paste(class(feature), collapse = "_"))
    
    if (any(!is.element(feature_types, c("numeric", "ordered_factor", "Surv"))))
        stop("data columns must be either of numeric, ordered factor or Surv type")
    
    .Object@feature_names <- colnames(data)
    .Object@feature_types <- unlist(lapply(feature_types, switch, "Surv" = c(2, 3), "ordered_factor" = 1, 0))
    .Object@data <- do.call(cbind, lapply(seq(feature_types), function(i) switch(feature_types[[i]],
                                "Surv" = cbind(event = data[, i][, "status"], time = data[, i][, "time"]),
                                "ordered_factor" = as.numeric(as.integer(data[, i]) - 1),
                                as.numeric(data[, i]))))
    
    rownames(.Object@data) <- rownames(data)
    colnames(.Object@data)[!.Object@feature_types %in% c(2, 3)] <- colnames(data)[feature_types != "Surv"]
    colnames(.Object@data)[.Object@feature_types %in% c(2, 3)] <- paste(rep(colnames(data)[feature_types == "Surv"],
                    each = 2), rep(c("event", "time"), sum(feature_types == "Surv")), sep = "@@@")
    
    if (missing(weights)) 
        .Object@weights <- rep(1, nrow(data))
    else if (length(weights) == nrow(data))
        .Object@weights <- as.numeric(weights)
    else
        stop("data and weight must contain the same number of samples")

    if (missing(strata)) 
        .Object@strata <- rep.int(0, nrow(data))
    else if (is.factor(strata))
        .Object@strata <- as.integer(strata) - 1
    else if (length(strata) != nrow(data))
        stop("data and strata must contain the same number of samples")
    else
        stop("strata must be provided as factors")
    
    return(.Object)
})





