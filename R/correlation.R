`build.mim` <- function(data, priors, prior_weights, strata, weights, uses_ranks=TRUE, outX=TRUE, bootstrap_count=0)
{
    if (!is.data.frame(data))
        stop("data must be of type data frame")
    
    feature_types <- sapply(data, function(x) paste(class(x), collapse="_"))
    
    if (any(!is.element(feature_types, c("numeric", "ordered_factor", "Surv"))))
        stop("data columns must be either of numeric, ordered factor or Surv type")
    
    if (missing(weights)) 
        weights <- rep.int(1, nrow(data))
    
    if (missing(strata)) 
        strata <- rep.int(0, nrow(data))
    else if (is.factor(strata))
        strata <- as.integer(strata) - 1
    else
        stop("strata must be provided as factors")
               
    if (missing(priors))
    {
        priors <- vector()
        prior_weights <- 0
    }
    else if (missing(prior_weights))
        stop("prior_weights must be provided with priors")
    else if (prior_weights > 1 || prior_weights < 0)
        stop("prior_weights must be [0, 1]")
	else if (max(priors) > 1 || min(priors) < 0)
		stop("prior must be [0,1]")
    
    expansion <- mRMRe:::.expand.input(feature_types=feature_types, data=data, priors=priors)
    data <- expansion$data
    feature_types <- expansion$feature_types
    feature_names <- expansion$feature_names
    priors <- expansion$priors
    
    data <- as.matrix(data)
    mi_matrix <- .Call(mRMRe:::.C_build_mim, as.vector(data), as.vector(priors), as.numeric(lambda), as.vector(strata),
            as.vector(weights), as.vector(feature_types), as.integer(nrow(data)), as.integer(ncol(data)),
            as.integer(length(unique(strata))), as.integer(uses_ranks), as.integer(outX), as.integer(bootstrap_count))
    mi_matrix <- matrix(mi_matrix, nrow=ncol(data), ncol=ncol(data))
    
    compression <- mRMRe:::.compress.output(feature_types=feature_types, feature_names=feature_names, mi_matrix=mi_matrix)
    mi_matrix <- -0.5 * log(1 - (compression$mi_matrix^2))
    
    return(mi_matrix)
}

`correlate` <- function(x, y, strata, weights, method=c("cramer", "pearson", "spearman", "cindex", "kendall"),
        outX=TRUE, bootstrap_count=0)
{
    method <- match.arg(method)
    
    type_x <- paste(class(x), collapse="_")
    type_y <- paste(class(y), collapse="_")
    
    if (type_x == "ordered_factor")
        x <- as.integer(x) - 1
    if (type_y == "ordered_factor")
        y <- as.integer(y) - 1
    
    is_survival <- FALSE
    
    if (type_x == "Surv" || type_y == "Surv")
    {
        is_survival <- TRUE
        
        if (method != "cindex")
            stop("survival data can only be used with the concordance index")
        else
        {
            if (type_y == "Surv")
            {
                t <- x
                x <- y
                y <- t
                t <- type_x
                type_x <- type_y
                type_y <- t
            }
            
            if (type_y != "numeric" && type_y != "integer")
                stop("survival data must be compared to numerical data")
        }
        
        x[, "status"] <- as.numeric(x[, "status"])
        x[, "time"] <- as.numeric(x[, "time"])
        y <- as.numeric(y)
        
        if (nrow(x) != length(y))
            stop("both sample sets must have the same length")
    }
    else
    {
        x <- as.numeric(x)
        y <- as.numeric(y)
        
        if (length(x) != length(y))
            stop("both sample sets must have the same length")
    }
    
    if (missing(weights)) 
        weights <- rep.int(1, length(y))
    
    if (missing(strata)) 
        strata <- rep.int(0, length(y))
    else if (is.factor(strata))
        strata <- as.integer(strata) - 1
    else
        stop("strata must be provided as factors")
    
    weights <- as.vector(weights)
    stratum_count <- as.integer(length(unique(strata)))
    bootstrap_count <- as.integer(bootstrap_count)
    value <- NA
        
    if (method == "cramer")
        value <- .Call(mRMRe:::.C_compute_cramers_v, x, y, weights, strata, stratum_count, bootstrap_count)
    else if (method == "pearson")
        value <- .Call(mRMRe:::.C_compute_pearson_correlation, x, y, weights, strata, stratum_count, bootstrap_count)
    else if (method == "spearman")
        value <- .Call(mRMRe:::.C_compute_spearman_correlation, x, y, weights, strata, stratum_count, bootstrap_count)
    else if (method == "cindex")
    {
        if (!is_survival)
            value <- .Call(mRMRe:::.C_compute_concordance_index, x, y, weights, strata, stratum_count, outX)
        else
            value <- .Call(mRMRe:::.C_compute_concordance_index_with_time, x[, "status"], y, x[, "time"], weights, strata, stratum_count, outX)
    }
    else if (method == "kendall")
        value <- (.Call(mRMRe:::.C_compute_concordance_index, x, y, weights, strata, stratum_count, outX)$statistic - 0.5) * 2
    
    return(value)
}
