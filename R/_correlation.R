`build.mim` <- function(data, priors, prior_weight, strata, weights, uses_ranks, outX, bootstrap_count, as_mi=TRUE)
{
    feature_names <- colnames(data)
    expansion <- mRMRe:::.expand.input(data=data, priors=priors, prior_weight=prior_weight, strata=strata,
            weights=weights, uses_ranks=uses_ranks, outX=outX, bootstrap_count=bootstrap_count)
    data <- expansion$data
    priors <- expansion$priors
    prior_weight <- expansion$prior_weight
    strata <- expansion$strata
    weights <- expansion$weights
    feature_types <- expansion$feature_types
    uses_ranks <- expansion$uses_ranks
    outX <- expansion$outX
    bootstrap_count <- expansion$bootstrap_count
    
    data <- as.matrix(data)
    mi_matrix <- .Call(mRMRe:::.C_export_mim, as.vector(data), as.vector(priors), as.numeric(prior_weight),
            as.vector(strata), as.vector(weights), as.vector(feature_types), as.integer(nrow(data)),
            as.integer(ncol(data)), as.integer(length(unique(strata))), as.integer(uses_ranks), as.integer(outX),
            as.integer(bootstrap_count))
    mi_matrix <- matrix(mi_matrix, nrow=ncol(data), ncol=ncol(data))
    
    compression <- mRMRe:::.compress.output(feature_types=feature_types, feature_names=feature_names, mi_matrix=mi_matrix)
    mi_matrix <- compression$mi_matrix
    
    if (as_mi)
        mi_matrix <- -0.5 * log(1 - (compression$mi_matrix^2))

    return(mi_matrix)
}

`correlate` <- function(x, y, strata, weights, method=c("cramer", "pearson", "spearman", "cindex"),
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
    
    z <- vector()
    if (method == "cindex" && is_survival)
    {
        z <- x[, "time"]
        y <- y
        x <- x[, "status"]
    }
    
    value <- .Call(mRMRe:::.C_export_association, x, y, z, strata, weights, stratum_count, outX, bootstrap_count, method)
    
    return(value)
}
