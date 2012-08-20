# feature_type = 0 for continous
# feature_type = 1 for discrete
# feature_type = 2 for survival
`build.mim` <- function(
        data,
        strata=rep.int(0, nrow(data)),
        weights=rep.int(1, nrow(data)),
        uses_ranks=TRUE,
        outX=TRUE,
        bootstrap_count=0)
{
    if (is.data.frame(data))
        stop("data must be of type data frame")
    
    expansion <- mRMRe:::.expand.data(data)
    data <- expansion$data
    feature_types <- expansion$feature_types
    
    data <- as.matrix(data)
    mi_matrix <- .Call(C_build_mim, as.vector(data), as.vector(strata), as.vector(weights), as.vector(feature_types),
            as.integer(nrow(data)), as.integer(ncol(data)), as.integer(length(unique(strata))), as.integer(uses_ranks),
            as.integer(outX), as.integer(bootstrap_count))
    mi_matrix <- matrix(mi_matrix, nrow=ncol(data), ncol=ncol(data))
    rownames(mi_matrix) <- colnames(data)
    colnames(mi_matrix) <- colnames(data)
    
    return(mi_matrix)
}

# For c-index, it is assumed that X is discrete and that Y is continuous
`correlate` <- function(
        x,
        y,
        strata=rep.int(0, length(x)),
        weights=rep.int(1, length(x)),
        method=c("cramer", "pearson", "spearman", "cindex"),
        outX=TRUE,
        bootstrap_count=0)
{
    x <- as.vector(x)
    y <- as.vector(y)
    strata <- as.vector(strata)
    weights <- as.vector(weights)
    method <- match.arg(method)
    stratum_count <- as.integer(length(unique(strata)))
    bootstrap_count <- as.integer(bootstrap_count)
    value <- NA
    
    if (method == "cramer")
        value <- .Call(C_compute_cramers_v, x, y, weights, strata, stratum_count, bootstrap_count)
    else if (method == "pearson")
        value <- .Call(C_compute_pearson_correlation, x, y, weights, strata, stratum_count, bootstrap_count)
    else if (method == "spearman")
        value <- .Call(C_compute_spearman_correlation, x, y, weights, strata, stratum_count, bootstrap_count)
    else if (method == "cindex")
        value <- .Call(C_compute_concordance_index, x, y, weights, strata, stratum_count, outX)
    
    return(value)
}
