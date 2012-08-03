`build.mim` <- function(data, strata=NULL, weights=NULL, feature_types=NULL)
{
    data <- as.matrix(data)
    
	if (is.null(strata))
		strata <- rep.int(0, nrow(data))
		
	if (is.null(weights))
		weights <- rep.int(1, nrow(data))
		
	if (is.null(feature_types))
		feature_types <- rep.int(0, ncol(data))	
	
    mi_matrix <- .Call(C_build_mim, as.vector(data), as.vector(strata), as.vector(weights), as.vector(feature_types),
            as.integer(nrow(data)), as.integer(ncol(data)), as.integer(length(unique(strata))))
    mi_matrix <- matrix(mi_matrix, nrow=ncol(data), ncol=ncol(data))
    rownames(mi_matrix) <- colnames(data)
    colnames(mi_matrix) <- colnames(data)
    
    return(mi_matrix)
}

`correlate` <- function(x, y, strata=NULL, weights=NULL, method=c("c-index" ,"cramer", "pearson", "spearman"))
{
    x <- as.vector(x)
    y <- as.vector(y)
    
    if (is.null(strata))
        strata <- rep.int(0, length(x))
    
    if (is.null(weights))
        weights <- rep.int(1, length(x))
    
    strata <- as.vector(strata)
    weights <- as.vector(weights)
    method <- match.arg(method)
    
    value <- NA
    stratum_count <- as.integer(length(unique(strata)))
    
    if (method == "cramer")
        value <- .Call(C_compute_cramers_v, x, y, weights, strata, stratum_count)
    else if (method == "pearson")
        value <- .Call(C_compute_pearson_correlation, x, y, weights, strata, stratum_count)
    else if (method == "spearman")
        value <- .Call(C_compute_spearman_correlation, x, y, weights, strata, stratum_count)
	else if (method == "c-index")
		value <- .Call(C_compute_concordance_index, x, y, weights, strata, stratum_count)
	else
		message("Invalid method specified")
    
    return(value)
}