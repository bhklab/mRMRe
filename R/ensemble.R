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