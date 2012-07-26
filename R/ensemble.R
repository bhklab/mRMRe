`build.mim` <- function(data_matrix, strata=NULL, weights=NULL, feature_types=NULL)
{
	if (is.null(strata))
		strata <- rep.int(0, nrow(data_matrix))
		
	if (is.null(weights))
		weights <- rep.int(1, nrow(data_matrix))
		
	if (is.null(feature_type))
		feature_type <- rep.int(0, ncol(data_matrix))	
	
    data_matrix <- as.matrix(data_matrix)
    mi_matrix <- .Call("build_mim", as.vector(data_matrix), as.vector(strata),
            as.vector(weights), as.vector(feature_types), as.integer(nrow(data_matrix)),
            as.integer(ncol(data_matrix)), as.integer(length(unique(strata))), "ensemble");
    mi_matrix <- matrix(mi_matrix, nrow=ncol(data_matrix), ncol=ncol(data_matrix))
    rownames(mi_matrix) <- colnames(data_matrix)
    colnames(mi_matrix) <- colnames(data_matrix)
    return(mi_matrix)
}