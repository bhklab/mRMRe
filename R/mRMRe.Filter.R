setClass("mRMRe.Filter", representation(solutions = "matrix", mi_matrix = "matrix", prior_weight = "numeric",
                target_index = "integer", levels = "integer", uses_ranks = "logical", outX = "logical",
                bootstrap_count = "integer"))

setMethod("initialize", "mRMRe.Filter",
        function(.Object, data, prior_weight, target_index, levels, uses_ranks, outX, bootstrap_count)
{
    if (class(data) != "mRMRe.Data")
        stop("data must be of type mRMRe.Data")
    
    if (prior_weight < 0 || prior_weight > 1)
        stop("prior weight must be a value ranging from 0 to 1")
    
    if (target_index < 1 || target_index > getFeatureCount(data))
        stop("target_index must be a value ranging from 1 to the amount of features in data")
        
    # Convert target_index
            
    if (missing(levels))
        stop("levels must be provided")
    
    if (missing(uses_ranks))
        uses_ranks <- TRUE
    
    if (missing(outX))
        outX <- TRUE
    
    if (missing(bootstrap_count))
        bootstrap_count <- 0
})