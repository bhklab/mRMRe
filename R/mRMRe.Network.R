setClass("mRMRe.Network", representation(topologies = "list", prior_weight = "numeric", target_indices = "integer",
                levels = "integer"))

setMethod("initialize", "mRMRe.Network", function(.Object, data, prior_weight, target_indices, levels, layers, ...)
{
    if (missing(layers))
        layers <- 1L
    
    topologies <- list()
    length(topologies) <- getFeatureCount(data)
    
    lapply(seq(layers), function(layer)
    {
        target_indices <<- unlist(lapply(target_indices, function(target_index)
        {
            filter <- new("mRMRe.Filter", data = data, prior_weight = prior_weight, target_index = target_index,
                    levels = levels, ...)

            topologies[[target_index]] <<- getSolutions(filter)

            return(as.vector(getSolutions(filter)))
        }))

        target_indices <<- intersect(target_indices, which(sapply(topologies, is.null)))
    })

    .Object@topologies <- topologies

    return(.Object)
})