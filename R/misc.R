# feature_type = 0 for continous
# feature_type = 1 for discrete
# feature_type = 2 for survival
`expand.data` <- function(data, feature_types)
{
    offset <- 0
    survival_feature_indices <- which(feature_types == 2)
    
    lapply(survival_feature_indices, function(index)
    {
        real_index <- index + offset
        decomposed_survival_features <- data[, real_index][, c("status", "time")]
        colnames(decomposed_survival_features) <- paste(colnames(data)[[real_index]], c("event", "time"))
        
        if (real_index == 1)
            data <<- cbind(decomposed_survival_features, data[, (real_index + 1):ncol(data), drop=FALSE])
        else
            data <<- cbind(data[, 1:(real_index - 1), drop=FALSE], decomposed_survival_features, data[, (real_index + 1):ncol(data), drop=FALSE])
        
        feature_types <<- c(feature_types[1:real_index], 3, feature_types[(real_index + 1):length(feature_types)]) 
        offset <<- offset + 1
    })

    return(list(data=data, feature_types=feature_types))
}

`set.thread.count` <- function(thread_count)
{
    .Call(C_set_thread_count, as.integer(thread_count))
}