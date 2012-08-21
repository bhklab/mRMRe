`.compress.output` <- function(feature_types, feature_names=NULL, mi_matrix=NULL, paths=NULL)
{
    if (!is.null(mi_matrix))
    {
        sub <- feature_types == 3
        mi_matrix <- mi_matrix[!sub, !sub]
        rownames(mi_matrix) <- feature_names
        colnames(mi_matrix) <- feature_names
    }
    
    if (!is.null(paths))
    {
        offsets <- rep(0, length(paths))
        lapply(which(feature_types == 3), function(index)
        {
            subset <- which(paths > index)
            offsets[subset] <<- offsets[subset] + 1
        })
        paths <- paths - offsets
    }
    
    return(list("mi_matrix"=mi_matrix, "paths"=paths))
}

`.expand.input` <- function(feature_types=NULL, data) # mim
{
    if (is.null(feature_types))
        feature_types <- sapply(data, function(x) paste(class(x), collapse="_"))
    
    new_feature_types <- unlist(lapply(feature_types, function(type)
    {
        if (type == "Surv")
            return(c(2, 3))
        else if (type == "ordered_factor")
            return(1)
        else # type == "numeric", type == "integer" or type == "double"
            return(0)  
    }))

    new_data <- do.call(cbind, lapply(data, function(column)
    {
        type <- paste(class(column), collapse="_")
        
        if (type == "Surv")
            return(cbind(event=column[, "status"], time=column[, "time"]))
        else if (type == "ordered_factor")
            return(as.numeric(as.integer(column) - 1))
        else
            return(as.numeric(column))
    }))
    rownames(new_data) <- rownames(data)
    colnames(new_data)[!new_feature_types %in% c(2, 3)] <- colnames(data)[feature_types != "Surv"]
    colnames(new_data)[new_feature_types %in% c(2, 3)] <- paste(rep(colnames(data)[feature_types == "Surv"], each=2),
        rep(c("time", "event"), sum(feature_types == "Surv")), sep="@@@")

    # new_mim <- sapply(seq(feature_types), function(i)
    # {
    #     a <- sapply(seq(feature_types), function(j)
    #     {
    #         b <- mim[i, j]
    #         if (feature_types[[j]] == "Surv")
    #             return(c(b, NA))
    #         return(b)
    #     })
    #     if (feature_types[[i]] == "Surv")
    #         return(cbind(a, rep(NA, length(a))
    #     return(a)
    # })

    return(list("data"=new_data, "feature_types"=new_feature_types, "feature_names"=colnames(data))) # new_mim=new_mim
}

`set.thread.count` <- function(thread_count)
{
    .Call(mRMRe:::.C_set_thread_count, as.integer(thread_count))
}