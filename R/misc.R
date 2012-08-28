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
            subset <- which(paths >= index)
            offsets[subset] <<- offsets[subset] + 1
        })
        paths <- paths - offsets
    }
    
    return(list("mi_matrix"=mi_matrix, "paths"=paths))
}

`.expand.input` <- function(feature_types, data, priors) # mim
{
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
        rep(c("event", "time"), sum(feature_types == "Surv")), sep="@@@")

    new_priors <- do.call(cbind, lapply(seq(feature_types), function(i)
    {
        column <- do.call(rbind, lapply(seq(feature_types), function(j)
        {
            item <- priors[j, i]
            if (feature_types[[j]] == "Surv")
                return(rbind(item, item, deparse.level=0))
            else
                return(item)
        }))

        if (feature_types[[i]] == "Surv")
            return(cbind(column, column))
        else
            return(column)
    }))

    return(list("data"=new_data, "feature_types"=new_feature_types, "feature_names"=colnames(data),
                    "priors"=new_priors)) # new_mim=new_mim
}

`set.thread.count` <- function(thread_count)
{
    .Call(mRMRe:::.C_set_thread_count, as.integer(thread_count))
}
