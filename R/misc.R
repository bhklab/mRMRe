`.expand.data` <- function(data) # mim
{
    feature_types <- unlist(lapply(data, function(column) class(column)[[1]]))

    new_feature_types <- unlist(lapply(feature_types, function(type)
    {
        if (type == "Surv")
            return(c(2, 3))
        else if (type == "factor")
            return(1)
        else # type == "numeric" or type == "integer"
            return(0)
    }))

    new_data <- do.call(cbind, lapply(data, function(type)
    {
        if (class(column)[[1]] == "Surv")
            return(cbind(event=column[, "status"], time=column[, "time"]))
        return(as.numeric(column))
    }))

    rownames(new_data) <- rownames(data)
    colnames(new_data)[!new_feature_types %in% c(2, 3)] <- colnames(data)[feature_types != "Surv"]
    colnames(new_data)[new_feature_types %in% c(2, 3)] <- paste(rep(colnames(data)[feature_types == "Surv"], each=2),
        rep(c("time", "event"), sum(feature_types == "Surv")), sep="_")

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

    return(list(data=new_data, feature_types=new_feature_types)) # new_mim=new_mim
}



`set.thread.count` <- function(thread_count)
{
    .Call(C_set_thread_count, as.integer(thread_count))
}
