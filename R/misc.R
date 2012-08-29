`.expand.input` <- function(data, priors, prior_weights, strata, weights, uses_ranks, outX, bootstrap_count,
        target_indices)
{
    if (!is.data.frame(data))
        stop("data must be of type data frame")
    
    feature_types <- sapply(data, function(x) paste(class(x), collapse="_"))
    
    if (any(!is.element(feature_types, c("numeric", "ordered_factor", "Surv"))))
        stop("data columns must be either of numeric, ordered factor or Surv type")
    
    if (missing(weights)) 
        weights <- rep.int(1, nrow(data))
    
    if (missing(strata)) 
        strata <- rep.int(0, nrow(data))
    else if (is.factor(strata))
        strata <- as.integer(strata) - 1
    else
        stop("strata must be provided as factors")
    
    if (missing(priors))
    {
        priors <- vector()
        prior_weights <- 0
    }
    else if (missing(prior_weights))
        stop("prior_weights must be provided with priors")
    else if (prior_weights >= 1 || prior_weights <= 0)
        stop("prior_weights must be [0, 1]")
    else if (max(priors) > 1 || min(priors) < 0)
        stop("prior matrix elements must be [0, 1]")
    
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
    
    if (length(priors) != 0 && ncol(priors) == ncol(data) && nrow(priors) == ncol(data))
    {
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
    }
    else
        new_priors <- priors
    
    if (missing(uses_ranks))
        uses_ranks <- TRUE
    
    if (missing(outX))
        outX <- TRUE
    
    if (missing(bootstrap_count))
        bootstrap_count <- 0

    if (missing(target_indices))
        new_target_indices <- NULL
    else
        new_target_indices <- mRMRe:::.expand.feature.indices(feature_types, target_indices)
    
    return(list("data"=new_data, "priors"=new_priors, "prior_weights"=prior_weights, "strata"=strata, 
                    "weights"=weights, "feature_types"=new_feature_types, "uses_ranks"=uses_ranks, "outX"=outX,
                    "bootstrap_count"=bootstrap_count, "target_indices"=new_target_indices))
}

`.expand.feature.indices` <- function(feature_types, feature_indices)
{
    offsets <- rep(0, length(feature_indices))
    lapply(which(feature_types == "Surv"), function(index)
    {
        subset <- which(feature_indices > index)
        offsets[subset] <<- offsets[subset] + 1
    })
    feature_indices <- feature_indices + offsets
    return(feature_indices)
}

`.compress.feature.indices` <- function(feature_types, feature_indices)
{
    offsets <- rep(0, length(feature_indices))
    lapply(which(feature_types == 3), function(index)
    {
        subset <- which(feature_indices >= index)
        offsets[subset] <<- offsets[subset] + 1
    })
    feature_indices <- feature_indices - offsets
    return(feature_indices)
}

`.compress.output` <- function(feature_types, feature_names=NULL, mi_matrix=NULL, paths=NULL, target_indices=NULL)
{
    if (!is.null(mi_matrix))
    {
        sub <- feature_types == 3
        mi_matrix <- mi_matrix[!sub, !sub]
        rownames(mi_matrix) <- feature_names
        colnames(mi_matrix) <- feature_names
    }
    
    if (!is.null(paths))
        paths <- mRMRe:::.compress.feature.indices(feature_types, paths)
    
    if (!is.null(target_indices))
        target_indices <- mRMRe:::.compress.feature.indices(feature_types, target_indices)
    
    return(list("mi_matrix"=mi_matrix, "paths"=paths, "target_indices"=target_indices))
}

`set.thread.count` <- function(thread_count)
{
    .Call(mRMRe:::.C_set_thread_count, as.integer(thread_count))
}
