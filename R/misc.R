## example of data
#library(survival)
#dd <- data.frame("surv1"=Surv(runif(100), sample(0:1, 100, replace=T)), "cont1"=runif(100), "cat1"=factor(sample(1:5, 100, replace=T), ordered=T), "surv2"=Surv(runif(100), sample(0:1, 100, replace=T)), "cont2"=runif(100), "cont3"=runif(100), "surv3"=Surv(runif(100), sample(0:1, 100, replace=T)), "cat2"=factor(sample(1:5, 100, replace=T), ordered=T))
#head(dd)
`.expand.data` <- function(data)
{
    ## expand features
    feature_types <- sapply(data, class)
    if(is.list(feature_types)) { feature_types <- unlist(lapply(feature_types, function(x) { paste(x, collapse="_") })) }
    ff <- unlist(sapply(feature_types, function(x) { if(x == "Surv") { return(c("Surv_time", "Surv_event")) } else { return(x) } }))
    ## expand data
    dd <- sapply(data, function(x) {
        if(class(x)[1] != "Surv") { return(as.numeric(x)) } else {
            ## split survival data
            x <- cbind("time"=x[ ,"time"], "event"=x[ ,"status"])
            return(x)
        }
    })
    dd <- do.call(cbind, dd)
    rownames(dd) <- rownames(data)
    ## update colnames
    colnames(dd)[!is.element(ff, c("Surv_time", "Surv_event"))] <- colnames(data)[feature_types != "Surv"]
    colnames(dd)[is.element(ff, c("Surv_time", "Surv_event"))] <- paste(rep(colnames(data)[feature_types == "Surv"], each=2), rep(c("time", "event"), sum(feature_types == "Surv")), sep="_")
    return(list("data"=dd, "feature_types"=ff))
}

`.expand.data.old` <- function(data, feature_types)
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