library(ensemble)

load("~/Testbed/x03.RData")

methods <- c("SINGLEGENE", "RANKMULTIV",  "mRMR", "ENSEMBLEmRMR")
drug <- drug_map["Irinotecan", ]
modes <- c("COMMON", "NEW")
folds <- 10
levels <- c(2,2,2,2,2)

get_predictions_per_method <- function(training_set, test_set, tree)
    # Returns a list with format: obj[[method]] is a vector of predictions
{
    df_training_set <- as.data.frame(training_set)
    df_test_set <- as.data.frame(test_set)
    
    predictions_per_method <- lapply(methods, function(method)
    {
        if (method == "SINGLEGENE")
        {
            formula <- as.formula(paste(colnames(training_set)[[1]], "~",
                            colnames(training_set)[tree$paths[1, 1]], collapse=" + "))
            model <- lm(data=df_training_set, formula=formula)
            predictions <- predict(object=model, newdata=df_test_set, type="response")
        }
        else if (method == "RANKMULTIV")
        {
            unique_indices <- unique(tree$paths[ , 1])
            formula <- as.formula(paste(colnames(training_set)[[1]], "~",
                            paste(sapply(unique_indices, function(element) colnames(training_set)[element]), collapse=" + ")))
            model <- lm(data=df_training_set, formula=formula)
            predictions <- predict(object=model, newdata=df_test_set, type="response")
        }
        else if (method == "mRMR")
        {
            formula <- as.formula(paste(colnames(training_set)[[1]], "~",
                            paste(sapply(tree$paths[1, ], function(element) colnames(training_set)[element]), collapse=" + ")))
            model <- lm(data=df_training_set, formula=formula)
            predictions <- predict(object=model, newdata=df_test_set, type="response")
        }
        else if (method == "ENSEMBLEmRMR")
        {
            predictions <- apply(apply(tree$paths, 1, function(path)
            {
                formula <- as.formula(paste(colnames(training_set)[[1]], "~",
                                paste(sapply(path, function(element) colnames(training_set)[element]), collapse=" + ")))
                model <- lm(data=df_training_set, formula=formula)
                prediction <- predict(object=model, newdata=df_test_set, type="response")
            }), 1, mean)
        }
        
        message(method)
        return(predictions)

    })

    names(predictions_per_method) <- methods
    return(predictions_per_method)
}

##
## Mode-specific routines (NEW, COMMON)
##

run_mephisto <- function()
    # Returns a list with format: obj[[mode]][[method]] is a vector of predictions
{
    predictions_per_mode <- lapply(modes, function(mode)
    {
        common_indices <- intersect(rownames(data_cgp), rownames(data_ccle))
        if (mode == "COMMON")
        {
            indices_cgp <- indices_ccle <- common_indices
        }
        else if (mode == "NEW")
        {
            indices_cgp <- which(!(rownames(data_cgp) %in% common_indices))
            indices_ccle <- which(!(rownames(data_ccle) %in% common_indices))
        }
        
        training_labels <- ic50_cgp[indices_cgp, drug["CGP"], drop=FALSE]
        training_labels_complete_cases <- complete.cases(training_labels)
        training_set <- cbind(training_labels,
                data_cgp[indices_cgp, , drop=FALSE])[training_labels_complete_cases, , drop=FALSE]
        colnames(training_set)[1] <- drug["CGP"]
        test_labels <- ic50_ccle[indices_ccle, drug["CCLE"], drop=FALSE]
        test_set <- data_ccle[indices_ccle, , drop=FALSE]
        
        tree <- ensemble::filter.mRMR_tree(levels=levels, data=training_set, uses_ranks=TRUE, target_feature_index=1)

        message(mode)
        return(get_predictions_per_method(training_set=training_set, test_set=test_set, tree=tree))
    })
    names(predictions_per_mode) <- modes
    return(predictions_per_mode)
}

graph_mephisto <- function(obj_mephisto)
    # Returns a list with format: obj[[mode]][[method]] is a correlation
{
    scores_per_mode <- lapply(modes, function(mode)
    {
        common_indices <- intersect(rownames(data_cgp), rownames(data_ccle))
        if (mode == "COMMON")
        {
            indices_cgp <- indices_ccle <- common_indices
        }
        else if (mode == "NEW")
        {
            indices_cgp <- which(!(rownames(data_cgp) %in% common_indices))
            indices_ccle <- which(!(rownames(data_ccle) %in% common_indices))
        }
        
        test_labels <- ic50_ccle[indices_ccle, drug["CCLE"], drop=TRUE]
        score_per_method <- sapply(methods, function(method)
                    ensemble::correlate(test_labels, obj_mephisto[[mode]][[method]], method="cindex"))
        names(score_per_method) <- methods

        pdf(paste("~/Testbed/Inter_", mode, ".pdf", sep=""))
        col <- rainbow(length(score_per_method), s=0.5, v=0.9)
        barplot(score_per_method, col=col, space=c(0.25, 5), las=1, horiz=F,
                ylab="concordance index", names.arg=drug["CCLE"])
        legend("topright", legend=names(score_per_method), col=col, bty="n", pch=15)
        dev.off()
        
        return(score_per_method)
    })
    names(scores_per_mode) <- modes
    return(scores_per_mode)
}

##
## CV-specific routines
##

run_baal <- function()
    # Returns a list with format: obj[[partition]][[method]] is a vector of predictions
{
    data <- cbind(ic50_cgp[, drug["CGP"], drop=FALSE], data_cgp)
    colnames(data)[1] <- drug["CGP"]
    complete_cases <- complete.cases(data[, 1])
    partitions <- mapply(function(i, j) c(i, j), split(which(complete_cases), seq(folds)),
            split(which(!complete_cases), seq(folds)), SIMPLIFY=FALSE)
    
    predictions_per_partition <- lapply(seq(partitions), function(partition)
    {
        training_indices <- Reduce(x=partitions[-partition], f=c)
        training_set <- data[training_indices, , drop=FALSE]
        test_indices <- Reduce(x=partitions[partition], f=c)
        test_set <- data[test_indices, , drop=FALSE]
        
        tree <- ensemble::filter.mRMR_tree(levels=levels, data=training_set, uses_ranks=TRUE, target_feature_index=1)
        
        message(partition)
        return(get_predictions_per_method(training_set=training_set, test_set=test_set, tree=tree))
    })
    names(predictions_per_partition) <- seq(partitions)
    return(predictions_per_partition)
}

graph_baal <- function(obj_baal)
    # Returns a list with format: obj[["scores"]][[method]] is a correlation
    #           				   obj[["p_values"]] is a matrix of pairwise per-fold wilcox tests
{
    # Flat scores
    predictions_per_method <- Reduce(f=rbind, x=lapply(obj_baal, as.data.frame))
    predictions_per_method <- as.list(predictions_per_method[rownames(data_cgp), , drop=FALSE])

    scores_per_method_flat <- lapply(methods, function(method)
    {
        labels <- ic50_cgp[, drug["CGP"], drop=TRUE]
        r <- ensemble::correlate(labels, predictions_per_method[[method]], method="cindex")
        return(r)
    })
    names(scores_per_method_flat) <- methods
    
    # Inter scores
    partitions <- seq(obj_baal)
    scores_per_partition <- lapply(partitions, function(partition)
    {
        labels <- ic50_cgp[names(obj_baal[[partition]][[methods[[1]]]]), drug["CGP"], drop=TRUE]
        scores_per_method <- lapply(methods, function(method)
                    ensemble::correlate(labels, obj_baal[[partition]][[method]], method="cindex"))
        names(scores_per_method) <- methods
        return(scores_per_method)
    })
    names(scores_per_partition) <- partitions
    
    scores_per_method_inter <- lapply(methods, function(method)
    {
        scores_per_partition_inter <- lapply(partitions, function(partition)
                    scores_per_partition[[partition]][[method]])
        names(scores_per_partition_inter) <- partitions
        return(scores_per_partition_inter)
    })
    names(scores_per_method_inter) <- methods

    p_values <- as.data.frame(sapply(methods, function(m1) sapply(methods, function(m2)
                                    wilcox.test(as.numeric(scores_per_method_inter[[m1]]),
                                            as.numeric(scores_per_method_inter[[m2]]))[["p.value"]])))
    rownames(p_values) <- methods
    colnames(p_values) <- methods

    return(list(scores=scores_per_method_flat, p_values=p_values))
}

##
## Heatmap-specific routines
##

##
## C-INDEX vs Depth
##
