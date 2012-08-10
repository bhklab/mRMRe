library(ensemble)

load("~/Testbed/x03.RData")

methods <- c("SINGLEGENE", "RANKMULTIV",  "mRMR", "ENSEMBLEmRMR")
drug <- drug_map["Irinotecan", ]
modes <- c("COMMON", "NEW")
levels <- c(2,2,2,2,2)

get_predictions_per_method <- function(training_set, test_set, tree) # Returns a list with format: obj[[method]] is a vector of predictions
{
    df_training_set <- as.data.frame(training_set)
    df_test_set <- as.data.frame(test_set)
    ranking <- apply(training_set[,-1], 2, cor, training_set[,1], method="pearson", use="complete.obs")
	
	
    predictions_per_method <- lapply(methods, function(method)
    {
        if (method == "SINGLEGENE")
        {	
			gene <- names(ranking)[order(abs(ranking), decreasing=TRUE)[1]]
            formula <- as.formula(paste(colnames(training_set)[[1]], "~", colnames(training_set)[tree$paths[1, 1]], collapse=" + "))
            model <- lm(data=df_training_set, formula=formula)
            predictions <- predict(object=model, newdata=df_test_set, type="response")
        }
        else if (method == "RANKMULTIV")
        {
			genes <- names(ranking)[order(abs(ranking_scc), decreasing=TRUE)[1:length(levels)]]
            unique_indices <- unique(tree$paths[ , 1])
            formula <- as.formula(paste(colnames(training_set)[[1]], "~", paste(sapply(genes, function(element) colnames(training_set)[element]), collapse=" + ")))
            model <- lm(data=df_training_set, formula=formula)
            predictions <- predict(object=model, newdata=df_test_set, type="response")
        }
        else if (method == "mRMR")
        {
            formula <- as.formula(paste(colnames(training_set)[[1]], "~", paste(sapply(tree$paths[1, ], function(element) colnames(training_set)[element]), collapse=" + ")))
            model <- lm(data=df_training_set, formula=formula)
            predictions <- predict(object=model, newdata=df_test_set, type="response")
        }
        else if (method == "ENSEMBLEmRMR")
        {
            predictions <- apply(apply(tree$paths, 1, function(path)
            {
                formula <- as.formula(paste(colnames(training_set)[[1]], "~", paste(sapply(path, function(element) colnames(training_set)[element]), collapse=" + ")))
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

run_mephisto <- function() # Returns a list with format: obj[[mode]][[method]] is a vector of predictions
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
        training_set <- cbind(training_labels, data_cgp[indices_cgp, , drop=FALSE])[training_labels_complete_cases, , drop=FALSE]
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

graph_mephisto <- function(obj_mephisto) # Returns a list with format: obj[[mode]][[method]] is a correlation
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
        score_per_method <- sapply(methods, function(method) ensemble::correlate(test_labels, obj_mephisto[[mode]][[method]], method="cindex"))
        names(score_per_method) <- methods

        pdf(paste("~/Testbed/Inter_", mode, ".pdf", sep=""))
        col <- rainbow(length(score_per_method), s=0.5, v=0.9)
        barplot(score_per_method, col=col, space=c(0.25, 5), las=1, horiz=F, ylab="concordance index", names.arg=drug["CCLE"])
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

run_baal <- function(folds=10) # Returns a list with format: obj[[method]] is a vector of predictions
{
    data <- cbind(ic50_cgp[, drug["CGP"], drop=FALSE], data_cgp)
    colnames(data)[1] <- drug["CGP"]
    complete_cases <- complete.cases(data[, 1])
    partitions <- mapply(function(i, j) c(i, j), split(which(complete_cases), seq(folds)), split(which(!complete_cases), seq(folds)), SIMPLIFY=FALSE)
    
    predictions <- Reduce(f=rbind, x=lapply(seq(partitions), function(partition)
    {
        training_indices <- Reduce(x=partitions[-partition], f=c)
        training_set <- data[training_indices, , drop=FALSE]
        test_indices <- Reduce(x=partitions[partition], f=c)
        test_set <- data[test_indices, , drop=FALSE]
        
        tree <- ensemble::filter.mRMR_tree(levels=levels, data=training_set, uses_ranks=TRUE, target_feature_index=1)
		browser()
        
        message(partition)
        return(as.data.frame(get_predictions_per_method(training_set=training_set, test_set=test_set, tree=tree)))
    }))
    predictions <- predictions[rownames(data_cgp), , drop=FALSE]
    return(as.list(predictions))
}

graph_baal <- function(obj_baal) # Returns a list with format: obj[[method]] is a correlation
{
    scores_per_method <- lapply(methods, function(method)
    {
        labels <- ic50_cgp[, drug["CGP"], drop=TRUE]
        ensemble::correlate(labels, obj_baal[[method]], method="cindex")
    })
    names(scores_per_method) <- methods
    return(scores_per_method)
}

# TODO:
# - Perhaps repeat these CVs 10 times ? A routine should be in place for that
# - graph_demonname takes run_demonname's return value as argument (graph_baal(run_baal(...)) -- provides some sort of caching so that we do not regenerate data if we change graphical schemes
# - A method for gathering correlations & graphs from this... however way we want it
# - The below routines

##
## Heatmap-specific routines
##

##
## C-INDEX vs Depth
##
