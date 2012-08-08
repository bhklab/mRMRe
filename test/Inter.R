library(ensemble)
library(multicore)
library(doParallel)
library(foreach)
library(glmnet)
registerDoParallel(4)
load("~/Testbed/x03.RData")

methods <- c("SINGLEGENE", "RANKMULTIV",  "mRMR", "ENSEMBLEmRMR")
drug <- drug_map["Irinotecan", ]
modes <- c("COMMON", "NEW")
levels <- c(2,2,2,2,2)

get_predictions <- function()
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
        df_training_set <- as.data.frame(training_set)
        df_test_set <- as.data.frame(test_set)
        
        tree <- ensemble::filter.mRMR_tree(levels=levels, data=training_set, uses_ranks=TRUE, target_feature_index=1)

        predictions_per_method <- lapply(methods, function(method)
        {
            if (method == "SINGLEGENE")
            {
                formula <- as.formula(paste(colnames(training_labels)[[1]], "~", colnames(training_set)[tree$paths[1, 1]], collapse=" + "))
                model <- lm(data=df_training_set, formula=formula)
                predictions <- predict(object=model, newdata=df_test_set, type="response")
            }
            else if (method == "RANKMULTIV")
            {
                unique_indices <- unique(tree$paths[ , 1])
                formula <- as.formula(paste(colnames(training_labels)[[1]], "~", paste(sapply(unique_indices, function(element) colnames(training_set)[element]), collapse=" + ")))
                model <- lm(data=df_training_set, formula=formula)
                predictions <- predict(object=model, newdata=df_test_set, type="response")
            }
            else if (method == "mRMR")
            {
                formula <- as.formula(paste(colnames(training_labels)[[1]], "~", paste(sapply(tree$paths[1, ], function(element) colnames(training_set)[element]), collapse=" + ")))
                model <- lm(data=df_training_set, formula=formula)
                predictions <- predict(object=model, newdata=df_test_set, type="response")
            }
            else if (method == "ENSEMBLEmRMR")
            {
                predictions <- apply(apply(tree$paths, 1, function(path)
                {
                    formula <- as.formula(paste(colnames(training_labels)[[1]], "~", paste(sapply(path, function(element) colnames(training_set)[element]), collapse=" + ")))
                    model <- lm(data=df_training_set, formula=formula)
                    prediction <- predict(object=model, newdata=df_test_set, type="response")
                }), 1, mean)
            }
            message(paste(drug[["CCLE"]], mode, method, sep="\t"))
            return(predictions)
        })
        names(predictions_per_method) <- methods
        return(predictions_per_method)
    })
    names(predictions_per_mode) <- modes
    return(predictions_per_mode)
}

generate_graph <- function(get_predictions_return)
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
        
        test_labels <- ic50_ccle[indices_ccle, drug["CCLE"], drop=FALSE]
        score_per_method <- sapply(methods, function(method) ensemble::correlate(test_labels, get_predictions_return[[mode]][[method]], method="cindex"))
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
