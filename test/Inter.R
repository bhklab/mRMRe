library(mRMRe)

load("~/Testbed/x03.RData")

methods <- c("SINGLEGENE", "RANKMULTIV", "mRMR", "ENSEMBLEmRMR")
drug <- drug_map["Irinotecan", ]
modes <- c("COMMON", "NEW")
folds <- 10
repetitions <- 100
levels <- c(200, rep(1, 29))
prefix <- "~/Testbed/"
limits <- c(0.5, 0.7)

get_predictions_per_method <- function(training_set, test_set)
    # Returns a list with format: obj[[method]] is a vector of predictions
{
    classic_tree <- mRMRe::filter.mRMR_tree(levels=rep(1, length(levels)), data=training_set, uses_ranks=FALSE, target_feature_index=1)
    ensemble_tree <- mRMRe::filter.mRMR_tree(levels=levels, data=training_set, uses_ranks=FALSE, target_feature_index=1)
    
    df_training_set <- as.data.frame(training_set)
    df_test_set <- as.data.frame(test_set)

    ranking <- apply(training_set[, -1, drop=FALSE], 2, cor, training_set[, 1, drop=TRUE], method="pearson", use="complete.obs")
    ranking_order <- order(abs(ranking), decreasing=TRUE)

    predictions_per_method <- lapply(methods, function(method)
    {
        if (method == "SINGLEGENE")
        { 
            gene <- names(ranking)[ranking_order[1]]
            formula <- as.formula(paste(colnames(training_set)[[1]], "~", gene, collapse=" + "))
            model <- lm(data=df_training_set, formula=formula)
            predictions <- predict(object=model, newdata=df_test_set, type="response")
        }
        else if (method == "RANKMULTIV")
        {
            genes <- names(ranking)[ranking_order[1:length(levels)]]
            formula <- as.formula(paste(colnames(training_set)[[1]], "~", paste(genes, collapse=" + ")))
            model <- lm(data=df_training_set, formula=formula)
            predictions <- predict(object=model, newdata=df_test_set, type="response")
        }
        else if (method == "mRMR")
        {
            formula <- as.formula(paste(colnames(training_set)[[1]], "~",
                            paste(sapply(classic_tree$paths[1, ], function(element) colnames(training_set)[element]), collapse=" + ")))
            model <- lm(data=df_training_set, formula=formula)
            predictions <- predict(object=model, newdata=df_test_set, type="response")
        }
        else if (method == "ENSEMBLEmRMR")
        {
            predictions <- apply(apply(ensemble_tree$paths, 1, function(path)
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
        
        training_set <- cbind(ic50_cgp[indices_cgp, drug["CGP"], drop=FALSE],
                data_cgp[indices_cgp, , drop=FALSE])
        training_set <- training_set[complete.cases(training_set), , drop=FALSE]
        colnames(training_set)[1] <- drug["CGP"]
        test_labels <- ic50_ccle[indices_ccle, drug["CCLE"], drop=FALSE]
        test_set <- data_ccle[indices_ccle, , drop=FALSE]

        message(mode)
        return(get_predictions_per_method(training_set=training_set, test_set=test_set))
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
        indices_cgp <- rownames(rownames(data_cgp))
        if (mode == "COMMON")
        {
            indices_ccle <- common_indices
        }
        else if (mode == "NEW")
        {
            #indices_cgp <- which(!(rownames(data_cgp) %in% common_indices))
            indices_ccle <- which(!(rownames(data_ccle) %in% common_indices))
        }
        
        test_labels <- ic50_ccle[indices_ccle, drug["CCLE"], drop=TRUE]
        score_per_method <- sapply(methods, function(method)
                    mRMRe::correlate(test_labels, obj_mephisto[[mode]][[method]], method="cindex"))
        names(score_per_method) <- methods

        pdf(paste(prefix, paste(limits, collapse="-"), "_", mode, "_BARPLOT.pdf", sep=""))
        col <- rainbow(length(score_per_method), s=0.5, v=0.9)
        barplot(score_per_method, col=col, space=c(0.25, 5), las=1, horiz=F,
                ylab="Concordance index", names.arg=names(score_per_method), ylim=limits, xpd=FALSE)
        #legend("topright", legend=names(score_per_method), col=col, bty="n", pch=15)
        dev.off()
        
        return(score_per_method)
    })
    names(scores_per_mode) <- modes
    return(scores_per_mode)
}

##
## CV-specific routines
##

run_andariel <- function()
{
    data <- cbind(ic50_cgp[, drug["CGP"], drop=FALSE], data_cgp)
    colnames(data)[1] <- drug["CGP"]
    
    predictions_per_repetition <- lapply(seq(repetitions), function(repetition) 
    {
        set.seed(repetition)
        order <- sample(nrow(data))
        message(paste("repetition ", repetition, sep=""))
        data <- data[order, , drop=FALSE]
        run_baal(data) # graph_baal(run_baal(data))
    })
    names(predictions_per_repetition) <- seq(repetitions)
    return(predictions_per_repetition)
}

graph_andariel <- function(obj_andariel)
{
    scores_per_repetition <- lapply(seq(repetitions), function(repetition) graph_baal(obj_andariel[[repetition]]))
    names(scores_per_repetition) <- seq(repetitions)
    return(scores_per_repetition)
}

work_andariel <- function(all_folds_avg)
{
    pdf(paste(prefix, paste(limits, collapse="-"), "_CV_BOXPLOT.pdf", sep=""))
    col <- rainbow(ncol(all_folds_avg), s=0.5, v=0.9)
    boxplot(all_folds_avg, col=col, space=c(0.25, 5), las=1, horiz=F,
            ylab="Concordance index", names.arg=colnames(all_folds_avg), ylim=limits, xpd=FALSE)
    #legend("topright", legend=names(score_per_method), col=col, bty="n", pch=15)
    dev.off()
    return(NULL)
}

run_baal <- function(data)
    # Returns a list with format: obj[[partition]][[method]] is a vector of predictions
{
    complete_cases <- complete.cases(data[, 1])
    partitions <- mapply(function(i, j) c(i, j), split(which(complete_cases), seq(folds)),
            split(which(!complete_cases), seq(folds)), SIMPLIFY=FALSE)
    
    predictions_per_partition <- lapply(seq(partitions), function(partition)
    {
        training_indices <- Reduce(x=partitions[-partition], f=c)
        training_set <- data[training_indices, , drop=FALSE]
        test_indices <- Reduce(x=partitions[partition], f=c)
        test_set <- data[test_indices, , drop=FALSE]

        message(paste("partition ", partition, sep=""))
        return(get_predictions_per_method(training_set=training_set, test_set=test_set))
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
        r <- mRMRe::correlate(labels, predictions_per_method[[method]], method="cindex")
        return(r)
    })
    names(scores_per_method_flat) <- methods
    
    # Inter scores
    #partitions <- seq(obj_baal)
    #scores_per_partition <- lapply(partitions, function(partition)
    #{
    #    labels <- ic50_cgp[names(obj_baal[[partition]][[methods[[1]]]]), drug["CGP"], drop=TRUE]
    #    scores_per_method <- lapply(methods, function(method)
    #                mRMRe::correlate(labels, obj_baal[[partition]][[method]], method="cindex"))
    #    names(scores_per_method) <- methods
    #    return(scores_per_method)
    #})
    #names(scores_per_partition) <- partitions
    
    #scores_per_method_inter <- lapply(methods, function(method)
    #{
    #    scores_per_partition_inter <- lapply(partitions, function(partition)
    #                scores_per_partition[[partition]][[method]])
    #    names(scores_per_partition_inter) <- partitions
    #    return(scores_per_partition_inter)
    #})
    #names(scores_per_method_inter) <- methods

    #p_values <- as.data.frame(sapply(methods, function(m1) sapply(methods, function(m2)
    #                                wilcox.test(as.numeric(scores_per_method_inter[[m1]]),
    #                                        as.numeric(scores_per_method_inter[[m2]]))[["p.value"]])))
    #rownames(p_values) <- methods
    #colnames(p_values) <- methods

    return(scores_per_method_flat)#list(scores=scores_per_method_flat, fold_score_per_method=scores_per_method_inter, p_values=p_values))
}

##
## Heatmap-specific routines
##

##
## C-INDEX vs Depth
##
