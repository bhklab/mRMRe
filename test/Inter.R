library(ensemble)
library(multicore)
library(doParallel)
library(foreach)
library(glmnet)
registerDoParallel(4)
load("~/Testbed/x03.RData")

drugs <- c("Lapatinib", "PD-0325901", "Irinotecan")
methods <- c("SINGLEGENE", "RANKMULTIV", "RANKENSEMBLE", "mRMR", "ENSEMBLEmRMR", "ELASTICNET")
drug_map <- drug_map[drugs, , drop=FALSE]
common_indices <- intersect(rownames(data_cgp), rownames(data_ccle))
training_data <- data_cgp
test_data <- data_ccle
metric <- apply(drug_map, 1, function(drug)
{
    training_labels <- ic50_cgp[, drug[[2]], drop=FALSE]
    indices <- complete.cases(training_labels)
    test_labels <- ic50_ccle[common_indices, drug[[1]], drop=FALSE]
    training_set <- cbind(training_labels, training_data)[indices, , drop=FALSE]
    colnames(training_set)[1] <- drug[[2]]
    tree <- ensemble::filter.mRMR_tree(levels=c(2,2,2,2,2), data=training_set, uses_ranks=TRUE, target_feature_index=1)
    predictions <- mclapply(methods, mc.preschedule=TRUE, mc.cores=length(methods), mc.cleanup=TRUE, function(method)
    {
        p <- NULL
        if (method == "SINGLEGENE")
        {
            formula <- as.formula(paste(colnames(training_labels)[[1]], "~", colnames(training_data)[tree$paths[1, 1] - 1], collapse=" + "))
            model <- lm(data=as.data.frame(training_set), formula=formula)
            p <- as.vector(predict(object=model, newdata=as.data.frame(test_data), type="response"))
        }
        else if (method == "RANKMULTIV")
        {
            unique_indices <- unique(tree$paths[ , 1])
            formula <- as.formula(paste(colnames(training_labels)[[1]], "~", paste(sapply(unique_indices, function(element) colnames(training_data)[element - 1]), collapse=" + ")))
            model <- lm(data=as.data.frame(training_set), formula=formula)
            p <- as.vector(predict(object=model, newdata=as.data.frame(test_data), type="response"))
        }
        else if (method == "RANKENSEMBLE")
        {
            unique_indices <- unique(tree$paths[ , 1])
            p <- apply(sapply(unique_indices, function(element)
            {
                formula <- as.formula(paste(colnames(training_labels)[[1]], "~", colnames(training_data)[element - 1], collapse= " + "))
                model <- lm(data=as.data.frame(training_set), formula=formula)
                prediction <- as.vector(predict(object=model, newdata=as.data.frame(test_data), type="response"))
            }), 1, mean)
        }
        else if (method == "mRMR")
        {
            formula <- as.formula(paste(colnames(training_labels)[[1]], "~", paste(sapply(tree$paths[1, ], function(element) colnames(training_data)[element - 1]), collapse=" + ")))
            model <- lm(data=as.data.frame(training_set), formula=formula)
            p <- as.vector(predict(object=model, newdata=as.data.frame(test_data), type="response"))
        }
        else if (method == "ENSEMBLEmRMR")
        {
            p <- apply(apply(tree$paths, 1, function(path)
            {
                formula <- as.formula(paste(colnames(training_labels)[[1]], "~",
                    paste(sapply(path, function(element) colnames(training_data)[element - 1]), collapse=" + ")))
                    model <- lm(data=as.data.frame(training_set), formula=formula)
                    prediction <- as.vector(predict(object=model, newdata=as.data.frame(test_data), type="response"))
            }), 1, mean)
        }
        else if (method == "ELASTICNET")
        {
            alpha_to_test <- seq(0.2, 1.0, length.out=10)
            lambda_to_test <- exp(rev(seq(from=-6, to=5, length.out=250)))
            cvmses <- foreach(i=1:length(alpha_to_test), .combine="cbind", .packages="glmnet") %dopar%
            {
                rr <- glmnet::cv.glmnet(x=training_set[,-1], y=training_set[,1], alpha=alpha_to_test[i],lambda=lambda_to_test)$cvm
            }
            dimnames(cvmses) <- list(paste("lambda", 1:length(lambda_to_test), sep="_"), paste("alpha", 1:length(alpha_to_test), sep="_"))
            minix <- order(cvmses, decreasing=FALSE)[1]
            idx <- minix%%length(lambda_to_test)
            if (idx == 0)
            {
                idx <- length(lambda_to_test)
            }
            best_params <- c(alpha_to_test[ceiling(minix/length(lambda_to_test))], lambda_to_test[idx])
            model <- glmnet::glmnet(x=training_set[,-1], y=training_set[,1], alpha=best_params[1], lambda=best_params[2])
            p <- as.vector(predict(object=model, newx=test_data, type="response")[, 1])
         }
         r <- cor(test_labels, p[common_indices], use="complete.obs")#, method="spearman")
         message(paste(drug[[1]], method, r, sep="\t"))
         return(r)
     })
     names(predictions) <- methods
     return(predictions)
})
names(metric) <- drugs
