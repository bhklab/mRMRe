#load('irinotecan_cgp_ccle.RData')
#library(mRMRe)

`compute.press` <- function(X, Y, lambda=1e-16)
{
	indices <- complete.cases(Y)
	X <- X[indices, ]
	Y <- Y[indices]
	X <- cbind(rep(1, nrow(X)), X)
	transposed_X <- t(X)
	H1 <- MASS::ginv(transposed_X %*% X + diag(x=lambda, ncol=ncol(X), nrow=ncol(X)))
	multiplied <- H1 %*% transposed_X
	coefficients <- multiplied %*% Y
	H <- X %*% multiplied
	Y.hat <- X %*% coefficients
	
	residuals <- Y-Y.hat
	residuals.loo <- residuals / (1-diag(H))
	
	return(mean(abs(residuals.loo)))
}

`bootstrap.press` <- function(input, target, features, bootstrap_count=1000)
{
	result <- sapply(1:1000, function(i){
				indices <- sample(nrow(input), replace=TRUE)
				X <- input[indices, features]
				Y <- target[indices]
				compute.press(X,Y)
			})
	return(result)
}

mrmr_spearman <- mRMR.classic(data.frame(target=ic50_cgp, data_cgp), 1, 30)$paths
mrmr_pearson <- mRMR.classic(data.frame(target=ic50_cgp, data_cgp), 1, 30, uses_ranks=FALSE)$paths
random <- sample(ncol(data_cgp), 30)
ranking <- order(apply(data_cgp, 2, cor, ic50_cgp, use="complete.obs"), decreasing=TRUE)

mrmr_spearman_press <- bootstrap.press(data_cgp, ic50_cgp, mrmr_spearman)
mrmr_pearson_press <- bootstrap.press(data_cgp, ic50_cgp, mrmr_pearson)
random_press <- bootstrap.press(data_cgp, ic50_cgp, random)
single_gene_press <- bootstrap.press(data_cgp, ic50_cgp, ranking[1])
rank_press <- bootstrap.press(data_cgp, ic50_cgp, ranking[1:30])



pdf("press_dist_multi.pdf")
plot(density(random_press), xlim=c(0.5,1), ylim=c(0,20), col="darkgrey")
sapply(1:9, function(i){
			random <- sample(ncol(data_cgp), 30)
			random_press <- bootstrap.press(data_cgp, ic50_cgp, random)
			lines(density(random_press), col="darkgrey")
		})
lines(density(mrmr_pearson_press), col=2)
lines(density(mrmr_spearman_press), col=3)
#lines(density(single_gene_press), col=4)
lines(density(single_gene_press), col=4)
legend("topright", legend=c("RANDOM", "PEARSON", "SPEARMAN", "RANK"), col=c("darkgrey",2:4), bty="n", pch=15)
dev.off()


formula <- paste("target ~ ", paste(colnames(data_cgp)[mrmr_pearson], collapse=" + "))
model <- lm(formula=formula, data=data.frame(target=ic50_cgp, data_cgp))
predictions <- predict(model, newdata=data.frame(data_ccle_new))
errors <- abs(predictions - ic50_ccle_new)