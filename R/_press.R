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
