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