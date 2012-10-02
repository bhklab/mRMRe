library(Hmisc)
library(mRMRe)

##
## Tests
##

#
# correlate - Testing the correlate method against other known and proven methods
#

## Cont vs Cont

a <- runif(100)
b <- runif(100)

test_correlate(a,b, method = "pearson")
test_correlate(a,b, method = "spearman")
test_correlate(a,b, method = "kendall")
test_correlate(a,b, method = "cindex")
test_correlate(a,b, method = "frequency")

## Cat vs Cont

a <- sample(0:3, 100, T)
b <- runif(100)

test_correlate(a,b, method = "pearson")
test_correlate(a,b, method = "spearman")
test_correlate(a,b, method = "kendall")
test_correlate(a,b, method = "cindex")
test_correlate(a,b, method = "frequency")

## Cat vs Cat

a <- sample(0:3, 100, T)
b <- sample(0:5, 100, T)

test_correlate(a,b, method = "pearson")
test_correlate(a,b, method = "spearman")
test_correlate(a,b, method = "kendall")
test_correlate(a,b, method = "cindex")
test_correlate(a,b, method = "frequency")
test_correlate(a,b, method = "cramersv")

#
# mim - Testing the mim method against the correlate method
#

dd <- data.frame(
        "surv1" = Surv(runif(100), sample(0:1, 100, replace = TRUE)),
        "cont1" = runif(100),
        "cat1"  = factor(sample(1:5, 100, replace = TRUE), ordered = TRUE),
        "surv2" = Surv(runif(100), sample(0:1, 100, replace = TRUE)),
        "cont2" = runif(100),
        "cat2"  = factor(sample(1:5, 100, replace = TRUE), ordered = TRUE))

data <- mRMR.data(data = dd)
cors <- mim(data, method = "cor")
combinations <- combn(featureNames(data), 2)
# results <- apply(combinations, 2, function(i) cors[i[[1]], i[[2]]] == correlate(i, j))
# FIXME: Finish this test

##
## Methods
##

test_correlate <- function(a,b, method)
{
	run_correlate_test(a, b, method)
	a[sample(1:length(a), round(length(a) / 10))] <- NA
	b[sample(1:length(b), round(length(b) / 10))] <- NA
	message("Test with NA")
	run_correlate_test(a, b, method)
}

run_correlate_test <- function(a, b, method)
{
	message("Testing with 3rd party")
	if (method == "pearson" || method == "spearman" || method == "kendall")
		confirmation <- cor(a, b, method = method, use = "complete.obs")
	else if (method == "frequency")
		confirmation <- mean (a > b, na.rm = T)
	else if (method == "cindex")
		confirmation <- as.numeric(Hmisc::rcorr.cens(a,b)[1])
	else if (method == "cramersv")
		confirmation <- as.numeric(sqrt(chisq.test(a, b, correct=FALSE)$statistic /
								(length(a) * (min(length(unique(a)),length(unique(b))) - 1))))
	
	correlation <- correlate(a,b, method=method)$s
	if(abs(confirmation - correlation) > 1e5)
		stop("Correlation is different that confirmation", confirmation, correlation)
	else
		message("\t3rd Party OK")
	
	message("Testing for symmetry")
	if (abs(correlate(b,a, method=method)$s - correlation) > 1e10)
		stop("Correlation is not symmetric")
	else
		message("\tSymmetry OK")
	
	message("Testing for stratification")
	s <- as.factor(c(rep(0,length(a)),rep(1,length(a)), rep(2,length(a))))
	w <- c(rep(runif(1),length(a)),rep(runif(1),length(a)), rep(runif(1),length(a)))
	a <- c(a,a,a)
	b <- c(b,b,b)
	
	if (abs(correlate(a,b, method=method, strata=s)$s - correlation) > 1e10)
		stop("Stratification is wrong")
	else if (abs(correlate(a,b, method=method, strata=s)$s - correlation) > 1e10)
		stop("Stratification is not symmetric")
	else
		message("\tStratification OK")
	message("Testing for weights")
	if (abs(correlate(a,b, method=method, weights=w)$s - correlation) > 1e10)
		stop("Weighting is wrong")
	else if (abs(correlate(b,a, method=method, weights=w)$s - correlation) > 1e10)
		stop("Weighting is not symmetric")
	else
		message("\tWeighting OK")
	
	message("Test for strata+weights")
	if(abs(correlate(a,b, method=method, weights=w, strata=s)$s - correlation) > 1e10)
		stop("Stratification + Weighting is wrong")
	else if (abs(correlate(b,a, method=method, weights=w, strata=s)$s - correlation) > 1e10)
		stop("Stratification + Weighting is not symmetric")
	else
		message("\tStratification + Weighting OK")
}
