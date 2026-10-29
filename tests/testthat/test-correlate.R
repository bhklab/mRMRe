tol <- 1e-6

# correlate() against base R / survival reference implementations.
# Categorical inputs are ordered factors (mRMR.data rejects plain integers).
ofac <- function(n, k) factor(sample(0:(k - 1), n, replace = TRUE), ordered = TRUE)

reference <- function(a, b, method)
{
    x <- as.numeric(a)
    y <- as.numeric(b)
    switch(method,
        pearson = , spearman = , kendall = cor(x, y, method = method, use = "complete.obs"),
        frequency = mean(x > y, na.rm = TRUE),
        # cindex(a, b) matches concordance with a as the response and b as the
        # predictor, when b has no ties. (With a tied, cindex(b, a) is not equal to
        # cindex(a, b), and tie handling differs from survival::concordance when both
        # are categorical, so those cases are not compared.)
        cindex = unname(survival::concordance(Surv(x) ~ y)$concordance),
        cramersv = {
            xy <- na.omit(data.frame(a, b))
            unname(sqrt(chisq.test(xy$a, xy$b, correct = FALSE)$statistic /
                (nrow(xy) * (min(nlevels(droplevels(xy$a)), nlevels(droplevels(xy$b))) - 1))))
        })
}

value <- function(...) {
    r <- correlate(...)
    if (is.null(r$estimate)) r$statistic else r$estimate
}

# Check against reference, symmetry, strata and weights.
# ref = FALSE skips the reference comparison (estimator conventions differ, see callers).
check_correlate <- function(a, b, method, symmetric = TRUE, ref = TRUE)
{
    got <- value(a, b, method = method)
    if (ref) expect_equal(got, suppressWarnings(reference(a, b, method)), tolerance = tol)

    if (symmetric) {
        sym <- value(b, a, method = method)
        # frequency is P(a > b), so swapping the arguments gives its complement
        expect_equal(sym, if (method == "frequency") 1 - got else got, tolerance = tol)
    }

    # three copies of the data, each its own stratum / constant weight, must not change the estimate
    n <- length(a)
    idx <- rep(seq_len(n), 3)
    a3 <- a[idx]
    b3 <- b[idx]
    s <- factor(rep(0:2, each = n))
    w <- rep(runif(3), each = n)
    expect_equal(value(a3, b3, method = method, strata = s), got, tolerance = tol)
    # Spearman with weights alone re-ranks the pooled weighted sample, so it differs from the
    # unweighted value; with strata (below) it matches.
    if (method != "spearman")
        expect_equal(value(a3, b3, method = method, weights = w), got, tolerance = tol)
    expect_equal(value(a3, b3, method = method, strata = s, weights = w), got, tolerance = tol)
}

with_na <- function(x) { x[sample(length(x), round(length(x) / 10))] <- NA; x }

test_that("continuous vs continuous matches reference", {
    set.seed(1)
    a <- runif(100)
    b <- runif(100)
    for (method in c("pearson", "spearman", "kendall", "cindex", "frequency")) {
        check_correlate(a, b, method)
        check_correlate(with_na(a), with_na(b), method)
    }
})

test_that("categorical vs continuous matches reference", {
    set.seed(2)
    a <- ofac(100, 4)
    b <- runif(100)
    for (method in c("pearson", "frequency")) {
        check_correlate(a, b, method)
        check_correlate(with_na(a), with_na(b), method)
    }
    # With tied values mRMRe's Spearman / Kendall differ from cor() (no tie correction,
    # tau-a rather than tau-b), so only check strata/weights consistency (and they are not symmetric under ties).
    for (method in c("spearman", "kendall")) {
        check_correlate(a, b, method, ref = FALSE, symmetric = FALSE)
        check_correlate(with_na(a), with_na(b), method, ref = FALSE, symmetric = FALSE)
    }
    # cindex is directional when the first argument has ties: only (a, b) is comparable
    check_correlate(a, b, "cindex", symmetric = FALSE)
    check_correlate(with_na(a), with_na(b), "cindex", symmetric = FALSE)
})

test_that("categorical vs categorical matches reference", {
    set.seed(3)
    a <- ofac(100, 4)
    b <- ofac(100, 6)
    for (method in c("pearson", "cramersv")) {
        check_correlate(a, b, method)
        check_correlate(with_na(a), with_na(b), method)
    }
    for (method in c("spearman", "kendall")) # tie handling and asymmetric with ties, see above
        check_correlate(a, b, method, ref = FALSE, symmetric = FALSE)
    # frequency is P(a > b) and not complementary under ties, so skip the symmetry check
    check_correlate(a, b, "frequency", symmetric = FALSE)
    # cindex is not compared to survival::concordance here: tie handling differs.
})

test_that("survival data requires cindex", {
    set.seed(4)
    s <- Surv(runif(50), sample(0:1, 50, replace = TRUE))
    expect_error(correlate(s, runif(50), method = "pearson"), "cindex")
})
