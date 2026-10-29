value <- function(...) {
    r <- correlate(...)
    if (is.null(r$estimate)) r$statistic else r$estimate
}

# Bootstrap paths used to segfault (null pTotalWeight), call Cramer's V from the
# Pearson loop, and use a broken variance. The C++ RNG is wall-clock seeded, so
# only finiteness and ranges are asserted.
test_that("mim with bootstrap_count 5 returns a finite matrix", {
    data(cgps)
    d <- mRMR.data(data = data.frame(cgps.ge[, 1:10]))
    m <- mim(d, bootstrap_count = 5)
    expect_equal(dim(m), c(10L, 10L))
    expect_true(all(is.finite(m[upper.tri(m)])))
    expect_true(all(abs(m[upper.tri(m)]) <= 1))
})

test_that("correlate with bootstrap_count 5 returns finite estimates", {
    data(cgps)
    x <- cgps.ge[, 1]
    y <- cgps.ge[, 2]
    for (method in c("pearson", "spearman", "cindex")) {
        r <- value(x, y, method = method, bootstrap_count = 5)
        expect_true(is.finite(r), info = method)
        expect_true(abs(r) <= 1, info = method)
    }
})

test_that("pearson bootstrap estimates Pearson, not Cramer's V", {
    set.seed(1)
    x <- rnorm(200)
    y <- x + rnorm(200, sd = 0.5)
    r <- value(x, y, method = "pearson", bootstrap_count = 50)
    expect_equal(r, value(x, y, method = "pearson"), tolerance = 0.1)
})

test_that("cramersv with a stratum smaller than the level count is finite", {
    set.seed(2)
    x <- factor(sample(0:5, 40, replace = TRUE), levels = 0:5, ordered = TRUE)
    y <- factor(sample(0:5, 40, replace = TRUE), levels = 0:5, ordered = TRUE)
    strata <- factor(c(rep(0, 37), rep(1, 3)))
    expect_true(is.finite(value(x, y, method = "cramersv")))
    expect_true(is.finite(value(x, y, method = "cramersv", strata = strata)))
})
