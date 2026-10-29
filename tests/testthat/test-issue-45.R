test_that("set.thread.count validates its argument", {
    for (bad in list(NA, NA_integer_, 0, -1, 1.5, c(1, 2), integer(0), "2", NULL))
        expect_error(set.thread.count(bad), "single whole number")
})

test_that("set.thread.count returns the effective count invisibly", {
    before <- get.thread.count()
    on.exit(if (before >= 1) set.thread.count(before))
    res <- withVisible(set.thread.count(1))
    expect_false(res$visible)
    expect_true(is.integer(res$value))
    expect_equal(res$value, if (before >= 1) 1L else 0L)
    expect_identical(get.thread.count(), res$value)
})

test_that("get.thread.count returns a single integer", {
    expect_true(is.integer(get.thread.count()))
    expect_length(get.thread.count(), 1)
})

test_that("mim, correlate and mRMR.classic results are unchanged", {
    set.seed(1)
    x <- matrix(rnorm(200), 40, 5)
    x[, 2] <- x[, 2] + x[, 1]
    d <- mRMR.data(data = data.frame(x))
    m <- mim(d, method = "cor")
    expect_equal(m[1, 2], 0.773328678399, tolerance = 1e-8)
    expect_equal(m[1, 3], -0.2514937910336, tolerance = 1e-8)
    expect_equal(m[2, 3], -0.0278404025624, tolerance = 1e-8)
    expect_equal(correlate(x[, 1], x[, 2], method = "spearman")$statistic, 0.754409005629, tolerance = 1e-8)
    expect_equal(correlate(x[, 1], x[, 2], method = "kendall")$statistic, 0.564102564103, tolerance = 1e-8)
    ci <- correlate(x[, 1], x[, 2], method = "cindex")
    expect_equal(ci$estimate, 0.782051282051, tolerance = 1e-8)
    expect_equal(ci$se, 0.0268103804178, tolerance = 1e-8)
    f <- mRMR.classic(data = d, target_indices = 1, feature_count = 2)
    expect_equal(c(solutions(f)[[1]]), c(2, 3))
    expect_equal(c(scores(f)[[1]]), c(0.455697921749, 0.032281260368), tolerance = 1e-8)
    expect_equal(mim(f)[1, 2], 0.455697921749, tolerance = 1e-8)
})
