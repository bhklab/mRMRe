test_that("mim(filter) equals mim(data) on non-NA cells", {
    set.seed(5)
    n <- 100
    dd <- data.frame(
        surv1 = Surv(runif(n), sample(0:1, n, replace = TRUE)),
        cont1 = runif(n),
        cat1  = factor(sample(1:5, n, replace = TRUE), ordered = TRUE),
        surv2 = Surv(runif(n), sample(0:1, n, replace = TRUE)),
        cont2 = runif(n),
        cont3 = runif(n),
        surv3 = Surv(runif(n), sample(0:1, n, replace = TRUE)),
        cat2  = factor(sample(1:5, n, replace = TRUE), ordered = TRUE))
    data <- mRMR.data(data = dd)

    master <- mim(data)
    # target_indices has the same length as the number of Surv columns (3) on purpose:
    # .expandFeatureIndices recycles `indices >= adaptor` and warns otherwise.
    # TODO: separate issue, not covered here.
    filters <- list(
        mRMR.ensemble("mRMRe.Filter", data = data, target_indices = c(1, 5, 8),
            feature_count = 2, solution_count = 2),
        mRMR.ensemble("mRMRe.Filter", data = data, target_indices = c(1, 2, 3),
            feature_count = 2, solution_count = 3))

    for (filter in filters) {
        m <- mim(filter)
        keep <- !is.na(m)
        expect_true(any(keep))
        expect_equal(m[keep], master[keep], tolerance = 1e-5)
    }
})

test_that("mim of a correlation method is symmetric with unit diagonal for continuous data", {
    set.seed(6)
    data <- mRMR.data(data = data.frame(matrix(rnorm(500), 100, 5)))
    m <- mim(data, method = "cor")
    expect_equal(m, t(m), tolerance = 1e-6)
    expect_equal(unname(diag(m)), rep(1, 5), tolerance = 1e-6)
    expect_equal(unname(m[1, 2]), cor(data@data[, 1], data@data[, 2]), tolerance = 1e-6)
})
