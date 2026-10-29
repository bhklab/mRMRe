test_that("classic filter picks the correlated features", {
    set.seed(0)
    x <- rnorm(100)
    dd <- data.frame(
        cont1 = x,
        cont2 = x + rnorm(100, 0.1),
        cont3 = x + rnorm(100, 0.1),
        cont4 = x + rnorm(100, 0.1),
        cont5 = x + rnorm(100, 0.1),
        cont6 = x + rnorm(100, 0.01))
    data <- mRMR.data(data = dd)
    filter <- mRMR.classic("mRMRe.Filter", data = data, target_indices = 3:5, feature_count = 2)

    sol <- solutions(filter)
    expect_equal(names(sol), c("3", "4", "5"))
    for (target in names(sol)) {
        # one solution of 2 distinct features, never the target itself
        expect_equal(dim(sol[[target]]), c(2, 1))
        expect_false(anyDuplicated(sol[[target]]) > 0)
        expect_false(as.integer(target) %in% sol[[target]])
        expect_false(anyNA(scores(filter)[[target]]))
    }
})

test_that("network, filter and data mim agree on non-NA cells (cgps)", {
    data(cgps)
    dd <- mRMR.data(data = data.frame(cgps.ge[, 1:20]))
    netw <- new("mRMRe.Network", data = dd, target_indices = c(1, 2), levels = c(1, 1, 1), layers = 1)
    filter <- new("mRMRe.Filter", data = dd, target_indices = c(1, 2), levels = c(1, 1, 1))

    network_mim <- mim(netw)
    keep <- !is.na(network_mim)
    expect_true(any(keep))
    expect_equal(mim(filter)[keep], mim(dd)[keep], tolerance = 1e-5)
    expect_equal(network_mim[keep], mim(dd)[keep], tolerance = 1e-5)

    expect_length(causality(filter), 2)
})
