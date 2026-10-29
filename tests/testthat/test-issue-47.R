cgps_data <- function() {
    data(cgps, package = "mRMRe", envir = environment())
    mRMR.data(data = data.frame(cgps.ge[, 1:8]))
}

test_that("last-layer targets whose links are all known stay in the network", {
    dd <- cgps_data()
    net <- new("mRMRe.Network", data = dd, target_indices = 1, levels = c(1, 1), layers = 2)
    sol <- solutions(net)
    in_solutions <- sort(unique(c(1L, na.omit(as.vector(unlist(sol))))))
    expect_true(all(in_solutions %in% as.integer(names(sol))))
    adj <- adjacencyMatrix(net)
    expect_true(all(colSums(adj)[in_solutions] > 0))
})

test_that("correlate returns the same names on the early-return path", {
    few <- correlate(survival::Surv(c(1, 2, 3, 4, 5), c(1, 1, 0, 1, 0)), c(5, 4, 3, 2, 1), method = "cindex")
    set.seed(2)
    many <- correlate(survival::Surv(rexp(60), rbinom(60, 1, 0.5)), rnorm(60), method = "cindex")
    expect_true("estimate" %in% names(few))
    expect_equal(names(few), names(many))
})

test_that("scores() guards NA in solutions for Data and Filter-style input", {
    dd <- cgps_data()
    s <- scores(dd, list("1" = matrix(NA, 3, 1)))
    expect_equal(dim(s[["1"]]), c(3L, 1L))
    expect_true(all(is.na(s[["1"]])))

    filter <- new("mRMRe.Filter", data = dd, target_indices = 1:2, levels = c(1, 1, 1))
    expect_false(anyNA(scores(filter)[["2"]]))
    net <- new("mRMRe.Network", data = dd, target_indices = 1, levels = c(1, 1), layers = 2)
    expect_named(scores(net), names(solutions(net)))
})

test_that("Filter no longer carries a scores slot and class check uses is()", {
    expect_false("scores" %in% slotNames("mRMRe.Filter"))
    expect_error(new("mRMRe.Filter", data = 1, target_indices = 1, levels = 1), "mRMRe.Data")
})

test_that(".map.continuous.estimator maps and rejects", {
    f <- mRMRe:::.map.continuous.estimator
    expect_identical(sapply(c("pearson", "spearman", "kendall", "frequency"), f, USE.NAMES = FALSE), 0:3)
    expect_error(f("bogus"))
})

test_that("visualize() draws without error or warning", {
    dd <- cgps_data()
    net <- new("mRMRe.Network", data = dd, target_indices = 1, levels = c(1, 1), layers = 2)
    pdf(NULL)
    on.exit(dev.off())
    expect_silent(visualize(net))
})

test_that("mim() method for Network has the generic's method formal", {
    expect_true("method" %in% names(formals(selectMethod("mim", "mRMRe.Network")@.Data)))
})
