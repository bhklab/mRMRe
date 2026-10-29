mixed <- function(n = 60) {
    set.seed(46)
    data.frame(
        a = rnorm(n),
        f = factor(sample(1:3, n, TRUE), ordered = TRUE),
        s = Surv(rexp(n), rbinom(n, 1, .5)),
        b = rnorm(n),
        g = factor(sample(1:4, n, TRUE), ordered = TRUE))
}

test_that("subsetData preserves types, Surv columns and names", {
    d <- mRMR.data(data = mixed())
    s <- subsetData(d, 1:30, c(2, 3, 5))
    expect_equal(s@feature_types, c(1, 2, 3, 1))
    expect_equal(s@feature_names, c("f", "s", "g"))
    expect_equal(s@data, d@data[1:30, c(2, 3, 4, 6), drop = FALSE])
    expect_equal(s@sample_names, d@sample_names[1:30])
    expect_equal(featureData(s)$s, featureData(d)$s[1:30])
    expect_equal(subsetData(d, 1:2, 1:2)@feature_types, c(0, 1))
})

test_that("mim of a subset equals the submatrix of mim of the full data", {
    d <- mRMR.data(data = mixed())
    s <- subsetData(d, column_indices = c(1, 2, 3, 5))
    expect_equal(mim(s), mim(d)[c(1, 2, 3, 5), c(1, 2, 3, 5)], tolerance = 1e-6)
})

test_that("subsetData subsets strata, weights and priors", {
    d <- mRMR.data(data = mixed(), strata = factor(rep(1:2, 30)), weights = seq_len(60))
    p <- matrix(0.5, 5, 5)
    diag(p) <- 1
    priors(d) <- p
    s <- subsetData(d, 1:10, c(1, 3, 4))
    expect_equal(unname(sampleWeights(s)), 1:10)
    expect_equal(unname(sampleStrata(s)), rep(0:1, 5))
    expect_equal(dim(priors(s)), c(3, 3))
})

test_that("correlate cramersv works on unordered factors", {
    set.seed(1)
    x <- factor(sample(letters[1:3], 30, TRUE))
    y <- factor(sample(letters[1:3], 30, TRUE))
    expect_true(is.finite(correlate(x, y, method = "cramersv")$statistic))
})

test_that("priors<- compares against featureCount and checks symmetry", {
    d <- mRMR.data(data = data.frame(x = rnorm(50), y = rnorm(50), s = Surv(rexp(50), rbinom(50, 1, .5))))
    priors(d) <- matrix(1, 3, 3)
    expect_equal(dim(priors(d)), c(3, 3))
    expect_error(priors(d) <- matrix(1, 4, 4), "symmetric")
    asym <- matrix(1, 3, 3)
    asym[1, 2] <- 0
    expect_error(priors(d) <- asym, "symmetric")
})

test_that("new('mRMRe.Data') with no arguments works", {
    expect_s4_class(new("mRMRe.Data"), "mRMRe.Data")
})

test_that("integer columns are accepted", {
    d <- mRMR.data(data = data.frame(x = 1:50, y = rnorm(50)))
    expect_equal(d@feature_types, c(0, 0))
    d <- mRMR.data(data = data.frame(x = 1:50, f = factor(rep(1:5, 10), ordered = TRUE)))
    expect_equal(d@feature_types, c(0, 1))
})

test_that("non-right-censored Surv is rejected clearly", {
    n <- 20
    cp <- data.frame(x = rnorm(n), s = Surv(rep(0, n), rexp(n) + 1, rbinom(n, 1, .5)))
    expect_error(mRMR.data(data = cp), "right-censored")
})

test_that("mim validates prior_weight length", {
    d <- mRMR.data(data = data.frame(x = rnorm(50), y = rnorm(50)))
    priors(d) <- matrix(1, 2, 2)
    expect_error(mim(d, prior_weight = c(0.1, 0.2)), "prior weight")
    expect_error(mim(d, prior_weight = 2), "prior weight")
})

test_that("sampleWeights<- rejects negative or non-numeric weights", {
    d <- mRMR.data(data = data.frame(x = rnorm(10), y = rnorm(10)))
    expect_error(sampleWeights(d) <- c(-1, rep(1, 9)), "non-negative")
    expect_error(sampleWeights(d) <- rep("a", 10), "numeric")
})

test_that("mim of a one-feature object keeps matrix shape", {
    d <- mRMR.data(data = data.frame(x = rnorm(20)))
    expect_equal(dim(mim(d)), c(1, 1))
})

test_that("expandFeatureIndices handles indices after a Surv column", {
    d <- mRMR.data(data = mixed())
    expect_equal(mRMRe:::.expandFeatureIndices(d, 1:5), c(1L, 2L, 3L, 5L, 6L))
    expect_equal(mRMRe:::.expandFeatureIndices(d, 1:2), 1:2)
})
