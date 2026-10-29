test_that("a survival column without its partner is rejected, not read out of bounds", {
    set.seed(51)
    n <- 30
    d <- mRMR.data(data = data.frame(a = rnorm(n), s = Surv(rexp(n), rbinom(n, 1, .5)), b = rnorm(n)))
    expect_equal(d@feature_types, c(0, 2, 3, 0))
    expect_silent(mim(d))

    last_event <- d
    last_event@feature_types <- c(0, 0, 0, 2)
    expect_error(mim(last_event), "must be followed by a survival time column")

    lone_time <- d
    lone_time@feature_types <- c(3, 2, 3, 0)
    expect_error(mim(lone_time), "must be preceded by a survival event column")
})
