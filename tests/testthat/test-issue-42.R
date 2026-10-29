dat <- function(n = 40) {
    set.seed(42)
    data.frame(a = rnorm(n), b = rnorm(n), c = rnorm(n))
}
s <- rep(c("a", "c"), 20)

test_that("unused strata levels give dense codes and match the 2-level factor", {
    x <- dat()
    d <- mRMR.data(data = x, strata = factor(s, levels = c("a", "b", "c")))
    d2 <- mRMR.data(data = x, strata = factor(s))
    expect_equal(unname(sampleStrata(d)[1:4]), c(0, 1, 0, 1))
    m <- mim(d)
    expect_true(all(is.finite(m[upper.tri(m)])))
    expect_equal(m, mim(d2))
})

test_that("sampleStrata<- drops unused levels and accepts character/integer", {
    d <- mRMR.data(data = dat())
    sampleStrata(d) <- factor(s, levels = c("a", "b", "c"))
    expect_equal(max(d@strata), 1)
    sampleStrata(d) <- s
    expect_equal(unique(d@strata), c(0L, 1L))
    sampleStrata(d) <- rep(c(5L, 9L), 20)
    expect_equal(unique(d@strata), c(0L, 1L))
})

test_that("NA strata and length mismatch error clearly", {
    d <- mRMR.data(data = dat())
    expect_error(sampleStrata(d) <- factor(c(NA, s[-1])), "missing values")
    expect_error(sampleStrata(d) <- factor(s[-1]), "same number of samples")
})
