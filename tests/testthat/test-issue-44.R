# Race-free MI matrix / Spearman order cache / filter prefetch. The thread-count tests are
# no-ops without OpenMP, but on an OpenMP build results must not depend on the thread count.
data(cgps)
d <- mRMR.data(data = data.frame(cgps.ge[, 1:30]))
orig_threads <- get.thread.count()

with_threads <- function(n, expr) {
    set.thread.count(n)
    on.exit(set.thread.count(orig_threads))
    force(expr)
}

test_that("spearman mim is identical across thread counts and repeats", {
    ref <- with_threads(1, mim(d, continuous_estimator = "spearman"))
    for (k in 1:5)
        expect_identical(with_threads(4, mim(d, continuous_estimator = "spearman")), ref)
})

test_that("filters are identical across thread counts", {
    mk <- function() new("mRMRe.Filter", data = d, target_indices = 1:2, levels = c(3, 2, 1),
                         continuous_estimator = "spearman")
    ref <- with_threads(1, mk())
    for (k in 1:5) {
        f <- with_threads(4, mk())
        expect_identical(solutions(f), solutions(ref))
        expect_identical(scores(f), scores(ref))
    }
})

test_that("spearman mim after a bootstrap filter equals a fresh spearman mim", {
    fresh <- mim(d, continuous_estimator = "spearman")
    set.seed(1)
    new("mRMRe.Filter", data = d, target_indices = 1, levels = c(3, 1, 1),
        method = "bootstrap", continuous_estimator = "spearman")
    expect_identical(mim(d, continuous_estimator = "spearman"), fresh)
})

test_that("NaN cells are stable across repeated calls", {
    x <- data.frame(cgps.ge[, 1:4])
    x[, 2] <- NA_real_
    dn <- mRMR.data(data = x)
    a <- mim(dn)
    expect_true(all(is.na(a[2, -2])))
    expect_identical(mim(dn), a)
    mk <- function() solutions(new("mRMRe.Filter", data = dn, target_indices = 1, levels = c(1, 1)))
    expect_identical(mk(), mk())
})
