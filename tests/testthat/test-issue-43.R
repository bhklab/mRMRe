library(mRMRe)

make_data <- function() {
    set.seed(1)
    mRMR.data(data = data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30)))
}

test_that("fixed_feature_count == length(levels) errors instead of hanging", {
    expect_error(mRMR.classic(data = make_data(), target_indices = 1, feature_count = 2,
                              fixed_feature_count = 2), "smaller than length")
})

test_that("a single out-of-range target index is rejected", {
    expect_error(mRMR.classic(data = make_data(), target_indices = 999, feature_count = 2),
                 "target_indices")
    expect_error(mRMR.classic(data = make_data(), target_indices = 0, feature_count = 2),
                 "target_indices")
})

test_that("feature_count = 1 works for solutions and scores", {
    f <- mRMR.classic(data = make_data(), target_indices = 1, feature_count = 1)
    expect_equal(dim(solutions(f)[[1]]), c(1L, 1L))
    expect_length(scores(f)[[1]], 1)
})

test_that("unfilled solution slots are NA, not the target", {
    f <- mRMR.ensemble(data = make_data(), target_indices = 1, solution_count = 2,
                       feature_count = 2)
    sol <- solutions(f)[[1]]
    expect_true(anyNA(sol))
    expect_false(1 %in% sol)
})

test_that("bootstrap rejects levels[-1] != 1", {
    expect_error(new("mRMRe.Filter", data = make_data(), target_indices = 1,
                     levels = c(2, 2), method = "bootstrap"), "bootstrap")
})
