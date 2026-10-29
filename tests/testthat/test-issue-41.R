# Bootstrap draws must come from R's RNG: reproducible under set.seed, different
# for different seeds (they used to be wall-clock seeded).
strata <- factor(rep(0:1, 30))
xy <- local({ set.seed(7); x <- rnorm(60); list(x = x, y = x + rnorm(60)) })
boot_cor <- function(seed) {
    set.seed(seed)
    correlate(xy$x, xy$y, method = "pearson", strata = strata, bootstrap_count = 20)$statistic
}

test_that("correlate bootstrap is reproducible under set.seed", {
    a <- boot_cor(1)
    Sys.sleep(1.1)
    expect_identical(a, boot_cor(1))
    expect_false(identical(a, boot_cor(2)))
})

test_that("mim bootstrap is reproducible under set.seed", {
    data(cgps)
    d <- mRMR.data(data = data.frame(cgps.ge[, 1:6]),
                   strata = factor(rep(0:1, length.out = nrow(cgps.ge))))
    run <- function(seed) { set.seed(seed); mim(d, bootstrap_count = 20) }
    a <- run(1)
    Sys.sleep(1.1)
    expect_identical(a, run(1))
    expect_false(identical(a, run(2)))
})

test_that("bootstrap filter is reproducible under set.seed", {
    data(cgps)
    d <- mRMR.data(data = data.frame(cgps.ge[, 1:20]))
    run <- function(seed) {
        set.seed(seed)
        new("mRMRe.Filter", data = d, target_indices = 1, levels = c(5, 1, 1),
            method = "bootstrap")
    }
    a <- run(1)
    Sys.sleep(1.1)
    expect_identical(solutions(a), solutions(run(1)))
    expect_false(identical(solutions(a), solutions(run(2))))
})

test_that("bootstrap_count 0 leaves the RNG stream untouched", {
    set.seed(3); u <- runif(1)
    set.seed(3); correlate(xy$x, xy$y, method = "pearson", strata = strata)
    expect_identical(runif(1), u)
})
