library(mRMRe)

dd <- data.frame(
        "surv1"=Surv(runif(100), sample(0:1, 100, replace=TRUE)),
        "cont1"=runif(100),
        "cat1"=factor(sample(1:5, 100, replace=TRUE), ordered=TRUE),
        "surv2"=Surv(runif(100), sample(0:1, 100, replace=TRUE)),
        "cont2"=runif(100),
        "cont3"=runif(100),
        "surv3"=Surv(runif(100),
                sample(0:1, 100, replace=TRUE)),
        "cat2"=factor(sample(1:5, 100, replace=TRUE), ordered=TRUE)
)

data <- mRMR.data(data = dd)
filter_1 <- mRMR.ensemble("mRMRe.Filter", data = data, target_index = 1, feature_count = 3, solution_count = 1)
filter_2 <- mRMR.ensemble("mRMRe.Filter", data = data, target_index = 2, feature_count = 3, solution_count = 1)


network <- new("mRMRe.Network", data = data, target_indices = c(1, 2), levels = c(1, 1, 1), layers = 1)

mim(filter_1)
mim(data)
mim(filter_2)
mim(network)