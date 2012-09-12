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
mim(data)

data <- mRMR.data(data = dd,
        strata = factor(sample(1:5, 100, replace=TRUE), ordered=TRUE),
        weights = runif(100))
mim(data) # Gives MI matrix
mim(data, method = "cor") # Gives correlation matrix

filter_1 <- mRMR.ensemble("mRMRe.Filter", data = data, target_index = 1, feature_count = 3, solution_count = 1)
mim(data)
mim(filter_1)

filter_2 <- mRMR.ensemble("mRMRe.Filter", data = data, target_index = 2, feature_count = 3, solution_count = 1)
mim(data)
mim(filter_2)

# No wrapper just yet
network <- new("mRMRe.Network", data = data, target_indices = c(1, 2), levels = c(1, 1, 1), layers = 1)
mim(data)
mim(network)

causality(network)[1, , ]
causality(filter_1)

causality(network)[2, , ]
causality(filter_2)

# I have no idea how to get igraph to print out vertex names

visualize(network)