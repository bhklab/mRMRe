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

data <- new("mRMRe.Data", data = dd)
filter <- new("mRMRe.Filter", data = data, target_index = 1, levels = c(1, 1))


#causality(filter)
#network <- new("mRMRe.Network", data = DATA, target_indices = c(1, 2), levels = c(3, 1, 1), layers=1)
#visualize(network)
mim(DATA)
mim(filter)