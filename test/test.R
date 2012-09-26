## Size

library(mRMRe)
set.seed(0)
data(cgps)
data <- mRMR.data(data = as.data.frame(cgps_ge))
system.time(filter <- mRMR.ensemble("mRMRe.Filter", data = data, target_indices = c(1, 2, 3, 4, 5, 6, 7, 8), feature_count = 500, solution_count = 10))
print(object.size(filter), units = "Mb")


## Simple Test

library(mRMRe)
set.seed(0)
dd <- data.frame(
        "surv1" = Surv(runif(100), sample(0:1, 100, replace = TRUE)),
        "cont1" = runif(100),
        "cat1"  = factor(sample(1:5, 100, replace = TRUE), ordered = TRUE),
        "surv2" = Surv(runif(100), sample(0:1, 100, replace = TRUE)),
        "cont2" = runif(100),
        "cont3" = runif(100),
        "surv3" = Surv(runif(100), sample(0:1, 100, replace = TRUE)),
        "cat2"=factor(sample(1:5, 100, replace = TRUE), ordered = TRUE))

data <- mRMR.data(data = dd)
filter <- mRMR.ensemble("mRMRe.Filter", data = data, target_indices = 3:5, feature_count = 2, solution_count = 2)


## NETWORK TEST

library(mRMRe)
set.seed(0)
dd <- data.frame(
        "surv1" = Surv(runif(100), sample(0:1, 100, replace = TRUE)),
        "cont1" = runif(100),
        "cat1"  = factor(sample(1:5, 100, replace = TRUE), ordered = TRUE),
        "surv2" = Surv(runif(100), sample(0:1, 100, replace = TRUE)),
        "cont2" = runif(100),
        "cont3" = runif(100),
        "surv3" = Surv(runif(100), sample(0:1, 100, replace = TRUE)),
        "cat2"=factor(sample(1:5, 100, replace = TRUE), ordered = TRUE))

data <- mRMR.data(data = dd)
network <- new("mRMRe.Network", data = data, target_indices = c(1, 2), levels = c(2, 1), layers = 1)
network@topologies
adjacencyMatrix(network)
visualize(network)


data <- mRMR.data(data = dd,
        strata = factor(sample(1:5, 100, replace=TRUE), ordered=TRUE),
        weights = runif(100))
mim(data) # Gives MI matrix
mim(data, method = "cor") # Gives correlation matrix



mim(filter_1)

filter_2 <- mRMR.ensemble("mRMRe.Filter", data = data, target_index = 2, feature_count = 2, solution_count = 1)
mim(filter_2)

mim(data)

# No wrapper just yet
network <- new("mRMRe.Network", data = data, target_indices = c(1, 2), levels = c(1, 1, 1), layers = 1)
mim(data)
mim(network)

adjacencyMatrix(network)

causality(network)[1, , ]
causality(filter_1)

causality(network)[2, , ]
causality(filter_2)

# I have no idea how to get igraph to print out vertex names

visualize(network)

## test for large network (1000 genes)
## install new version of the package if needed
library(devtools)
install_github("mRMRe", username="bhaibeka", branch="master")
system("chmod -R 775 /stockage/Laboratoires/HAI/Rlib")


library(mRMRe)
## set the number of threads
set.thread.count(8)
## run the network inference
data(cgps)
ge <- mRMR.data(data = data.frame(cgps_ge[ ,1:1000]))
#Rprof(filename = "Rprof.out", append = FALSE, interval = 0.02, memory.profiling=TRUE)
exect <- system.time(netw <- new("mRMRe.Network", data = ge, target_indices = 1:10, levels = c(8, 1, 1, 1, 1), layers = 2))
print(exect)
#summaryRprof(filename = "Rprof.out", chunksize = 5000, memory=c("both"), index=2, diff=TRUE, exclude=NULL)


print(table(adjacencyMatrix(netw)))

pdf("temp.pdf", width=14, height=14)
visualize(netw)
dev.off()

