## Size

library(mRMRe)
data(cgps)
data <- mRMR.data(data = as.data.frame(cgps.ge))
system.time(filter <- mRMR.ensemble("mRMRe.Filter", data = data, target_indices = c(1), feature_count = 200, solution_count = 160))
print(object.size(filter), units = "Mb")


library(mRMRe)
load('~/Downloads/irinotecan_cgp_ccle.RData')
data <- mRMR.data(data = data.frame(data_cgp[, 1:10000]))
system.time(mim(data))

## Simple Test

library(mRMRe)
set.seed(0)

x <- rnorm(100, 0)
dd <- data.frame(
        "cont1" = x,
        "cont2" = x + rnorm(100, 0.1),
		"cont3" = x + rnorm(100, 0.1),
		"cont4" = x + rnorm(100, 0.1),
		"cont5" = x + rnorm(100, 0.1),
		"cont6" = x + rnorm(100, 0.01))
		
data <- mRMR.data(data = dd)
filter <- mRMR.classic("mRMRe.Filter", data = data, target_indices = 3:5, feature_count = 2)
scores(filter)
solutions(filter)

mim(filter)
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
## install new version of the package if needed on gen01
library(devtools)
install_github("mRMRe", username="bhaibeka", ref="master")
system("chmod -R 775 /stockage/Laboratoires/HAI/Rlib")


library(mRMRe)
## set the number of threads
set.thread.count(8)
## run the network inference
data(cgps)
ge <- mRMR.data(data = data.frame(cgps.ge[ ,1:1000]))
#Rprof(filename = "Rprof.out", append = FALSE, interval = 0.02, memory.profiling=TRUE)
exect <- system.time(netw <- new("mRMRe.Network", data = ge, target_indices = 1:10, levels = c(8, 1, 1, 1, 1), layers = 2))
print(exect)
#summaryRprof(filename = "Rprof.out", chunksize = 5000, memory=c("both"), index=2, diff=TRUE, exclude=NULL)


print(table(adjacencyMatrix(netw)))

pdf("temp.pdf", width=14, height=14)
visualize(netw)
dev.off()


# Generate the basic stuff
data(cgps)
dd <- mRMR.data(data=data.frame(cgps.ge))
netw <- new("mRMRe.Network", data = dd, target_indices = c(1, 2), levels = c(1, 1, 1), layers = 1)
filter <- new("mRMRe.Filter", data = dd, target_indices = c(1, 2), levels = c(1,1,1))

# Create all 3 mim and remove NA's
network_mim <- mim(netw)
indices <- !is.na(network_mim)
network_mim <- network_mim[indices]
filter_mim <- mim(filter)[indices]
data_mim <- mim(dd)[indices]

# Compare
quantile(network_mim)
quantile(filter_mim)
quantile(data_mim)

# Causality
causality(filter) # a list of causality scores for each solution
causality(filter, merge=T) # a single vector of causality scores that is the min value between all solutions
causality(netw) # network-wide minimal causality scores

