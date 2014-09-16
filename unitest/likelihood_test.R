# Given a data_set and network, computes barcode for each edge
# Returns a matrix of 3 columns where col 1 is the higher value
# col 2 is the lower value, and col 3 is the frequency of such a
# relation.
load('~/Testbed/irinotecan_cgp_ccle.RData')
`compute.barcode` <- function(data_set, network)
{
	edges <- which(network == 1)
	barcode <- do.call(rbind, lapply(edges, function(edge) {
						j <- edge%%ncol(network)
						if(j == 0)
							j <- ncol(network)           
						i <- ceiling(edge/ncol(network))
						c(i, j, mean(data_set[, i] > data_set[, j], use="complete.obs"))
					}))
	return(barcode)
}

# Returns the product of likelihoods for a given barcode and sample
`compute.likelihood` <- function(barcode, sample)
{
	likelihood <- apply(barcode, 1, function(feature) {
				if (sample[feature[1]] > sample[feature[2]])
					return(feature[3])
				else
					return(1 - feature[3])
			})
	return(prod(likelihood))
}

`discretize.labels` <- function(labels, probabilities)
{
	quantiles <- quantile(labels, probs=probabilities, na.rm=TRUE)
	classes <- Hmisc::cut2(x=labels, cuts=quantiles)
	levels(classes) <- c(0, NA, 1)
	names(classes) <- names(labels)
	as.integer(classes) - 1
}



training_set <- data_cgp
training_labels <- ic50_cgp
validating_set <- data_cgp


# Select the top 10 genes most correlated with phenotype
genes <- order(apply(training_set, 2, cor, training_labels, use="complete.obs"))[1:10]

# Create empty network
network <- matrix(0, ncol(training_set), ncol(training_set))

# Connect the 10 selected genes (every pair of nodes)
sapply(1:(length(genes) - 1), function(i) {
			sapply((i + 1):length(genes), function(j) {
						network[genes[i], genes[j]] <<- 1
					})
		})

# Discretize CGP into senstive (0) and resistant (1)
discrete_labels <- discretize.labels(training_labels, c(0.25, 0.75))

# Create barcodes for resistant and senstive cgp samples
resistant_barcode <- compute.barcode(training_set[which(discrete_labels == 1), ], network)
sensitive_barcode <- compute.barcode(training_set[which(discrete_labels == 0), ], network)

# Compute the likelihoods
likelihoods <- t(apply(validating_set, 1, function(sample) {
			resistant_like <- compute.likelihood(resistant_barcode, sample)
			sensitive_like <- compute.likelihood(sensitive_barcode, sample)
			return(c(resistant_like, sensitive_like))
		}))

hist(likelihoods[, 1] / ( likelihoods[, 1] + likelihoods[, 2]))

