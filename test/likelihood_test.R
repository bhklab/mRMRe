# Given a data_set and network, computes barcode for each edge
# Returns a matrix of 3 columns where col 1 is the higher value
# col 2 is the lower value, and col 3 is the frequency of such a
# relation.
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

# Select the top 10 genes most correlated with phenotype
genes <- order(apply(data_cgp, 2, cor, ic50_cgp, use="complete.obs"))[1:10]

# Create empty network
network <- matrix(0, ncol(data_cgp), ncol(data_cgp))

# Connect the 10 selected genes (every pair of nodes)
sapply(1:length(genes), function(i) {
			sapply(i:length(genes), function(j) {
						network[genes[i], genes[j]] <<- 1
					})
		})

# Discretize CGP into senstive (0) and resistant (1)
cgp_discrete <- discretize.labels(ic50_cgp, c(0.25, 0.75))

# Create barcodes for resistant and senstive cgp samples
resistant_barcode <- compute.barcode(data_cgp[which(cgp_discrete == 1), ], network)
sensitive_barcode <- compute.barcode(data_cgp[which(cgp_discrete == 0), ], network)

# Compute the likelihoods
ccle_likelihoods <- t(apply(data_ccle_common, 1, function(sample) {
			resistant_like <- compute.likelihood(resistant_barcode, sample)
			sensitive_like <- compute.likelihood(sensitive_barcode, sample)
			return(c(resistant_like, sensitive_like))
		}))

# Compute concordance-index for likelihoods
cor(ccle_likelihoods[,1], ic50_ccle_common, use="complete.obs", method="spearman")
cor(ccle_likelihoods[,2], ic50_ccle_common, use="complete.obs", method="spearman")


