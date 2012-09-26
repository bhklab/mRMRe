library(mRMRe)

##
## Tests
##

dd <- data.frame("surv1" = Surv(runif(100), sample(0:1, 100, replace = TRUE)),
                 "cont1" = runif(100),
                 "cat1"  = factor(sample(1:5, 100, replace = TRUE), ordered = TRUE),
                 "surv2" = Surv(runif(100), sample(0:1, 100, replace = TRUE)),
                 "cont2" = runif(100),
                 "cont3" = runif(100),
                 "surv3" = Surv(runif(100), sample(0:1, 100, replace = TRUE)),
                 "cat2"  = factor(sample(1:5, 100, replace = TRUE), ordered = TRUE))

data <- mRMR.data(data = dd)
filter_1 <- mRMR.ensemble("mRMRe.Filter", data = data, target_indices = sample(1:8, 4, replace = FALSE),
        feature_count = 2, solution_count = 2)
filter_2 <- mRMR.ensemble("mRMRe.Filter", data = data, target_indices = c(1, 2), feature_count = 2, solution_count = 3)

mi_master <- mim(data)
diag(mi_master) <- 0

mi_slave_list <- list()
mi_slave_list[[1]] <- mim(filter_1)
mi_slave_list[[2]] <- mim(filter_2)

for (slave_index in seq(length(mi_slave_list)))
{
    diag(mi_slave_list[[slave_index]]) <- 0
    mi_slave <- as.vector(mi_slave_list[[slave_index]])
    
    mi_slave_list[[slave_index]] <- 0
    
    for (feature in seq(mi_slave))
    {
        if (!is.na(mi_slave[[feature]]))
            mi_slave_list[[slave_index]] <- mi_slave_list[[slave_index]] +
                    (abs(mi_slave[[feature]] - mi_master[[feature]]) < 1e-5)
    }
    
    mi_slave_list[[slave_index]] <- length(mi_slave) - mi_slave_list[[slave_index]]
}

error_rate <- mean(as.integer(mi_slave_list)) / length(as.numeric(mi_master))

print(paste("Error rate: ", error_rate))

# FIXME: Add tests on network object