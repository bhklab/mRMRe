`compute.causality` <- function(
        data,
        target_index,
        mim,
        solutions,
        estimator=c("pearson", "spearman", "kendall"))
{
    if (class(data) == "mRMReObject")
        return(mRMRe:::.compute.causality.mRMReObject(data=data))

    if (missing(mim))
        allcor <- as.matrix(mim)
    else
        allcor <- matrix(ncol=ncol(data), nrow=ncol(data))
        
    estimator <- match.arg(estimator)
    
    causality_coefficients <- matrix(ncol=ncol(allcor), nrow=ncol(allcor))
    apply(solutions, 1, function(row)
    {
        pairs <- combn(row, 2)
        apply(pairs, 2, function(pair)
        {
            i <- pair[1]
            j <- pair[2]
            
            if (is.na(causality_coefficients[i, j]))
            {
                if (is.na(allcor[i, target_index]))
                {
                    r <- cor(data[,i], data[,target_index], method=estimator)
                    allcor[i, target_index] <<- r
                    allcor[target_index, i] <<- r
                }
                
                if (is.na(allcor[i, j]))
                {
                    r <- cor(data[,i], data[,j], method=estimator)
                    allcor[i, j] <<- r
                    allcor[j, i] <<- r
                }
                
                if (is.na(allcor[j, target_index]))
                {
                    r <- cor(data[,j], data[,target_index], method=estimator)
                    allcor[target_index, j] <<- r
                    allcor[j, target_index] <<- r
                }
                
                causality <- -1/2 * log(((1 - allcor[i, j]^2) * (1 - allcor[i, target_index]^2)
                                    * (1 - allcor[j, target_index]^2)) / (1 + 2 * allcor[i, j] * allcor[i, target_index] 
                                    * allcor[j, target_index] - allcor[i, j]^2 - allcor[i, target_index]^2 - 
                                    allcor[j, target_index]^2))
                
                causality_coefficients[i, j] <<- causality
                causality_coefficients[j, i] <<- causality
            }
        })
    })
    return(causality_coefficients)
}

`.compute.causality.mRMReObject` <- function(data)
{
    tree <- data
    target_index <- 1
    causality_coefficients <- matrix(ncol=ncol(tree$mim), nrow=ncol(tree$mim))
    apply(tree$paths, 1, function(row)
    {
        pairs <- combn(row, 2)
        apply(pairs, 2, function(pair)
        {
            i <- pair[1]
            j <- pair[2]
 
            if (is.na(causality_coefficients[i, j]))
            {
                causality <- -1/2 * log(((1 - tree$mim[i, j]^2) * (1 - tree$mim[i, target_index]^2)
                                    * (1 - tree$mim[j, target_index]^2)) / (1 + 2 * tree$mim[i, j] * tree$mim[i, target_index] 
                                    * tree$mim[j, target_index] - tree$mim[i, j]^2 - tree$mim[i, target_index]^2 - 
                                    tree$mim[j, target_index]^2))
                causality_coefficients[i, j] <<- causality
                causality_coefficients[j, i] <<- causality
            }
        })
    })
    return(causality_coefficients)
}