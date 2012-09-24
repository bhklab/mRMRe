

library(mRMRe)
data(cgps)
ge <- mRMR.data(data = data.frame(cgps_ge))

zabou <- function(f, t)
{

## set the number of threads
set.thread.count(t)
## run the network inference
a <- mRMR.ensemble(data = ge, target_index = 1, solution_count = 10, feature_count = f)
t(solutions(a)) -> Z


nnn <- apply(Z, 2, function(i)   apply(Z, 2, function(j) {  length(intersect(i, j))   } )   )
diag(nnn) <- NA


list(glove=Z, stats=nnn)
}
