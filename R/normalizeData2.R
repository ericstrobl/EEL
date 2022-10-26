normalizeData2 <- function(X){
  
  X = as.matrix(X)
  ms = c()
  sds = c() # record standard deviation for adjustment of total effects and Shapley values
  for (i in seq_len(ncol(X))){
    ms = c(ms, mean(X[,i]))
    sds = c(sds, sd(X[,i]))
    if (sds[i] == 0){
      X[,i] = X[,i] - ms[i]
    } else{
      X[,i] = (X[,i] - ms[i])/sds[i]
    }
  }
  X = as.matrix(X)
  
  return(list(X=X,ms=ms,sds=sds))
}
