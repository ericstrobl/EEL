normalizeData <- function(X){
  
  X = as.matrix(X)
  ms = c()
  sds = c() # record standard deviation for adjustment of total effects and Shapley values
  for (i in seq_len(ncol(X))){
    sd = sd(X[,i])
    if (sd == 0){
      X[,i] = X[,i] - mean(X[,i])
    } else{
      X[,i] = (X[,i] - mean(X[,i]))/sd
    }
  }
  X = as.matrix(X)
  
  return(X)
}