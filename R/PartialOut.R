PartialOut <- function(X,Y,w=NULL){
  
  Y = as.matrix(Y)
  X = as.matrix(X)
  
  Res = lm.fit(X,Y)$residuals
  
  return(Res)
}
