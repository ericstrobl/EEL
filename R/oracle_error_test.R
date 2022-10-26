oracle_error_test <- function(xis,yi,oracle){
  
  Ge_tot = oracle$Ge_tot
  errors = oracle$errors
  
  Ex = matrix(0,nrow(errors),length(xis)) # one error for each in xis
  ixs = c()
  for (xi in seq_len(length(xis))){
    ix = which(abs(Ge_tot[,xis[xi]])>1E-8) # find parental errors
    ixs = c(ixs, ix) # take union of parental error indices
    # print( Ge_tot[ix,xis[xi],drop=FALSE]  )
    Ex[,xi] = errors[,ix,drop=FALSE] %*% Ge_tot[ix,xis[xi],drop=FALSE] # compute Ex
  }

  
  ixs = unique(ixs)
  iy = which(abs(Ge_tot[,yi])>1E-15) # find parental errors
  
  iy = intersect(iy,ixs) # only use overlapping errors, so mse is basically zero
  # print(iy)
  # print(Ge_tot[iy,yi,drop=FALSE])
  
  Ey = errors[,iy,drop=FALSE] %*% Ge_tot[iy,yi,drop=FALSE]
  
  metric = mean(lm.fit(Ex,Ey)$residuals^2)
  # if (yi == 9){
  #   print(xis)
  #   print(metric)
  #   # print(lm.fit(Ex,Ey)$coefficients)
  # }
  
  if (metric <1E-15){
    return(p=1)
  } else{
    return(p=0)
  }
  
  
}