compute_Shapley <- function(E,G,beta){
  
  require(earth)
  
  S = E
  for (i in 1:ncol(E) ){
    nn = sum(G[i,])
    if (nn==0){
      S[,i] = E[,i]*beta[i]
    } else{
      if (nn < 10){ # exhaustive Shapley ####
        S[,i] = E[,i]*beta[i]
        for (k in 1:nn){ #k = 0 not needed because E(E_i) = 0
          nn_k = 1/choose(nn,k)
          # print(c("ab",k))
          # print(sum_all(E[,i],E[,G[i,],drop=FALSE],k)[1:5])
          S[,i] = S[,i] - beta[i]*(nn_k/(nn+1))*sum_all(E[,i],E[,G[i,],drop=FALSE],k) # average over all expectations with conditioning set size k
        }
      } else{
        S[,i] = stochastic_Shapley(E[,i],E[,G[i,],drop=FALSE],beta[i])
      }
    }
  }
  
  return(S)
}

sum_all <- function(E1,E1s,k,regressor){
  E1s = as.matrix(E1s) # make sure predictors in a matrix
  L = ncol(E1s)
  
  W = seq_len(k) # first W
  Ef = rep(0,length(E1)) # final expectation
  # d = 0
  repeat{
    # d=d+1
    # print(d)
    Ef = Ef + earth(x=E1s[,W,drop=FALSE],y=E1)$fitted.values # sum over expectations
    
    nextSet = getNextSet(L, k, W) #
    if (nextSet$wasLast){  #
      break #
    }
    W <- nextSet$nextSet #
  }
  
  return(Ef)
}


stochastic_Shapley <- function(Ei, Eg, betai, nsamp = 1000){
  
  Si = Ei * betai
  nn = ncol(Eg)
  Ef = rep(0,length(Si))
  for (j in 1:nsamp){
    k = sample(0:nn,1) # sample k uniformly with p = 1/(|C_i|+1)
    if (k>0){
      W = sample(1:nn,k,replace=FALSE) # sample random permutation of variables with k
      Ef = Ef + earth(x=Eg[,W,drop=FALSE],y=Ei)$fitted.values
    }
  }
  Ef = (betai*Ef)/nsamp
  
  return(Si-Ef)
}
