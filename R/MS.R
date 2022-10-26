MS <- function(X,Y){
  # Causal structure-based root cause analysis of outliers, ICML 2022
  
  require(glmnet)
  X = as.matrix(X)
  
  listL = DirectLiNGAM_regressCO(X) # Direct LiNGAM
  B = adaptive_log_lasso(listL$E,Y) # obtain coefficients for Y
  
  d = ncol(X)
  S = matrix(0,sum(Y==1),d)
  for (i in 1:d){
    if (d<0){
      S[,i] = get_Shap_exact_MS(listL$E,Y,B,i)
    } else{
      S[,i] = get_Shap_approx_MS(listL$E,Y,B,i)
    }
  }
  
  return(list(scores=S,order=1:d))
}

get_Shap_exact_MS <- function(E,Y,B,i){
  
  E = as.matrix(E)
  E1 = E[Y==1,,drop=FALSE]
  E0 = E[Y==0,,drop=FALSE]
  
  n = nrow(E1)
  rnd = sample(1:n,nrow(E0),replace=TRUE) # random sample
  # rnd = nrow(E1):1
  d = ncol(E1) 
  idx = setdiff(1:d,i) # exclude ith variable
  B = as.matrix(B)
  
  # k = 0
  E0t = E0
  S = -mean(E0t %*% B) # Example 2
  E0t[,i] = E1[rnd,i,drop=FALSE] # must draw from diseased population
  S = S + mean(E0t %*% B) # Example 2
  
  # k > 0
  for (k in 1:(d-1)){
    
    W = 1:k
    Sn = 0
    repeat{
      E0t = E0
      E0t[,idx[W]] = E1[rnd,idx[W],drop=FALSE] # randomly sample from diseased population
      Sn = Sn - mean(E0t %*% B) # Example 2
      E0t[,i] = E1[rnd,i,drop=FALSE] # sample ith index from diseased population
      Sn = Sn + mean(E0t %*% B) # Example 2
      
      nextSet = getNextSet(d-1, k, W) #
      if (nextSet$wasLast){  #
        break #
      }
      W <- nextSet$nextSet #
    }
    
    S = S + (1/choose(d-1,k))*Sn
    
  }
  
  return(S/d)
  
}


get_Shap_approx_MS <- function(E,Y,B,i,nsamp=200){
  
  E = as.matrix(E)
  E1 = E[Y==1,,drop=FALSE]
  E0 = E[Y==0,,drop=FALSE]
  
  n = nrow(E1)
  rnd = sample(1:n,nrow(E0),replace=TRUE) # random sample
  # rnd = nrow(E1):1
  d = ncol(E1) 
  idx = setdiff(1:d,i) # exclude ith variable
  B = as.matrix(B)
  
  S = 0
  for (j in 1:nsamp){
    k = sample(0:(d-1),1) # sample k uniformly with p = 1/(|C_i|+1)
    E0t = E0
    
    if (k>0){
      W = sample(idx,k,replace=FALSE) # sample random permutation of variables with k
      E0t[,W] = E1[rnd,W,drop=FALSE] # randomly sample from diseased population
    }
    
    S = S - mean(E0t %*% B) # equation right before Equation (10) in paper
    E0t[,i] = E1[rnd,i,drop=FALSE] # sample ith index from diseased population
    S = S + mean(E0t %*% B) # equation right before Equation (10) in paper
    
  }
  
  return(S/nsamp)
}

adaptive_log_lasso <- function(X,Y){
  
  B_ols = glm.fit(cbind(X,1),Y,family=binomial())$coefficients[1:ncol(X)]
  B = ada_log_lasso(X,Y,B_ols);
  
  return(B)
  
}

ada_log_lasso <- function(X,Y,B_ols){
  
  w = 1 / (abs(B_ols)^1)
  if (length(w)>1){
    X = X %*% diag(1/w)
    s = 2
  } else{
    X = X * (1/w)
    X = cbind(999,X)
    s = 3
  }
  cvfit = cv.glmnet(X,Y,family="binomial",type.measure = "class")
  B_lasso = c(as.matrix(coef(cvfit, s = "lambda.min")))
  return(B_lasso[s:length(B_lasso)])
}
