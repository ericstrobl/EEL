DirectLiNGAM_regressCO <- function(X){
  
  time = proc.time()
  
  K = c() #1 
  U = 1:ncol(X) #2
  X = normalizeData(X) #3
  Cov = cov(X) #4
  
  X0 = X
  Cov0 = Cov
  
  repeat{ #11
    
    root = FindRoot(X,U,Cov) #6
    
    K = c(K,root) #7
    U = U[-which(U==root)] #8
    
    if (length(U)==0){
      break
    }
    X = UpdateData(X,U,Cov,root) #9
    Cov = UpdateCovMat(U,Cov,root) #10
  }
  
  time = (proc.time()-time)[3]

  return(list(E = X, time = time)) #output
}

FindRoot <- function(X,U,Cov){
  
  r = length(U) #4
  
  if (r==1){ #1
    return(U) #2
  }
  
  p = ncol(Cov)
  # M = matrix(0,p,p)
  S = rep(0,r)
  for (i in 1:(r-1)){
    # V = U[U>i]
    for (j in (i+1):r){
      
      # M[i,j] = Compare2(X,i,j,Cov)
      # M[j,i] = -M[i,j]
      
      score = Compare2(X,U[i],U[j],Cov)
      S[i] = S[i] + min(0,score)^2
      S[j] = S[j] + min(0,-score)^2
    }
  }
  
  # S = rowSums(pmin(M[U,U],0)^2)
  
  # print(S)
  
  root = U[S==min(S)][1]
  
  return(root) #output
  
}

Compare2 <- function(X,i,j,Cov){
  
  if (i==j){  ####
    return(0) ####
  }
  
  rij = X[,i] - Cov[i,j]/Cov[j,j] * X[,j] #9
  rji = X[,j] - Cov[j,i]/Cov[i,i] * X[,i] #10
  
  rsd = sqrt(1-Cov[i,j]^2)
  rij = rij / rsd #11
  rji = rji / rsd #12
  
  score = LRT(X[,i],X[,j],rij,rji)
  
  return(score) #13  ####
}

LRT <- function(xi,xj,rij,rji){
  
  return(Hu(xj)+Hu(rij)-Hu(xi)-Hu(rji))
}

Hu <- function(u){
  
  k1 = 79.047
  k2 = 7.4129
  beta = 0.37457
  
  H = 0.5*(1+log(2*pi))-k1*(mean(log(cosh(u)))-beta)^2 - k2*mean(u*exp(-(u^2)/2))^2
  
  return(H)
}

UpdateData <- function(X,U,Cov,root){
  
  for (j in setdiff(U,root)){
    X[,j] = (X[,j] - Cov[j,root] / Cov[root,root] * X[,root,drop=FALSE])/
      sqrt(1-Cov[j,root]^2)
  }
  
  return(X)
  
}

UpdateCovMat <- function(U,Cov,root){
  
  for (j in setdiff(seq_len(ncol(Cov)),root)){
    Cov[U,j] = (Cov[U,j] - Cov[U,root] * Cov[root,j])/
      (sqrt(1-Cov[U,root]^2)*sqrt(1-Cov[j,root]^2)) 
  }
  
  return(Cov)
  
}

