#' Root Causal Inference (RCI) algorithm.
#'
#' @param X a matrix of continuous variables, rows are samples and columns are variables
#' @param Y a vector of a binary target
#' @param alpha the alpha value for the t-test, default is 0.2
#' @return A list containing Shapley values \code{scores}, the errors \code{E}, the variable ordering \code{order}, and time of Local Plus \code{time}
#' @export

RCI <- function(X,Y,alpha=0.2){
  time = proc.time()
  outL = DirectLiNGAM_fast_Y(X,Y) # Local Plus
  # outL = DirectLiNGAM(X)
  time = (proc.time() - time)[3]
  K = outL$K;
  if (length(K)==0){
    delta = matrix(0,nrow(X),ncol(X))
    return(list(delta=delta,K=K))
  }
  E = outL$X[,order(K),drop=FALSE]
  K = sort(K)
  
  ## patient-specific root causes
  beta = glm.fit(cbind(E,1),Y,family=binomial())$coefficients[1:length(K)]
  # print(beta)
  if (length(beta)==1){
    scores = E * beta
  } else{
    scores = E %*% diag(beta)
  }
  
  colnames(scores) = K
  colnames(E) = K
  
  return(list(scores=scores, order=K, E=E, time = time))
}

DirectLiNGAM_fast_Y <- function(X,Y,alpha=0.2){
  X = as.matrix(X)
  
  K = c() #1 
  X = normalizeData(X) #3
  Cov = cov(X) #4
  
  U = 1:ncol(X) #2
  
  repeat{ #11
    
    U = U[ttest_fast(X[,U,drop=FALSE],Y,alpha=alpha)] ###
    if (length(U)==0){ ###
      break ###
    } ###
    
    root = FindRoot_fast(X,U,Cov) #6
    K = c(K,root) #7
    U = U[-which(U==root)] #8
    
    if (length(U)==0){ ###
      break ###
    } ###
    X = UpdateData(X,U,Cov,root) #9
    Cov = UpdateCovMat(U,Cov,root) #10
  }
  
  return(list(K=K,X=X[,K,drop=FALSE])) #output
}

UpdateCovMat <- function(U,Cov,root){
  
  for (j in setdiff(seq_len(ncol(Cov)),root)){
    Cov[U,j] = (Cov[U,j] - Cov[U,root] * Cov[root,j])/
      (sqrt(1-Cov[U,root]^2)*sqrt(1-Cov[j,root]^2)) 
  }
  
  return(Cov)
  
}

UpdateData <- function(X,U,Cov,root){
  
  for (j in setdiff(U,root)){
    X[,j] = (X[,j] - Cov[j,root] / Cov[root,root] * X[,root,drop=FALSE])/
      sqrt(1-Cov[j,root]^2)
  }
  
  return(X)
  
}

FindRoot_fast <- function(X,U,Cov){
  
  r = length(U) #4
  
  if (r==1){ #1
    return(U) #2
  }
  
  O = rep(1,r)
  S = rep(0,r)
  M = matrix(TRUE,r,r)
  
  repeat{
    I = (S == min(S))
    
    if ( sum(O[I] > r) > 0){
      break
    }
    
    for (i in which(I)){
      if (M[i,O[i]]){
        score = Compare2(X,U[i],U[O[i]],Cov)
        S[i] = S[i] + min(0,score)^2
        S[O[i]] = S[O[i]] + min(0,-score)^2
        
        M[i,O[i]]=FALSE
        M[O[i],i]=FALSE
      }
      
      O[i] = O[i]+1
    }
  }
  
  root = U[S==min(S)][1]
  
  return(root) #output
  
}

ttest_fast <- function(X,Y,suffStat=NULL,alpha=0.2){
  
  X0 = X[Y==0,,drop=FALSE]
  X1 = X[Y==1,,drop=FALSE]
  
  ms = colMeans(X0) - colMeans(X1)
  
  n0 = nrow(X0)
  n1 = nrow(X1)
  
  var0 = apply(X0,2,var)
  var1 = apply(X1,2,var)
  
  t = ms/ sqrt(var0/n0 + var1/n1)
  
  df = ((var0 + var1)^2) / ((var0^2)/(n0-1) + (var1^2)/(n1-1))
  
  ps = (1-pt(abs(t),df=df))*2
  
  return( which(ps<alpha) )
  
  
  
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
