compute_ground_truth <- function(DAG,E,Y,sds){
  
  
  # construct directed error graph, *error terms* connected to variables by directed inducing path
  p = ncol(DAG$graph)
  dG = matrix(0,p,p)
  for (i in 1:p){
    Ancs = isAnc(DAG$graph,1:p,i)
    dG[er_IP_to_b(DAG$graph,DAG$L,Ancs,i),i ] = 1
  }
  # print(dG)
  # diag(dG)=0

  # print(dG)
  # get coefficients
  p = ncol(DAG$graph)
  A = diag(sds)%*%solve(diag(p) - DAG$weights)

  # compute E_star values with directed dependence graph dG
  E_star = E
  for (i in 1:p){
    # get Estar
    ix = which(dG[,i]==1)
    E_star[,i] = E[,ix,drop=FALSE] %*% A[ix,i,drop=FALSE]
  }
  E_star = E_star[,-c(DAG$L,DAG$Y)] # remove latents and target

  # convert directed dependence graph dG to dependence graph G
  G = dG # find parents
  # now find children of parents
  for (i in 1:ncol(G)){
    pa = which(dG[,i]>0)
    for (p in pa){
      ch = which(dG[p,]>0)
      G[ch,i] = 1
    }
  }
  
  # compute Shapley values with dependence graph G
  
  betas = lm.fit(cbind(E_star,1),Y)$coefficients[1:ncol(E_star)]
  betas[which(abs(betas)<1E-8)] = 0 # set basically zero to zero
  
  # print(dG)
  # compute Shapley with linear trees
  require(Cubist)
  S = matrix(NaN,nrow(E_star),ncol(E_star))
  diag(G)=0
  G = G[-c(DAG$L,DAG$Y),-c(DAG$L,DAG$Y)] # remove latents and target
  for (i in 1:ncol(E_star) ){
    # if (i==DAG$Y) next
    nn = sum(G[,i])
    if (nn==0){
      S[,i] = E_star[,i] * betas[i] # A[i,DAG$Y] is beta_{iY}
    } else{
      S[,i] = E_star[,i] * betas[i]
      if (nn < 10){ # exhaustive Shapley ####
        for (k in 1:nn){ #k = 0 not needed because E(E_i) = 0
          nn_k = 1/choose(nn,k)
          # print(c("ab",k))
          # print(sum_all(E[,i],E[,G[i,],drop=FALSE],k)[1:5])
          S[,i] = S[,i] - betas[i]*(nn_k/(nn+1))*sum_all_cubist(E_star[,i],E_star[,G[,i]>0,drop=FALSE],k) # average over all expectations with conditioning set size k
        }
      } else{
        S[,i] = stochastic_Shapley_cubist(E_star[,i],E_star[,G[,i]>0,drop=FALSE],betas[i])
      }
    }
  }
    
  return(S)
}

er_IP_to_b <- function(graph,Ls,Ancb,b,visited=rep(FALSE,nrow(graph))){
  
  visited[b] = TRUE
  
  if (b %in% Ls){ #if last vertex was not observed
    adj = setdiff( intersect(which(graph[b,]>0),Ancb), which(visited) ) # then next vertex must be observed
  } else{ #if last vertex was observed
    adj = setdiff( intersect(which(graph[,b]>0),Ls), which(visited) ) # then next vertex must be latent
  }
  
  for (j in adj){
    visited = IP_to_b(graph,Ls,Ancb,j,visited)
  }
  
  return(visited)
}


isAnc <- function(graph,candidates,target){
  
  Anc = c()
  for (c in candidates){
    if (isAnc_fast_LE(graph,c,target)){
      Anc = c(Anc, c)
    }
  }
  return(Anc)
}

isAnc_fast_LE <- function(graph,a,b,visited=rep(FALSE,nrow(graph)))
{
  
  if (a %in% b){
    return(TRUE)
  }
  
  visited[a] = TRUE;
  
  adj = which(graph[a,] & !visited);
  
  out=FALSE;
  for (j in adj){
    out=isAnc_fast_LE(graph, j, b, visited);
    if(out==TRUE){
      break;
    }
  }
  
  return(out)
  
}


sum_all_cubist <- function(E1,E1s,k){
  E1s = as.matrix(E1s) # make sure predictors in a matrix
  L = ncol(E1s)
  
  W = seq_len(k) # first W
  Ef = rep(0,length(E1)) # final expectation
  # d = 0
  repeat{
    # d=d+1
    # print(d)
    mod = cubist(x=as.data.frame(E1s[,W,drop=FALSE]),y=E1, committees=10)
    Ef = Ef + Cubist:::predict.cubist(mod,as.data.frame(E1s[,W,drop=FALSE])) # sum over expectations
    
    nextSet = getNextSet(L, k, W) #
    if (nextSet$wasLast){  #
      break #
    }
    W <- nextSet$nextSet #
  }
  
  return(Ef)
}

stochastic_Shapley_cubist <- function(Ei, Eg, betai, nsamp = 1000){
  
  Si = Ei * betai
  nn = ncol(Eg)
  Ef = rep(0,length(Si))
  for (j in 1:nsamp){
    k = sample(0:nn,1) # sample k uniformly with p = 1/(|C_i|+1)
    if (k>0){
      W = sample(1:nn,k,replace=FALSE) # sample random permutation of variables with k
      mod = cubist(x=as.data.frame(Eg[,W,drop=FALSE]),y=Ei, committees=10)
      Ef = Ef + Cubist:::predict.cubist(mod,as.data.frame(Eg[,W,drop=FALSE]))
    }
  }
  Ef = (betai*Ef)/nsamp
  
  return(Si-Ef)
}
