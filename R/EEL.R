EEL <- function (X, D, indepTest=hoeffding.D.test, alpha=0.05, oracle=NULL)
{
  
  require(Rfast)
  require(independence)
  
  # X = normalizeData(X) #
  
  p = ncol(X) 
  G <- matrix(TRUE,p,p)
  diag(G) <- FALSE
  
  iTest = lapply(seq_len(p), function(.) vector("list", p))
  
  out = prune_step(X,G,iTest,alpha,indepTest,oracle=oracle) # neighbor prune
  
  beta = glm.fit(cbind(out$X,1),D,family=binomial())$coefficients[1:ncol(X)]
  S = compute_Shapley(out$X,out$G,beta)
  # S = 0 ###
  
  return(list(E=out$X,order=1:p,scores=S,beta=beta,G=out$G,oracle=out$oracle))
}

prune_step <-function(X,G,iTest,alpha,indepTest,Gt=NULL,oracle=NULL){
  
  p = ncol(X)
  
  l = 0
  repeat{ # first loop
    
    l = l + 1 
    
    inds <- which(rowSums(G)>=l)
    
    if (!is.null(Gt)){
      inds <- intersect(which(rowSums(Gt)>=1),inds) # has at least one TRUE in Gt
    }
    
    Yf = c(); Vf = c()
    while(length(inds)>0){ # can be performed in parallel?
      ind = inds[1] # pop inds
      inds = inds[-1] # pop inds
      
      nbrs = which(G[ind,])
      length_nbrs <- length(nbrs)# for each neighbor of ind
      if (is.null(Gt)){
        ni = nbrs
      } else{
        ni = which(Gt[ind,])
        if (length(intersect(ni,nbrs))==0) next
      }
      W = getNextW(length_nbrs, l, nbrs, ni, W=NULL)$nextSet # find W containing at least one member of ni
      repeat{ # third loop
        
        Y = nbrs[W]
        
        IV = getIndepVars(X[,Y,drop=FALSE],X[,ind,drop=FALSE],Y,ind,indepTest,alpha,iTest,oracle) # residuals in rem independent of Y
        pmax = IV$p
        V = ind[IV$V]
        iTest = IV$iTest
        
        if (length(V)>0){
          Yf = Y; Vf = V
          break
        }
        
        nextSet = getNextW(length_nbrs, l, nbrs, ni, W) # find W containing at least one member of ni
        if (nextSet$wasLast) break 
        W <- nextSet$nextSet 
      }
      
      if (length(Vf)>0){
        # print(c(oracle$ix[Yf],oracle$ix[Vf]))
        
        X[,Vf] = partial_out(X[,Yf,drop=FALSE],X[,Vf,drop=FALSE]) #
        
        if (!is.null(oracle)){
          ixx = c() ########
          for (yy in Yf){ ########
            ixt = which(oracle$Ge_tot[,oracle$ix[yy]]!=0) ########
            ixx = c(ixx, ixt) ########
          } ########
          oracle$Ge_tot[ixx,oracle$ix[Vf]] = 0 ###### remove those errors from Vf
        }
        
        iTest = reset_iTest(Vf,G,iTest) # must reset before removing adjacencies, is it this one?
        
        l = 0
        if (is.null(Gt)){ # if fixed set of variables to condition on
          G[Vf,Yf] = FALSE; G[Yf,Vf] = FALSE # remove adjacencies
        } else{
          Gt[Vf,Yf] = FALSE; # conditioning sets in G dont change
        }
        break
      }
      
      
    }
    if (is.null(Gt)){
      if ( sum(rowSums(G)<l)==p ){ # all vertices have less than or equal to l adjacencies
        break
      }
    } else{
      ll = intersect(which(rowSums(G)>=l), which(rowSums(Gt)>=1))
      if ( length(ll)==0 ){ # all vertices have less than l adjacencies or no adjacencies in Gt (or both)
        break
      }
    }
  }
  
  if (!is.null(Gt)) G = Gt
  
  return(list(G=G,X=X,iTest = iTest, oracle=oracle))
}

getNextW <- function(length_nbrs, l, nbrs, ni, W=NULL){
  
  if( is.null(W)){
    nextSet = list()
    nextSet$nextSet = seq_len(l)
    nextSet$wasLast = FALSE
    if( length(intersect(ni,nbrs[nextSet$nextSet]))>0 ){ #if set contains at least one member of ni
      return(nextSet)
    } else{
      W = nextSet$nextSet
    }
  }
  
  repeat{ # find set containing at least one member of ni
    nextSet = getNextSet(length_nbrs, l, W)
    if ( nextSet$wasLast ){  
      break 
    } else if( length(intersect(ni,nbrs[nextSet$nextSet]))>0 ){
      break
    }
    W <- nextSet$nextSet
  }
  
  return(nextSet)
  
}

getNextSet <- function (n, k, set) 
{
  chInd <- k - (zeros <- sum((seq(n - k + 1, n) - set) == 0))
  wasLast <- (chInd == 0)
  if (!wasLast) {
    set[chInd] <- s.ch <- set[chInd] + 1
    if (chInd < k) 
      set[(chInd + 1):k] <- seq(s.ch + 1L, s.ch + zeros)
  }
  list(nextSet = set, wasLast = wasLast)
}


ttest_fast <- function(X,Y,alpha=0.2){
  
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


partial_out <- function(X,Y){
  
  Y = Y - X %*% (spdinv(t(X) %*% X +1E-10*diag(ncol(X))) %*% t(X) %*% Y)
  
  return(as.matrix(Y))
  
}

getIndepVars <- function(Y, Rem, Yi, Remi, indepTest, alpha, iTest, oracle=NULL){
  V = c()
  R = partial_out(Y,Rem)
  pmaxs = c()
  for (j in 1:length(Remi)){ # for each residual
    p = 1
    iy = intersect(Remi,Yi)
    for (k in seq_len(length(Yi))){ # for each predictor
      pn = check_iTest(Yi[k],Remi[j],setdiff(Yi,Yi[k]),iTest)
      if (length(pn)==0){
        if (!is.null(oracle)){
          pn = oracle_error_test(oracle$ix[Yi],oracle$ix[Remi[j]],oracle) ###*****
        } else{
          pn = indepTest(Y[,k],R[,j],precision=(alpha/2))$p #precision = 2e-3
        }
        # if ( oracle$GE[oracle$ix[Yi[k]],oracle$ix[Remi[j]]]==1 ){
        #   pn = 0 
        # }###***** refinement
        # pn = oracle_error_test(oracle$ix[Yi],oracle$ix[Remi[j]],oracle) ###*****
        iTest[[Yi[k]]][[Remi[j]]] = append(iTest[[Yi[k]]][[Remi[j]]],list(c(pn-1,setdiff(Yi,Yi[k]))) ) # subtract one from new p-value so no confusion with indices
        p = min(p,pn) # find residual independent of predictor
      } else{
        p = min(p,pn)
      }
    }
    if (p > alpha){ # if residual independent of all predictors
      pmaxs = c(pmaxs, p)
      V = c(V, j)
    }
  }
  
  if (length(V)>0){
    return(list(V=V[which(pmaxs==max(pmaxs))[1]],pval=max(pmaxs),iTest = iTest)) ## just one V with max p-value
  } else{
    return(list(V=V,pval=0,iTest = iTest))
  }
  
}

check_iTest <- function(k,j,oth,iTest){
  p = c()
  for (i in seq_len(length(iTest[[k]][[j]]))){ # for each conditioning set
    L = length(iTest[[k]][[j]][[i]])
    if ( L>0 ){
      iSec = intersect(oth,iTest[[k]][[j]][[i]][2:L]) # first index is p-value
      if (length(iSec)==length(oth)){ # if oth and iTest[[k]][[j]][[i]] are the same
        p = iTest[[k]][[j]][[i]][1] + 1 # then test already run, so get the p-value, adding back one
        # print(p)
        return(p)
      }
    }
  }
  return(p)
}

reset_iTest <- function(Vf,G,iTest){
  
  for (i in Vf){ # predictor
    iTest[[i]] = vector("list", length(iTest)) #
  }
  
  for (i in seq_len(length(iTest))){ # target
    for (j in Vf){
      iTest[[i]][j] = list(NULL) #
    }
  }
  
  for (i in which(G[Vf,])){ # in set of other predictors
    for (j in seq_len(length(iTest))){
      L = length(iTest[[i]][[j]])
      for (k in seq_len(L)){
        if (length(intersect(Vf,iTest[[i]][[j]][[k]]))>0){
          iTest[[i]][[j]][k] = list(NULL) 
        }
      }
      if (L>0){
        # iTest[[i]][[j]] <- iTest[[i]][[j]][!sapply(iTest[[i]][[j]],is.null)]
        iTest[[i]][[j]][sapply(iTest[[i]][[j]], is.null)] <- NULL
      }
    }
  }
  
  return(iTest)
  
}
