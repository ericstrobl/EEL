GRCI <- function(X,Y,outL=NULL){
  require(RANN)
  require(gtools)
  
  p = ncol(X)
  if (is.null(outL)){
    X = normalizeData(X)
    
    ## SKELETON DISCOVERY
    # print('skeleton')
    # require(pcalg)
    # suffStat = list()
    # suffStat$data = X
    # G = pc(suffStat, earth_wrap, alpha=0.10, p=ncol(X))
    # G = as(G@graph, "matrix")
    # G = ((G + t(G))>0)
    G = matrix(TRUE,p,p)
    diag(G) = FALSE
    
    ## EXTRACT ERRORS
    # print('errors')
    
    outL = DirectHNM_fast_Y(X,Y,G)
  }
  
  ## LOGISTIC REGRESSION
  ## patient-specific root causes
  beta = glm.fit(cbind(outL$X,1),Y,family=binomial())$coefficients[1:ncol(X)]
  # print(beta)
  if (length(beta)==1){
    scores = outL$X * beta
  } else{
    scores = outL$X %*% diag(beta)
  }
  
  return(list(E=outL$X, order=1:ncol(X), scores=scores, G=G))
}

CompareG <- function(X,i,js,w=NULL){
  
  if (identical(i,js)){  ####
    return(0) ####
  }
  
  # rXY = PartialOut(X[,js],X[,i],w[js,,drop=FALSE])
  rXY = lm.fit(X[,js,drop=FALSE],X[,i])$residuals
  
  score = c()
  for (k in 1:length(js)){
    score = c(score, hoeffding.D.test(X[,js[k]],rXY,precision=1)$scale)
  }
  score = max(score)
  
  return(score) #13  ####
}

DirectHNM_fast_Y <- function(X,Y,G,alpha=0.2){
  X = as.matrix(X)
  # E = X
  
  K = c()
  S = rep(-Inf,ncol(X))
  U = 1:ncol(X)
  update = U
  
  # penalty = 1
  # w = rdirichlet(200,rep(1,ncol(X))*penalty) # random projections
  # # w = normalize_rep(w)
  # w = t(w)
  
  repeat{ 
    s_out = FindSink(X,U,S,G,update)
    sink = s_out$sink
    S = s_out$S
    
    K = c(K,sink)
    U = U[-which(U==sink)]
    
    if (length(U)==0){ ###
      break ###
    } ###
    
    if (sum(G[sink,])){
      X[,sink] = PartialOut(X[,G[sink,]],X[,sink],w[G[sink,],,drop=FALSE]) # partial out neighbors from sink
    }
    update = intersect(U,which(G[sink,]))
    G[sink,] = FALSE; G[,sink]=FALSE # remove node
    
  }
  
  return(list(K=K,X=X)) #output
}


FindSink <- function(X,U,S,G,update,w=NULL){
  r = length(update)
  
  ## BASELINE SCORES
  for (i in seq_len(r)){
    V = which(G[update[i],])
    if (length(V)==0){ # if no neighbors, then break because score is -Inf
      next
    }
    
    S[update[i]] = CompareG(X,update[i],V,w) # baseline
    
  }
  
  sink = U[S[U]==min(S[U])][1] 
  
  return(list(sink=sink,S=S)) #output
  
}

