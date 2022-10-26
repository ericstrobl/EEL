generate_DAG_lat <- function(p,en,perc_lat){
  
  N = p*p - p;
  
  DAG = list()
  
  nL = p*perc_lat/100
  if (nL != round(nL)){
    sw = sample(0:1,1)
    if (sw == 0){
      nL = floor(nL)
    } else{
      nL = ceiling(nL)
    }
  }
  repeat{
    samplesB = rbinom(N/2,1, en/(p-1) ); # sample edges
    graph = matrix(0,p,p)
    graph[upper.tri(graph, diag=FALSE)] <- samplesB; # put in edges in upper triangular
  
    ord = sample(1:p,p,replace=FALSE) # permute order
    DAG$graph = graph[ord,ord]
    
    Ys = which( (rowSums(DAG$graph)==0) & (colSums(DAG$graph)>0) ) # no children, some observed parents
    # print((colSums(DAG$graph[DAG$L,,drop=FALSE])==0))
    DAG$Y = Ys[sample(length(Ys),1)]
    
    pL = which( (rowSums(DAG$graph)>=2) & (colSums(DAG$graph)==0) & (DAG$graph[,DAG$Y]==0) ); #variables with >=2 children and no parents, target is not a child
    AncY = isAnc(DAG$graph,(1:p)[-DAG$Y],DAG$Y)
    pL = intersect(pL,AncY) # possible latent must also be an ancestor of target
    # print(length(pL))
    DAG$L = pL[sample(length(pL), min(length(pL),nL))]
    
    if ((length(pL)<nL) | length(Ys)==0){
      next
    } else{
      # DAG$L = pL[sample(length(pL), nL, replace=FALSE)]
      break
    }
  }
  
  weights = matrix((0.75*runif(p^2)+0.25)*sample(c(-1,1),p^2,replace=TRUE),p,p)
  DAG$weights = weights*DAG$graph
  
  DAG$errors = sample(3,p,replace=TRUE)
  
  DAG$idx = setdiff(1:ncol(DAG$graph),DAG$L)
  
  return(DAG)
}
