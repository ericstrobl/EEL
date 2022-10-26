ground_truth_dependence_graph <- function(DAG){

  # construct directed error graph, *error terms* connected to variables by directed inducing path
  p = ncol(DAG$graph)
  dG = matrix(0,p,p)
  for (i in 1:p){
    Ancs = isAnc(DAG$graph,1:p,i)
    dG[er_IP_to_b(DAG$graph,DAG$L,Ancs,i),i ] = 1 #er_IP_to_b is only for *error term* on directed inducing path
  }
  # print(dG)
  # diag(dG)=0
  
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

  G = G>0
  diag(G) = FALSE
  return(G)
}