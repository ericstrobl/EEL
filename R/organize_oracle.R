organize_oracle <- function(G,X,p,L){
  
  oracle=list();
  oracle$Ge_tot = solve(diag(p) - G$weights)
  oracle$Ge_tot[which(abs(oracle$Ge_tot)<1E-8)]=0
  oracle$errors = X$E
  oracle$ix = (1:p)[-c(G$Y,G$L)]
  oracle$yi = G$Y
  
  return(oracle)
}
