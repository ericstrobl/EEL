data = read.csv('diabetes.csv')
data = as.matrix(data[complete.cases(data),])

Xo = data[,1:8]
X = normalizeData(Xo)
Y = data[,9]

Gt = matrix(0,8,8)
Gt[8,c(3,1,2)]=1
Gt[7,c(6,2,5)]=1
Gt[6,c(3,4)] = 1
Gt[2,5] = 1

G = list()
G$graph = Gt
G$weights = Gt
err = X
for (c in 1:8){
  pa = which(Gt[,c]>0)
  if (length(pa)>0){
    mod = lm.fit(as.matrix(cbind(data[,pa],1)),data[,c])
    
    G$weights[pa,c] = mod$coefficients[1:length(pa)]
    err[,c] = mod$residuals # regress out intervention nodes
  }
}

reps = 200

EEL_res = vector("list",reps)
RCI_res = EEL_res
GRCI_res = EEL_res
ICA_res = EEL_res
RCAO_res = EEL_res
MS_res = EEL_res

r = nrow(err)
L1 = which(colSums(G$graph)==0)
G$Y = 9

for (ii in 1:reps){
  print(ii)
  
  ######
  G$L = sample(L1,sample(1:length(L1),1)) # sample target
  
  n = nrow(X)
  ib = sample(1:n,n,replace=TRUE)
  Xn = X[ib,]
  errn = err[ib,]
  Yn = Y[ib]
  
  nD = normalizeData2(cbind(Xn,Yn)) # normalize to prevent gaming of variances
  
  # compute ground truth
  
  ms = rep(0,8); ms[-c(G$Y,G$L)] = nD$ms[-c(G$Y,G$L)]
  sds = rep(1,8); sds[-c(G$Y,G$L)] = nD$sds[-c(G$Y,G$L)]
  S = compute_ground_truth_diabetes(G,sweep(errn,2,ms)%*% diag(1/sds),Yn,sds) # ground truth, adjusted for normalization
  
  Xn = nD$X[,-c(G$Y,G$L)]
  #######
  
  ix = which(Yn==1)
  
  # EEL
  ptm <- proc.time()
  out = EEL(Xn+1E-8*rnorm(length(Xn)),Yn)
  out$scores = out$scores[ix,]
  EEL_res[[ii]]$time = (proc.time() - ptm)[3]
  EEL_res[[ii]]$rank_overlap = eval_scores(S[ix,],out)
  EEL_res[[ii]]$MSE = eval_L2(S[ix,],out)
  
  # RCI
  ptm <- proc.time()
  out = RCI(Xn,Yn)
  out$scores = out$scores[ix,]
  RCI_res[[ii]]$time = (proc.time() - ptm)[3]
  RCI_res[[ii]]$rank_overlap = eval_scores(S[ix,],out)
  RCI_res[[ii]]$MSE = eval_L2(S[ix,],out)
  
  # GRCI
  ptm <- proc.time()
  out = GRCI(Xn+1E-8*rnorm(length(Xn)),Yn)
  out$scores = out$scores[ix,]
  GRCI_res[[ii]]$time = (proc.time() - ptm)[3]
  GRCI_res[[ii]]$rank_overlap = eval_scores(S[ix,],out)
  GRCI_res[[ii]]$MSE = eval_L2(S[ix,],out)
  
  # ICA
  ptm <- proc.time()
  out = ICA_predict(Xn,Yn)
  out$scores = out$scores[ix,]
  ICA_res[[ii]]$time = (proc.time() - ptm)[3]
  ICA_res[[ii]]$rank_overlap = eval_scores(S[ix,],out)
  ICA_res[[ii]]$MSE = eval_L2(S[ix,],out)
  
  # RCAO
  ptm <- proc.time()
  out = RCAO(Xn,Yn)
  RCAO_res[[ii]]$time = (proc.time() - ptm)[3]
  RCAO_res[[ii]]$rank_overlap = eval_scores(S[ix,],out)
  RCAO_res[[ii]]$MSE = eval_L2(S[ix,],out)
  
  # MS
  ptm <- proc.time()
  out = MS(Xn,Yn)
  MS_res[[ii]]$time = (proc.time() - ptm)[3]
  MS_res[[ii]]$rank_overlap = eval_scores(S[ix,],out)
  MS_res[[ii]]$MSE = eval_L2(S[ix,],out)
  
  
  save(file="Res_diabetes.RData", EEL_res, RCI_res, GRCI_res, 
       ICA_res, RCAO_res, MS_res)
  
}

### PBC RESULTS
RBO = matrix(0,reps1,6)
for (i in 1:reps){
  RBO[i,1] = EEL_res[[i]]$rank_overlap
  RBO[i,2] = RCI_res[[i]]$rank_overlap
  RBO[i,3] = GRCI_res[[i]]$rank_overlap
  RBO[i,4] = ICA_res[[i]]$rank_overlap
  RBO[i,5] = RCAO_res[[i]]$rank_overlap
  RBO[i,6] = MS_res[[i]]$rank_overlap
}
print(colMeans(RBO))


### PBC RESULTS
time = matrix(0,reps1,6)
for (i in 1:reps){
  time[i,1] = EEL_res[[i]]$time
  time[i,2] = RCI_res[[i]]$time
  time[i,3] = GRCI_res[[i]]$time
  time[i,4] = ICA_res[[i]]$time
  time[i,5] = RCAO_res[[i]]$time
  time[i,6] = MS_res[[i]]$time
}
print(colMeans(time))

