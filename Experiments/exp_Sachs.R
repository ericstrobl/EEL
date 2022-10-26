data = read.csv('Sachs_data.csv') 
for (c in 12:20){
  ix1 = which(data[,c]>=0.5)
  data[ix1,c] = 1
  ix0 = which(data[,c]<0.5)
  data[ix0,c] = 0
  
  for (c1 in 1:11){
    data[ix1,c1] = (data[ix1,c1] - mean(data[ix1,c1]))/sd(data[ix1,c1])
    data[ix0,c1] = (data[ix0,c1] - mean(data[ix0,c1]))/sd(data[ix0,c1])
    # data[,c1] = lm.fit(cbind(1,data[,c]),data[,c1])$residuals # regress out intervention nodes
  }
}
err = as.matrix(data[,1:11])

Gt = read.csv('Sachs_graph.csv'); colnames(Gt)[1] = "raf"
Gt = as.matrix(Gt)
rownames(Gt) = colnames(Gt)

G$graph = Gt
G$weights = Gt
for (c in 1:11){
  pa = which(Gt[,c]>0)
  if (length(pa)>0){
    mod = lm.fit(as.matrix(cbind(data[,pa],1)),data[,c])

    G$weights[pa,c] = mod$coefficients[1:length(pa)]
    err[,c] = mod$residuals # regress out intervention nodes
  }
}

ix = which(abs(G$weights)<0.01)
G$weights[ix] = 0
G$graph[ix] = 0
plot(as(G$graph,"graphNEL"))

reps = 1000

EEL_res = vector("list",reps)
RCI_res = EEL_res
GRCI_res = EEL_res
ICA_res = EEL_res
RCAO_res = EEL_res
MS_res = EEL_res

r = nrow(err)
L1 = which(colSums(G$graph)==0)
Y1 = (1:11)[-L1]
for (ii in 20:reps){
  print(ii)
  
  G$Y = sample(Y1,sample(1:3,1)) # potential latents
  G$L = L1
  
  ######
  G$L = sample(G$L,sample(1:length(G$L),1)) # sample target
  
  samps = sample_DAG_lat_err(r,G,as.matrix(err))
  G$Y = 12
  
  nD = normalizeData2(samps$data) # normalize to prevent gaming of variances

  # compute ground truth
  
  ms = rep(0,11); ms[-c(G$Y,G$L)] = nD$ms[-c(G$Y,G$L)]
  sds = rep(1,11); sds[-c(G$Y,G$L)] = nD$sds[-c(G$Y,G$L)]
  S = compute_ground_truth(G,sweep(samps$E,2,ms)%*% diag(1/sds),samps$Y0,sds) # ground truth, adjusted for normalization
  
  X = samps$data[,-c(G$Y,G$L)]
  Y = samps$data[,G$Y]
  #######
  
  ix = which(Y==1)
  
  # EEL
  ptm <- proc.time()
  out = EEL(X+1E-8*rnorm(length(X)),Y)
  out$scores = out$scores[ix,]
  EEL_res[[ii]]$time = (proc.time() - ptm)[3]
  EEL_res[[ii]]$rank_overlap = eval_scores(S[ix,],out)
  EEL_res[[ii]]$MSE = eval_L2(S[ix,],out)
  
  # RCI
  ptm <- proc.time()
  out = RCI(X,Y)
  out$scores = out$scores[ix,]
  RCI_res[[ii]]$time = (proc.time() - ptm)[3]
  RCI_res[[ii]]$rank_overlap = eval_scores(S[ix,],out)
  RCI_res[[ii]]$MSE = eval_L2(S[ix,],out)
  
  # GRCI
  ptm <- proc.time()
  out = GRCI(X+1E-8*rnorm(length(X)),Y)
  out$scores = out$scores[ix,]
  GRCI_res[[ii]]$time = (proc.time() - ptm)[3]
  GRCI_res[[ii]]$rank_overlap = eval_scores(S[ix,],out)
  GRCI_res[[ii]]$MSE = eval_L2(S[ix,],out)
  
  # ICA
  ptm <- proc.time()
  out = ICA_predict(X,Y)
  out$scores = out$scores[ix,]
  ICA_res[[ii]]$time = (proc.time() - ptm)[3]
  ICA_res[[ii]]$rank_overlap = eval_scores(S[ix,],out)
  ICA_res[[ii]]$MSE = eval_L2(S[ix,],out)
  
  # RCAO
  ptm <- proc.time()
  out = RCAO(X,Y)
  RCAO_res[[ii]]$time = (proc.time() - ptm)[3]
  RCAO_res[[ii]]$rank_overlap = eval_scores(S[ix,],out)
  RCAO_res[[ii]]$MSE = eval_L2(S[ix,],out)
  
  # MS
  ptm <- proc.time()
  out = MS(X,Y)
  MS_res[[ii]]$time = (proc.time() - ptm)[3]
  MS_res[[ii]]$rank_overlap = eval_scores(S[ix,],out)
  MS_res[[ii]]$MSE = eval_L2(S[ix,],out)
  
  
  save(file="Res_Sachs2.RData", EEL_res, RCI_res, GRCI_res, 
       ICA_res, RCAO_res, MS_res)
  
}

### PBC RESULTS
RBO = matrix(0,reps,6)
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
MSE = matrix(0,reps,6)
for (i in 1:reps){
  MSE[i,1] = EEL_res[[i]]$MSE
  MSE[i,2] = RCI_res[[i]]$MSE
  MSE[i,3] = GRCI_res[[i]]$MSE
  MSE[i,4] = ICA_res[[i]]$MSE
  MSE[i,5] = RCAO_res[[i]]$MSE
  MSE[i,6] = MS_res[[i]]$MSE
}

print(colMeans(MSE))


### PBC TIMING
time = matrix(0,reps,6)
for (i in 1:reps){
  time[i,1] = EEL_res[[i]]$time
  time[i,2] = RCI_res[[i]]$time
  time[i,3] = GRCI_res[[i]]$time
  time[i,4] = ICA_res[[i]]$time
  time[i,5] = RCAO_res[[i]]$time
  time[i,6] = MS_res[[i]]$time
}
print(colMeans(time))
