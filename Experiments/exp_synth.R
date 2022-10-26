library(xgboost)
library(pcalg)
library(DirichletReg)
library(Rfast)
library(RANN)
library(treeshap)
library(RCI)
library(independence)

ns = c(1000, 10000, 100000)
ps = 15
Ls = c(0,10,20)
reps = 500

Gs = lapply(1:reps, function (.) lapply(1:length(ps), function(.) lapply(1:length(Ls),function(.) vector("list",length(ns)))))
EEL_res = Gs
RCI_res = Gs
GRCI_res = Gs
ICA_res = Gs
RCAO_res = Gs
MS_res = Gs

for (i in 1:200){
  print(i)
  for (p in 1:length(ps)){
    for (L in 1:length(Ls)){
      
      # generate DAG and sample data
      G = generate_DAG_lat(p=ps[p],en=2,perc_lat = Ls[L])
      plot(as(G$graph,"graphNEL"))
      # print(G$L); print(G$Y)
      Gs[[i]][[p]][[L]][[1]]$G = G
      
      X = sample_DAG_lat(nsamps = 100000,G)
      nD = normalizeData2(X$data) # normalize to prevent gaming of variances
      X$data[,-c(G$Y,G$L)] = nD$X[,-c(G$Y,G$L)]
      
      # compute ground truth
      
      ms = rep(0,ps[p]); ms[-c(G$L,G$Y)] = nD$ms[-c(G$L,G$Y)]
      sds = rep(1,ps[p]); sds[-c(G$L,G$Y)] = nD$sds[-c(G$L,G$Y)]
      S = compute_ground_truth(G,sweep(X$E,2,ms)%*% diag(1/sds),X$Y0,sds) # ground truth, adjusted for normalization

      for (n in 1:length(ns)){
        
          ix = intersect(1:ns[n],which(X$data[,G$Y]==1))
          cx = (1:ps[p])[-c(G$Y,G$L)]
          
          # EEL
          ptm <- proc.time()
          out = EEL(X$data[1:ns[n],-c(G$Y,G$L)],X$data[1:ns[n],G$Y])
          out$scores = out$scores[ix,]
          # colnames(out$G) <- cx[out$order]; rownames(out$G) <- cx[out$order]
          # print(out$G)
          # colnames(S) <- cx
          EEL_res[[i]][[p]][[L]][[n]]$time = (proc.time() - ptm)[3]
          EEL_res[[i]][[p]][[L]][[n]]$rank_overlap = eval_scores(S[ix,],out)
          EEL_res[[i]][[p]][[L]][[n]]$MSE = eval_L2(S[ix,],out)
          
          # RCI
          ptm <- proc.time()
          out = RCI(X$data[1:ns[n],-c(G$Y,G$L)],X$data[1:ns[n],G$Y])
          out$scores = out$scores[ix,]
          RCI_res[[i]][[p]][[L]][[n]]$time = (proc.time() - ptm)[3]
          RCI_res[[i]][[p]][[L]][[n]]$rank_overlap = eval_scores(S[ix,],out)
          RCI_res[[i]][[p]][[L]][[n]]$MSE = eval_L2(S[ix,],out)
          
          # GRCI
          ptm <- proc.time()
          out = GRCI(X$data[1:ns[n],-c(G$Y,G$L)],X$data[1:ns[n],G$Y])
          out$scores = out$scores[ix,]
          GRCI_res[[i]][[p]][[L]][[n]]$time = (proc.time() - ptm)[3]
          GRCI_res[[i]][[p]][[L]][[n]]$rank_overlap = eval_scores(S[ix,],out)
          GRCI_res[[i]][[p]][[L]][[n]]$MSE = eval_L2(S[ix,],out)
          
          # ICA
          ptm <- proc.time()
          out = ICA_predict(X$data[1:ns[n],-c(G$Y,G$L)],X$data[1:ns[n],G$Y])
          out$scores = out$scores[ix,]
          ICA_res[[i]][[p]][[L]][[n]]$time = (proc.time() - ptm)[3]
          ICA_res[[i]][[p]][[L]][[n]]$rank_overlap = eval_scores(S[ix,],out)
          ICA_res[[i]][[p]][[L]][[n]]$MSE = eval_L2(S[ix,],out)
          
          # RCAO
          ptm <- proc.time()
          out = RCAO(X$data[1:ns[n],-c(G$Y,G$L)],X$data[1:ns[n],G$Y])
          RCAO_res[[i]][[p]][[L]][[n]]$time = (proc.time() - ptm)[3]
          RCAO_res[[i]][[p]][[L]][[n]]$rank_overlap = eval_scores(S[ix,],out)
          RCAO_res[[i]][[p]][[L]][[n]]$MSE = eval_L2(S[ix,],out)
          
          # MS
          ptm <- proc.time()
          out = MS(X$data[1:ns[n],-c(G$Y,G$L)],X$data[1:ns[n],G$Y])
          MS_res[[i]][[p]][[L]][[n]]$time = (proc.time() - ptm)[3]
          MS_res[[i]][[p]][[L]][[n]]$rank_overlap = eval_scores(S[ix,],out)
          MS_res[[i]][[p]][[L]][[n]]$MSE = eval_L2(S[ix,],out)
          
          save(file="Results_synth15.RData", Gs, EEL_res, RCI_res, GRCI_res,
               ICA_res, RCAO_res, MS_res)
      }
    }
  }
}

## RBO
RBO_EEL = lapply(1:length(Ls),function(.) vector("list",length(ns)))
RBO_RCI = RBO_EEL
RBO_GRCI = RBO_EEL
RBO_ICA = RBO_EEL
RBO_RCAO = RBO_EEL
RBO_MS = RBO_EEL

RBO_EELm = matrix(0,length(Ls),length(ns))
RBO_RCIm = RBO_EELm
RBO_GRCIm = RBO_EELm
RBO_ICAm = RBO_EELm
RBO_RCAOm = RBO_EELm
RBO_MSm = RBO_EELm

imax = 120
for (i in 1:imax){
  for (L in 1:3){
    for (n in 1:3){
      RBO_EEL[[L]][[n]] = c(RBO_EEL[[L]][[n]], EEL_res[[i]][[p]][[L]][[n]]$rank_overlap)
      RBO_RCI[[L]][[n]] = c(RBO_RCI[[L]][[n]], RCI_res[[i]][[p]][[L]][[n]]$rank_overlap)
      RBO_GRCI[[L]][[n]] = c(RBO_GRCI[[L]][[n]], GRCI_res[[i]][[p]][[L]][[n]]$rank_overlap)
      RBO_ICA[[L]][[n]] = c(RBO_ICA[[L]][[n]], ICA_res[[i]][[p]][[L]][[n]]$rank_overlap)
      RBO_RCAO[[L]][[n]] = c(RBO_RCAO[[L]][[n]], RCAO_res[[i]][[p]][[L]][[n]]$rank_overlap)
      RBO_MS[[L]][[n]] = c(RBO_MS[[L]][[n]], MS_res[[i]][[p]][[L]][[n]]$rank_overlap)
      
      RBO_EELm[L,n] = RBO_EELm[L,n] + EEL_res[[i]][[p]][[L]][[n]]$rank_overlap
      RBO_RCIm[L,n] = RBO_RCIm[L,n] + RCI_res[[i]][[p]][[L]][[n]]$rank_overlap
      RBO_GRCIm[L,n] = RBO_GRCIm[L,n] + GRCI_res[[i]][[p]][[L]][[n]]$rank_overlap
      RBO_ICAm[L,n] = RBO_ICAm[L,n] + ICA_res[[i]][[p]][[L]][[n]]$rank_overlap
      RBO_RCAOm[L,n] = RBO_RCAOm[L,n] + RCAO_res[[i]][[p]][[L]][[n]]$rank_overlap
      RBO_MSm[L,n] = RBO_MSm[L,n] + MS_res[[i]][[p]][[L]][[n]]$rank_overlap
    }
  }
}

RBO_EELm/imax
RBO_RCIm/imax
RBO_GRCIm/imax
RBO_ICAm/imax
RBO_RCAOm/imax
RBO_MSm/imax


## MSE
MSE_EEL = lapply(1:length(Ls),function(.) vector("list",length(ns)))
MSE_RCI = MSE_EEL
MSE_GRCI = MSE_EEL
MSE_ICA = MSE_EEL
MSE_RCAO = MSE_EEL
MSE_MS = MSE_EEL

MSE_EELm = matrix(0,length(Ls),length(ns))
MSE_RCIm = MSE_EELm
MSE_GRCIm = MSE_EELm
MSE_ICAm = MSE_EELm
MSE_RCAOm = MSE_EELm
MSE_MSm = MSE_EELm

imax = 120
for (i in 1:imax){
  for (L in 1:3){
    for (n in 1:3){
      MSE_EEL[[L]][[n]] = c(MSE_EEL[[L]][[n]], EEL_res[[i]][[p]][[L]][[n]]$MSE)
      MSE_RCI[[L]][[n]] = c(MSE_RCI[[L]][[n]], RCI_res[[i]][[p]][[L]][[n]]$MSE)
      MSE_GRCI[[L]][[n]] = c(MSE_GRCI[[L]][[n]], GRCI_res[[i]][[p]][[L]][[n]]$MSE)
      MSE_ICA[[L]][[n]] = c(MSE_ICA[[L]][[n]], ICA_res[[i]][[p]][[L]][[n]]$MSE)
      MSE_RCAO[[L]][[n]] = c(MSE_RCAO[[L]][[n]], RCAO_res[[i]][[p]][[L]][[n]]$MSE)
      MSE_MS[[L]][[n]] = c(MSE_MS[[L]][[n]], MS_res[[i]][[p]][[L]][[n]]$MSE)
      
      MSE_EELm[L,n] = MSE_EELm[L,n] + EEL_res[[i]][[p]][[L]][[n]]$MSE
      MSE_RCIm[L,n] = MSE_RCIm[L,n] + RCI_res[[i]][[p]][[L]][[n]]$MSE
      MSE_GRCIm[L,n] = MSE_GRCIm[L,n] + GRCI_res[[i]][[p]][[L]][[n]]$MSE
      MSE_ICAm[L,n] = MSE_ICAm[L,n] + ICA_res[[i]][[p]][[L]][[n]]$MSE
      MSE_RCAOm[L,n] = MSE_RCAOm[L,n] + RCAO_res[[i]][[p]][[L]][[n]]$MSE
      MSE_MSm[L,n] = MSE_MSm[L,n] + MS_res[[i]][[p]][[L]][[n]]$MSE
    }
  }
}

MSE_EELm/imax
MSE_RCIm/imax
MSE_GRCIm/imax
MSE_ICAm/imax
MSE_RCAOm/imax
MSE_MSm/imax

## time
time_EEL = lapply(1:length(Ls),function(.) vector("list",length(ns)))
time_RCI = time_EEL
time_GRCI = time_EEL
time_ICA = time_EEL
time_RCAO = time_EEL
time_MS = time_EEL

time_EELm = matrix(0,length(Ls),length(ns))
time_RCIm = time_EELm
time_GRCIm = time_EELm
time_ICAm = time_EELm
time_RCAOm = time_EELm
time_MSm = time_EELm

imax = 120
for (i in 1:imax){
  for (L in 1:3){
    for (n in 1:3){
      time_EEL[[L]][[n]] = c(time_EEL[[L]][[n]], EEL_res[[i]][[p]][[L]][[n]]$time)
      time_RCI[[L]][[n]] = c(time_RCI[[L]][[n]], RCI_res[[i]][[p]][[L]][[n]]$time)
      time_GRCI[[L]][[n]] = c(time_GRCI[[L]][[n]], GRCI_res[[i]][[p]][[L]][[n]]$time)
      time_ICA[[L]][[n]] = c(time_ICA[[L]][[n]], ICA_res[[i]][[p]][[L]][[n]]$time)
      time_RCAO[[L]][[n]] = c(time_RCAO[[L]][[n]], RCAO_res[[i]][[p]][[L]][[n]]$time)
      time_MS[[L]][[n]] = c(time_MS[[L]][[n]], MS_res[[i]][[p]][[L]][[n]]$time)
      
      time_EELm[L,n] = time_EELm[L,n] + EEL_res[[i]][[p]][[L]][[n]]$time
      time_RCIm[L,n] = time_RCIm[L,n] + RCI_res[[i]][[p]][[L]][[n]]$time
      time_GRCIm[L,n] = time_GRCIm[L,n] + GRCI_res[[i]][[p]][[L]][[n]]$time
      time_ICAm[L,n] = time_ICAm[L,n] + ICA_res[[i]][[p]][[L]][[n]]$time
      time_RCAOm[L,n] = time_RCAOm[L,n] + RCAO_res[[i]][[p]][[L]][[n]]$time
      time_MSm[L,n] = time_MSm[L,n] + MS_res[[i]][[p]][[L]][[n]]$time
    }
  }
}

time_EELm/imax
time_RCIm/imax
time_GRCIm/imax
time_ICAm/imax
time_RCAOm/imax
time_MSm/imax
