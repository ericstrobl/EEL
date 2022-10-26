# Extract Errors with Latents (EEL)

EEL is an algorithm that discovers sample-specific root causes even when confounding exists among the predictors. The algorithm also identifies the sample-specific total effects of each root cause. EEL is designed to handle a noisy binary target -- like a diagnostic label -- and is therefore particularly suitable for identifying root causes of disease.

# Installation

> library(devtools)

> install_local("Directory of EEL-main.zip")

> library(RCI)

# Run the Sample Version

# Run the Oracle Version

Instantiate number variables p, number samples n, number latents L:
> p = 25; n = 200; L = runif(1)*20

Generate the grouth truth DAG:
> G = generate_DAG_lat(p=p,en=2,perc_lat = L)
 
> X = sample_DAG_lat(nsamps = n,G); nD = normalizeData2(X$data) # samle DAG and normalize

> X$data[,-c(G$Y,G$L)] = nD$X[,-c(G$Y,G$L)] # remove target and latents
 
> oracle = organize_oracle(G,X,p,L) # get all the info required by oracle

> out = EEL(X$data[,-c(G$Y,G$L)],X$data[,G$Y],oracle=oracle) # run EEL 
 
> G_theory = ground_truth_dependence_graph(G) # compute ground truth dependence graph of inducing terms

> print(sum(out$G != G_theory[-c(G$Y,G$L),-c(G$Y,G$L)])) # print number of errors, should be zero
