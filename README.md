# Extract Errors with Latents (EEL)

EEL is an algorithm that discovers sample-specific root causes even when confounding exists among the predictors. The algorithm does not require knowledge about the underlying causal graph. EEL also identifies the sample-specific total effects of each root cause when theory allows; otherwise, it outputs a measure of predictivity, and these two cases are unified under Shapley values. EEL is designed to handle a noisy binary target -- like a diagnostic label -- and is therefore particularly suitable for identifying root causes of disease.

# Installation

> library(devtools)

> install_local("Directory of EEL-main.zip")

> library(RCI)

# Run the Sample Version

Instantiate number of variables p, number of samples n, number of latents L:
> p = 15; n = 10000; L = runif(1)*20

Generate the grouth truth DAG:
> G = generate_DAG_lat(p=p,en=2,perc_lat = L)

Sample DAG and normalize to prevent gaming of variances:
> X = sample_DAG_lat(nsamps = n,G)

> nD = normalizeData2(X$data) 

Remove target and latents from data:
> X$data[,-c(G$Y,G$L)] = nD$X[,-c(G$Y,G$L)]

Run EEL:
> out = EEL(X$data[1:n,-c(G$Y,G$L)],X$data[1:n,G$Y])

Print first 5 Shapley values:
> ix = intersect(1:n,which(X$data[,G$Y]==1))

> print(out$scores[ix[1:5],])


# Run the Oracle Version

Instantiate number of variables p, number of samples n, number of latents L:
> p = 25; n = 200; L = runif(1)*20

Generate the grouth truth DAG:
> G = generate_DAG_lat(p=p,en=2,perc_lat = L)
 
Sample DAG and normalize. The oracle needs some samples to compute regression residuals:
> X = sample_DAG_lat(nsamps = n,G)

> nD = normalizeData2(X$data) 

Remove target and latents from data:
> X$data[,-c(G$Y,G$L)] = nD$X[,-c(G$Y,G$L)]
 
Get all the info required by the oracle including ground truth total effects, ground truth errors, index of target, index of latents:
> oracle = organize_oracle(G,X,p,L)

Run EEL:
> out = EEL(X$data[,-c(G$Y,G$L)],X$data[,G$Y],oracle=oracle) 
 
Compute ground truth dependence graph of inducing terms:
> G_theory = ground_truth_dependence_graph(G)

Print number of errors, should be zero:
> print(sum(out$G != G_theory[-c(G$Y,G$L),-c(G$Y,G$L)])) 
