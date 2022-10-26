# Extract Errors with Latents (EEL)

EEL is an algorithm that discovers sample-specific root causes even when confounding exists among the predictors. The algorithm also identifies the sample-specific total effects of each root cause. EEL is designed to handle a noisy binary target -- like a diagnostic label -- and is therefore particularly suitable for identifying root causes of disease.

# Installation

> library(devtools)

> install_local("Directory of EEL-main.zip")

> library(RCI)

# Run the Sample Version

# Run the Oracle Version

Instantiate number of variables p, number of samples n, number of latents L:
> p = 25; n = 200; L = runif(1)*20

Generate the grouth truth DAG:
> G = generate_DAG_lat(p=p,en=2,perc_lat = L)
 
Sample DAG and normalize. The oracle needs some samples to compute regression residuals:
> X = sample_DAG_lat(nsamps = n,G)

> nD = normalizeData2(X$data) 

Remove target and latents from data
> X$data[,-c(G$Y,G$L)] = nD$X[,-c(G$Y,G$L)]
 
Get all the info required by the oracle including ground truth total effects, ground truth errors, index of target, index of latents:
> oracle = organize_oracle(G,X,p,L)

Run EEL:
> out = EEL(X$data[,-c(G$Y,G$L)],X$data[,G$Y],oracle=oracle) 
 
Compute ground truth dependence graph of inducing terms:
> G_theory = ground_truth_dependence_graph(G)

Print number of errors, should be zero
> print(sum(out$G != G_theory[-c(G$Y,G$L),-c(G$Y,G$L)])) 
