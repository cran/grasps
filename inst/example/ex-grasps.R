library(grasps)

## reproducibility for everything
set.seed(1234)

## block-structured precision matrix based on SBM
sim <- gen_prec_sbm(d = 30, K = 3,
                    within.prob = 0.25, between.prob = 0.05,
                    weight.dists = list("gamma", "unif"),
                    weight.paras = list(c(shape = 20, rate = 10),
                                        c(min = 0, max = 5)),
                    cond.target = 100)
## visualization
plot(sim)

## n-by-d data matrix
library(MASS)
X <- mvrnorm(n = 20, mu = rep(0, 30), Sigma = sim$Sigma)

## adapt, HBIC
res <- grasps(X = X, membership = sim$membership, penalty = "adapt", crit = "HBIC")

## visualization
plot(res)

## performance
performance(hatOmega = res$hatOmega, Omega = sim$Omega)
