library(grasps)

## reproducibility for everything
set.seed(1234)

## block-structured precision matrix based on SBM
sim <- gen_prec_sbm(d = 100, K = 10,
                    within.prob = 0.2, between.prob = 0.05,
                    weight.dists = list("gamma", "unif"),
                    weight.paras = list(c(shape = 100, scale = 10),
                                        c(min = 0, max = 5)),
                    cond.target = 100)

## visualization
plot(sim)
