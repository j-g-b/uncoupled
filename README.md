# uncoupled
R package implementing minimum-Wasserstein isotonic regression estimator from Rigollet and Weed (2019)

# Implementation details

To implement [Rigollet and Weed](https://arxiv.org/abs/1806.10648)'s relaxed minimum-Wasserstein estimator, we use projected subgradient descent on the Wasserstein-2 objective. This amounts to a two-step procedure where evaluation of the objective at each iteration also yields subgradients that are combined to update the true measure weights. The deconvolved measure weights are then updated with a subgradient step and projected onto the probability simplex.