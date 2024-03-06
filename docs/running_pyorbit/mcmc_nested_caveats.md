(mcmc_nested_caveats)=

# Caveats regarding MCMC and NS

`PyORBIT` is one of the few Exoplanetary codes that supports both Monte Carlo Markov Chain (MCMC) and Nested Sampling (NS), thus allowing a direct comparison of the two approaches over the same datasets and models. In my experience, when the model is appropriate, the two algorithms deliver the *same* results, i.e., the posteriors are perfectly superimposed. Although it does not guarantee your analysis's correctness, it is indeed a useful check. It is essential to understand the different implementations of the two algorithms, as they affect how priors 

## Priors on parameters

In MCMC, a prior acts as a penalisation to the log-likelihood: the farther a parameter from the prior, the more negative (with the other parameters frozen) the log-likelihood will be. This is true for both *sampler* parameters and *derived* parameters. When using MCMC, **it is possible to set a prior on a derived parameter**. For example, it is possible to set a prior on the eccentricity when using the Eastman parametrisation, even if eccentricity is not a sampler parameter.

In NS, the prior is actually a reparametrisation of the unitary cube, in such a way that a uniform sampling of the unitary parameter will result into a distribution of the sampler parameter that follows the given prior. As a consequence, **priors can be imposed only on sampler parameters**. For example, you cannot use the Eastman parametrization and a prior on the eccentricity simultaneously.  

TO BE COMPLETED