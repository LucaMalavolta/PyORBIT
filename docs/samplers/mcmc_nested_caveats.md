(mcmc_nested_caveats)=

# Caveats regarding MCMC and NS

`PyORBIT` is one of the few Exoplanetary codes that supports both Monte Carlo Markov Chain (MCMC) and Nested Sampling (NS), thus allowing a direct comparison of the two approaches over the same datasets and models. In my experience, when the model is appropriate the two algorithms deliver the *same* results, in the sense that the posteriors are perfectly superimposed. Although it does not guarantee you the correctness of your analysis, it is indeed an useful checks to do. 

It is important then to understand the different implementations of the two algorithms, as they affect how priors 

## Priors on parameters

In MCMC, a prior acts as a penalization to the log-likelihood: the further is a parameter from the prior, the more negative (with the other parameters frozen) the log-likelihood will be. This is true for both *sampler* parameters and *derived* parameters. When using **it is possible to set a prior on a derived parameters**. For example, it is possible to set a prior on th eeccentricty when using the Eastman parametrization, even if eccentricity is not a sampler parameter.

In NS, the prior is actually a reparametrization of the unitary cube, in such a way that an uniform sampling of the unitary parameter will result into a distribution of the sampler parameter that follows the given prior. As a consequence, **priors can be imposed only on sampler parameters**. For example, you cannot use the Eastman parametrization and a prior on the eccentricity simultaneously.  

TO BE COMPLETED