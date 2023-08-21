(samplers)=

# Samplers

`PyORBIT` supports both Markov Chain Monte Carlo (MCMC) and Nesgted Sampling (NS). 
There are some important differences in how the priors are dealt with, sometimes directly affecting the involved parametrizations.
Please refer to [Caveats regarding MCMC and NS](mcmc_nested_caveats) for more details. 
It must be noted that in the majority of my analysis MCMC and NS delivered perfeclty supoerimposed posteriors.

## Supported samplers 
This is the list of samplers supported by `PyORBIT`. 
The keyword to be used as first argument of `pyorbit_run`/`pyorbit_results` is reported in boldface.
The same configuration file can be used with different sampler, as the results will be stored in different folders

Fully testes and regularly updated samplers:

- **emcee**: it is actually the combination of two algorithms: the *Global optimization using differential evolution in Python* [PyDE](https://github.com/hpparvi/PyDE) by Hannu Parviainen, and the *The Python ensemble sampling toolkit for affine-invariant MCMC* [emcee](https://emcee.readthedocs.io/en/stable/) by Dan Foreman-Mackey.
- **dynesty**: the *Dynamic Nested Sampling Package* [dynesty](https://dynesty.readthedocs.io/en/stable/) by John Speagle.
- **ultranest**: [ultranest](https://johannesbuchner.github.io/UltraNest/index.html)

The samplers listed above are those that I use more frequently, and I have verified that they provide consistent results among each other. Other algorithms have been implemented as well, although their performances have not been evaluated thoroughly:

- **zeus**: *Lightning Fast MCMC* [zeus](https://zeus-mcmc.readthedocs.io/en/latest/).
- **nestle**: *pure-Python implementation of nested sampling algorithms* [nestle](http://kylebarbary.com/nestle/).
- **polychord**: *next-generation nested sampling* [PolyChordLite](https://github.com/PolyChord/PolyChordLite)
- **multinest**: *Multimodal nested sampling* [MultiNest](https://github.com/farhanferoz/MultiNest) through [PyMultiNest](https://johannesbuchner.github.io/PyMultiNest/)
- **optimize**: Scipy function to be used as alternative to PyDE


```{toctree}
:maxdepth: 1
samplers/mcmc_nested_caveats.md
```