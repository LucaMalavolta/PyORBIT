(running_pyorbit)=

# Running PyORBIT

To launch the analysis, just run this command on your terminal:

```bash
pyorbit_run sampler configuration_file.yaml
```

Where `sampler` is one of the boldface keywords listed in [Samplers](samplers).

The files produced by the analysis will be stored in the folder `configuration_file/sampler/`, automatically created when a given combination *sampler* + *configuration file* is run for the first time.

```{admonition} Analysis preservation
Running again a *sampler* + *configuration file* combo will not overwrite previous results. If you change your configuration file, most likely you will need to either delete the *configuration file* folder or (better) rename the configuration file.
 ```
<!--  Include exceptions for changes in configuration files -->

To save your log on a text file, just redirect the output to a file (see [these examples](https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file)):


```bash
pyorbit_run sampler configuration_file.yaml > configuration_file_sampler_run.log
```

When the analysis is completed, you can retrieve the results using the
`pyorbit_results` command. The syntax is similar to `pyorbit_run`, with the exception that now the command accepts several flags depending on the desired output.

```bash
pyorbit_results sampler_name configuration_file.yaml -flags
```

For a list of possible flags, just run `pyorbit_results` with the `-h` (for *help*) flag. Not using any flag will just produce the standard terminal output. The `-all` flag will produce any conceivable plots and files. Flags are described on this [dedicated](results_flags) page](results_flags).

The output of `pyorbit_results` will be saved into the folder `configuration_file/sampler_plot/` - notice the `_plot` at the end of the name.

```{warning} Analysis preservation
Differently from `pyorbit_run`, all the files in the `configuration_file/sampler_plot/` folder will be overwritten at each execution of the command.
 ```


## Samplers

`PyORBIT` supports both Markov Chain Monte Carlo (MCMC) and Nested Sampling (NS). 
There are some important differences in how the priors are dealt with, sometimes directly affecting the involved parametrizations.
Please refer to [Caveats regarding MCMC and NS](mcmc_nested_caveats) for more details. 
It must be noted that in the majority of my analysis, MCMC and NS delivered perfectly superimposed posteriors.

### Supported samplers 
This is the list of samplers supported by `PyORBIT`. 
The keyword to be used as the first argument of `pyorbit_run`/`pyorbit_results` is reported in boldface.
The same configuration file can be used with different samplers, as the results will be stored in different folders

Fully tested and regularly updated samplers:

- **emcee**: it is actually the combination of two algorithms: the *Global optimization using differential evolution in Python* [PyDE](https://github.com/hpparvi/PyDE) by Hannu Parviainen, and the *The Python ensemble sampling toolkit for affine-invariant MCMC* [emcee](https://emcee.readthedocs.io/en/stable/) by Dan Foreman-Mackey.
- **dynesty**: the *Dynamic Nested Sampling Package* [dynesty](https://dynesty.readthedocs.io/en/stable/) by John Speagle.
- **ultranest**: [ultranest](https://johannesbuchner.github.io/UltraNest/index.html)

The samplers listed above are those that I use more frequently, and I have verified that they provide consistent results among each other. Other algorithms have been implemented as well, although their performances have not been evaluated thoroughly:

- **zeus**: *Lightning Fast MCMC* [zeus](https://zeus-mcmc.readthedocs.io/en/latest/).
- **nestle**: *pure-Python implementation of nested sampling algorithms* [nestle](http://kylebarbary.com/nestle/).
- **polychord**: *next-generation nested sampling* [PolyChordLite](https://github.com/PolyChord/PolyChordLite).
- **multinest**: *Multimodal nested sampling* [MultiNest](https://github.com/farhanferoz/MultiNest) through [PyMultiNest](https://johannesbuchner.github.io/PyMultiNest/).
- **optimize**: Scipy function to be used as an alternative to PyDE.


```{toctree}
:maxdepth: 1
running_pyorbit/mcmc_nested_caveats
running_pyorbit/tinygp_caveats
running_pyorbit/terminal_output
running_pyorbit/results_flags
```
