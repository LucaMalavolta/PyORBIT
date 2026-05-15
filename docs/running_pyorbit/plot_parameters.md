(plot_parameters)=

# Plot parameters

Starting forom version `11.2.3` it is possible to provide additional plotting parameters to speed up the calculation of the plot model when running `pyorbit_results`.

The keywords can be listed in the `parameters` section or in a new `plot_parameters` section:

```{code} yaml
parameters:
  Tref: 2456000.00
  cpu_threads: 128
  safe_reload: False
plot_parameters:
  low_ram_plot: True
  plot_split_threshold: 1000
  ...
```

The default keyword (i.e., the value used by `PyORBIT` if not provided in the `yaml` configuration file) is highlighted in boldface.

**save_pdf**
* accepted values: `True` | **`False`** 
* If `True`, files are saved to PDF instead of PNG format

**plot_config_parameters**
* accepted values: any integer | **9**
* Default font size for the majority of the plots
Plots har egenerated with matplotlib using `figure.figsize=(9,9)` and Arial font family. Some codes may override these settings. 

**rv_lnlike_samplings**
* accepted values: any integer | **`n_sample`**
* aliases: `rv_loglike_samplings', `rv_lnlike_samplings`, `rv_like_samplings`
* Number of samplings to be used to compute the RV log-likehood. The computation is not parallelized, so computing the RV log-likehood for all the samples could take too much time. If a number is provided, the RV log-likelihood will be computed on *rv_lnlike_samplings* random selection of the sample  

**force_full_correlation_plot**
* accepted values: `True` | **`False`**
* When the model has more than 30 parameters, the full correlation plot is not generated anymore as the crowding would render it useless. Set this flag to `True` to disable this behaviour

## Plot step size

It is a common practice to generate a model with a much denser sampling than the original data, so that smooth variations of the model can be appreciated. 

`PyORBIT` will automatically decide the step size according to this scheme:

- The shortest orbital period among all the planets in the sample  is taken as reference as `minimum_planet_period`, even if that planet is not used in the modelling of a specific dataset.
- The minimum non-zero step size of the sorted data is computed, and divided by two to compute the `minimum_step_size`
- For `photometry` data, the step size is automatically increased to 5 minutes if the range of the dataset is larger than three times 'minimum_planet_period', otherwise `minimum_step_size` is used
- For all the other kinds of dataset, the code will use the minimum value among the `minimum_planet_period` divided by 20, the `minimum_step_size`, and the range of the dataset divided by ten times the number of points


