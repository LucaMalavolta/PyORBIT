(plot_parameters)=

# Plot parameters

Starting with version `11.2.3`, it is possible to provide additional plotting parameters to speed up the plot model calculation when running `pyorbit_results`.

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
* If `True`, files are saved to PDF instead of PNG format.

**plot_config_parameters**
* accepted values: any integer | **9**
* Default font size for the majority of the plots
Plots have been generated with matplotlib using `figure.figsize=(9,9)` and the Arial font family. Some codes may override these settings. 

**rv_lnlike_samplings**
* accepted values: any integer | **`n_sample`**
* aliases: `rv_loglike_samplings', `rv_lnlike_samplings`, `rv_like_samplings`
* Number of samplings to be used to compute the RV log-likelihood. The computation is not parallelised, so computing the RV log-likelihood for all the samples could take too much time. If a number is provided, the RV log-likelihood will be computed on a random selection of *rv_lnlike_samplings* samples. 

**correlation_plot**
* accepted values: **`pygtc`** | `corner` | `getdist`
* Default code to compute the corner plot, among [pygtc](https://pygtc.readthedocs.io/en/latest/), [corner](https://corner.readthedocs.io/en/latest/), [getdist](https://getdist.readthedocs.io/en/latest/).

**force_full_correlation_plot**
* accepted values: `True` | **`False`**
* When the model has more than 30 parameters, the full correlation plot is not generated anymore as the crowding would render it useless. Set this flag to `True` to disable this behaviour.


## Plot step size

It is common practice to generate a model with much denser sampling than the original data, so that smooth variations in the model can be appreciated. 

`PyORBIT` will automatically decide the step size according to this scheme:

- The shortest orbital period among all the planets in the sample is taken as reference as `minimum_planet_period`, even if that planet is not used in the modelling of a specific dataset.
- The minimum non-zero step size of the sorted data is computed, and divided by two to compute the `minimum_step_size`
- For `photometry` data, the step size is automatically increased to 5 minutes if the range of the dataset is larger than three times 'minimum_planet_period', otherwise `minimum_step_size` is used
- For all the other kinds of datasets, the code will use the minimum value among the `minimum_planet_period` divided by 20, the `minimum_step_size`, and the range of the dataset divided by ten times the number of points




**use_shared_axis_for_rv**
* accepted values: **`True`** | `False`
* All radial velocity datasets will share the same axis 

**use_shared_axis_for_activity**
* accepted values: **`True`** | `False`
* All activity index datasets will share the same axis. Photometry is never considered an activity index even if it is modelled with an activity model.




## GP regression 

Computing the model with denser sampling than the original data can be computationally expensive for datasets with a large temporal baseline, as the computational time can increase with the cube of the covariance matrix size. These options should reduce runtime and RAM usage.

```{tip}
If you get a *segmentation fault* error or the code stops abruptly while generating plots, you most likely exceeded the maximum RAM usage.
 ```


**plot_split_threshold** 
* accepted values: any integer | **10000** 
* The calculation of models involving the covariance matrix (e.g., GP regression) will be split into smaller chunks of $N=$`plot_split_threshold` data points. Helpful to speed up the calculation of GP regression for large datasets or multivariate GP analysis.


**low_ram_plot** 
* accepted values: `True` | **`False`** 
* If `True`, the `plot_split_threshold` is further enforced for larger datasets.

**compute_gp_variance**
* accepted values: `True` | **`False`** 
* If `True`, the variance of the GP regression is computed together with the model.


**progress_bar**
* accepted values: **`True`** | `False`
* It will show a progress bar when computing the model using a split threshold. **Note**_ the progress bar is not entirely accurate, so it may not reach 100% even for a successful job.




                #if dataset.kind in activity_datatype or dataset.kind == 'radial_velocity':
                #    input_step_size = plot_config_parameters['model_step_size'].get('activity', input_step_size)                
                #if dataset.kind == 'radial_velocity':
                #    input_step_size = plot_config_parameters['model_step_size'].get('radial_velocity', input_step_size)                
                #if dataset.kind == 'photometry':
                #    input_step_size = plot_config_parameters['model_step_size'].get('photometry', input_step_size)     
                #input_step_size = plot_config_parameters['model_step_size'].get(dataset_name, input_step_size)

