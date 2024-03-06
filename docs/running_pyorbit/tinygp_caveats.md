(tinygp_caveats)=




# Caveats on the use of `tinyGP`

The [`tinyGP`](https://tinygp.readthedocs.io/en/stable/) package is a lightweight alternative to the `george` package or the internally implemented Gaussian Process regression.

Models implemented using tinyGP are:

* `tinygp_quasiperiodic`, as an alternative to `gp_quasiperiodic`
* `tinygp_multidimensional_quasiperiodic`, as an alternative to `gp_multidimensional_quasiperiodic`

The `tinyGP` alternatives deliver the same results as their traditional counterparts (as far as I know). The main advantages over the traditional are faster loglikelihood computation and the support of GPU acceleration.

Likely due to the use of [`jax`](https://jax.readthedocs.io/en/latest/) and a slightly different approach compared to the other models, at the moment `PyORBIT` is not able to start the `emcee` analysis after the `pyDE` step, or reload an intermediate `emcee` analysis, unless the code is stopped and launch again. The code will not throw any error, it will simply be stuck without performing any further analysis.

To circumvent this problem, I introduced a flag called `safe_reload` in the `parameter` section of the configuration file.

```yaml
parameters:
  Tref: 59200.00
  safe_reload: True
```

If the `safe_reload` flag is set to `True`, the code will exit after completing the global search with `pyDE`. It is then necessary to launch `PyORBIT` again, but this task can be accomplished by writing a simple `bash` script where `PyORBIT` is just executed twice on the same configuration file.

Following the example in the [`PyORBIT_examples/gaussian_processes`](https://github.com/LucaMalavolta/PyORBIT_examples/tree/main/gaussian_processes) repository, we prepare a simple text file code `exec_all_fits.source` including these lines:

```{code} bash
pyorbit_run emcee RV_tinyGPframework_2p.yaml > RV_tinyGPframework_2p_emcee_run.log
pyorbit_run emcee RV_tinyGPframework_2p.yaml >> RV_tinyGPframework_2p_emcee_run.log
```

```{code} bash
source exec_all_fits.source
```

The first line will execute the `PyDE` step, the second line will launch the `emcee` analysis. Note that in the second line, we have replaced `>` with `>>`, so that the terminal output of the second run will not overwrite the previous terminal output file, but the output will be appended at the end to the existing file.


If you are using the `nsave` option to save the temporary emcee results after a given number of steps, you will need to launch `PyORBIT` a number of time equal to the total number of steps `nsteps` divided by the number of temporary steps `nsave`, plus an additional run for  `PyDE`. Consider the example below:

In this case, you would need to launch `PyORBIT` for `nsteps`/`nsave` = 3 times, plus one for `PyDE`. The resulting script would look like this:

```{code} bash
pyorbit_run emcee RV_tinyGPframework_2p.yaml > RV_tinyGPframework_2p_emcee_run.log
pyorbit_run emcee RV_tinyGPframework_2p.yaml >> RV_tinyGPframework_2p_emcee_run.log
pyorbit_run emcee RV_tinyGPframework_2p.yaml >> RV_tinyGPframework_2p_emcee_run.log
pyorbit_run emcee RV_tinyGPframework_2p.yaml >> RV_tinyGPframework_2p_emcee_run.log
```

You can add as a last line a `pyorbit_results` call to automatically generate the results files.
