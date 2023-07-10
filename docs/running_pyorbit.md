(running_pyorbit)=

# Running PyORBIT

To launch the analysis, just run this command on your terminal:

```bash
pyorbit_run sampler configuration_file.yaml
```

Where `sampler` is one of the boldface keywords listed in [Samplers](samplers).

The files produced by the analysis will be stored in the folder `configuration_file/sampler/`, automatically created when specific combination *sampler* + *configuration file* is run for the first time.

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

For a list of possible flags, just run `pyorbit_results` with the `-h` (for *help*) flag. Not using any flag will just produce the standard terminal output. The `-all` flag will produce any conceivable plots and files. Flags are described in this [dedicated page](results_flags).

The output of `pyorbit_results` will be saved into the folder `configuration_file/sampler_plot/` - notice the `_plot` at the end of the name.

```{warning} Analysis preservation
Differently from `pyorbit_run`, all the files in the `configuration_file/sampler_plot/` folder will be overwritten at each execution of the command.
 ```



```{toctree}
:maxdepth: 1
running_pyorbit/terminal_output
running_pyorbit/results_flags
```
