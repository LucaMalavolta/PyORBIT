(results_interpretation)=

# Interpreting the output

When the analysis is completed, you can retrieve the results using the 
`pyorbit_results` command. 

The syntax is similar to `pyorbit_run`, with the exception that now the command
accepts several flags depending on the desired output. 

```bash
pyorbit_results sampler_name configuration_file.yaml -flags
```

For a list of possible flags 

blah bla

```bash

(base) [~]$ pyorbit_results -h

usage: PyORBIT_GetResults.py [-h] [-p [P]] [-w [W]] [-wp [WP]] [-ws [WS]] [-c [C]] [-ln [LN]] [-t [T]] [-fc [FC]] [-cc [CC]] [-dc [DC]] [-v [V]] [-all_corners [ALL_CORNERS]] [-all [ALL]] [-dfm_corner [DFM_CORNER]]
                             [-getdist_corner [GETDIST_CORNER]]
                             sampler config_file

positional arguments:
  sampler               sampler (emcee or polychord)
  config_file           config file

optional arguments:
  -h, --help            show this help message and exit
  -p [P]                Plot model files
  -w [W]                Write model files
  -wp [WP]              Write samples for orbital parameters
  -ws [WS]              Write all samples
  -c [C]                Save chains plots
  -ln [LN]              Save ln_prob chain plot
  -t [T]                Compute and save Gelman-Rubin traces
  -fc [FC]              Save full correlation plot - it may be slow!
  -cc [CC]              Save corner plots of common variables
  -dc [DC]              Save individual corner plots of reach dataset
  -v [V]                Write Veusz files for corner plot
  -all_corners [ALL_CORNERS]
                        Do all the corner plots
  -all [ALL]            Active all flags
  -dfm_corner [DFM_CORNER]
                        Use DFM corner script for corner plots
  -getdist_corner [GETDIST_CORNER]
                        Use getdist script for corner plots

```


```{toctree}
:maxdepth: 1
results_interpretation/terminal_output
```
