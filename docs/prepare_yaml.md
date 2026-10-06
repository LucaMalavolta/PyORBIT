(prepare_yaml)=

# Prepare a configuration file

The configuration file is the most important part of your analysis because it will contain the specifics of your model,
the priors you want to assign to the parameters, and the way you want to analyze your data.
Thanks to the readability of the YAML language, it will also serve as a useful memento at a glance for the set-up of your model.

## File structure

The file is divided in several sections, each one with a clear purpose:

- ``inputs``: here you list the files containing the datasets.
- ``common``: the physical parameters of your planetary system that are going to be shared between the models (and are independent from them): orbital parameters of the planets, mass and radius of the star, rotational period and decay timescale of the active regions...
- ``models``: this is where you specify how the physical elements of the system in the *common* section are going to influence your data: which planets should be included in the radial velocity computation, what kind of model/kernel you want to use to model the activity...
- ``parameters``: values that are not system parameters per se, but are required to properly perform the analysis.
- ``plot_parameters`` (optional): settings used by `pyorbit_results` to control plot output and model calculations, such as the output format and memory usage. If this section is absent, plotting settings can be placed under ``parameters`` instead.
- ``solver``: the parameters required by the samplers are all listed here.

<!---
To have a glance at how a configuration file looks like, check the .. _documentation_example.yaml: http://cnn.com/:
--->

```{toctree}
:maxdepth: 1
configuration_file/prepare_input
configuration_file/prepare_common
configuration_file/plot_parameters
```
