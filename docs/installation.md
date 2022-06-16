# Installing PyORBIT

## Setting up an environment

Before proceeding with the installation, I suggest to create an environment dedicated to `PyORBIT` using python\<=3.9 .  At the moment of writing I hve received a few complaints (unrelated to PyORBIT) about Python 3.10, so you may use it at your own risk. 
With conda/anaconda:

```{code} bash
conda create --name pyorbit python=3.9
```

To list the available environments do:

```{code} bash
conda env list
```

The active environment will be marked with a \*

To activate the `pyorbit` environment:

```{code} bash
WINDOWS: activate pyorbit
LINUX, macOS: conda activate pyorbit
```

## Install using pip

You can then install `PyORBIT` by using `pip` inside the code repository:

```{code} bash
 pip install pyorbit-package
```

Note that the name is  `pyorbit-package` and not  `pyorbit`, as the former was already taken by another package in PyPI (altough not installable). The name for package importation will still be `pyorbit`:

```{code} bash
 python -c "import pyorbit"
```

## Install from the repository

Download the latest version from the GitHub repository:

```{code} bash
 git clone https://github.com/LucaMalavolta/PyORBIT.git
```

You can then install `PyORBIT` by using `pip` inside the code repository:

```{code} bash 
 cd PyORBIT
 pip install .
```

Alternatively, you can install `PyORBIT` using the `setup.py` file:

```{code} bash
 cd PyORBIT
 python setup.py install
```

## Requirements

```{admonition} Give people credit for their work

If you are using any of those packages listed above, *please be sure to cite the proper references*, as stated in the relative web page. 
```

These packages are installed automatically when using pip.

- `numpy`, `scipy`, `matplotlib`: pretty standard
- `numba`: open source JIT compiler, actually required as undeclared dependency by some packages ([numba home page])
- `argparse`: Parser for command-line options, arguments and sub-commands, required to pass terminal keywords ([argpares home page])
- `pyyaml`: a full-featured YAML framework for the Python programming language.  YAML is the language used for the configuration file ([pyyaml home page], [yaml home page])
- `h5py`: HDF5 for Python ([h5py home page])
- `pygtc`: Make a publication-ready giant-triangle-confusogram (GTC) ([pygtc home page])
- `tqdm`: A Fast, Extensible Progress Bar for Python and CLI ([tqdm home page])

Basic analysis can be performed using the `scipy.optimize` package, however to fully unwind the power of `PyORBIT` these two packages should be installed:

- `pyDE`: global optimization package ([PyDE home page])
- `emcee`: ensemble sampling toolkit for affine-invariant MCMC ([emcee home page]). 

`emcee` is already included in the requirements, `pyDE` needs to be installed separately as the GitHub version supports multiprocessing: 

```{code} bash
 pip install git+https://github.com/hpparvi/PyDE.git
```

THe alternative ensemble slice sampler `zeus` is supported as well ([zeus home page]).

````{tip}
 `pyDE` and all the additional requirements can be installed by locating the `extra_requirements.txt` file in the [PyORBIT repository] or by downloading it directly [from here](https://github.com/LucaMalavolta/PyORBIT/blob/main/extra_requirements.txt) and then running from a terminal:

 ```{code} bash
 pip install -r extra_requirements.txt
 ```

````

[extra_requirements.txt]: https://github.com/LucaMalavolta/PyORBIT/blob/main/extra_requirements.txt

## Additional requirements

<!---
Simply speaking, `PyDE` searches for the best global solution and passes it to `emcee`, ensuring that the MCMC will not be stuck around a local minimum of the chi-square. The `PyDE` + `emcee` combination is the easiest to install and set up, but it is possible to specify the starting point of `emcee` instead of using the outcome of `PyDE`.
It is possible to use other samplers as well, such as:

- `MultiNEST` ([MultiNest home page] and [PyMultiNest home page])
- `PolyChordLite`, previously known as just `PolyChord` ([PolyChordLite home page])
- `dynesty` ([dynesty home page])

Additional packages may be required to perform certain types of analysis:

- `batman`: Bad-Ass Transit Model cAlculatioN ([BATMAN home page])
- `george` : Fast and flexible Gaussian Process regression in Python ([george home page])
- `celerite` : scalable 1D Gaussian Processes ([celerite home page])
- `TRADES` : dynamical simulation of exoplanetary systems ([TRADES home page])
- `TTVfast` : transit times for and radial velocities for n-planet systems ([TTVfast home page])
- `cython` : C extension for Python ([Cython home page])
- `getdist`: For the analysis of MCMC chains ([getdist home page])
-->
Depending on the kind of analysis you want to perform, you may need the following packages:

Transit modelling:

- `batman`: Bad-Ass Transit Model cAlculatioN ([BATMAN home page])
- `PyTransit`: a package for exoplanet transit light curve modelling ([Pttransit home page])

Gaussian Process Regression:

- `george`: a fast and flexible Python library for Gaussian process regression ([george home page])
- `celerite` or `celerite2`: an algorithm for fast and scalable Gaussian process regression ([celerite home page], [celerite2 home page])

Nersted Sampling:

- `dynesty`: a pure Python Dynamic Nested Sampling package for estimating Bayesian posteriors and evidences ([dynesty home page])
- `UltraNest`: fit and compare complex models reliably and rapidly with advanced sampling techniques ([UltraNest home page])

Other models:

- `PyAstronomy`: a collection of astronomy related packages ([PyAstronomy home page])
- `starry`: a suite of tools for mapping stars and exoplanets based on timeseries data ([starry home page])
- `spiderman`: A fast code to simulate secondary transits and phase curves ([spiderman home page])

````{warning}
 `spiderman` installation with recent version of matplotlib will fail, a working version can be found [at this repository] and can be installed directly from terminal with this comand:

```{code} bash
pip install git+https://github.com/LucaMalavolta/SPIDERMAN.git
```
````

[PyORBIT repository]: https://github.com/LucaMalavolta/PyORBIT

[numba home page]: https://numba.pydata.org/
[tqdm home page]: https://tqdm.github.io/
[pygtc home page]: https://pygtc.readthedocs.io/
[argpares home page]: https://docs.python.org/3/library/argparse.html
[pyyaml home page]: https://pyyaml.org/
[yaml home page]: https://yaml.org/
[emcee home page]: https://emcee.readthedocs.io/
[h5py home page]: http://docs.h5py.org/
[zeus home page]: https://zeus-mcmc.readthedocs.io/

[BATMAN home page]: https://github.com/lkreidberg/batman
[celerite home page]: https://celerite.readthedocs.io/
[celerite2 home page]: https://celerite2.readthedocs.io/
[PyTransit home page]: https://pytransit.readthedocs.io/
[george home page]: https://george.readthedocs.io/
[dynesty home page]: https://dynesty.readthedocs.io/
[UltraNest home page]: https://johannesbuchner.github.io/UltraNest/index.html
[PyAstronomy home page]: https://pyastronomy.readthedocs.io/
[starry home page]: https://starry.readthedocs.io/
[spiderman home page]: https://spiderman.readthedocs.io/en/latest/
[at this repository]: https://github.com/LucaMalavolta/SPIDERMAN
