(installation)=

# Installing PyORBIT

## Using PyORBIT on Windows

Short version: don't do it.
Long version: I never tested `PyORBIT` on Windows, so I cannot guarantee that the outcome of the analysis is what it is supposed to be, or if it will work at all.  

Rather than trying to fix `PyORBIT` or Python on Windows, I usually suggest [setting up an Ubuntu machine with Windows Subsystem for Linux 2](https://canonical-ubuntu-wsl.readthedocs-hosted.com/en/latest/). With respect to installing a virtual machine, WSL2 will use fewer computer resources while giving direct access from Windows to Linux files and vice versa. 


## Setting up an environment

Before proceeding with the installation, I suggest creating an environment dedicated to `PyORBIT` using **Python 3.10**. This is the version I'm currently use to test the code.

With conda/anaconda:

```{code} bash
conda create --name pyorbit python=3.10
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

You can then install `PyORBIT` by using `pip`:

```{code} bash
 pip install pyorbit-package
```

Note that the name is `pyorbit-package` and not `pyorbit`, as the former was already taken by another package in PyPI (although not installable). The name for package importation will still be `pyorbit`:

```{code} bash
 python -c "import pyorbit"
```

`PyORBIT` only comes with the basic packages; to exploit its full potential, follow the instructions in the section [Requirements](Requirements). For the impatient, just do:
```{code} bash
wget https://raw.githubusercontent.com/LucaMalavolta/PyORBIT/main/extra_requirements.txt
pip install -r extra_requirements.txt
```

```{tip}
The requirement file has weaker constraints on package versioning (see below), so remember to install `PyORBIT` first  to avoid version incompatibilities
```

## `starry` support **updated**

The [`starry`](https://starry.readthedocs.io/en/latest/) code package is a suite of tools for mapping stars and exoplanets based on time series data, and it has been implemented in several models within `PyORBIT`. Due to the use of discontinued libraries as [`Theano`](https://github.com/Theano), its installation requires some different steps. 

### `g++` installation
In order to work, starry requires `g++` (available through `gcc`), while `blas` libraries are optional.
It may be possible that you already have `g++` installed, in this cause you should get a  `fatal error` when trying to run iot without an input file 
```{code} bash
g++
   g++: fatal error: no input files
   compilation terminated.
```

On my Ubuntu computer, I managed to install `g++` and `blas` libraries only after installing `aptitude`. 
```{code} bash
sudo apt install build-essential manpages-dev software-properties-common
sudo apt install aptitude
sudo aptitude install g++
sudo aptitude install libopenblas-dev
```
In the case of `g++`, I had to take a step further by checking the proposed solutions to solve conflicts.  

Installation on Fedora was much more straightforward. 

```{code} bash
sudo yum install gcc-c++
sudo yum install blas blas-devel
```

I'm just reporting these issues for your convenience, if you are running into any trouble, please check with your IT crowd / best friend.  

### installing `starry` 

```{note} 
Thefollowing instructions will work only with PyORBIT version 10.3.0 or above
```


First of all, create a dedicated environment:
```{code} bash
conda create --name starry python=3.9
```
This step is **strongly** suggested as the installation of `starry` will downgrade several packages, causing dependency issues.

The best way to make sure that our packages are compatible with `starry` is to install it as the first package in the newly created environment, together with some extra packages.

```{code} bash
conda activate starry

```

We all the packages required by `starry`, after I painfully checked all the required versions to avoid dependency errors:
```{code} bash
wget https://raw.githubusercontent.com/LucaMalavolta/PyORBIT/main/starry_requirements.txt
pip install -r starry_requirements.txt
```
This command will install `starry` as well.

We finally install `PyORBIT` using `pip`, but without checking for dependencies
```{code} bash
pip install --no-dependencies pyorbit-package
```

```{warning} 
`tinygp` and all the packages relying on `jax` will not work with this installation.
```

You should see something like this:
```{code} bash
(starry) [~]$ python
Python 3.9.17 (main, Jul  5 2023, 20:41:20)
[GCC 11.2.0] :: Anaconda, Inc. on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import starry
>>> starry.__version__
'1.2.0'
>>>
```

If you get an error relative to `g++`, it means that something went wrong during its installation and you will not be able to use `starry` until you fix it. If you get an error about the `blas` library, you may still use `starry` but with possibly slower performances.



## Install from the repository

Download the latest version from the GitHub repository:

```{code} bash
 git clone https://github.com/LucaMalavolta/PyORBIT.git
```

If you are downloading the code from GitHub, most likely it's because you want to try the *development* version

```{code} bash
 cd PyORBIT
 git checkout development
```

```{note} 
The development version may not be available at all times.
```

You can then install `PyORBIT` by using `pip` inside the code repository:

```{code} bash
 pip install .
```

Alternatively, you can install `PyORBIT` using the `setup.py` file:

```{code} bash
 python setup.py install
```

Keep in mind that you can still run PyORBIT by specifying the full path of the code. For example, if you cloned the repository in the folder ``~/CODE/PyORBIT``, you can run the analysis and explore the results with these scripts:

```{code} bash
python ~/CODE/PyORBIT/PyORBIT_Run.py emcee configuration_file.yaml
python ~/CODE/PyORBIT/PyORBIT_Results.py emcee configuration_file.yaml -all
```

Again, I suggest to install the extra requirements following the instructions given above

If you repent and you want to go back to a previous version of  `PyORBIT`, just install the desired version with pip:

```{code} bash
 pip install pyorbit-package==10.0.0
 ```

## Requirements

```{admonition} Give people credit for their work

If you are using any of those packages listed above, *please be sure to cite the proper references*, as stated in the relative web page.
```

These packages are installed automatically when using pip.

- `numpy`, `scipy`, `matplotlib`: pretty standard
- `numba`: open source JIT compiler, actually required as an undeclared dependency by some packages ([numba home page])
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
 `pyDE` and all the additional requirements can be installed by locating the `extra_requirements.txt` file in the [PyORBIT repository] or by downloading it directly [from here](https://raw.githubusercontent.com/LucaMalavolta/PyORBIT/main/extra_requirements.txt) and then running from a terminal:

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

- `dynesty`: a pure Python Dynamic Nested Sampling package for estimating Bayesian posteriors and evidence ([dynesty home page])
- `UltraNest`: fit and compare complex models reliably and rapidly with advanced sampling techniques ([UltraNest home page])

Other models:

- `PyAstronomy`: a collection of astronomy-related packages ([PyAstronomy home page])
- `starry`: a suite of tools for mapping stars and exoplanets based on timeseries data ([starry home page])
- `spiderman`: A fast code to simulate secondary transits and phase curves ([spiderman home page])

````{warning}
 `spiderman` installation with recent version of matplotlib will fail, a working version can be found [at this repository] and can be installed directly from terminal with this comand:

```{code} bash
pip install git+https://github.com/LucaMalavolta/SPIDERMAN.git
```

as the installation of this package may result problematic for its many dependencies, it has been excluded from the `extra_requirements.txt` file.
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
