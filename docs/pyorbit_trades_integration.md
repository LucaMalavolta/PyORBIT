(pyorbit_trades_integration)=

# Integration with TRADES

`TRADES` (TRAnsits and Dynamics of Exoplanetary Systems) is a Fortran 90 code developed by [Borsato et al. 2014](https://ui.adsabs.harvard.edu/abs/2014A%2526A...571A..38B/abstract) to perform dynamical modelling of exoplanetary signals. The code has been expanded with a Python interface ([Borsato et al. 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.3233B/abstract)) and to perform photodynamical modelling ([Borsato et al. 2024](https://ui.adsabs.harvard.edu/abs/2024A%2526A...689A..52B/abstract)).

## Installing PyORBIT on a dedicated environment

`PyORBIT` and `TRADES` rely on different versions of `numpy` due to some requirements set by the Fortran compiler.

I strongly suggest preparing a new Python environment specifically dedicated to running PyORBIT with the TRADES integration.

First, install PyORBIT using the instructions in <project:./installation.md> minding to use a different name for the environment, such as **pyorbit_trades**

```{code} bash
conda create --name pyorbit_trades python=3.10 
conda activate pyorbit_trades
pip install pyorbit-package
```

Remember to install the extra dependencies from the `PyORBIT` repository:

```{code} bash
wget https://raw.githubusercontent.com/LucaMalavolta/PyORBIT/main/extra_requirements.txt
pip install -r extra_requirements.txt
```

## Installing TRADES

`TRADES` is freely available on [this GitHub repository](https://github.com/lucaborsato/trades) by Luca Borsato; however, it is not pip-installable due to the Fortran compilation.

The following installation instructions come from the official `TRADES` repository, but they have been adapted and shortened to take into account the packages already installed by `PyORBIT`.

First of all, be sure you have a Fortran 90 compiler installed on your Linux box, for example, by executing the command `gfortran --version`. If it isn't already preinstalled on your Linux distribution, it can usually be installed with a single command via your distribution's repository manager.



```{code} bash
gfortran --version
> GNU Fortran (Ubuntu 11.4.0-1ubuntu1\~22.04.2) 11.4.0
> Copyright (C) 2021 Free Software Foundation, Inc.
```

In your `pyorbit_trades` environment, install the following packages with the corresponding version specifiers:

```{code} bash
pip install  pybind11  numpy==1.23.5 setuptools==65.6.3 numba==0.60.0 jax==0.4.30 jaxlib==0.4.30
``` 

You will get the following warning:

```{code} bash 
ERROR: pip's dependency resolver does not currently take into account all the packages that are installed. This behaviour is the source of the following dependency conflicts.
pyorbit-package 11.2.0 requires jax==0.4.31, but you have jax 0.4.30 which is incompatible.
pyorbit-package 11.2.0 requires jaxlib==0.4.31, but you have jaxlib 0.4.30 which is incompatible.
pyorbit-package 11.2.0 requires numpy==1.24.3, but you have numpy 1.23.5 which is incompatible.
```

However, `PyORBIT` will work fine even with these older versions of the codes.

```{tip}
Remember to run the command above every time you update `PyORBIT` in the `pyorbit_trades` environment.
```

Download the `TRADES` repository in a folder of your choice:

```{code} bash
git clone https://github.com/lucaborsato/trades.git
```

Enter the `src` directory inside the repository and execute the installation command to

```{code} bash
cd trades/src
make cleanall
make full_parallel_release
```

Ignore the following warning if you get it:

```{code} bash
WARN: Could not locate executable armflang
```

Finally, you need to install the Python wrapper of the code:

```{code} bash
cd ..
pip install .
```

That's all!

