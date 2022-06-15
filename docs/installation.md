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

Note that the name is 


To download `PyORBIT` from the repository:

```{code} bash
git clone https://github.com/LucaMalavolta/PyORBIT.git
```

You can then install `PyORBIT` by using `pip` inside the code repository:

```{code} bash cd PyORBIT pip install .
```

(requirements-label)=

## Requirements

This is the list of packages required by PyORBIT to work out-of-the-box:

- `numpy`, `scipy`, `matplotlib`: pretty standard
- `argparse`: required to pass terminal keywords
- `pyyaml`: YAML is the language used for the configuration file
- `corner`: to make corner plots ([corner.py home page])
- `h5py`: HDF5 for Python ([h5py home page])

Basic analysis can be performed using the `scipy.optimize` package, however to fully unwind the power of `PyORBIT` I strongly recommend these two packages (already included in the requirements):
\- `pyDE`: global optimization package ([PyDE home page])
\- `emcee`: ensemble sampling toolkit for affine-invariant MCMC ([emcee home page])

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

If you are using any of those packages listed above, please be sure to cite the proper references, as stated in their web page

## Instructions

### `pyDE`

Installing from pip will results in an error, so you have to install the most up-to-date version from source using the following commands:

```{code} bash
git clone https://github.com/hpparvi/PyDE.git
cd PyDE
python setup.py install
```

From the [pyDE source repository]

[batman home page]: https://www.cfa.harvard.edu/~lkreidberg/batman/
[brew]: https://brew.sh
[celerite home page]: https://github.com/dfm/celerite
[celerite installation page]: http://celerite.readthedocs.io/en/stable/python/install/
[corner.py home page]: https://github.com/dfm/corner.py
[cython]: http://cython.org/
[cython home page]: http://cython.org/
[distutils]: https://docs.python.org/2/extending/building.html
[dynesty home page]: https://github.com/joshspeagle/dynesty
[emcee home page]: https://github.com/dfm/emcee
[fixing missing headers for homebrew in mac os x mojave (from the caffeinated engineer)]: https://silvae86.github.io/sysadmin/mac/osx/mojave/beta/libxml2/2018/07/05/fixing-missing-headers-for-homebrew-in-mac-osx-mojave/
[george home page]: https://github.com/dfm/george
[george installation page]: http://george.readthedocs.io/en/latest/user/quickstart/#installation
[getdist home page]: https://github.com/cmbant/getdist
[h5py home page]: http://docs.h5py.org/en/stable
[hyperthreading]: https://superuser.com/questions/96001/why-does-my-intel-i7-920-display-8-cores-instead-of-4-cores
[multinest home page]: https://github.com/farhanferoz/MultiNest
[openmpi]: https://www.open-mpi.org/
[polychordlite home page]: https://github.com/PolyChord/PolyChordLite
[pyde home page]: https://github.com/hpparvi/PyDE
[pyde source repository]: https://github.com/hpparvi/PyDE
[pymultinest]: https://github.com/JohannesBuchner/PyMultiNest
[pymultinest documentation]: http://johannesbuchner.github.io/PyMultiNest/
[pymultinest home page]: https://github.com/JohannesBuchner/PyMultiNest
[trades home page]: https://github.com/lucaborsato/trades
[ttvfast home page]: https://github.com/kdeck/TTVFast

