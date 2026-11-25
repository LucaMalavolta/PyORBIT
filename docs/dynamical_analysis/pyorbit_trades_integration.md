(pyorbit_trades_integration)=

# Integration with TRADES

Install TRADES from the [LucaMalavolta/trades repository](https://github.com/LucaMalavolta/trades) in the same Python 3.14 environment as PyORBIT. Follow the {ref}`PyORBIT installation guide <installation>` to set up the environment first.

## Installing TRADES

On Ubuntu or WSL, install the tools needed to build the Fortran extension.

```{code} bash
sudo apt install gcc gfortran make git
```

Activate the PyORBIT environment, clone the TRADES repository and install it from the repository root. Pip resolves the Python dependencies and compiles the Fortran extension for the active interpreter.

```{code} bash
conda activate pyorbit
git clone https://github.com/LucaMalavolta/trades.git
cd trades
python -m pip install .
```

Run the import check from outside the TRADES checkout to load the installed extension.

```{code} bash
cd ..
python -c "from pytrades import pytrades_lib; print(pytrades_lib.__file__)"
```

## Updating TRADES or PyORBIT

Run the following commands from the TRADES checkout in the PyORBIT environment to pull changes and rebuild the extension.

```{code} bash
git pull --ff-only
python -m pip install .
```

To update PyORBIT, pull the changes in its development checkout and repeat the {ref}`checkout installation <full-checkout-installation>` from its root.
