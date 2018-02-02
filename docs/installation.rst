.. _installation:

Installation
============

Requirements
++++++++++++

This is the list of packages required by PyORBIT to work out-of-the-box:

- ``numpy``, ``scipy``, ``matplotlib``: pretty standard
- ``argparse``:
- ``pyyaml``: YAML is the language used for the configuration file
- ``pyDE``: global optimization package (`PyDE home page`_)
- ``emcee``: ensemble sampling toolkit for affine-invariant MCMC (`emcee home page`_)

The ``PyDE`` + ``emcee`` combination is the easiest to install and set up, but it is possible to specifiy the starting point of ``emcee`` instead of using the outcome of ``PyDE``.
It is possible to use other sapler as well,

- ``MultiNEST``
- ``PolyChord``


Additional packages may be required to perform certain types of analysis:

- ``corner``: to make corner plots (`corner.py home page`_)
- ``cython`` : C extension for Python (`Cython home page`_)
- ``george`` : Fast and flexible Gaussian Process regression in Python (`george home page`_)
- ``celerite`` : scalable 1D Gaussian Processes (`celerite home page`_)
- ``TRADES`` : dynamical simulation of exoplanetary systems (`TRADES home page`_)
- ``TTVfast`` : transit times for and radial velocities for n-planet systems (`TTVfast home page`_)

.. _Cython home page: http://cython.org/
.. _george home page: https://github.com/dfm/george
.. _celerite home page: https://github.com/dfm/celerite
.. _TRADES home page: https://github.com/lucaborsato/trades
.. _TTVfast home page: https://github.com/kdeck/TTVFast
.. _PyDE home page: https://github.com/hpparvi/PyDE
.. _emcee home page: https://github.com/dfm/emcee
.. _corner.py home page: https://github.com/dfm/corner.py

Instructions
++++++++++++

``pyDE``
--------

Installing from pip results in an error, so you have to install from source using the following commands:

.. code:: bash

  git clone https://github.com/hpparvi/PyDE.git
  cd PyDE
  python setup.py install

From the `pyDE source repository`_

.. _pyDE source repository: https://github.com/hpparvi/PyDE

``emcee``
---------

I’m currently using the latest version of emcee (Version 3.0 at the moment of writing), installed from the source repository:

.. code:: bash

  git clone https://github.com/dfm/emcee.git
  cd emcee
  python setup.py install

(from here: http://emcee.readthedocs.io/en/stable/user/install.html#from-source)



Additional codes
++++++++++++++++

When running PyORBIT you may get one of the following warnings:

.. code:: bash

  WARNING! Imported dummy celerite, models relying on this package will not work
  WARNING! Imported dummy PyPolyChord, models relying on this package will not work
  WARNING! Imported dummy TRADES, models relying on this package will not work
  WARNING! Imported dummy TTVFAST, models relying on this package will not work
  WARNING! Imported dummy george, models relying on this package will not work

*Simple* RV fit and analysis will still work, but if you want to use one of these packages and you are getting one of these error, the code will fail miserabily. You will still have some of these warnings because the relative module is loaded anyway even if you are not actually using it.

The following codes may be required to do some specific kind of analysis.

``george``
----------

Please refer to the `george installation page`_ for detailed instruction on your preferred method of installation.
At the time of writing this guide, using conda installed version 0.3.1 of the package.

.. code:: bash

  conda install -c conda-forge george

Please check that your installed version is equal or higher than ``3.0``:

::

  import george
  george.__version__


.. _george installation page: http://george.readthedocs.io/en/latest/user/quickstart/#installation

``celerite``
------------

Follow the instructions at `celerite installation page`_

.. _celerite installation page: http://celerite.readthedocs.io/en/stable/python/install/


``PolyChord``
-------------

Download the code at `PolyChord home page`_ .
I tested the code with version 1.12, earlier version have a different Python interface and may not work.

.. code:: bash

  tar -xvf PolyChord_v1.12.tar.gz
  cd PolyChord/

Change the Makefile appropriately if you are using weird C/Fortran compilers or Linux architecture. For ``anaconda2`` and ``Ubuntu 16.04 LTS`` I didn't have to change any setting. Then run ``make`` to build the code.

The next step is to configure your ``LD_LIBRARY_PATH`` to point to your PolyChord installation, and your ``LD_PRELOAD`` to point to your mpi installation. PolyChord will tell you the exact line to be added to your ``~\.bashrc`` file by executing:

.. code:: bash

  python run_PyPolyChord.py

Remeber to load the modified ``~\.bashrc`` file by running ``source ~\.bashrc`` in a terminal.


Finally, to use the MPI functionalities, prepend the MPI command before the python one, specyfying the number of processor you want to use after ``-np`` (20 in the example).

.. code:: bash

  mpirun -np 20 python run_PyPolyChord.py

If you already ran the command without the MPI instruction or with a different number of CPU, remember to delete the ``chains`` directory or the execution will fail.

PoliChord+MPI troubleshooting
+++++++++++++++++++++++++++++

Here I report the three errors I encountered so far when I try to install or run PolyChord in MPI mode. For other errors, please refer to the README that comes with the source code.


*Broken  MPI*

If you get the following errors when executing ``run_PyPolyChord.py`` , your MPI/OpenMPI installation is likely broken and you have to re-install it. You need to have a working MPI installation even when you are using PolyChord in single-CPU mode!

.. code:: bash

  [ghoul:30446] [[INVALID],INVALID] ORTE_ERROR_LOG: A system-required executable either could not be found or was not executable by this user in file ess_singleton_module.c at line 231
  [ghoul:30446] [[INVALID],INVALID] ORTE_ERROR_LOG: A system-required executable either could not be found or was not executable by this user in file ess_singleton_module.c at line 140
  [ghoul:30446] [[INVALID],INVALID] ORTE_ERROR_LOG: A system-required executable either could not be found or was not executable by this user in file runtime/orte_init.c at line 128

In my case, I decided to re-build `OpenMPI`_ by following these `instructions <https://www.open-mpi.org/faq/?category=building>`_. Be sure to modify the ``LD_PRELOAD`` in your ``~\.bashrc`` accordingly.
If you are not able to fix the problem, you can still run PolyChord without using the MPI/OpenMPI support (but be ready to wait a lot of time when executing a program...). Open the ``Makefile`` file end switch the MPI flag to zero:

.. code:: bash

  # Whether to use MPI
  MPI=1

then run:

.. code:: bash

  make veryclean
  make

*MPI non starting*

If you get the following error when executing ``mpirun -np 20 python run_PyPolyChord.py`` :

.. code:: bash

  -----------------------------------------------------------------------------
  It seems that there is no lamd running on the host.

  This indicates that the LAM/MPI runtime environment is not operating.
  The LAM/MPI runtime environment is necessary for the "mpirun" command.

  Please run the "lamboot" command the start the LAM/MPI runtime
  environment.  See the LAM/MPI documentation for how to invoke
  "lamboot" across multiple machines.
  -----------------------------------------------------------------------------

Then check if the mpirun executable belongs to the same installation of the library that have been used to compile PolyChord.
For example, in my case I re-installed OpenMPI in the directory ``/home/malavolta/CODE/others/openmpi_dir`` . This is how ```LD_PRELOAD`` is configured in my ``~\.bashrc`` file:

.. code:: bash

  export LD_PRELOAD=/home/malavolta/CODE/others/openmpi_dir/lib/libmpi.so:$LD_PRELOAD
  export LD_LIBRARY_PATH=/home/malavolta/CODE/others/PolyChord/lib:$LD_LIBRARY_PATH

I have to add the path of the binaries of my OpenMPI installation
The correct ``mpirun`` is:

.. code:: bash

  $ which mpirun
  /home/malavolta/CODE/others/openmpi_dir/bin/mpirun

If your ``mpirun`` is not coming from the same installation directory of your MPI libraries, add to the ``PATH`` environment variable the ``bin`` directory of the MPI distribution you are crrently using, at the end of your ``~\.bashrc`` file:

.. code:: bash

  export PATH=/home/malavolta/CODE/others/openmpi_dir/bin:$PATH

*Crash after a few iterations*

If you have an error similar to this one:

.. code:: bash

  -------------------------------------------------------
  Primary job  terminated normally, but 1 process returned
  a non-zero exit code. Per user-direction, the job has been aborted.
  -------------------------------------------------------

  --------------------------------------------------------------------------
  mpirun noticed that process rank 0 with PID 0 on node ghoul exited on signal 11 (Segmentation fault).
  --------------------------------------------------------------------------

You are experiencing a problem already reported in the README file of th ePolyChord source:

Try increasing the stack size:
Linux:    ulimit -s unlimited
OSX:      ulimit -s hard
and resume your job.
The slice sampling & clustering steps use a recursive procedure. The default memory allocated to recursive procedures is embarassingly small (to guard against memory leaks).


.. _OpenMPI: https://www.open-mpi.org/
.. _PolyChord home page: https://ccpforge.cse.rl.ac.uk/gf/project/polychord/
