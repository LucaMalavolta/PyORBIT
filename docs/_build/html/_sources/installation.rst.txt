.. _installation:

Installation
============

.. _requirements-label:

Requirements
++++++++++++

This is the list of packages required by PyORBIT to work out-of-the-box:

- ``numpy``, ``scipy``, ``matplotlib``: pretty standard
- ``argparse``: required to pass terminal keywords
- ``pyyaml``: YAML is the language used for the configuration file
- ``corner``: to make corner plots (`corner.py home page`_)
- ``h5py``: HDF5 for Python (`h5py home page`_)

Basic analysis can be performed using the ``scipy.optimize`` package, however to fully unwind the power of ``PyORBIT`` I strongly recommend these two packages:
- ``pyDE``: global optimization package (`PyDE home page`_)
- ``emcee``: ensemble sampling toolkit for affine-invariant MCMC (`emcee home page`_)

Simply speaking, ``PyDE`` searches for the best global solution and passes it to ``emcee``, ensuring that the MCMC will not be stuck around a local minimum of the chi-square. The ``PyDE`` + ``emcee`` combination is the easiest to install and set up, but it is possible to specify the starting point of ``emcee`` instead of using the outcome of ``PyDE``.
It is possible to use other samplers as well, such as:

- ``MultiNEST`` (`MultiNest home page`_)
- ``PolyChordLite``, previously known as just ``PolyChord`` (`PolyChordLite home page`_)
- ``dynesty`` (`dynesty home page`_)

Additional packages may be required to perform certain types of analysis:

- ``batman``: Bad-Ass Transit Model cAlculatioN (`BATMAN home page`_)
- ``george`` : Fast and flexible Gaussian Process regression in Python (`george home page`_)
- ``celerite`` : scalable 1D Gaussian Processes (`celerite home page`_)
- ``TRADES`` : dynamical simulation of exoplanetary systems (`TRADES home page`_)
- ``TTVfast`` : transit times for and radial velocities for n-planet systems (`TTVfast home page`_)
- ``cython`` : C extension for Python (`Cython home page`_)
- ``getdist``: For the analysis of MCMC chains (`getdist home page`_)


.. _BATMAN home page: https://www.cfa.harvard.edu/~lkreidberg/batman/
.. _Cython home page: http://cython.org/
.. _george home page: https://github.com/dfm/george
.. _celerite home page: https://github.com/dfm/celerite
.. _TRADES home page: https://github.com/lucaborsato/trades
.. _TTVfast home page: https://github.com/kdeck/TTVFast
.. _PyDE home page: https://github.com/hpparvi/PyDE
.. _emcee home page: https://github.com/dfm/emcee
.. _corner.py home page: https://github.com/dfm/corner.py
.. _h5py home page: http://docs.h5py.org/en/stable
.. _getdist home page: https://github.com/cmbant/getdist
.. _MultiNest home page: https://github.com/farhanferoz/MultiNest
.. _PolyChordLite home page: https://github.com/PolyChord/PolyChordLite
.. _dynesty home page: https://github.com/joshspeagle/dynesty

If you are using any of those packages listed above, please be sure to cite the proper references, as stated in their web page

Instructions
++++++++++++

``pyDE``
--------

Installing from pip will results in an error, so you have to install the most up-to-date version from source using the following commands:

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

In principle PyORBIT_GetResults should be able to recognize if the output files have been produced by version ``2.x`` or ``3.x``. To save you some trouble, however, I suggest you to check that you have actually installed version ``3.x``:

::

  import emcee
  emcee.__version__


Additional codes
++++++++++++++++

When running PyORBIT you may get one of the following warnings:

.. code:: bash

  WARNING! Imported dummy celerite, models relying on this package will not work
  WARNING: Could not preload libmpi.so.If you are running with MPI, this may cause segfaults
  WARNING! Imported dummy PyPolyChord, models relying on this package will not work
  WARNING! Imported dummy TRADES, models relying on this package will not work
  WARNING! Imported dummy TTVFAST, models relying on this package will not work
  WARNING! Imported dummy george, models relying on this package will not work
  WARNING! Imported pyorbit.classes.dummy batman, models relying on this package will not work

*Simple* RV fit and analysis will still work, but if you want to use one of these packages and you
are getting one of these error, the code will fail miserably. You will still have some of these
warnings because ``PyORBIT`` will try to load the relative module even if you are not actually using it.
So be worried only if you want to do some of the things listed here :ref:`requirements-label` and the appropriate package is not installed (the code will crash anyway).

The following codes may be required to do some specific kind of analysis.

george
------

Please refer to the `george installation page`_ for detailed instruction on your preferred method of installation.
At the time of writing this guide, using conda installed version 0.3.1 of the package.

.. code:: bash

  conda install -c conda-forge george

Please check that your installed version is equal or higher than ``0.3``:

::

  import george
  george.__version__


.. _george installation page: http://george.readthedocs.io/en/latest/user/quickstart/#installation

celerite
--------

On some systems, importing ``george`` and ``celerite`` during the same Python session may cause a segmentation fault. The only workaround I found is to install ``celerite`` using conda-forge instead of pip.
If you are not planning to use celerite, you may proceed with the standard installation through conda-forge:

.. code:: bash

  conda install -c conda-forge celerite


If you plan to use celerite, you may be interested in compiling from source in order to enable improved automatic differentiation. Please refer to the `celerite installation page`_ to check the best option for your installation.

.. _celerite installation page: http://celerite.readthedocs.io/en/stable/python/install/



PolyChordLite
-------------

Download the code at `PolyChordLite home page`_ .
``pypolychord``, the Python interface of ``PolyChord``, has been revamped starting from version ``1.12`` and then renamed after its transformation to ``PolyChordLite``. Earlier versions will likely not work with ``PyORBIT``.

.. code:: bash

  git clone https://github.com/PolyChord/PolyChordLite.git
  cd PolyChordLite/

Change the Makefile appropriately if you are using weird C/Fortran compilers or Linux distributions. With ``anaconda`` on ``Ubuntu 16.04 LTS`` and ``Ubuntu 18.04 LTS`` I didn't have to change any setting.
In the past, ``MPI`` was disabled by default when installing on ``macOS``, I didn't manage to make MPI and PolyChord work together on my laptop so decided to leave it that way. Right now I'm not sure what is the situation.

When you have finished modifying the Makefile, to build the code run

.. code:: bash

  make pypolychord
  python setup.py install --user

The next step is to configure your ``LD_LIBRARY_PATH`` to point to your PolyChord installation, and your ``LD_PRELOAD`` to point to your mpi installation. PolyChord will tell you the exact line to be added to your ``~\.bashrc`` file by executing:

.. code:: bash

  python run_pypolychord.py

Remeber to load the modified ``~\.bashrc`` file by running ``source ~\.bashrc`` in a terminal.


Finally, to use the MPI functionalities, prepend the MPI command before the python one, specifying the number of processor you want to use after ``-np`` (20 in the example).

.. code:: bash

  mpirun -np 20 python run_pypolychord.py

If you already ran the command without the MPI instruction or with a different number of CPU, remember to delete the ``chains`` directory or the execution will fail.

** Compilation tricks for Mac users **


Cythonizing your code
+++++++++++++++++++++

You can improve the performance of the code by compiling it with ``Cython`` and ``distutils``. To compile the code, just execute

.. code:: bash

  ./compile.sh

in the main directory of the source code of ``PyORBIT``. Note that you have to run the command every time you change a file in the code,
 otherwise the compiled version will stay behind.

.. code:: bash

  ./compile.sh

To clean the repository fro the compiled version, .i.e. if frequent changes are made to the code and you want to avoid recompiling each time, simply run:

.. code:: bash

  ./clean_compile.sh

Note that in order to allow cythonization, the ``.py`` files in the ``pyorbit/classes`` and ``pyorbit/models``
directory are actually symbolic links to the ``.pyx`` files in the same directory.

More information on `Cython`_ and `distutils`_ can be found at their respective web pages.

.. _Cython: http://cython.org/
.. _distutils: https://docs.python.org/2/extending/building.html


Mac Troubleshooting
+++++++++++++++++++

I run my code a Linux Box, but if I need to do a quick test or debug and I’m not in the office I do it on my Mac. Unfortunately some things are not as straightforward as they should be.
Until now the most problematic external code has been PolyChord, particularly if you try to install it with the MPI support.
Below you can find a collection of errors I found along the way and how I fix them.

In the following, I assume you have installed the Command Line Tools with the command ``xcode-select --install`` and the package manager for macOS `brew`_. If you are using macOS 10.14 (Mojave), follow this additional instructions here: `Fixing missing headers for homebrew in Mac OS X Mojave (from The caffeinated engineer)`_. Note that I had to download again the Command Line Tools from the Apple Developer website in order to have the  macOS SDK headers appearing in the correct folder.

**I'm not a IT expert, use these advices at your own risk!**

.. _Fixing missing headers for homebrew in Mac OS X Mojave (from The caffeinated engineer): https://silvae86.github.io/sysadmin/mac/osx/mojave/beta/libxml2/2018/07/05/fixing-missing-headers-for-homebrew-in-mac-osx-mojave/
.. _brew: https://brew.sh


PolyChord - check gcc/gfortran/g++ versions
-------------------------------------------

**Note:** This error may show up if PolyChord is being compiled without MPI support, which is disabled by default.

I have ``gfortran`` installed through ``brew`` on my ``macOS 10.14``, but when I run ``make pypolychord`` it keeps asking for ``gfortran-8`` when installing ``PolyChord 1.16``. The offending lines are from 11 to 13 of the ``Makefile_gnu`` file, in the main directory:

.. code:: bash
  FC = gfortran-8
  CC = gcc-8
  CXX = g++-8

To fix this, first check the version of your fortran compiler with ``gfortran -v``:

.. code:: bash
  $ gfortran -v
    Using built-in specs.
    COLLECT_GCC=gfortran
    COLLECT_LTO_WRAPPER=/usr/local/Cellar/gcc/9.1.0/libexec/gcc/x86_64-apple-darwin18/9.1.0/lto-wrapper
    Target: x86_64-apple-darwin18
    Configured with: .... [cut]
    Thread model: posix
    gcc version 9.1.0 (Homebrew GCC 9.1.0)

From the last line I can see that my ``gfortran`` is part of version 9 of the ``gcc`` compiler provided by ``brew``. However, a version check of ``gcc`` gives a different answer:

.. code:: bash
  $ gcc -v
    Configured with: ...[cut]
    Apple clang version 11.0.0 (clang-1100.0.20.17)
    Target: x86_64-apple-darwin18.6.0
    Thread model: posix
    InstalledDir: /Library/Developer/CommandLineTools/usr/bin

In other words, the command ``gcc`` will call the version provided by Apple, while ``gfortran`` comes with the ``brew`` version of ``gcc`` (and apparently it's not provided by Apple at all). To avoid conflicts with libraries, be sure to use to identify the correct commands to call ``gcc``, ``gfortran`` and ``g++`` from the same installation. Most of the time, you just have to append the version number at the end, i.e. ``gcc-9``, ``gfortran-9``, and ``g++-9``.

Finally, modify the ``Makefile_gnu`` accordingly:

.. code:: bash
  FC = gfortran-9
  CC = gcc-9
  CXX = g++-9

Run ``make pypolychord``, ignore the warnings, and then execute the command suggested at the end (if compilation was successful), in my case ``CC=gcc-9 CXX=g++-9 python setup.py install --user``

PolyChord + MPI - check gcc/gfortran/g++ versions
-------------------------------------------------

``PolyChord`` installation with ``MPI`` activated will make use of ``mpicc``, ``mpicxx`` and  ``mpif90``, instead of ``gcc`` and co.. Has done for ``gcc``


``./configure CC=/usr/local/bin/gcc-9 FC=/usr/local/bin/gfortran-9 CXX=/usr/local/bin/g++-9``


PolyChord + MPI - symbol(s) not found for architecture x86_64
-------------------------------------------------------------

**Note:** This error is showing up only when PolyChord is being compiled with MPI support.
When installing ``PolyChord`` with ``MPI`` support, you may get this long list of error at the time of compiling the library:

.. code:: bash

  gfortran -shared abort.o array_utils.o calculate.o chordal_sampling.o clustering.o feedback.o generate.o ini.o interfaces.o mpi_utils.o nested_sampling.o params.o priors.o random_utils.o read_write.o run_time_info.o settings.o utils.o c_interface.o -o /Users/malavolta/Astro/CODE/others/PolyChord/lib/libchord.so
  Undefined symbols for architecture x86_64:
    "std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >::_M_create(unsigned long&, unsigned long)", referenced from:
        run_polychord(double (*)(double*, int, double*, int), void (*)(int, int, int, double*, double*, double*, double, double), Settings) in c_interface.o
        run_polychord(double (*)(double*, int, double*, int), void (*)(double*, double*, int), Settings) in c_interface.o
        run_polychord(double (*)(double*, int, double*, int), Settings) in c_interface.o
    ... [cut] ...
    "___gxx_personality_v0", referenced from:
        Dwarf Exception Unwind Info (__eh_frame) in c_interface.o
  ld: symbol(s) not found for architecture x86_64
  collect2: error: ld returned 1 exit status
  make[1]: *** [/Users/malavolta/Astro/CODE/others/PolyChord/lib/libchord.so] Error 1
  make: *** [/Users/malavolta/Astro/CODE/others/PolyChord/lib/libchord.so] Error 2

Change directory to ``src/polychord/``, copy the full command starting with ``gfortran -shared .. `` and add the end ``-lstdc++ -lc++``

.. code:: bash

  gfortran -shared abort.o array_utils.o calculate.o chordal_sampling.o clustering.o feedback.o generate.o ini.o interfaces.o mpi_utils.o nested_sampling.o params.o priors.o random_utils.o read_write.o run_time_info.o settings.o utils.o c_interface.o -o /Users/malavolta/Astro/CODE/others/PolyChord/lib/libchord.so -lstdc++ -lc++

Go back to the main directory and execute again ``make pypolychord``.



PolyChord - Segmentation fault
------------------------------

If you get this error using Conda/Anaconda and running ``python run_pypolychord.py``:

.. code:: bash

  *** Process received signal ***
  Signal: Segmentation fault: 11 (11)
  Signal code: Address not mapped (1)
  Failing at address: 0x2000000020
  [ 0] 0   libsystem_platform.dylib            0x00007fff7991cf5a _sigtramp + 26
  [ 1] 0   ???                                 0x000000005a21bf38 0x0 + 1512161080
  [ 2] 0   libsystem_c.dylib                   0x00007fff7972fc3d __vfprintf + 4711
  [ 3] 0   libsystem_c.dylib                   0x00007fff79757091 __v2printf + 473
  [ 4] 0   libsystem_c.dylib                   0x00007fff7973c4af _vsnprintf + 415
  [ 5] 0   libsystem_c.dylib                   0x00007fff7973c562 vsnprintf + 80
  [ 6] 0   libgfortran.3.dylib                 0x000000010e8b5d9b _gfortran_convert_char4_to_char1 + 3963
  *** End of error message ***

My guess is that ``lib/libchord.so`` has been compiled with different system libraries than those called by Conda. I don't have a solution for this problem, but using the system python seems the easiest workaround:

.. code:: bash

  /usr/bin/python pypolychord.py


MPI - Crash after a few iterations
----------------------------------

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
The slice sampling & clustering steps use a recursive procedure. The default memory allocated to recursive procedures is embarrassingly small (to guard against memory leaks).

MPI - No available slots
------------------------

The solution to this error:

.. code:: bash

  mpirun -np 8 python run_PyPolyChord.py

  --------------------------------------------------------------------------
  There are not enough slots available in the system to satisfy the 8 slots
  that were requested by the application:
    /usr/bin/python

  Either request fewer slots for your application, or make more slots available
  for use.
  --------------------------------------------------------------------------

Is quite simple: use a lower number after ``-np``. If `HyperThreading`_ is activated, the number of cores you see in your favorite task manager (or just ``htop``) is the number of _logical_ processor, while MPI cannot go further than the real number of cores in your machine.


Magically fixed problems
------------------------

Here I list some problems that I encountered in the past while installing some code, but that dind't appear anymore when a tried a new installation on more recent computers.

*ldd: command not found*

This error seems to be fixed in ``PolyChord v1.14``, but I'll leave it here for reference.

.. code:: bash

  /bin/sh: ldd: command not found

Open the ``Makefile`` in the main directory and substitute ``ldd`` with ``otool -L``. In version 1.12 this is the only line you have to change, from this:

.. code:: bash

  $(shell touch PyPolyChord/.ld_preload.sh; ldd $(LIB_DIR)/libchord.so | grep -o '/.*libmpi.so[^/]* ' | awk '{print "export LD_PRELOAD="$$1":$$LD_PRELOAD"}' > PyPolyChord/.ld_preload.sh)

to this:

.. code:: bash

    $(shell touch PyPolyChord/.ld_preload.sh; otool -L $(LIB_DIR)/libchord.so | grep -o '/.*libmpi.so[^/]* ' | awk '{print "export LD_PRELOAD="$$1":$$LD_PRELOAD"}' > PyPolyChord/.ld_preload.sh)

Executing ``make clean`` will not delete the library files created in the ``lib`` folder, so you have to delete them manually:

.. code:: bash

  make clean
  rm lib/polychord*.*
  make


Magically fixed MPI problems
----------------------------

Here I report errors I encountered so far when I try to install or run PolyChord in MPI mode. I had all these problems using ``PolyChord 1.12`` on ``Ubuntu 16.04 LTS``. Intalling and running ``PolyChord 1.14`` on ``Ubuntu 18.04 LTS`` didn't result in any of these errors. MAGIC!
For other errors, please refer to the README that comes with the source code.

*Broken MPI*

If you get the following errors when executing ``run_PyPolyChord.py`` , your MPI/OpenMPI installation is likely broken and you have to re-install it. You need to have a working MPI installation even when you are using PolyChord in single-CPU mode!

.. code:: bash

  [[INVALID],INVALID] ORTE_ERROR_LOG: A system-required executable either could not be found or was not executable by this user in file ess_singleton_module.c at line 231
  [[INVALID],INVALID] ORTE_ERROR_LOG: A system-required executable either could not be found or was not executable by this user in file ess_singleton_module.c at line 140
  [[INVALID],INVALID] ORTE_ERROR_LOG: A system-required executable either could not be found or was not executable by this user in file runtime/orte_init.c at line 128

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



.. _OpenMPI: https://www.open-mpi.org/
.. _Hyperthreading: https://superuser.com/questions/96001/why-does-my-intel-i7-920-display-8-cores-instead-of-4-cores
