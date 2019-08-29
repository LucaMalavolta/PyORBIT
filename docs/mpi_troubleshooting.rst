.. _scratch_page:

Scratch page
============

PolyChord + MPI - check mpicc/mpicxx/gmpif90 versions
-----------------------------------------------------

``PolyChord`` installation with ``MPI`` activated will make use of ``mpicc``, ``mpicxx`` and  ``mpif90`` from ``openMPI``, instead of ``gcc``, ``g++`` and ``gfortran``. Therefore the same check that has been done for ``gcc`` & co has to be performed for those other compilers.
In my actual configuration (``macOS 10.14 (Mojave)`` + ``brew``) there is mix situation again, with ``mpif90`` coming from the ``brew`` installation and the other compilers relying on system libraries, with the difference however that no ``brew`` counterparts of ``mpicc`` and ``mpicxx`` seem to be available.
After a few failed attempts to install ``open MPI`` from ``brew``, I decided to go directly to for a manual installation from from the source code. The procedure is quite simple:

- Download the source code from the `Open MPI download page`_ (stable version is 4.0.1 at the moment of writing this guide)
- Untar the package and and move inside the source main folder:

.. code:: bash

  $ tar -xvf openmpi-4.0.1.tar.bz2
  $ cd openmpi-4.0.1

- Run this configure script, with the appropriate compiler options as in :ref:`gcc/gfortran/g++ versions-label`. You can select a different installation folder with ``--prefix``, if you don't want to overwrite the current installation, as suggested by the `Building Open MPI FAQ page`_

.. code:: bash

  $ ./configure --prefix=/usr/local/openmpi-4/ CC=/usr/local/bin/gcc-9 FC=/usr/local/bin/gfortran-9 CXX=/usr/local/bin/g++-9

- Run ``make install all`` to install ``Open MPI`` in the desired folder

More information can be found at the `Open MPI FAQ page`_

Before installin


Once the installation is complete, open the ``Makefile_gnu`` file in the main directory of ``PolyChord`` and change the first three lines after `ifdef MPI` accordingly, keeping in mind that the path of the executables will depend on your choice for the `--prefix` option:

.. code:: bash

  ifdef MPI
  FC = /usr/local/openmpi-4/bin/mpif90
  CC = /usr/local/openmpi-4/bin/mpicc
  CXX = /usr/local/openmpi-4/bin/mpicxx


Now unistall and re-install ``mpi4py``

Before compiling ``PolyChord``, open ``Makefile`` to check if the ``MPI`` flag is activated also for non-Linux installations:

.. code:: bash

  # Whether to use MPI
  ifeq "$(shell uname)" "Linux"
  MPI=1
  else
  MPI=1 # it was empy by default
  endif


.. _Open MPI download page: https://www.open-mpi.org/software
.. _Open MPI FAQ page: https://www.open-mpi.org/faq
.. _Building Open MPI FAQ page: https://www.open-mpi.org/faq/?category=building



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



**Note:** This error may show up if PolyChord is being compiled without MPI support, which is disabled by default.
