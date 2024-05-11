### MultiNest support

Link: [MultiNest repository by Johannes Buchner](https://github.com/JohannesBuchner/MultiNest)

```bash



Could NOT find LAPACK (missing: LAPACK_LIBRARIES)

Fedora:
sudo yum install lapack-devel
sudo yum install blas-devel

MPI not found, only non-MPI MultiNest libraries will be built

MPI toolkit: https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html
https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html#gs.wefwpr

source /opt/intel/oneapi/mpi/2021.12/env/vars.sh

pip install mpi4py

(original repository): https://github.com/farhanferoz/MultiNest

git clone https://github.com/JohannesBuchner/MultiNest.git
cd MultiNest/build

export LD_LIBRARY_PATH=/home/malavolta/CODE/others/MultiNest/lib/:$LD_LIBRARY_PATH


https://johannesbuchner.github.io/pymultinest-tutorial/install.html

git clone https://github.com/JohannesBuchner/PyMultiNest.git
pip install .