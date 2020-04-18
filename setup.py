from distutils.core import setup
from Cython.Build import cythonize
import sys

setup(
    name="pyorbit",
    author="Luca Malavolta",
    ext_modules = cythonize("./pyorbit/*/*.pyx", compiler_directives={'language_level' : sys.version_info[0]}),
)
