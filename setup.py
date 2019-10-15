from distutils.core import setup
from Cython.Build import cythonize

setup(
    name="pyorbit",
    author="Luca Malavolta",
    ext_modules = cythonize("./pyorbit/*/*.pyx"),
    language_level="3"
)
