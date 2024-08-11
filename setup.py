from setuptools import setup

# Inspired from here:
# https://hynek.me/articles/sharing-your-labor-of-love-pypi-quick-and-dirty/
# https://realpython.com/pypi-publish-python-package/#prepare-your-package-for-publication

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="pyorbit-package",
    version='10.6.0',
    author="Luca Malavolta",
	author_email = 'luca.malavolta@unipd.it',
	url = 'https://github.com/LucaMalavolta/PyORBIT',
	packages =['pyorbit', 'pyorbit.common', 'pyorbit.classes', 'pyorbit.models', 'pyorbit.samplers', 'pyorbit.subroutines'],
	license = 'MIT License',
	description ='PyORBIT: a code for exoplanet orbital parameters and stellar activity',
    long_description=long_description,
    long_description_content_type='text/markdown',
	classifiers = [
		'Development Status :: 5 - Production/Stable',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering',
        'Operating System :: OS Independent',
		'Programming Language :: Python :: 3'
		],
    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'pyorbit_run=pyorbit.pyorbit_run:pyorbit_run',
            'pyorbit_results=pyorbit.pyorbit_results:pyorbit_results'
        ]
    },
    zip_safe=False,
    install_requires=[
        'numpy==1.23.5',
        'numba>=0.55.2',
        'scipy>=1.8.1,<1.12.0',
        'matplotlib>=3.5.2,<3.8',
        'astropy>=5.1',
        'astroquery>=0.4',
        'pyerfa>=2.0',
        'argparse',
        'emcee>=3.1.2',
        'pyyaml',
        'h5py',
        'tqdm',
        'pygtc',
        'jaxlib>=0.4.28',
        'jax >=0.4.28',
        'tinygp>=0.3.0',
        'ordered-set'
        #'schwimmbad'
    ],
    setup_requires=['setuptools']
)
