from setuptools import setup
import numpy as np

setup(
    name="pyorbit",
    version='9.0',
    author="Luca Malavolta",
	author_email = 'luca.malavolta@unipd.it',
	url = 'https://github.com/LucaMalavolta/PyORBIT',
	packages =['pyorbit', 'pyorbit.common', 'pyorbit.classes', 'pyorbit.models', 'pyorbit.samplers', 'pyorbit.subroutines'],
	#license = 'GNU GPLv3',
	#description ='',
	classifiers = [
		'Development Status :: 5 - Production/Stable',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering',
		'Programming Language :: Python'
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
     zip_safe=False
)
