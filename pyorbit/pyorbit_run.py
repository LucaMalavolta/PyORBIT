from __future__ import print_function
import pyorbit
import argparse
import os
import sys
import warnings

def pyorbit_run():

    print()
    print('PyORBIT v{0}'.format(pyorbit.__version__))
    print()
    print('Python version in use:')
    print(sys.version)
    #if sys.version_info[0] == 3 and sys.version_info[1] > 7:
    #    print('WARNING MESSAGES SUPPRESSED!')
    #print()

    #from packaging.version import Version
    #Version("2.3.1") < Version(tinygp.__version__)

    parser = argparse.ArgumentParser(prog='PyORBIT_run.py', description='PyORBIT runner')
    parser.add_argument('sampler', type=str, nargs=1, help='sampler (emcee or polychord)')
    parser.add_argument('config_file', type=str, nargs=1, help='config file')

    args = parser.parse_args()
    sampler = args.sampler[0]
    file_conf = args.config_file[0]

    config_in = pyorbit.subroutines.input_parser.yaml_parser(file_conf)

    sampler_keyword = {
        'pyde': ['pyde', 'PyDE', 'pyDE', 'global', 'global_solution', 'global_optimization'],
        'multinest':['multinest', 'MultiNest', 'multi'],
        'polychord':['polychord', 'PolyChord', 'polychrod', 'poly'],
        'emcee': ['emcee', 'MCMC', 'Emcee'],
        'emcee_legacy': ['emcee_legacy', 'MCMC_legacy', 'Emcee_legacy'],
        'emcee_mpi': ['emcee_MPI', 'MCMC_MPI', 'Emcee_MPI','emcee_mpi', 'MCMC_mpi', 'Emcee_mpi'],
        'zeus_legacy': ['zeus', 'ZEUS', 'Zeus', 'zeus-mcmc', 'zeus_legacy', 'ZEUS_legacy', 'Zeus_legacy', 'zeus-mcmc_legacy'],
        'dynesty': ['dynesty', 'DyNesty', 'Dynesty', 'DYNESTY'],
        'dynesty_legacy': ['dynesty_legacy', 'DyNesty_legacy', 'Dynesty_legacy', 'DYNESTY_legacy'],
        'dynesty_dryrun': ['dynesty_dryrun', 'DyNesty_dryrun', 'Dynesty_dryrun', 'DYNESTY_dryrun'],
        'dynesty_static': ['dynesty_static', 'DyNesty_static', 'Dynesty_static', 'DYNESTY_static'],
        'nestle': ['nestle', 'Nestle', 'NESTLE', 'nelste'],
        'ultranest': ['ultranest', 'UltraNest', 'Ultranest', 'ULTRANEST', 'ultra','Unest'],
        'ultranest_stepsampler': ['ultranest_stepsampler', 'UltraNest_StepSampler', 'Ultranest_Stepsampler', 'ULTRANEST_STEPSAMPLER', 'ultra_step','Unest_StepS'],
        'ultranest_warmstart': ['ultranest_warmstart', 'UltraNest_WarmStart', 'Ultranest_Warmstart', 'ULTRANEST_WARMSTART', 'ultra_warm','Unest_WarmS', 'UltraWarm', 'ultrawarm'],
        'optimize': ['optimize', 'scipy', 'Optimize', 'OPTIMIZE'],
    }

    if sampler in sampler_keyword['pyde']:
        pyorbit.pyorbit_pyde(config_in)

    if sampler in sampler_keyword['emcee']:
        pyorbit.pyorbit_emcee(config_in)

    if sampler in sampler_keyword['emcee_legacy']:
        pyorbit.pyorbit_emcee_legacy(config_in)

    if sampler in sampler_keyword['emcee_mpi']:
        pyorbit.pyorbit_emcee_mpi(config_in)

    if sampler in sampler_keyword['zeus_legacy']:
        pyorbit.pyorbit_zeus_legacy(config_in)

    if sampler in sampler_keyword['multinest']:
        config_in = pyorbit.subroutines.input_parser.yaml_fix_nested(config_in)
        pyorbit.pyorbit_multinest(config_in)

    if sampler in sampler_keyword['polychord']:
        config_in = pyorbit.subroutines.input_parser.yaml_fix_nested(config_in)
        pyorbit.pyorbit_polychord(config_in)

    if sampler in sampler_keyword['dynesty']:
        config_in = pyorbit.subroutines.input_parser.yaml_fix_nested(config_in)
        pyorbit.pyorbit_dynesty(config_in)

    if sampler in sampler_keyword['dynesty_dryrun']:
        config_in = pyorbit.subroutines.input_parser.yaml_fix_nested(config_in)
        pyorbit.pyorbit_dynesty(config_in, run_nested=False)

    if sampler in sampler_keyword['dynesty_static']:
        config_in = pyorbit.subroutines.input_parser.yaml_fix_nested(config_in)
        pyorbit.pyorbit_dynesty_static(config_in)

    if sampler in sampler_keyword['dynesty_legacy']:
        config_in = pyorbit.subroutines.input_parser.yaml_fix_nested(config_in)
        pyorbit.pyorbit_dynesty_legacy(config_in)

    if sampler in sampler_keyword['nestle']:
        config_in = pyorbit.subroutines.input_parser.yaml_fix_nested(config_in)
        pyorbit.pyorbit_nestle(config_in)

    if sampler in sampler_keyword['ultranest']:
        config_in = pyorbit.subroutines.input_parser.yaml_fix_nested(config_in)
        pyorbit.pyorbit_ultranest(config_in)

    if sampler in sampler_keyword['ultranest_stepsampler']:
        config_in = pyorbit.subroutines.input_parser.yaml_fix_nested(config_in)
        pyorbit.pyorbit_ultranest_stepsampler(config_in)

    if sampler in sampler_keyword['ultranest_warmstart']:
        config_in = pyorbit.subroutines.input_parser.yaml_fix_nested(config_in)
        pyorbit.pyorbit_ultranest_warmstart(config_in)

    if sampler in sampler_keyword['optimize']:
        pyorbit.pyorbit_optimize(config_in)
