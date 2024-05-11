import pyorbit
import argparse
import sys

def pyorbit_results():
    # print 'This program is being run by itself'


    parser = argparse.ArgumentParser(prog='PyORBIT_GetResults.py', description='PyDE+emcee runner')
    parser.add_argument('sampler', type=str, nargs=1, help='sampler (emcee or polychord)')
    parser.add_argument('config_file', type=str, nargs=1, help='config file')
    parser.add_argument('-p', type=str, nargs='?', default=False, help='Plot model files')
    parser.add_argument('-w', type=str, nargs='?', default=False, help='Write model files')
    parser.add_argument('-wp', type=str, nargs='?', default=False, help='Write samples for orbital parameters')
    parser.add_argument('-ws', type=str, nargs='?', default=False, help='Write all samples')
    parser.add_argument('-c', type=str, nargs='?', default=False, help='Save chains plots')
    parser.add_argument('-a', type=str, nargs='?', default=False, help='Save ACF plots')
    parser.add_argument('-ln', type=str, nargs='?', default=False, help='Save ln_prob chain plot')
    parser.add_argument('-t', type=str, nargs='?', default=False, help='Compute and save Gelman-Rubin traces')
    parser.add_argument('-fc', type=str, nargs='?', default=False, help='Save full correlation plot - it may be slow!')
    parser.add_argument('-cc', type=str, nargs='?', default=False, help='Save corner plots of common parameters')
    parser.add_argument('-dc', type=str, nargs='?', default=False, help='Save individual corner plots of reach dataset')
    parser.add_argument('-v', type=str, nargs='?', default=False, help='Write Veusz files for corner plot')
    parser.add_argument('-vv', type=str, nargs='?', default=False, help='Write Veusz-compatible files')
    parser.add_argument('-ddp', type=str, nargs='?', default=False, help='Save dynesty default plots')
    parser.add_argument('-all_corners', type=str, nargs='?', default=False, help='Do all the corner plots')
    parser.add_argument('-all', type=str, nargs='?', default=False, help='Active all flags')
    parser.add_argument('-noacf', type=str, nargs='?', default=False, help='Skip ACF, even when all the flags are active')
    parser.add_argument('-dfm_corner', type=str, nargs='?', default=False, help='Use DFM corner script for corner plots')
    parser.add_argument('-getdist_corner', type=str, nargs='?', default=False, help='Use getdist script for corner plots')

    plot_dictionary = {
        'plot_models': False,
        'write_models': False,
        'write_planet_samples': False,
        'write_all_samples': False,
        'print_acf': True,
        'plot_acf': False,
        'chains': False,
        'traces': False,
        'lnprob_chain': False,
        'full_correlation': False,
        'dataset_corner': False,
        'common_corner': False,
        'use_getdist': False,
        'use_corner': False,
        'veuz_corner_files': False,
        'veuz_compatibility': False,
        'P_versus_lnprob': False,
        'dynesty_default_plots': False,
        'oversampled_models': ['transit',
                               'eclipse',
                               'phasecurve',
                               'eclipse_phasecurve',
                               'transit_eclipse_phasecurve',
                               'rv_keplerian']
    }

    # Moved here from pyorbit_getresults.py for consistency with PyORBIT_run.py
    sampler_keyword = {
        'multinest':['multinest', 'MultiNest', 'multi'],
        'polychord':['polychord', 'PolyChord', 'polychrod', 'poly'],
        'emcee': ['emcee', 'MCMC', 'Emcee'],
        'emcee_legacy': ['emcee_legacy', 'MCMC_legacy', 'Emcee_legacy'],
        'emcee_mpi': ['emcee_MPI', 'MCMC_MPI', 'Emcee_MPI','emcee_mpi', 'MCMC_mpi', 'Emcee_mpi'],
        'zeus_legacy': ['zeus', 'ZEUS', 'Zeus', 'zeus-mcmc', 'zeus_legacy', 'ZEUS_legacy', 'Zeus_legacy', 'zeus-mcmc_legacy'],
        'dynesty': ['dynesty', 'DyNesty', 'Dynesty', 'DYNESTY'],
        'dynesty_legacy': ['dynesty_legacy', 'DyNesty_legacy', 'Dynesty_legacy', 'DYNESTY_legacy'],
        'dynesty_static': ['dynesty_static', 'DyNesty_static', 'Dynesty_static', 'DYNESTY_static'],
        #'nestle': ['nestle', 'Nestle', 'NESTLE', 'nelste'],
        'ultranest': ['ultranest', 'UltraNest', 'Ultranest', 'ULTRANEST', 'ultra','Unest'],
        'emcee_warmstart': ['emcee_warmstart', 'emcee_warmup'],
        'ultranest_warmstart': ['ultranest_warmstart', 'UltraNest_WarmStart', 'Ultranest_Warmstart', 'ULTRANEST_WARMSTART', 'ultra_warm','Unest_WarmS', 'UltraWarm', 'ultrawarm'],

        #'optimize': ['optimize', 'scipy', 'Optimize', 'OPTIMIZE'],
    }

    unchained_samplers = ['ultranest_warmstart','polychord', 'multinest', 'dynesty', 'dynesty_legacy','dynesty_static', 'nestle', 'ultranest']

    args = parser.parse_args()
    sampler = args.sampler[0]
    file_conf = args.config_file[0]

    sampler_name = None
    for sampler_key, sampler_value in sampler_keyword.items():
        if sampler in sampler_value:
            sampler_name = sampler_key

    if not sampler_name:
        print(' *** Sampler not supported by GetResults, exiting *** ')
        quit()

    if args.p is not False :
        plot_dictionary['plot_models'] = True
    if args.w is not False:
        plot_dictionary['write_models'] = True
    if args.wp is not False:
        plot_dictionary['write_planet_samples'] = True
    if args.ws is not False:
        plot_dictionary['write_all_samples'] = True
    if args.a is not False :
        plot_dictionary['print_acf'] = True
        plot_dictionary['plot_acf'] = True
    if args.c is not False :
        plot_dictionary['chains'] = True
    if args.t is not False:
        plot_dictionary['traces'] = True
    if args.ln is not False:
        plot_dictionary['lnprob_chain'] = True
    if args.fc is not False:
        plot_dictionary['full_correlation'] = True
    if args.dc is not False:
        plot_dictionary['dataset_corner'] = True
    if args.cc is not False:
        plot_dictionary['common_corner'] = True
    if args.ddp is not False:
        plot_dictionary['dynesty_default_plots'] = True
    if args.v is not False:
        plot_dictionary['veuz_corner_files'] = True
    if args.vv is not False:
        plot_dictionary['veuz_compatibility'] = True
    if args.getdist_corner is not False :
        plot_dictionary['use_getdist'] = True
        plot_dictionary['full_correlation'] = True
    if args.dfm_corner is not False :
        plot_dictionary['use_getdist'] = False
        plot_dictionary['use_corner'] = True
        plot_dictionary['full_correlation'] = True

    if args.all is not False:
        plot_dictionary['plot_models'] = True
        plot_dictionary['write_models'] = True
        plot_dictionary['lnprob_chain'] = True
        plot_dictionary['print_acf'] = True
        plot_dictionary['plot_acf'] = True
        plot_dictionary['chains'] = True
        plot_dictionary['traces'] = True
        plot_dictionary['full_correlation'] = True
        plot_dictionary['common_corner'] = True
        plot_dictionary['dataset_corner'] = True
        plot_dictionary['dynesty_default_plots'] = True
        #plot_dictionary['veuz_corner_files'] = True
        plot_dictionary['write_planet_samples'] = True
        plot_dictionary['write_all_samples'] = True

    if args.all_corners is not False:
        plot_dictionary['full_correlation'] = True
        plot_dictionary['common_corner'] = True
        plot_dictionary['dataset_corner'] = True

    if args.noacf is not False:
        plot_dictionary['print_acf'] = False
        plot_dictionary['plot_acf'] = False

    if sampler_name in unchained_samplers:
        plot_dictionary['lnprob_chain'] = False
        plot_dictionary['print_acf'] = False
        plot_dictionary['plot_acf'] = False
        plot_dictionary['chains'] = False
        plot_dictionary['traces'] = False
        plot_dictionary['P_versus_lnprob'] = True


    print()
    print('PyORBIT v{0}'.format(pyorbit.__version__))
    print()
    print('Python version in use:')
    print(sys.version)

    config_in = pyorbit.subroutines.input_parser.yaml_parser(file_conf)

    pyorbit.pyorbit_getresults(config_in, sampler_name, plot_dictionary)
    #pyorbit.pyorbit_getresults_getdist(config_in, sampler, plot_dictionary)
