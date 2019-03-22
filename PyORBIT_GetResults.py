import pyorbit
import argparse

if __name__ == '__main__':
    # print 'This program is being run by itself'

    parser = argparse.ArgumentParser(prog='PyORBIT_GetResults.py', description='PyDE+emcee runner')
    parser.add_argument('sampler', type=str, nargs=1, help='sampler (emcee or polychord)')
    parser.add_argument('config_file', type=str, nargs=1, help='config file')
    parser.add_argument('-p', type=str, nargs='?', default=False, help='Plot model files')
    parser.add_argument('-w', type=str, nargs='?', default=False, help='Write model files')
    parser.add_argument('-c', type=str, nargs='?', default=False, help='Save chains plots')
    parser.add_argument('-ln', type=str, nargs='?', default=False, help='Save ln_prob chain plot')
    parser.add_argument('-t', type=str, nargs='?', default=False, help='Compute and save Gelman-Rubin traces')
    parser.add_argument('-fc', type=str, nargs='?', default=False, help='Save full corellation plot - it may be slow!')
    parser.add_argument('-cc', type=str, nargs='?', default=False, help='Save corner plots of common variables')
    parser.add_argument('-dc', type=str, nargs='?', default=False, help='Save individual corner plots of reach dataset')
    parser.add_argument('-v', type=str, nargs='?', default=False, help='Write Veusz files for corner plot')
    parser.add_argument('-all_corners', type=str, nargs='?', default=False, help='Do all the corner plots')
    parser.add_argument('-all', type=str, nargs='?', default=False, help='Active all flags')
    parser.add_argument('-dfm_corner', type=str, nargs='?', default=False, help='Use DFM corner script for corner plots')

    plot_dictionary = {
        'plot_models': False,
        'write_models': False,
        'chains': False,
        'traces': False,
        'lnprob_chain': False,
        'full_correlation': False,
        'dataset_corner': False,
        'common_corner': False,
        'use_getdist': True,
        'veuz_corner_files': True
    }

    args = parser.parse_args()
    sampler = args.sampler[0]
    file_conf = args.config_file[0]

    if args.p is not False :
        plot_dictionary['plot_models'] = True
    if args.w is not False:
        plot_dictionary['write_models'] = True
    if args.c is not False :
        plot_dictionary['chains'] = True
    if args.t is not False:
        plot_dictionary['traces'] = True
    if args.fc is not False:
        plot_dictionary['lnprob_chain'] = True
    if args.fc is not False:
        plot_dictionary['full_correlation'] = True
    if args.dc is not False:
        plot_dictionary['dataset_corner'] = True
    if args.cc is not False:
        plot_dictionary['common_corner'] = True
    if args.v is not False:
        plot_dictionary['veuz_corner_files'] = True
    if args.dfm_corner is not False :
        plot_dictionary['use_getdist'] = False

    if args.all is not False:
        plot_dictionary['plot_models'] = True
        plot_dictionary['write_models'] = True
        plot_dictionary['lnprob_chain'] = True
        plot_dictionary['chains'] = True
        plot_dictionary['traces'] = True
        plot_dictionary['full_correlation'] = True
        plot_dictionary['common_corner'] = True
        plot_dictionary['dataset_corner'] = True
        plot_dictionary['veuz_corner_files'] = True

    if args.all_corners is not False:
        plot_dictionary['full_correlation'] = True
        plot_dictionary['common_corner'] = True
        plot_dictionary['dataset_corner'] = True

    config_in = pyorbit.yaml_parser(file_conf)

    pyorbit.pyorbit_getresults(config_in, sampler, plot_dictionary)
    #pyorbit.pyorbit_getresults_getdist(config_in, sampler, plot_dictionary)
