import pyorbit
import argparse

if __name__ == '__main__':
    print 'This program is being run by itself'

    parser = argparse.ArgumentParser(prog='PyORBIT_GetResults.py', description='PyDE+emcee runner')
    parser.add_argument('sampler', type=str, nargs=1, help='sampler (emcee or polychord)')
    parser.add_argument('config_file', type=str, nargs=1, help='config file')
    parser.add_argument('-p', type=str, nargs='?', default=False, help='Create plot files')
    parser.add_argument('-c', type=str, nargs='?', default=False, help='Create chains plots')
    parser.add_argument('-t', type=str, nargs='?', default=False, help='Create Gelman-Rubin traces')
    parser.add_argument('-f', type=str, nargs='?', default=False, help='Create model files using median parameters')
    parser.add_argument('-fc', type=str, nargs='?', default=False, help='Create full corellation plot - it may be slow!')
    parser.add_argument('-cc', type=str, nargs='?', default=False, help='Create corner plots of common variables')
    parser.add_argument('-dc', type=str, nargs='?', default=False, help='Create individual corner plots of reach dataset')
    parser.add_argument('-map', type=str, nargs='?', default=False, help='Create model files using MAP parameters')
    parser.add_argument('-all_corners', type=str, nargs='?', default=False, help='Do all the corner plots')
    parser.add_argument('-all', type=str, nargs='?', default=False, help='Active all flags')

    plot_dictionary = {
        'chains': False,
        'plot': False,
        'traces': False,
        'full_correlation': False,
        'common_corner': False,
        'dataset_corner': False,
        'model_files': False,
        'MAP_model_files': False
    }

    args = parser.parse_args()
    sampler = args.sampler[0]
    file_conf = args.config_file[0]

    if args.p is not False :
        plot_dictionary['plot'] = True
    if args.c is not False :
        plot_dictionary['chains'] = True
    if args.t is not False:
        plot_dictionary['traces'] = True
    if args.fc is not False:
        plot_dictionary['full_correlation'] = True
    if args.dc is not False:
        plot_dictionary['dataset_corner'] = True
    if args.cc is not False:
        plot_dictionary['common_corner'] = True
    if args.f is not False:
        plot_dictionary['model_files'] = True
    if args.map is not False:
        plot_dictionary['MAP_model_files'] = True
    if args.all is not False:
        plot_dictionary['plot'] = True
        plot_dictionary['chains'] = True
        plot_dictionary['traces'] = True
        plot_dictionary['common_corner'] = True
        plot_dictionary['dataset_corner'] = True
        plot_dictionary['full_correlation'] = True
        plot_dictionary['model_files'] = True
        plot_dictionary['MAP_model_files'] = True

    if args.all_corners is not False:
        plot_dictionary['common_corner'] = True
        plot_dictionary['dataset_corner'] = True
        plot_dictionary['full_correlation'] = True

    print plot_dictionary
    config_in = pyorbit.yaml_parser(file_conf)

    pyorbit.pyorbit_getresults(config_in, sampler, plot_dictionary)