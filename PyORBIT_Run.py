import pyorbit
import argparse

if __name__ == '__main__':
    print 'This program is being run by itself'

    parser = argparse.ArgumentParser(prog='PyORBIT_run.py', description='PyORBIT runner')
    parser.add_argument('sampler', type=str, nargs=1, help='sampler (emcee or polychord)')
    parser.add_argument('config_file', type=str, nargs=1, help='config file')

    args = parser.parse_args()
    sampler = args.sampler[0]
    file_conf = args.config_file[0]

    config_in = pyorbit.yaml_parser(file_conf)

    sampler_keyword = {
        'polychord':['polychord', 'PolyChord', 'polychrod', 'poly'],
        'polychord_v1.9': ['polychord_1.9', 'PolyChord_1.9', 'polychrod_1.9', 'poly_1.9',
                           'polychord_1_9', 'PolyChord_1_9', 'polychrod_1_9', 'poly_1:9',
                           'polychord_v1.9', 'PolyChord_v1.9', 'polychrod_v1.9', 'poly_v1.9',
                           'polychord_v1_9', 'PolyChord_v1_9', 'polychrod_v1_9', 'poly_v1:9',
                           ],
        'emcee': ['emcee', 'MCMC', 'Emcee']
    }

    if sampler in sampler_keyword['emcee']:
        pyorbit.pyorbit_emcee(config_in)

    if sampler in sampler_keyword['polychord']:
        pyorbit.pyorbit_polychord(config_in)

    if sampler in sampler_keyword['polychord_v1.9']:
        pyorbit.pyorbit_polychord(config_in, pl_version='v1.9')

# This line was used to check if imprtation was working
# else:
#     print 'I am being imported from another module'
