import pyorbit
import argparse

if __name__ == '__main__':
    print 'This program is being run by itself'

    parser = argparse.ArgumentParser(prog='PyORBIT_run.py', description='PyORBIT runner')
    parser.add_argument('sampler', type=str, nargs=1, help='sample (emcee or polychord)')
    parser.add_argument('config_file', type=str, nargs=1, help='config file')

    args = parser.parse_args()
    sampler = args.sampler[0]
    file_conf = args.config_file[0]

    config_in = pyorbit.yaml_parser(file_conf)

    sampler_keyword = {
        'polychord':['polychord', 'PolyChord', 'polychrod', 'poly'],
        'emcee': ['emcee', 'MCMC', 'Emcee']
    }

    if sampler in sampler_keyword['emcee']:
        pyorbit.pyorbit_emcee(config_in)

    if sampler in sampler_keyword['polychord']:
        pyorbit.pyorbit_polychord(config_in)
# This line was used to check if imprtation was working
# else:
#     print 'I am being imported from another module'
