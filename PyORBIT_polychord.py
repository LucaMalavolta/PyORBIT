import pyorbit
import argparse

if __name__ == '__main__':
    print 'This program is being run by itself'

    parser = argparse.ArgumentParser(prog='PyORBIT_emcee.py', description='PyDE+emcee runner')
    # parser.add_argument('-l', type=str, nargs='+', help='line identificator')
    parser.add_argument('config_file', type=str, nargs=1, help='config file')

    args = parser.parse_args()
    file_conf = args.config_file[0]
    config_in = pyorbit.yaml_parser(file_conf)

    pyorbit.pyorbit_emcee(config_in)

# This line was used to check if imprtation was working
# else:
#     print 'I am being imported from another module'
