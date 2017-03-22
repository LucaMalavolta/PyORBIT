from PyORBIT_V3_Classes import *
import numpy as np
import h5py
import cPickle as pickle
import scipy.optimize
import csv
import os
import corner
import argparse
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(prog='PyORBIT_V3_PyDEcheck.py', description='Extract results from PyDE output')
# parser.add_argument('-l', type=str, nargs='+', help='line identificator')
parser.add_argument('-i', type=str, nargs='+', required=True, help='config file')

args = parser.parse_args()

file_conf = args.i[0]

# file_conf = raw_input()

mc = ModelContainer()
yaml_parser(file_conf, mc)
#get_pyorbit_input(file_conf, mc)

mc.create_bounds()
print 'Dimensions = ', mc.ndim


dir_output = './' + mc.planet_name + '/'

# mc.variable_list = pickle.load(open(dir_output + 'vlist.pick', 'rb'))
# mc.scv.use_offset = pickle.load(open(dir_output + 'scv_offset.pick', 'rb'))

population = pickle.load(open(dir_output + 'pyde_pops.pick', 'rb'))
pyde_mean = pickle.load(open(dir_output + 'pyde_mean.pick', 'rb'))


mc.results_resumen(pyde_mean)

