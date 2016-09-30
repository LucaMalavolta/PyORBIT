from PyORBIT_V2_Classes import *
import numpy as np
import emcee
from pyde.de import DiffEvol
import h5py
import cPickle as pickle
import os
import argparse
#import json

parser = argparse.ArgumentParser(prog='PyORBIT_V2_GetResults.py', description='Extract results from output MCMC')
# parser.add_argument('-l', type=str, nargs='+', help='line identificator')
parser.add_argument('config_file', type=str, nargs=1, help='config file')


args = parser.parse_args()

file_conf = args.config_file[0]

mc = ModelContainer()

yaml_parser(file_conf, mc)
#get_pyorbit_input(file_conf, mc)

dir_output = './' + mc.planet_name + '/'
if not os.path.exists(dir_output):
    os.makedirs(dir_output)

mc.create_bounds()
print 'Dimensions = ', mc.ndim

mc.nwalkers = mc.ndim * mc.npop_mult
if mc.nwalkers%2 == 1: mc.nwalkers += 1

print 'Nwalkers = ', mc.nwalkers

if os.path.isfile(dir_output + 'pyde_pops.pick'):
        population = pickle.load(open(dir_output + 'pyde_pops_orig.pick', 'rb'))
        pyde_mean = pickle.load(open(dir_output + 'pyde_mean.pick', 'rb'))
        bounds = pickle.load(open(dir_output + 'bounds_orig.pick', 'rb'))
print

mc.recenter_bounds(pyde_mean, population)

for ii in xrange(0,19):
    #print population
    sel = (population[:,ii]>bounds[ii,0]) &  (population[:,ii]<bounds[ii,1])
    print ii, np.sum(sel), pyde_mean[ii], bounds[ii]

e1 = population[:,6]**2 + population[:,7]**2
e2 = population[:,10]**2 + population[:,11]**2
e3 = population[:,15]**2 + population[:,16]**2
for ii in xrange(0,18):
    print e1[ii], e2[ii], e3[ii]


#print population[:, 15]
#print population[:, 16]

print
print
if os.path.isfile(dir_output + 'pyde_pops.pick'):
        population = pickle.load(open(dir_output + 'pyde_pops.pick', 'rb'))
        pyde_mean = pickle.load(open(dir_output + 'pyde_mean.pick', 'rb'))
        bounds = pickle.load(open(dir_output + 'bounds.pick', 'rb'))

for ii in xrange(0,19):
    sel = (population[:,ii]>bounds[ii,0]) &  (population[:,ii]<bounds[ii,1])
    print ii, np.sum(sel), pyde_mean[ii], bounds[ii]
print
e1 = population[:,6]**2 + population[:,7]**2
e2 = population[:,10]**2 + population[:,11]**2
e3 = population[:,15]**2 + population[:,16]**2
for ii in xrange(0,18):
    print e1[ii], e2[ii], e3[ii]
print population[:, 15]
print population[:, 16]
