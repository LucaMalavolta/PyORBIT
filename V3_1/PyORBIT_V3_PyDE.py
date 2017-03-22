from classes.model_container import ModelContainer
from classes.input_parser import yaml_parser
import numpy as np
import emcee
from pyde.de import DiffEvol
import h5py
import cPickle as pickle
import os
import argparse
#import json

parser = argparse.ArgumentParser(prog='PyORBIT_V3_PyDE.py', description='Run PyDE for frequentist analysis')
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
print '   '
print 'Variable list:', mc.variable_list
print
print 'Variable bounds:', mc.bounds
print
mc.nwalkers = mc.ndim * mc.npop_mult
if mc.nwalkers%2 == 1: mc.nwalkers += 1

print 'Nwalkers = ', mc.nwalkers


print 'PyDE'
de = DiffEvol(mc, mc.bounds, mc.nwalkers, maximize=True)
de.optimize(mc.ngen)
print 'PyDE completed'

population = de.population
pyde_mean = np.mean(population, axis=0)
pickle.dump(pyde_mean, open(dir_output + 'pyde_mean.pick', 'wb'))

#np.savetxt(dir_output + 'pyDEout_original_bounds.dat', mc.bounds)
#np.savetxt(dir_output + 'pyDEout_original_pops.dat', population)

# bounds redefinition and fix for PyDE anomalous results
if mc.recenter_bounds_flag:
    pickle.dump(mc.bounds, open(dir_output + 'bounds_orig.pick', 'wb'))
    pickle.dump(population, open(dir_output + 'pyde_pops_orig.pick', 'wb'))
    mc.recenter_bounds(pyde_mean, population)
    pickle.dump(mc.bounds, open(dir_output + 'bounds.pick', 'wb'))
    pickle.dump(population, open(dir_output + 'pyde_pops.pick', 'wb'))

    #np.savetxt(dir_output + 'pyDEout_redefined_bounds.dat', mc.bounds)
    #np.savetxt(dir_output + 'pyDEout_redefined_pops.dat', de.population)
    print 'REDEFINED BOUNDS'

else:
    pickle.dump(mc.bounds, open(dir_output + 'bounds.pick', 'wb'))
    pickle.dump(population, open(dir_output + 'pyde_pops.pick', 'wb'))

mc.results_resumen(pyde_mean)

#json.dump(mc.variable_list, open('output/' + mc.planet_name + '_vlist.json', 'wb'))
pickle.dump(mc.variable_list, open(dir_output + 'vlist.pick', 'wb'))
pickle.dump(mc.scv.use_offset,  open(dir_output + 'scv_offset.pick', 'wb'))
