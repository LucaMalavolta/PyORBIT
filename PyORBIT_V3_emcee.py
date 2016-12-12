from PyORBIT_V3_Classes import *
import numpy as np
import emcee
from pyde.de import DiffEvol
import h5py
import cPickle as pickle
import os
import argparse
#import json

parser = argparse.ArgumentParser(prog='PyORBIT_V3_emcee.py', description='PyDE+emcee runner')
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



if bool(mc.pcv.dynamical):
    mc.pcv.prepare_dynamical(mc)

print 'Dimensions = ', mc.ndim
print '   '
print 'Variable list:', mc.variable_list
print
print 'Variable bounds:', mc.bounds
print
mc.nwalkers = mc.ndim * mc.npop_mult
if mc.nwalkers%2 == 1: mc.nwalkers += 1

print 'Nwalkers = ', mc.nwalkers

print mc.starting_point

if mc.starting_point_flag:
    mc.create_starting_point()
    starting_point = mc.starting_point
    population = np.zeros([mc.nwalkers, mc.ndim], dtype=np.double)
    for ii in xrange(0, mc.nwalkers):
        population[ii, :] = np.random.normal(starting_point, 0.0000001)

else:
    if os.path.isfile(dir_output + 'pyde_pops.pick'):
        print os.path.isfile(dir_output + 'pyde_pops.pick')
        population = pickle.load(open(dir_output + 'pyde_pops.pick', 'rb'))
        starting_point = np.median(population, axis=0)
        #pyde_mean = pickle.load(open(dir_output + 'pyde_mean.pick', 'rb'))
        mc.recenter_bounds(starting_point, population)

    else:
        print 'PyDE'
        de = DiffEvol(mc, mc.bounds, mc.nwalkers, maximize=True)
        de.optimize(mc.ngen)
        print 'PyDE completed'

        population = de.population
        starting_point = np.median(population, axis=0)
        pickle.dump(starting_point, open(dir_output + 'pyde_mean.pick', 'wb'))

        #np.savetxt(dir_output + 'pyDEout_original_bounds.dat', mc.bounds)
        #np.savetxt(dir_output + 'pyDEout_original_pops.dat', population)

        # bounds redefinition and fix for PyDE anomalous results
        if mc.recenter_bounds_flag:
            pickle.dump(mc.bounds, open(dir_output + 'bounds_orig.pick', 'wb'))
            pickle.dump(population, open(dir_output + 'pyde_pops_orig.pick', 'wb'))
            mc.recenter_bounds(starting_point, population)
            pickle.dump(mc.bounds, open(dir_output + 'bounds.pick', 'wb'))
            pickle.dump(population, open(dir_output + 'pyde_pops.pick', 'wb'))

            #np.savetxt(dir_output + 'pyDEout_redefined_bounds.dat', mc.bounds)
            #np.savetxt(dir_output + 'pyDEout_redefined_pops.dat', de.population)
            print 'REDEFINED BOUNDS'

        else:
            pickle.dump(mc.bounds, open(dir_output + 'bounds.pick', 'wb'))
            pickle.dump(population, open(dir_output + 'pyde_pops.pick', 'wb'))

print 'PyDE completed'
mc.results_resumen(starting_point)

    #json.dump(mc.variable_list, open('output/' + mc.planet_name + '_vlist.json', 'wb'))
pickle.dump(mc.variable_list, open(dir_output + 'vlist.pick', 'wb'))
pickle.dump(mc.scv.use_offset,  open(dir_output + 'scv_offset.pick', 'wb'))

print 'emcee'
sampler = emcee.EnsembleSampler(mc.nwalkers, mc.ndim, mc, threads=mc.nwalkers)
sampler.run_mcmc(population, mc.nsteps, thin=mc.thin)

print 'emcee completed'

# json.dump(mc.variable_list, open('output/' + mc.planet_name + '_vlist.json', 'wb'))
# pickle.dump(mc.variable_list, open(dir_output + 'vlist.pick', 'wb'))
# pickle.dump(mc.scv.use_offset,  open(dir_output + 'scv_offset.pick', 'wb'))

h5f = h5py.File(dir_output + mc.planet_name + '.hdf5', "w")

data_grp = h5f.create_group("data")
data_grp.attrs.create('file_conf',data=file_conf)

data_grp.create_dataset("starting_point", data=starting_point, compression="gzip")
data_grp.create_dataset("starting_population", data=population, compression="gzip")

emcee_grp = h5f.create_group("emcee")
emcee_grp.attrs.create("nwalkers", data=mc.nwalkers)
emcee_grp.attrs.create("ndim", data=mc.ndim)


emcee_grp.create_dataset("bound", data=mc.bounds, compression="gzip")
emcee_grp.create_dataset("chain", data=sampler.chain, compression="gzip")

emcee_grp.create_dataset("lnprobability", data=sampler.lnprobability, compression="gzip")
emcee_grp.create_dataset("acceptance_fraction", data=sampler.acceptance_fraction, compression="gzip")
emcee_grp.create_dataset("acor", data=sampler.acor, compression="gzip")

h5f.close()

