from PyORBIT_V2_Classes import *
import numpy as np
import emcee
from pyde.de import DiffEvol


file_conf="OC102_test.config"


mc = Model_Container()

Get_PyOrbit_Input(file_conf,mc)

mc.Model_setup()
mc.create_bounds()
print 'Dimensions = ', mc.ndim

mc.nwalkers = mc.ndim * mc.npop_mult
if (mc.nwalkers%2==1): mc.nwalkers+=1

print 'Nwalkers = ', mc.nwalkers

print 'PyDE'
de = DiffEvol(mc, mc.bounds, mc.nwalkers, maximize=True)
de.optimize(mc.ngen)
print 'PyDE completed'

pyde_mean = np.mean(de.population, axis=0)

# fix for PyDE anomalous results
for ii in xrange(0,mc.ndim):
    if np.amax(de.population[:,ii])-np.amin(de.population[:,ii]) < 10e-7 :
        range_restricted = (mc.bounds[ii,1]-mc.bounds[ii,0])/1000.
        min_bound = np.maximum((pyde_mean[ii]-range_restricted/2.0),mc.bounds[ii,0])
        max_bound = np.minimum((pyde_mean[ii]+range_restricted/2.0),mc.bounds[ii,1])
        de.population[:,ii] =  np.random.uniform(min_bound,max_bound,mc.nwalkers)

if mc.recenter_bounds_flag:
    de.population = mc.recenter_bounds(pyde_mean,de.population)
    np.savetxt('output/' + mc.planet_name + '_pyDEout_redefined_bounds.dat',mc.bounds)
    np.savetxt('output/' + mc.planet_name + '_pyDEout_redefined_pops.dat',de.population)


print pyde_mean
print mc.results_resumen(pyde_mean)
print 'REDEFINED BOUNDS'
print mc.bounds


print 'emcee'
sampler = emcee.EnsembleSampler(mc.nwalkers, mc.ndim, mc, threads=24)
sampler.run_mcmc(de.population, mc.nsteps, thin=mc.thin)

print 'emcee completed'
h5f = h5py.File('output/' + mc.planet_name + '.hdf5', "w")

data_grp = h5f.create_group("data")
data_grp.attrs.create('file_conf',data=file_conf)
data_grp.create_dataset("variable_list", data=mc.variable_list, compression="gzip")

emcee_grp = h5f.create_group("emcee")
emcee_grp.attrs.create("nwalkers", data=mc.nwalkers)
emcee_grp.attrs.create("ndim", data=mc.ndim)

emcee_grp.create_dataset("bound", data=mc.bounds, compression="gzip")
emcee_grp.create_dataset("chain", data=sampler.chain, compression="gzip")

emcee_grp.create_dataset("lnprobability", data=sampler.lnprobability, compression="gzip")
emcee_grp.create_dataset("acceptance_fraction", data=sampler.acceptance_fraction, compression="gzip")
emcee_grp.create_dataset("acor", data=sampler.acor, compression="gzip")

h5f.close()

