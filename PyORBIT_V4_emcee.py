from classes.model_container import ModelContainer
from classes.input_parser import yaml_parser
import numpy as np
import emcee
from pyde.de import DiffEvol
import h5py
import cPickle as pickle
import os
import argparse


def save_to_hdf5(samples):

    h5f = h5py.File(emcee_dir_output + mc.planet_name + '.hdf5', "w")

    data_grp = h5f.create_group("data")
    data_grp.attrs.create('file_conf', data=file_conf)

    data_grp.create_dataset("starting_point", data=starting_point, compression="gzip")
    data_grp.create_dataset("starting_population", data=population, compression="gzip")

    emcee_grp = h5f.create_group("emcee")
    emcee_grp.attrs.create("nwalkers", data=mc.emcee_parameters['nwalkers'])
    emcee_grp.attrs.create("ndim", data=mc.ndim)
    emcee_grp.attrs.create("ndof", data=mc.ndof)

    emcee_grp.attrs.create("nsave", data=mc.emcee_parameters['nsave'])
    emcee_grp.attrs.create("nsample", data=samples)

    emcee_grp.create_dataset("bound", data=mc.bounds, compression="gzip")
    emcee_grp.create_dataset("chain", data=sampler.chain, compression="gzip")

    emcee_grp.create_dataset("lnprobability", data=sampler.lnprobability, compression="gzip")
    emcee_grp.create_dataset("acceptance_fraction", data=sampler.acceptance_fraction, compression="gzip")
    # emcee_grp.create_dataset("acor", data=sampler.acor, compression="gzip")

    h5f.close()

parser = argparse.ArgumentParser(prog='PyORBIT_V4_emcee.py', description='PyDE+emcee runner')
# parser.add_argument('-l', type=str, nargs='+', help='line identificator')
parser.add_argument('config_file', type=str, nargs=1, help='config file')

args = parser.parse_args()
file_conf = args.config_file[0]

mc = ModelContainer()

yaml_parser(file_conf, mc)
mc.initialize_model()

pyde_dir_output = './' + mc.planet_name + '/pyde/'
emcee_dir_output = './' + mc.planet_name + '/emcee/'

if not os.path.exists(emcee_dir_output):
    os.makedirs(emcee_dir_output)

if bool(mc.pcv.dynamical):
        mc.dynamical_model.prepare(mc, mc.pcv)

print
print 'Reference Time Tref: ', mc.Tref
print
print '*************************************************************'
print
print 'Dimensions = ', mc.ndim
print '   '
print 'Variable list:', mc.variable_list
print
print 'Variable bounds:', mc.bounds
print
print '*************************************************************'
print

mc.emcee_parameters['nwalkers'] = mc.ndim * mc.emcee_parameters['npop_mult']
if mc.emcee_parameters['nwalkers']%2 == 1: mc.emcee_parameters['nwalkers'] += 1

print 'Nwalkers = ', mc.emcee_parameters['nwalkers']

print mc.starting_point

if mc.starting_point_flag:
    mc.create_starting_point()
    starting_point = mc.starting_point
    print starting_point
    population = np.zeros([mc.emcee_parameters['nwalkers'], mc.ndim], dtype=np.double)
    for ii in xrange(0, mc.emcee_parameters['nwalkers']):
        population[ii, :] = np.random.normal(starting_point, 0.0000001)

else:
    if not os.path.exists(pyde_dir_output):
        os.makedirs(pyde_dir_output)

    if os.path.isfile(pyde_dir_output + 'pyde_pops.pick'):
        print os.path.isfile(pyde_dir_output + 'pyde_pops.pick')
        population = pickle.load(open(pyde_dir_output + 'pyde_pops.pick', 'rb'))
        starting_point = np.median(population, axis=0)
        mc.recenter_bounds(starting_point, population)

    else:
        print 'PyDE'
        de = DiffEvol(mc, mc.bounds, mc.emcee_parameters['nwalkers'], maximize=True)
        de.optimize(mc.pyde_parameters['ngen'])
        print 'PyDE completed'

        population = de.population
        starting_point = np.median(population, axis=0)
        pickle.dump(starting_point, open(pyde_dir_output + 'pyde_mean.pick', 'wb'))

        #np.savetxt(pyde_dir_output + 'pyDEout_original_bounds.dat', mc.bounds)
        #np.savetxt(pyde_dir_output + 'pyDEout_original_pops.dat', population)

        # bounds redefinition and fix for PyDE anomalous results
        if mc.recenter_bounds_flag:
            pickle.dump(mc.bounds, open(pyde_dir_output + 'bounds_orig.pick', 'wb'))
            pickle.dump(population, open(pyde_dir_output + 'pyde_pops_orig.pick', 'wb'))
            mc.recenter_bounds(starting_point, population)
            pickle.dump(mc.bounds, open(pyde_dir_output + 'bounds.pick', 'wb'))
            pickle.dump(population, open(pyde_dir_output + 'pyde_pops.pick', 'wb'))

            #np.savetxt(pyde_dir_output + 'pyDEout_redefined_bounds.dat', mc.bounds)
            #np.savetxt(pyde_dir_output + 'pyDEout_redefined_pops.dat', de.population)
            print 'REDEFINED BOUNDS'

        else:
            pickle.dump(mc.bounds, open(pyde_dir_output + 'bounds.pick', 'wb'))
            pickle.dump(population, open(pyde_dir_output + 'pyde_pops.pick', 'wb'))

print 'PyDE completed'
mc.results_resumen(starting_point)

    #json.dump(mc.variable_list, open('output/' + mc.planet_name + '_vlist.json', 'wb'))
pickle.dump(mc.variable_list, open(emcee_dir_output + 'vlist.pick', 'wb'))
pickle.dump(mc.scv.use_offset,  open(emcee_dir_output + 'scv_offset.pick', 'wb'))

if mc.emcee_parameters['MultiRun'] is not None:

    if os.path.isfile(emcee_dir_output + 'emcee_MR_pops.pick'):
        print os.path.isfile(emcee_dir_output + 'pyde_pops.pick')
        meds = pickle.load(open(emcee_dir_output + 'emcee_MR_meds.pick', 'rb'))
        population = pickle.load(open(emcee_dir_output + 'emcee_MR_pops.pick', 'rb'))
        print 'output from emcee exploratory runs retrieved'
    else:
        for ii in xrange(0, mc.emcee_parameters['MultiRun_iter']):
            print 'emcee exploratory run #', ii, ' of ', mc.emcee_parameters['MultiRun_iter']
            sampler = emcee.EnsembleSampler(mc.emcee_parameters['nwalkers'], mc.ndim, mc,
                                            threads=mc.emcee_parameters['nwalkers'])
            population, prob, state = sampler.run_mcmc(population, mc.emcee_parameters['MultiRun'])
            max_ind = np.argmax(prob)
            meds = population[max_ind, :]
            population = np.asarray([meds + 1e-4*np.random.randn(mc.ndim) for i in range(mc.emcee_parameters['nwalkers'])])
            sampler.reset()

        mc.results_resumen(meds)
        pickle.dump(meds, open(emcee_dir_output + 'emcee_MR_meds.pick', 'wb'))
        pickle.dump(population, open(emcee_dir_output + 'emcee_MR_pops.pick', 'wb'))
        print 'emcee exploratory runs completed'


print 'emcee'
state = None
sampler = emcee.EnsembleSampler(mc.emcee_parameters['nwalkers'], mc.ndim, mc, threads=mc.emcee_parameters['nwalkers'])

if mc.emcee_parameters['nsave'] > 0:
    print ' Saving temporary steps'
    niter = int(mc.emcee_parameters['nsteps']/mc.emcee_parameters['nsave'])
    sampled = 0
    for i in xrange(0, niter):
        population, prob, state = sampler.run_mcmc(population, mc.emcee_parameters['nsave'], thin=mc.emcee_parameters['thin'], rstate0=state)
        sampled += mc.emcee_parameters['nsave']
        save_to_hdf5(sampled)
        print sampled, '  steps completed, average lnprob:, ', np.median(prob)

else:
    population, prob, state = sampler.run_mcmc(population, mc.emcee_parameters['nsteps'], thin=mc.emcee_parameters['thin'])
    save_to_hdf5(mc.emcee_parameters['nsteps'])

print 'emcee completed'

