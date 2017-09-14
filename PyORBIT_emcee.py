from classes.model_container import ModelContainer
from classes.input_parser import yaml_parser, pars_input
import numpy as np
import emcee
from pyde.de import DiffEvol
import h5py
import cPickle as pickle
import os
import argparse


def pyde_save_to_pickle(mc, population, starting_point, prefix=''):

    add_prefix = (prefix + '_' if prefix else '')
    pickle.dump(mc, open(mc.pyde_dir_output + add_prefix + "model_container.p", "wb"))
    pickle.dump(population, open(mc.pyde_dir_output + add_prefix + "population.p", "wb"))
    pickle.dump(starting_point, open(mc.pyde_dir_output + add_prefix + "starting_point.p", "wb"))


def pyde_load_from_cpickle(pyde_dir_output, prefix=''):

    add_prefix = (prefix + '_' if prefix else '')

    mc = pickle.load(open(pyde_dir_output + add_prefix + "model_container.p", "rb"))
    population = pickle.load(open(pyde_dir_output + add_prefix + "population.p", "rb"))
    starting_point = pickle.load(open(pyde_dir_output + add_prefix + "starting_point.p", "rb"))

    return mc, population, starting_point


def emcee_save_to_cpickle(mc, starting_point, population, prob, state, sampler, samples=None, prefix=None):

    if samples:
        mc.emcee_parameters['nsteps'] = samples
    add_prefix = (prefix + '_' if prefix else '')

    pickle.dump(mc, open(mc.emcee_dir_output + add_prefix + "model_container.p", "wb"))
    pickle.dump(starting_point, open(mc.emcee_dir_output + add_prefix + "starting_point.p", "wb"))
    pickle.dump(population, open(mc.emcee_dir_output + add_prefix + "starting_population.p", "wb"))
    pickle.dump(prob, open(mc.emcee_dir_output + add_prefix + "prob.p", "wb"))
    pickle.dump(state, open(mc.emcee_dir_output + add_prefix + "state.p", "wb"))
    pickle.dump(sampler.chain, open(mc.emcee_dir_output + add_prefix + "sampler_chain.p", "wb"))
    pickle.dump(sampler.lnprobability, open(mc.emcee_dir_output + add_prefix + "sampler_lnprobability.p", "wb"))
    pickle.dump(sampler.acceptance_fraction, open(mc.emcee_dir_output + add_prefix + "sampler_acceptance_fraction.p", "wb"))


def emcee_load_from_cpickle(emcee_dir_output, prefix=''):

    add_prefix = (prefix + '_' if prefix else '')

    mc = pickle.load(open(emcee_dir_output + add_prefix + "model_container.p", "rb"))
    starting_point = pickle.load(open(emcee_dir_output + add_prefix + "starting_point.p", "rb"))
    population = pickle.load(open(emcee_dir_output + add_prefix + "starting_population.p", "rb"))
    prob = pickle.load(open(mc.emcee_dir_output + add_prefix + "prob.p", "rb"))
    state = pickle.load(open(mc.emcee_dir_output + add_prefix + "state.p", "rb"))
    sampler_chain = pickle.load(open(emcee_dir_output + add_prefix + "sampler_chain.p", "rb"))
    sampler_lnprobability = pickle.load(open(emcee_dir_output + add_prefix + "sampler_lnprobability.p", "rb"))
    sampler_acceptance_fraction = pickle.load(open(emcee_dir_output + add_prefix + "sampler_acceptance_fraction.p", "rb"))

    return mc, starting_point, population, prob, state, \
            sampler_chain, sampler_lnprobability, sampler_acceptance_fraction


def emcee_flatchain(chain, nburnin, nthin):
    """flattening of the emcee chains with removal of burn-in"""
    _, d, _ = np.shape(chain)
    nburn = nburnin / nthin
    if nburn > d:
        nburn = d/4

    s = chain[:, nburn:, :].shape
    return chain[:, nburn:, :].reshape(s[0] * s[1], s[2])



def pyorbit_emcee(config_in):

    pyde_dir_output = './' + config_in['output'] + '/pyde/'
    emcee_dir_output = './' + config_in['output'] + '/emcee/'

    reloaded_pyde = False
    reloaded_emcee_multirun = False
    reloaded_emcee = False

    try:
        mc, population, starting_point = pyde_load_from_cpickle(pyde_dir_output, prefix='')
        reloaded_pyde = True
    except:
        pass

    try:
        mc, starting_point, population, _, _, sampler_chain, _, _ = \
            emcee_load_from_cpickle(emcee_dir_output, prefix='MR')
        reloaded_emcee_multirun = True
    except:
        pass

    try:
        mc, starting_point, population, _, _, sampler_chain, _, _ = emcee_load_from_cpickle(emcee_dir_output)
        reloaded_emcee = True
    except:
        pass

    print
    print 'reloaded_pyde: ', reloaded_pyde
    print 'reloaded_emcee_multirun: ', reloaded_emcee_multirun
    print 'reloaded_emcee: ', reloaded_emcee
    print

    if reloaded_emcee:
        """ There's no need to do anything"""
        flatchain = emcee_flatchain(sampler_chain, mc.emcee_parameters['nburn'], mc.emcee_parameters['thin'])
        mc.results_resumen(flatchain)
        return

    reloaded_mc = reloaded_pyde or reloaded_emcee_multirun or reloaded_emcee_multirun

    if not reloaded_mc:
        mc = ModelContainer()
        pars_input(config_in, mc)
        mc.initialize_model()

        mc.pyde_dir_output = pyde_dir_output
        mc.emcee_dir_output = emcee_dir_output

        mc.emcee_parameters['nwalkers'] = mc.ndim * mc.emcee_parameters['npop_mult']
        if mc.emcee_parameters['nwalkers']%2 == 1: mc.emcee_parameters['nwalkers'] += 1

        if mc.dynamical_model is not None:
            mc.dynamical_model.prepare(mc)

    if not os.path.exists(mc.emcee_dir_output):
        os.makedirs(mc.emcee_dir_output)

    print
    print 'Reference Time Tref: ', mc.Tref
    print
    print 'Dimensions = ', mc.ndim
    print 'Nwalkers = ', mc.emcee_parameters['nwalkers']
    print
    print '*************************************************************'
    print

    if mc.starting_point_flag:
        mc.create_starting_point()
        starting_point = mc.starting_point
        population = np.zeros([mc.emcee_parameters['nwalkers'], mc.ndim], dtype=np.double)
        for ii in xrange(0, mc.emcee_parameters['nwalkers']):
            population[ii, :] = np.random.normal(starting_point, 0.0000001)
        reloaded_pyde = True

    if not reloaded_pyde:
        if not os.path.exists(mc.pyde_dir_output):
            os.makedirs(mc.pyde_dir_output)

        print 'PyDE'
        de = DiffEvol(mc, mc.bounds, mc.emcee_parameters['nwalkers'], maximize=True)
        de.optimize(mc.pyde_parameters['ngen'])
        print 'PyDE completed'

        population = de.population
        starting_point = np.median(population, axis=0)

        #np.savetxt(mc.pyde_dir_output + 'pyDEout_original_bounds.dat', mc.bounds)
        #np.savetxt(mc.pyde_dir_output + 'pyDEout_original_pops.dat', population)

        """ bounds redefinition and fix for PyDE anomalous results """
        if mc.recenter_bounds_flag:
            pyde_save_to_pickle(mc, population, starting_point, prefix='orig')

            mc.recenter_bounds(starting_point)
            population = mc.fix_population(starting_point, population)
            starting_point = np.median(population, axis=0)

            print 'REDEFINED BOUNDS'

        pyde_save_to_pickle(mc, population, starting_point)

        print 'PyDE completed'

    mc.results_resumen(starting_point)

    if mc.emcee_parameters['multirun'] and not reloaded_emcee_multirun:

        for ii in xrange(0, mc.emcee_parameters['multirun_iter']):
            print 'emcee exploratory run #', ii, ' of ', mc.emcee_parameters['multirun_iter']
            sampler = emcee.EnsembleSampler(mc.emcee_parameters['nwalkers'], mc.ndim, mc,
                                            threads=mc.emcee_parameters['nwalkers'])
            population, prob, state = sampler.run_mcmc(population, mc.emcee_parameters['multirun'])
            flatchain = emcee_flatchain(sampler.chain, mc.emcee_parameters['nburn'], mc.emcee_parameters['thin'])
            mc.results_resumen(flatchain)

            max_ind = np.argmax(prob)
            starting_point = population[max_ind, :]
            population = np.asarray([starting_point + 1e-4*np.random.randn(mc.ndim) for i in range(mc.emcee_parameters['nwalkers'])])
            sampler.reset()

            emcee_save_to_cpickle(mc, starting_point, population, prob, state, sampler, prefix='MR_'+repr(ii))

        emcee_save_to_cpickle(mc, starting_point, population, prob, state, sampler, prefix='MR')

        flatchain = emcee_flatchain(sampler.chain, mc.emcee_parameters['nburn'], mc.emcee_parameters['thin'])
        mc.results_resumen(flatchain)

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
            emcee_save_to_cpickle(mc, starting_point, population, prob, state, sampler, samples=sampled)

            print sampled, '  steps completed, average lnprob:, ', np.median(prob)

            flatchain = emcee_flatchain(sampler.chain, mc.emcee_parameters['nburn'], mc.emcee_parameters['thin'])
            mc.results_resumen(flatchain)

    else:
        population, prob, state = sampler.run_mcmc(population, mc.emcee_parameters['nsteps'], thin=mc.emcee_parameters['thin'])

        emcee_save_to_cpickle(mc, starting_point, population, prob, state, sampler)

        flatchain = emcee_flatchain(sampler.chain, mc.emcee_parameters['nburn'], mc.emcee_parameters['thin'])
        mc.results_resumen(flatchain)

    print 'emcee completed'


if __name__ == '__main__':
    print 'This program is being run by itself'

    parser = argparse.ArgumentParser(prog='PyORBIT_emcee.py', description='PyDE+emcee runner')
    # parser.add_argument('-l', type=str, nargs='+', help='line identificator')
    parser.add_argument('config_file', type=str, nargs=1, help='config file')

    args = parser.parse_args()
    file_conf = args.config_file[0]
    config_in = yaml_parser(file_conf)

    pyorbit_emcee(config_in)

else:
    print 'I am being imported from another module'