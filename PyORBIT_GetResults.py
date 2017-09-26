from classes.model_container import ModelContainer
from classes.input_parser import yaml_parser, pars_input
from classes.io_subroutines import pyde_save_to_pickle, pyde_load_from_cpickle, \
    emcee_save_to_cpickle, emcee_load_from_cpickle, emcee_flatchain, emcee_flatlnprob
import numpy as np
import emcee
from pyde.de import DiffEvol
import os
import argparse
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
import sys
mpl.use('Agg')
from matplotlib import pyplot as plt
import corner
import classes.constants as constants



def pyorbit_getresults(config_in, sampler):

    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    G_grav = constants.Gsi # Gravitational Constants in SI system [m^3/kg/s^2]
    G_ttvfast = constants.Giau  # G [AU^3/Msun/d^2]
    M_SJratio = constants.Msjup
    M_SEratio = constants.Msear
    M_JEratio = constants.Mjear

    R_SJratio = constants.Rsjup
    R_JEratio = constants.Rjear
    R_SEratio = constants.Rsjup * constants.Rjear

    M_sun = constants.Msun
    M_jup = constants.Mjup

    Mu_sun = constants.Gsi * constants.Msun
    seconds_in_day = constants.d2s
    AU_km = constants.AU
    AUday2ms = AU_km / seconds_in_day * 1000.0


    sample_keyword = {
        'polychord':['polychord', 'PolyChord', 'polychrod', 'poly'],
        'emcee': ['emcee', 'MCMC', 'Emcee']
    }


    if sampler in sample_keyword['emcee']:
        print 'Wowow'

        dir_input = './' + config_in['output'] + '/emcee/'
        dir_output = './' + config_in['output'] + '/emcee_plot/'
        os.system('mkdir -p ' + dir_output)

        mc, starting_point, population, prob, state, \
        sampler_chain, sampler_lnprobability, sampler_acceptance_fraction = \
            emcee_load_from_cpickle(dir_input)

        theta_dictionary = mc.get_theta_dictionary()

        nburnin = mc.emcee_parameters['nburn']
        nthin = mc.emcee_parameters['thin']

        flat_chain = emcee_flatchain(sampler_chain, nburnin, nthin)
        flat_lnprob = emcee_flatlnprob(sampler_lnprobability, nburnin, nthin)

        lnprob_median = np.median(flat_lnprob)

        fig = plt.figure(figsize=(12, 12))
        plt.xlabel('$\ln \mathcal{L}$')
        plt.plot(flat_lnprob.T, '-', alpha=0.5)
        plt.axhline(lnprob_median)
        plt.axvline(nburnin/nthin, c='r')
        plt.savefig(dir_output + 'LNprob_chain.png', bbox_inches='tight', dpi=300)
        plt.close(fig)
        print 'LNprob median = ', lnprob_median

        #for ii in xrange(0, mc.ndim):
        for theta_name, ii in  theta_dictionary.iteritems():
            #print theta_name, chain_med[ii, 0], ' +\sigma ', chain_med[ii, 1], ' -\sigma ', chain_med[ii, 2]
            file_name = dir_output + 'chain_' + repr(ii) + '_' + theta_name + '.png'
            fig = plt.figure(figsize=(12, 12))
            plt.plot(sampler_chain[:, :, ii].T, '-', alpha=0.5)
            plt.axvline(nburnin/nthin, c='r')
            #plt.axhline(chain_med[ii, 0], c='k')
            plt.savefig(file_name, bbox_inches='tight', dpi=300)
            plt.close(fig)


if __name__ == '__main__':
    print 'This program is being run by itself'

    parser = argparse.ArgumentParser(prog='PyORBIT_GetResults.py', description='PyDE+emcee runner')
    parser.add_argument('sample', type=str, nargs=1, help='sample (emcee or polychord)')
    parser.add_argument('config_file', type=str, nargs=1, help='config file')

    args = parser.parse_args()
    sampler = args.sample[0]
    file_conf = args.config_file[0]

    config_in = yaml_parser(file_conf)

    pyorbit_getresults(config_in, sampler)

else:
    print 'I am being imported from another module'
