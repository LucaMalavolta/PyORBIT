from __future__ import print_function

from pyorbit.subroutines.common import *
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
import pygtc
import pyorbit.subroutines.kepler_exo as kepler_exo
import pyorbit.subroutines.common as common
import pyorbit.subroutines.results_analysis as results_analysis
import h5py
import re
from pyorbit.classes.model_container_multinest import ModelContainerMultiNest
from pyorbit.classes.model_container_polychord import ModelContainerPolyChord
from pyorbit.classes.model_container_emcee import ModelContainerEmcee
from pyorbit.classes.model_container_zeus import ModelContainerZeus

from pyorbit.subroutines.input_parser import pars_input
from pyorbit.subroutines.io_subroutines import *
import numpy as np
import os
import matplotlib as mpl
import sys

mpl.use('Agg')


__all__ = ["pyorbit_getresults"]


def pyorbit_getresults(config_in, sampler_name, plot_dictionary):

    print(' LaTeX disabled by default')
    use_tex = False
    plt.rc('text', usetex=use_tex)


    if plot_dictionary['use_corner']:
        import corner

    if plot_dictionary['use_getdist']:
        from getdist import plots, MCSamples

    # plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
    plt.rcParams["font.family"] = "Times New Roman"

    oversampled_models = plot_dictionary['oversampled_models']

    if sampler_name == 'emcee':

        dir_input = './' + config_in['output'] + '/emcee/'
        dir_output = './' + config_in['output'] + '/emcee_plot/'
        os.system('mkdir -p ' + dir_output)

        mc, starting_point, population, prob, state, \
            sampler_chain, sampler_lnprobability, sampler_acceptance_fraction, _, _ = \
            emcee_load_from_cpickle(dir_input)

        if hasattr(mc.emcee_parameters, 'version'):
            emcee_version = mc.emcee_parameters['version'][0]
        else:
            import emcee
            emcee_version = emcee.__version__[0]

        pars_input(config_in, mc, reload_emcee=True)

        if hasattr(mc.emcee_parameters, 'version'):
            emcee_version = mc.emcee_parameters['version'][0]
        else:
            import emcee
            emcee_version = emcee.__version__[0]

        mc.model_setup()
        theta_dictionary = results_analysis.get_theta_dictionary(mc)

        nburnin = int(mc.emcee_parameters['nburn'])
        nthin = int(mc.emcee_parameters['thin'])
        nsteps = int(sampler_chain.shape[1] * nthin)
        nwalkers = mc.emcee_parameters['nwalkers']

        """ Computing a new burn-in if the computation has been interrupted suddenly"""
        nburn, modified = emcee_burnin_check(sampler_chain, nburnin, nthin)

        if modified:
            print()
            print('WARNING: burn-in value is larger than the length of the chains, resized to 1/4 of the chain length')
            print('new burn-in will be used for statistical analysis, but kept in the plots as a reminder of your mistake')

        flat_chain = emcee_flatchain(sampler_chain, nburnin, nthin)
        flat_lnprob, sampler_lnprob = emcee_flatlnprob(
            sampler_lnprobability, nburnin, nthin, population, nwalkers)

        flat_BiC = -2 * flat_lnprob + mc.ndim * np.log(mc.ndata)

        lnprob_med = common.compute_value_sigma(flat_lnprob)
        chain_med = common.compute_value_sigma(flat_chain)
        chain_MAP, lnprob_MAP = common.pick_MAP_parameters(
            flat_chain, flat_lnprob)

        n_samplings, n_pams = np.shape(flat_chain)

        print()
        print('emcee version: ', emcee.__version__)
        if mc.emcee_parameters['version'] == '2':
            print('WARNING: upgrading to version 3 is strongly advised')
        print()
        print(' Reference Time Tref: {}'.format(mc.Tref))
        print()
        print(' Dimensions = {}'.format(mc.ndim))
        print(' Nwalkers = {}'.format(mc.emcee_parameters['nwalkers']))
        print()
        print(' Steps: {}'.format(nsteps))

        results_analysis.print_integrated_ACF(
            sampler_chain, theta_dictionary, nthin)

    if sampler_name == 'zeus':

        dir_input = './' + config_in['output'] + '/zeus/'
        dir_output = './' + config_in['output'] + '/zeus_plot/'
        os.system('mkdir -p ' + dir_output)

        mc, starting_point, population, prob, state, \
            sampler_chain, sampler_lnprobability, sampler_acceptance_fraction, _, _ = \
            zeus_load_from_cpickle(dir_input)

        pars_input(config_in, mc, reload_zeus=True)

        if hasattr(mc.zeus_parameters, 'version'):
            zeus_version = mc.zeus_parameters['version'][0]
        else:
            import zeus
            zeus = zeus.__version__[0]

        mc.model_setup()
        theta_dictionary = results_analysis.get_theta_dictionary(mc)

        nburnin = int(mc.zeus_parameters['nburn'])
        nthin = int(mc.zeus_parameters['thin'])
        nsteps = int(sampler_chain.shape[1] * nthin)
        nwalkers = mc.zeus_parameters['nwalkers']

        """ Computing a new burn-in if the computation has been interrupted suddenly"""
        nburn, modified = emcee_burnin_check(sampler_chain, nburnin, nthin)

        if modified:
            print()
            print('WARNING: burn-in value is larger than the length of the chains, resized to 1/4 of the chain length')
            print('new burn-in will be used for statistical analysis, but kept in the plots as a reminder of your mistake')

        flat_chain = emcee_flatchain(sampler_chain, nburnin, nthin)
        flat_lnprob, sampler_lnprob = emcee_flatlnprob(
            sampler_lnprobability, nburnin, nthin, population, nwalkers)

        flat_BiC = -2 * flat_lnprob + mc.ndim * np.log(mc.ndata)

        lnprob_med = common.compute_value_sigma(flat_lnprob)
        chain_med = common.compute_value_sigma(flat_chain)
        chain_MAP, lnprob_MAP = common.pick_MAP_parameters(
            flat_chain, flat_lnprob)

        n_samplings, n_pams = np.shape(flat_chain)

        print()
        print('zeus version: ', zeus.__version__)
        print()
        print(' Reference Time Tref: {}'.format(mc.Tref))
        print()
        print(' Dimensions = {}'.format(mc.ndim))
        print(' Nwalkers = {}'.format(mc.zeus_parameters['nwalkers']))
        print()
        print(' Steps: {}'.format(nsteps))

        results_analysis.print_integrated_ACF(
            sampler_chain, theta_dictionary, nthin)


    if sampler_name == 'multinest':

        dir_input = './' + config_in['output'] + '/multinest/'
        dir_output = './' + config_in['output'] + '/multinest_plot/'
        os.system('mkdir -p ' + dir_output)

        mc = nested_sampling_load_from_cpickle(dir_input)

        mc.model_setup()
        mc.initialize_logchi2()
        results_analysis.results_resumen(mc, None, skip_theta=True)

        theta_dictionary = results_analysis.get_theta_dictionary(mc)

        data_in = np.genfromtxt(dir_input + 'post_equal_weights.dat')
        flat_lnprob = data_in[:, -1]
        flat_chain = data_in[:, :-1]
        # nsample = np.size(flat_lnprob)
        n_samplings, n_pams = np.shape(flat_chain)

        lnprob_med = common.compute_value_sigma(flat_lnprob)
        chain_med = common.compute_value_sigma(flat_chain)
        chain_MAP, lnprob_MAP = common.pick_MAP_parameters(
            flat_chain, flat_lnprob)

        print()
        print(' Reference Time Tref: {}'.format(mc.Tref))
        print()
        print(' Dimensions: {}'.format(mc.ndim))
        print()
        print(' Samples: {}'.format(n_samplings))

    if sampler_name == 'polychord':

        dir_input = './' + config_in['output'] + '/polychord/'
        dir_output = './' + config_in['output'] + '/polychord_plot/'
        os.system('mkdir -p ' + dir_output)

        mc = nested_sampling_load_from_cpickle(dir_input)

        # pars_input(config_in, mc)

        mc.model_setup()
        mc.initialize_logchi2()
        results_analysis.results_resumen(mc, None, skip_theta=True)

        theta_dictionary = results_analysis.get_theta_dictionary(mc)

        data_in = np.genfromtxt(dir_input + 'pyorbit_equal_weights.txt')
        flat_lnprob = data_in[:, 1]
        flat_chain = data_in[:, 2:]
        # nsample = np.size(flat_lnprob)

        n_samplings, n_pams = np.shape(flat_chain)

        lnprob_med = common.compute_value_sigma(flat_lnprob)
        chain_med = common.compute_value_sigma(flat_chain)

        chain_MAP, lnprob_MAP = common.pick_MAP_parameters(
            flat_chain, flat_lnprob)

        print()
        print(' Reference Time Tref: {}'.format(mc.Tref))
        print()
        print(' Dimensions: {}'.format(mc.ndim))
        print()
        print(' Samples: {}'.format(n_samplings))

    if sampler_name == 'dynesty':

        from dynesty import utils as dyfunc
        from dynesty import plotting as dyplot

        dir_input = './' + config_in['output'] + '/dynesty/'
        dir_output = './' + config_in['output'] + '/dynesty_plot/'
        os.system('mkdir -p ' + dir_output)

        mc = nested_sampling_load_from_cpickle(dir_input)

        mc.model_setup()
        mc.initialize_logchi2()
        results_analysis.results_resumen(mc, None, skip_theta=True)

        theta_dictionary = results_analysis.get_theta_dictionary(mc)

        pfrac = mc.nested_sampling_parameters['pfrac']

        try:
            results = dynesty_results_maxevidence_load_from_cpickle(dir_input)
            pfrac = 0.00
            print('Model evidence from dynesty run with posterior/evidence split = {0:4.3f}'.format(pfrac))

            labels_array = [None] * len(theta_dictionary)
            for key_name, key_value in theta_dictionary.items():
                labels_array[key_value] = re.sub('_', '-', key_name)
            print(results)

            print()
            # Plot a summary of the run.
            try:
                print('Plot a summary of the run.')
                rfig, raxes = dyplot.runplot(results)
                rfig.savefig(dir_output + 'dynesty_results_maxevidence_summary.pdf', bbox_inches='tight', dpi=300)
                plt.close(rfig)
            except:
                print('Unable to plot a summary of the run using the internal dynesty routine - skipped')

            # Plot traces and 1-D marginalized posteriors.
            try:
                print('Plot traces and 1-D marginalized posteriors.')
                tfig, taxes = dyplot.traceplot(results, labels=labels_array)
                tfig.savefig(dir_output + 'dynesty_results_maxevidence_traceplot.pdf', bbox_inches='tight', dpi=300)
                plt.close(tfig)
            except:
                print('Unable to plot traces and 1-D marginalized posteriors using the internal dynesty routine - skipped')

            # Plot the 2-D marginalized posteriors.
            try:
                print('Plot the 2-D marginalized posteriors.')
                cfig, caxes = dyplot.cornerplot(results, labels=labels_array)
                cfig.savefig(dir_output + 'dynesty_results_maxevidence_cornerplot.pdf', bbox_inches='tight', dpi=300)
                plt.close(cfig)
            except:
                print('Unable to plot the 2-D marginalized posteriors using the internal dynesty routine - skipped')

            pfrac = 1.00

        except:
            results = dynesty_results_load_from_cpickle(dir_input)


        #taken from dynesty/dynesty/results.py  but without the nlive point causing an error
        res = ("niter: {:d}\n"
                "ncall: {:d}\n"
                "eff(%): {:6.3f}\n"
                "logz: {:6.3f} +/- {:6.3f}"
                .format(results.niter, sum(results.ncall),
                        results.eff, results.logz[-1], results.logzerr[-1]))

        print()
        print('Summary - \n=======\n'+res)

        try:
            # Generate a new set of results with statistical+sampling uncertainties.
            results_sim = dyfunc.simulate_run(results)
            res = ("niter: {:d}\n"
                    "ncall: {:d}\n"
                    "eff(%): {:6.3f}\n"
                    "logz: {:6.3f} +/- {:6.3f}"
                    .format(results_sim.niter, sum(results_sim.ncall),
                            results_sim.eff, results_sim.logz[-1], results_sim.logzerr[-1]))

            print()
            print('Summary - statistical+sampling errors - \n=======\n'+res)
        except AttributeError:
            print('Computation of statistical+sampling errors skipped - workaround after dynesty>1.2 update')


        results = dynesty_results_load_from_cpickle(dir_input)
        print('Posteriors analysis from dynesty run with posterior/evidence split = {0:4.3f}'.format(pfrac))

        labels_array = [None] * len(theta_dictionary)
        for key_name, key_value in theta_dictionary.items():
            labels_array[key_value] = re.sub('_', '-', key_name)

        print()

        # Plot a summary of the run.
        try:
            print('Plot a summary of the run.')
            #span= [(0.0, 2100.0), (0.0, 1.05), (0.0, 0.10444691225380567), (0.0, 1000)]
            #rfig, raxes = dyplot.runplot(results, span=span)
            rfig, raxes = dyplot.runplot(results)
            rfig.savefig(dir_output + 'dynesty_results_summary.pdf', bbox_inches='tight', dpi=300)
            plt.close(rfig)
        except:
            print('Unable to plot a summary of the run using the internal dynesty routine - skipped')

        # Plot traces and 1-D marginalized posteriors.
        try:
            print('Plot traces and 1-D marginalized posteriors.')
            tfig, taxes = dyplot.traceplot(results, labels=labels_array)
            tfig.savefig(dir_output + 'dynesty_results_traceplot.pdf', bbox_inches='tight', dpi=300)
            plt.close(tfig)
        except:
            print('Unable to plot traces and 1-D marginalized posteriors using the internal dynesty routine - skipped')

        # Plot the 2-D marginalized posteriors.
        try:
            print('Plot the 2-D marginalized posteriors.')
            cfig, caxes = dyplot.cornerplot(results, labels=labels_array)
            cfig.savefig(dir_output + 'dynesty_results_cornerplot.pdf', bbox_inches='tight', dpi=300)
            plt.close(cfig)
        except:
            print('Unable to plot the 2-D marginalized posteriors using the internal dynesty routine - skipped')

        # Extract sampling results.
        samples = results.samples  # samples

        # normalized weights
        weights = np.exp(results.logwt - results.logz[-1])

        # Compute 5%-95% quantiles.
        quantiles = dyfunc.quantile(samples, [0.05, 0.95])

        # Compute weighted mean and covariance.
        mean, cov = dyfunc.mean_and_cov(samples, weights)
        print()
        print('Weighted mean and convariance from original samplings')
        for key_name, key_value in theta_dictionary.items():
            print('  {0:s}  {1:15.6f} +- {2:15.6f}'.format(key_name, mean[key_value], cov[key_value, key_value]))
        print('From now on, all results are from weighted samples')


        # Resample weighted samples.
        flat_chain = dyfunc.resample_equal(samples, weights)
        flat_lnprob = dyfunc.resample_equal(results.logl, weights)

        n_samplings, n_pams = np.shape(flat_chain)

        lnprob_med = common.compute_value_sigma(flat_lnprob)
        chain_med = common.compute_value_sigma(flat_chain)

        chain_MAP, lnprob_MAP = common.pick_MAP_parameters(
            flat_chain, flat_lnprob)

        #data_in = np.genfromtxt(dir_input + 'post_equal_weights.dat')
        #flat_lnprob = data_in[:, -1]
        #flat_chain = data_in[:, :-1]
        # nsample = np.size(flat_lnprob)
        #n_samplings, n_pams = np.shape(flat_chain)

        #lnprob_med = common.compute_value_sigma(flat_lnprob)
        #chain_med = common.compute_value_sigma(flat_chain)
        # chain_MAP, lnprob_MAP = common.pick_MAP_parameters(
        #    flat_chain, flat_lnprob)

        print()
        print(' Reference Time Tref: {}'.format(mc.Tref))
        print()
        print(' Dimensions: {}'.format(mc.ndim))
        print()
        print(' Samples: {}'.format(n_samplings))


    if sampler_name == 'ultranest':

        import json


        dir_input = './' + config_in['output'] + '/ultranest/'
        dir_output = './' + config_in['output'] + '/ultranest_plot/'
        os.system('mkdir -p ' + dir_output)



        mc = nested_sampling_load_from_cpickle(dir_input)

        mc.model_setup()
        mc.initialize_logchi2()
        results_analysis.results_resumen(mc, None, skip_theta=True)

        with open(dir_input + 'info/results.json') as f:
            results = json.load(f)

        theta_dictionary = results_analysis.get_theta_dictionary(mc)

        res = ("niter: {:d}\n"
                "ncall: {:d}\n"
                "logz: {:6.3f} +/- {:6.3f}"
                .format(results['niter'], results['ncall'],
                        results['logz'], results['logzerr']))

        print()
        print('Summary - \n=======\n'+res)

        """ Copy the plots from default output directory of ultranest in
        ultranest_plots"""
        os.system(' cp '+ dir_input + 'plots/corner.pdf ' + dir_output + 'ultranest_corner.pdf')
        os.system(' cp '+ dir_input + 'plots/run.pdf ' + dir_output + 'ultranest_run.pdf')
        os.system(' cp '+ dir_input + 'plots/trace.pdf ' + dir_output + 'ultranest_trace.pdf')

        #labels_array = [None] * len(theta_dictionary)
        #for key_name, key_value in theta_dictionary.items():
        #    labels_array[key_value] = re.sub('_', '-', key_name)


        flat_chain = np.genfromtxt(dir_input + 'chains/equal_weighted_post.txt', skip_header=1)
        n_samplings, n_pams = np.shape(flat_chain)

        """ Filling the lnprob array the hard way """
        flat_lnprob = np.empty(n_samplings)
        for ii in range(0,n_samplings):
            flat_lnprob[ii] = mc.ultranest_call(flat_chain[ii,:])


        lnprob_med = common.compute_value_sigma(flat_lnprob)
        chain_med = common.compute_value_sigma(flat_chain)

        chain_MAP, lnprob_MAP = common.pick_MAP_parameters(
            flat_chain, flat_lnprob)

        un_lnprob_MAP = results['maximum_likelihood']['logl']
        un_chain_MAP = results['maximum_likelihood']['point']

        """ Weighted mean and standard deviation from original distribution """
        print()
        print('Weighted median, mean and convariance from original samplings, internal MAP and ultranest MAP')
        for key_name, key_value in theta_dictionary.items():
            print('  {0:15s} {1:15.6f}, {2:15.6f} +- {3:15.6f}, {4:15.6f} {5:15.6f} '.format(key_name,
                                                           results['posterior']['median'][key_value],
                                                           results['posterior']['mean'][key_value],
                                                           results['posterior']['stdev'][key_value],
                                                           chain_MAP[key_value],
                                                           un_chain_MAP[key_value]))

        print('From now on, all results are from weighted samples')


        #data_in = np.genfromtxt(dir_input + 'post_equal_weights.dat')
        #flat_lnprob = data_in[:, -1]
        #flat_chain = data_in[:, :-1]
        # nsample = np.size(flat_lnprob)
        #n_samplings, n_pams = np.shape(flat_chain)

        #lnprob_med = common.compute_value_sigma(flat_lnprob)
        #chain_med = common.compute_value_sigma(flat_chain)
        # chain_MAP, lnprob_MAP = common.pick_MAP_parameters(
        #    flat_chain, flat_lnprob)

        print()
        print(' Reference Time Tref: {}'.format(mc.Tref))
        print()
        print(' Dimensions: {}'.format(mc.ndim))
        print()
        print(' Samples: {}'.format(n_samplings))


    print()
    print(' LN posterior: {0:12f}   {1:12f} {2:12f} (15-84 p) '.format(
        lnprob_med[0], lnprob_med[2], lnprob_med[1]))

    med_log_priors, med_log_likelihood = mc.log_priors_likelihood(
        chain_med[:, 0])
    BIC = -2.0 * med_log_likelihood + np.log(mc.ndata) * mc.ndim
    AIC = -2.0 * med_log_likelihood + 2.0 * mc.ndim
    AICc = AIC + (2.0 + 2.0 * mc.ndim) * mc.ndim / (mc.ndata - mc.ndim - 1.0)
    # AICc for small sample

    print()
    print(' Median log_priors     = {}'.format(med_log_priors))
    print(' Median log_likelihood = {}'.format(med_log_likelihood))
    print(' Median BIC  (using likelihood) = {}'.format(BIC))
    print(' Median AIC  (using likelihood) = {}'.format(AIC))
    print(' Median AICc (using likelihood) = {}'.format(AICc))

    med_log_posterior = med_log_likelihood + med_log_priors
    BIC = -2.0 * med_log_posterior + np.log(mc.ndata) * mc.ndim
    AIC = -2.0 * med_log_posterior + 2.0 * mc.ndim
    AICc = AIC + (2.0 + 2.0 * mc.ndim) * mc.ndim / (mc.ndata - mc.ndim - 1.0)

    print()
    print(' Median BIC  (using posterior)  = {}'.format(BIC))
    print(' Median AIC  (using posterior)  = {}'.format(AIC))
    print(' Median AICc (using posterior)  = {}'.format(AICc))

    MAP_log_priors, MAP_log_likelihood = mc.log_priors_likelihood(chain_MAP)
    BIC = -2.0 * MAP_log_likelihood + np.log(mc.ndata) * mc.ndim
    AIC = -2.0 * MAP_log_likelihood + 2.0 * mc.ndim
    AICc = AIC + (2.0 + 2.0 * mc.ndim) * mc.ndim / (mc.ndata - mc.ndim - 1.0)
    # AICc for small sample

    print()
    print(' MAP log_priors     = {}'.format(MAP_log_priors))
    print(' MAP log_likelihood = {}'.format(MAP_log_likelihood))
    print(' MAP BIC  (using likelihood) = {}'.format(BIC))
    print(' MAP AIC  (using likelihood) = {}'.format(AIC))
    print(' MAP AICc (using likelihood) = {}'.format(AICc))

    MAP_log_posterior = MAP_log_likelihood + MAP_log_priors
    BIC = -2.0 * MAP_log_posterior + np.log(mc.ndata) * mc.ndim
    AIC = -2.0 * MAP_log_posterior + 2.0 * mc.ndim
    AICc = AIC + (2.0 + 2.0 * mc.ndim) * mc.ndim / (mc.ndata - mc.ndim - 1.0)

    print()
    print(' MAP BIC  (using posterior)  = {}'.format(BIC))
    print(' MAP AIC  (using posterior)  = {}'.format(AIC))
    print(' MAP AICc (using posterior)  = {}'.format(AICc))

    if mc.ndata < 40 * mc.ndim:
        print()
        print(
            ' AICc suggested over AIC because NDATA ( {0:12f} ) < 40 * NDIM ( {1:12f} )'.format(mc.ndata, mc.ndim))
    else:
        print()
        print(
            ' AIC suggested over AICs because NDATA ( {0:12f} ) > 40 * NDIM ( {1:12f} )'.format(mc.ndata, mc.ndim))

    print()
    print('****************************************************************************************************')
    print('****************************************************************************************************')
    print()
    print(' Confidence intervals (median value, 34.135th percentile from the median on the left and right side)')

    planet_variables = results_analysis.results_resumen(
        mc, flat_chain, chain_med=chain_MAP, return_samples=True)

    print()
    print('****************************************************************************************************')
    print()
    print(' Parameters corresponding to the Maximum a Posteriori probability ( {} )'.format(lnprob_MAP))
    print()

    results_analysis.results_resumen(mc, chain_MAP)

    print()
    print('****************************************************************************************************')
    print()

    sys.stdout.flush()

    # Computation of all the planetary variables
    planet_variables_med = results_analysis.get_planet_variables(
        mc, chain_med[:, 0])
    star_variables = results_analysis.get_stellar_parameters(
        mc, chain_med[:, 0])

    planet_variables_MAP = results_analysis.get_planet_variables(mc, chain_MAP)
    star_variables_MAP = results_analysis.get_stellar_parameters(mc, chain_MAP, warnings=False)

    if plot_dictionary['lnprob_chain'] or plot_dictionary['chains']:

        print(' Plot FLAT chain ')

        fig = plt.figure(figsize=(12, 12))
        plt.xlabel('$\ln \mathcal{L}$')
        plt.plot(sampler_lnprob, '-', alpha=0.5)
        plt.axhline(lnprob_med[0])
        plt.axvline(nburnin / nthin, c='r')
        plt.savefig(dir_output + 'LNprob_chain.png',
                    bbox_inches='tight', dpi=300)
        plt.close(fig)

        print()
        print('****************************************************************************************************')
        print()
        sys.stdout.flush()

    if plot_dictionary['full_correlation']:

        corner_plot = {
            'samples': np.zeros([np.size(flat_chain, axis=0), np.size(flat_chain, axis=1) + 1]),
            'labels': [],
            'truths': []
        }

        i_corner = 0
        for var, var_dict in theta_dictionary.items():
            corner_plot['samples'][:, i_corner] = flat_chain[:, var_dict]
            if len(theta_dictionary) > 10:
                corner_plot['labels'].append(repr(var_dict))
            else:
                corner_plot['labels'].append(re.sub('_', '-', var))
            corner_plot['truths'].append(chain_med[var_dict, 0])
            i_corner += 1

        corner_plot['samples'][:, -1] = flat_lnprob[:]
        corner_plot['labels'].append('ln-prob')
        corner_plot['truths'].append(lnprob_med[0])

        if plot_dictionary['use_getdist']:
            try:
                print(' Plotting full_correlation plot with GetDist')
                print()
                print(' Ignore the no burn in error warning from getdist')
                print(' since burn in has been already removed from the chains')

                samples = MCSamples(samples=corner_plot['samples'], names=corner_plot['labels'],
                                    labels=corner_plot['labels'])

                g = plots.getSubplotPlotter()
                g.settings.num_plot_contours = 6
                g.triangle_plot(samples, filled=True)
                g.export(dir_output + "all_internal_variables_corner_getdist.pdf")

            except AttributeError:
                print(' Something went wrong when plotting the coner plot with GetDist')
                print(' Please Run PyORBIT_GetResults.py with without corner plot flat to get an alternative corner plot')

            print()

        elif plot_dictionary['use_corner']:
            # plotting mega-corner plot
            print('Plotting full_correlation plot with Corner')
            fig = corner.corner(
                corner_plot['samples'], labels=corner_plot['labels'], truths=corner_plot['truths'])
            fig.savefig(dir_output + "all_internal_variables_corner_dfm.pdf",
                        bbox_inches='tight', dpi=300)
            plt.close(fig)

        else:
            print('Plotting full_correlation plot with pygtc')

            GTC = pygtc.plotGTC(chains=corner_plot['samples'],
                                paramNames=corner_plot['labels'],
                                truths=corner_plot['truths'],
                                plotName=dir_output + "all_internal_variables_corner_pygtc.pdf")
            GTC = None
        print()
        print('****************************************************************************************************')
        print()
        sys.stdout.flush()
        corner_plot = None

    if plot_dictionary['chains']:

        print(' Plotting the chains... ')

        os.system('mkdir -p ' + dir_output + 'chains')
        for theta_name, ii in theta_dictionary.items():
            file_name = dir_output + 'chains/' + \
                repr(ii) + '_' + theta_name + '.png'
            fig = plt.figure(figsize=(12, 12))
            plt.plot(sampler_chain[:, :, ii].T, '-', alpha=0.5)
            plt.axvline(nburnin / nthin, c='r')
            plt.savefig(file_name, bbox_inches='tight', dpi=300)
            plt.close(fig)

        print()
        print('****************************************************************************************************')
        print()
        sys.stdout.flush()

    if plot_dictionary['traces']:

        print(' Plotting the Gelman-Rubin traces... ')
        print()
        """
        Gelman-Rubin traces are stored in the dedicated folder iniside the _plot folder
        Note that the GR statistics is not robust because the wlakers are not independent
        """
        os.system('mkdir -p ' + dir_output + 'gr_traces')

        step_sampling = np.arange(nburn, nsteps / nthin, 1, dtype=int)

        for theta_name, th in theta_dictionary.items():
            rhat = np.array([GelmanRubin_v2(sampler_chain[:, :steps, th])
                             for steps in step_sampling])
            print(
                '     Gelman-Rubin: {0:5d} {1:12f} {2:s} '.format(th, rhat[-1], theta_name))
            file_name = dir_output + 'gr_traces/v2_' + \
                repr(th) + '_' + theta_name + '.png'
            fig = plt.figure(figsize=(12, 12))
            plt.plot(step_sampling, rhat[:], '-', color='k')
            plt.axhline(1.01, c='C0')
            plt.savefig(file_name, bbox_inches='tight', dpi=300)
            plt.close(fig)

        print()
        print('****************************************************************************************************')
        print()
        sys.stdout.flush()

    if plot_dictionary['common_corner']:

        print(' Plotting the common models corner plots')

        for common_name, common_model in mc.common_models.items():

            print('     Common model: ', common_name)

            corner_plot = {
                'var_list': [],
                'samples': [],
                'labels': [],
                'truths': []
            }
            variable_values = common_model.convert(flat_chain)
            variable_median = common_model.convert(chain_med[:, 0])

            if len(variable_median) < 1.:
                continue

            """
            Check if the eccentricity and argument of pericenter were set as free parameters or fixed by simply
            checking the size of their distribution
            """
            for var in variable_values.keys():
                if np.size(variable_values[var]) == 1:
                    variable_values[var] = variable_values[var] * \
                        np.ones(n_samplings)
                else:
                    corner_plot['var_list'].append(var)

            corner_plot['samples'] = []
            corner_plot['labels'] = []
            corner_plot['truths'] = []
            for var_i, var in enumerate(corner_plot['var_list']):
                corner_plot['samples'].extend([variable_values[var]])
                corner_plot['labels'].append(var)
                corner_plot['truths'].append(variable_median[var])

            """ Check if the semi-amplitude K is among the parameters that
                have been fitted. If so, it computes the corresponding
                planetary mass with uncertainty """

            try:
                if plot_dictionary['use_corner']:
                    fig = corner.corner(np.asarray(corner_plot['samples']).T,
                                        labels=corner_plot['labels'],
                                        truths=corner_plot['truths'])
                    fig.savefig(dir_output + common_name + "_corners.pdf",
                                bbox_inches='tight', dpi=300)
                    plt.close(fig)
                else:
                    GTC = pygtc.plotGTC(chains=np.asarray(corner_plot['samples']).T,
                                        paramNames=corner_plot['labels'],
                                        truths=corner_plot['truths'],
                                        #figureSize='MNRAS_page',
                                        plotName=dir_output + common_name + "_corners.pdf")
                    GTC = None

            except AssertionError:
                print('     Something went wrong, plot skipped ')
                print()

            corner_plot = None
            sys.stdout.flush()

        print()
        print('****************************************************************************************************')
        print()
        sys.stdout.flush()

    if plot_dictionary['dataset_corner']:

        print(' Dataset + models corner plots ')
        print()

        for dataset_name, dataset in mc.dataset_dict.items():

            for model_name in dataset.models:

                corner_plot = {
                    'samples': [],
                    'labels': [],
                    'truths': []
                }


                print('     Dataset: ', dataset_name, '    model: ',
                      model_name, ' corner plot  starting ')

                variable_values = dataset.convert(flat_chain)
                variable_median = dataset.convert(chain_med[:, 0])

                for common_ref in mc.models[model_name].common_ref:
                    variable_values.update(
                        mc.common_models[common_ref].convert(flat_chain))
                    variable_median.update(
                        mc.common_models[common_ref].convert(chain_med[:, 0]))

                variable_values.update(
                    mc.models[model_name].convert(flat_chain, dataset_name))
                variable_median.update(mc.models[model_name].convert(
                    chain_med[:, 0], dataset_name))

                for var_i, var in enumerate(variable_values):
                    if np.size(variable_values[var]) <= 1:
                        continue
                    corner_plot['samples'].extend([variable_values[var]])
                    corner_plot['labels'].append(var)
                    corner_plot['truths'].append(variable_median[var])


                try:
                    if plot_dictionary['use_corner']:
                        fig = corner.corner(np.asarray(corner_plot['samples']).T,
                                            labels=corner_plot['labels'], truths=corner_plot['truths'])
                        fig.savefig(dir_output + dataset_name + '_' + model_name +
                                    "_corners.pdf", bbox_inches='tight', dpi=300)
                        plt.close(fig)
                        fig = None
                    else:
                        GTC = pygtc.plotGTC(chains=np.asarray(corner_plot['samples']).T,
                                            paramNames=corner_plot['labels'],
                                            truths=corner_plot['truths'],
                                            #figureSize='MNRAS_page',
                                            plotName=dir_output + dataset_name + '_' + model_name + "_corners.pdf")
                        GTC = None
                except AssertionError:
                    print('     Something went wrong, plot skipped ')
                    print()

                print('     Dataset: ', dataset_name, '    model: ',
                      model_name, ' corner plot  done ')
                corner_plot = None
                sys.stdout.flush()

        print()
        print('****************************************************************************************************')
        print()
        sys.stdout.flush()

    if plot_dictionary['write_planet_samples']:

        print(' Saving the planet variable samplings to files (with plots)')

        samples_dir = dir_output + '/planet_samples/'
        os.system('mkdir -p ' + samples_dir)
        sys.stdout.flush()

        for common_ref, variable_values in planet_variables.items():
            for variable_name, variable in variable_values.items():

                rad_filename = samples_dir + common_ref + '_' + variable_name
                fileout = open(rad_filename + '.dat', 'w')
                for val in variable:
                    fileout.write('{0:.12f} \n'.format(val))
                fileout.close()

                try:
                    fig = plt.figure(figsize=(10, 10))
                    plt.hist(variable, bins=50, color='C0',
                             alpha=0.75, zorder=0)

                    perc0, perc1, perc2 = np.percentile(
                        variable, [15.865, 50, 84.135], axis=0)

                    plt.axvline(planet_variables_med[common_ref][variable_name], color='C1', zorder=1,
                                label='Median-corresponding value')
                    plt.axvline(planet_variables_MAP[common_ref][variable_name], color='C2', zorder=1,
                                label='MAP-corresponding value')
                    plt.axvline(perc1, color='C3', zorder=2,
                                label='Median of the distribution')
                    plt.axvline(perc0, color='C4', zorder=2,
                                label='15.865th and 84.135th percentile')
                    plt.axvline(perc2, color='C4', zorder=2)
                    plt.xlabel(
                        re.sub('_', '-', variable_name + '_' + common_ref))
                    plt.legend()
                    plt.ticklabel_format(useOffset=False)
                    plt.savefig(rad_filename + '.png',
                                bbox_inches='tight', dpi=300)
                    plt.close(fig)
                except:
                    print(
                        '   Error while producing the histogram plot for variable ', variable_name)

        print()
        print('****************************************************************************************************')
        print()
        sys.stdout.flush()

    if plot_dictionary['write_all_samples']:

        print(' Saving all the variable samplings to files (with plots)')

        samples_dir = dir_output + '/all_samples/'
        os.system('mkdir -p ' + samples_dir)

        rad_filename = samples_dir + 'log_likelihood'
        fileout = open(rad_filename + '.dat', 'w')
        for val in flat_lnprob:
            fileout.write('{0:.12f} \n'.format(val))
        fileout.close()

        for theta_name, th in theta_dictionary.items():

            rad_filename = samples_dir + repr(th) + '_' + theta_name
            fileout = open(rad_filename + '.dat', 'w')
            for val in flat_chain[:, th]:
                fileout.write('{0:.12f} \n'.format(val))
            fileout.close()

            try:
                fig = plt.figure(figsize=(10, 10))
                plt.hist(flat_chain[:, th], bins=50,
                         color='C0', alpha=0.75, zorder=0)

                perc0, perc1, perc2 = np.percentile(
                    flat_chain[:, th], [15.865, 50, 84.135], axis=0)

                plt.axvline(chain_med[th, 0], color='C1', zorder=1,
                            label='Median-corresponding value')
                plt.axvline(chain_MAP[th], color='C2', zorder=1,
                            label='MAP-corresponding value')
                plt.axvline(perc1, color='C3', zorder=2,
                            label='Median of the distribution')
                plt.axvline(perc0, color='C4', zorder=2,
                            label='15.865th and 84.135th percentile')
                plt.axvline(perc2, color='C4', zorder=2)
                plt.xlabel(re.sub('_', '-', 'theta: ' +
                                  repr(th) + ' variable: ' + theta_name))
                plt.legend()
                plt.ticklabel_format(useOffset=False)
                plt.savefig(rad_filename + '.png',
                            bbox_inches='tight', dpi=300)
                plt.close(fig)
            except:
                print('   Error while producing the histogram plot for sample variable {0:s} (id: {1:5.0f})'.format(
                    theta_name, th))

        print()
        print('****************************************************************************************************')
        print()
        sys.stdout.flush()

    if plot_dictionary['plot_models'] or plot_dictionary['write_models']:

        print(' Computing the models for plot/data writing ')

        bjd_plot = {
            'full': {
                'start': None, 'end': None, 'range': None
            }
        }

        kinds = {}

        P_minimum = 2.0  # this temporal range will be divided in 20 subsets
        for key_name, key_val in planet_variables_med.items():
            P_minimum = min(key_val.get('P', 2.0), P_minimum)

        for dataset_name, dataset in mc.dataset_dict.items():
            if dataset.kind in kinds.keys():
                kinds[dataset.kind].extend([dataset_name])
            else:
                kinds[dataset.kind] = [dataset_name]

            bjd_plot[dataset_name] = {
                'start': np.amin(dataset.x),
                'end': np.amax(dataset.x),
                'range': np.amax(dataset.x) - np.amin(dataset.x),
            }

            if bjd_plot[dataset_name]['range'] < 0.1:
                bjd_plot[dataset_name]['range'] = 0.1

            bjd_plot[dataset_name]['start'] -= bjd_plot[dataset_name]['range'] * 0.10
            bjd_plot[dataset_name]['end'] += bjd_plot[dataset_name]['range'] * 0.10

            if dataset.kind == 'Phot':

                if bjd_plot[dataset_name]['range'] > P_minimum:
                    # more than one transit:
                    step_size = 5.  / (24 * 60) #five minute stepsize
                else:
                    step_size = np.min(
                        bjd_plot[dataset_name]['range'] / dataset.n / 10.)
            else:
                step_size = P_minimum / 20.

            bjd_plot[dataset_name]['x_plot'] = \
                np.arange(bjd_plot[dataset_name]['start'],
                          bjd_plot[dataset_name]['end'], step_size)
            bjd_plot[dataset_name]['x0_plot'] = bjd_plot[dataset_name]['x_plot'] - mc.Tref

            if bjd_plot['full']['range']:
                bjd_plot['full']['start'] = min(
                    bjd_plot['full']['start'], np.amin(dataset.x))
                bjd_plot['full']['end'] = max(
                    bjd_plot['full']['end'], np.amax(dataset.x))
                bjd_plot['full']['range'] = bjd_plot['full']['end'] - \
                    bjd_plot['full']['start']
            else:
                bjd_plot['full']['start'] = np.amin(dataset.x)
                bjd_plot['full']['end'] = np.amax(dataset.x)
                bjd_plot['full']['range'] = bjd_plot['full']['end'] - \
                    bjd_plot['full']['start']

        bjd_plot['full']['start'] -= bjd_plot['full']['range'] * 0.50
        bjd_plot['full']['end'] += bjd_plot['full']['range'] * 0.50
        bjd_plot['full']['x_plot'] = np.arange(
            bjd_plot['full']['start'], bjd_plot['full']['end'], P_minimum / 20.)
        bjd_plot['full']['x0_plot'] = bjd_plot['full']['x_plot'] - mc.Tref

        # Special cases
        for dataset_name, dataset in mc.dataset_dict.items():
            if dataset.kind == 'RV':
                bjd_plot[dataset_name] = bjd_plot['full']
            if dataset.kind == 'Tcent':
                bjd_plot[dataset_name]['x_plot'] = dataset.x
                bjd_plot[dataset_name]['x0_plot'] = dataset.x

        bjd_plot['model_out'], bjd_plot['model_x'] = results_analysis.get_model(
            mc, chain_med[:, 0], bjd_plot)
        bjd_plot['MAP_model_out'], bjd_plot['MAP_model_x'] = results_analysis.get_model(
            mc, chain_MAP, bjd_plot)

        if plot_dictionary['plot_models']:
            print(' Writing the plots ')

            for kind_name, kind in kinds.items():
                for dataset_name in kind:

                    try:
                        error_bars = np.sqrt(mc.dataset_dict[dataset_name].e**2
                                             + bjd_plot['model_out'][dataset_name]['jitter']**2)
                    except ValueError:
                        error_bars = mc.dataset_dict[dataset_name].e

                    fig = plt.figure(figsize=(12, 12))

                    # Partially taken from here:
                    # http://www.sc.eso.org/~bdias/pycoffee/codes/20160407/gridspec_demo.html
                    gs = gridspec.GridSpec(2, 1, height_ratios=[3.0, 1.0])
                    # Also make sure the margins and spacing are apropriate
                    gs.update(left=0.3, right=0.95, bottom=0.08,
                              top=0.93, wspace=0.15, hspace=0.05)

                    ax_0 = plt.subplot(gs[0])
                    ax_1 = plt.subplot(gs[1], sharex=ax_0)

                    # Adding minor ticks only to x axis
                    minorLocator = AutoMinorLocator()
                    ax_0.xaxis.set_minor_locator(minorLocator)
                    ax_1.xaxis.set_minor_locator(minorLocator)

                    # Disabling the offset on top of the plot
                    ax_0.ticklabel_format(useOffset=False)
                    ax_1.ticklabel_format(useOffset=False)

                    ax_0.scatter(mc.dataset_dict[dataset_name].x,
                                 mc.dataset_dict[dataset_name].y
                                 - bjd_plot['model_out'][dataset_name]['systematics']
                                 - bjd_plot['model_out'][dataset_name]['time_independent'],
                                 color='C0', zorder=4, s=16)
                    ax_0.errorbar(mc.dataset_dict[dataset_name].x,
                                  mc.dataset_dict[dataset_name].y
                                  - bjd_plot['model_out'][dataset_name]['systematics']
                                  - bjd_plot['model_out'][dataset_name]['time_independent'],
                                  yerr=error_bars,
                                  color='C0', fmt='o', ms=0, zorder=3, alpha=0.5)
                    ax_0.plot(bjd_plot[dataset_name]['x_plot'], bjd_plot['model_x'][dataset_name]['complete'],
                              label='Median-corresponding model',
                              color='C1', zorder=2)
                    ax_0.plot(bjd_plot[dataset_name]['x_plot'], bjd_plot['MAP_model_x'][dataset_name]['complete'],
                              label='MAP-corresponding model',
                              color='C2', zorder=1)

                    ax_0.set_ylabel('Same as input data')
                    ax_0.legend()

                    ax_1.scatter(mc.dataset_dict[dataset_name].x,
                                 mc.dataset_dict[dataset_name].y -
                                 bjd_plot['model_out'][dataset_name]['complete'],
                                 color='C0', zorder=4, s=16)
                    ax_1.errorbar(mc.dataset_dict[dataset_name].x,
                                  mc.dataset_dict[dataset_name].y -
                                  bjd_plot['model_out'][dataset_name]['complete'],
                                  yerr=error_bars,
                                  color='C0', fmt='o', ms=0, zorder=3, alpha=0.5)
                    ax_1.axhline(0.0, color='k', alpha=0.5, zorder=0)

                    ax_1.set_xlabel('Time [d] (offset as the input data)')
                    ax_1.set_ylabel('Residuals (wrt median model)')

                    plt.savefig(dir_output + 'model_' + kind_name + '_' + dataset_name + '.png', bbox_inches='tight',
                                dpi=300)
                    plt.close(fig)

        if plot_dictionary['write_models']:

            for prepend_keyword in ['', 'MAP_']:

                print(' Writing the ', prepend_keyword, 'data files ')

                plot_out_keyword = prepend_keyword + 'model_out'
                plot_x_keyword = prepend_keyword + 'model_x'
                file_keyword = prepend_keyword + 'model_files'

                if prepend_keyword == '':
                    planet_vars = planet_variables_med
                    # star_vars = star_variables # leaving here, it could be useful for the future
                    chain_ref = chain_med[:, 0]
                elif prepend_keyword == 'MAP_':
                    planet_vars = planet_variables_MAP
                    # star_vars = star_variables_MAP
                    chain_ref = chain_MAP

                dir_models = dir_output + file_keyword + '/'
                os.system('mkdir -p ' + dir_models)

                for dataset_name, dataset in mc.dataset_dict.items():
                    for model_name in dataset.models:

                        if getattr(mc.models[model_name], 'systematic_model', False):
                            continue

                        fileout = open(dir_models + dataset_name +
                                       '_' + model_name + '.dat', 'w')

                        phase = np.zeros(dataset.n)
                        tc_folded = np.zeros(dataset.n)
                        phase_plot = np.zeros(
                            np.size(bjd_plot[dataset_name]['x_plot']))
                        tc_folded_plot = np.zeros(
                            np.size(bjd_plot[dataset_name]['x_plot']))

                        for common_ref in mc.models[model_name].common_ref:
                            if common_ref in planet_vars:
                                if 'P' in planet_vars[common_ref]:
                                    phase = (dataset.x0 /
                                             planet_vars[common_ref]['P']) % 1
                                    phase_plot = ((bjd_plot[dataset_name]['x_plot'] - mc.Tref) /
                                                  planet_vars[common_ref]['P']) % 1
                                    if 'Tc' in planet_vars[common_ref]:
                                        tc_folded = (dataset.x - planet_vars[common_ref]['Tc']
                                                     + planet_vars[common_ref]['P'] / 2.) \
                                            % planet_vars[common_ref]['P'] \
                                            - planet_vars[common_ref]['P'] / 2.
                                        tc_folded_plot = (bjd_plot[dataset_name]['x_plot'] - planet_vars[common_ref][
                                            'Tc']
                                            + planet_vars[common_ref]['P'] / 2.) \
                                            % planet_vars[common_ref]['P'] \
                                            - planet_vars[common_ref]['P'] / 2.
                                    else:
                                        tc_folded = dataset.x0 % planet_vars[common_ref]['P']
                                        tc_folded_plot = (bjd_plot[dataset_name]['x_plot'] - mc.Tref) % \
                                            planet_vars[common_ref]['P']

                        fileout.write(
                            'descriptor BJD Tc_folded pha val,+- sys mod full val_compare,+- res,+- jit \n')

                        try:
                            len(bjd_plot[plot_out_keyword]
                                [dataset_name][model_name])
                        except:
                            bjd_plot[plot_out_keyword][dataset_name][model_name] = \
                                bjd_plot[plot_out_keyword][dataset_name][model_name] * \
                                np.ones(dataset.n)

                            bjd_plot[plot_x_keyword][dataset_name][model_name] = \
                                bjd_plot[plot_x_keyword][dataset_name][model_name] * \
                                np.ones(dataset.n)

                        for x, tcf, pha, y, e, sys_val, mod, com, obs_mod, res, jit in zip(
                                dataset.x, tc_folded, phase, dataset.y, dataset.e,
                                bjd_plot[plot_out_keyword][dataset_name]['systematics'],
                                bjd_plot[plot_out_keyword][dataset_name][model_name],
                                bjd_plot[plot_out_keyword][dataset_name]['complete'],
                                dataset.y - bjd_plot[plot_out_keyword][dataset_name]['complete'] +
                                bjd_plot[plot_out_keyword][dataset_name][model_name],
                                dataset.y -
                            bjd_plot[plot_out_keyword][dataset_name]['complete'],
                                bjd_plot[plot_out_keyword][dataset_name]['jitter']):
                            fileout.write('{0:f} {1:f} {2:f} {3:f} {4:f} {5:f} {6:1f} {7:f} {8:f} {9:f} {10:f} {11:f} {12:f}'
                                          '\n'.format(x, tcf, pha, y, e, sys_val, mod, com, obs_mod, e, res, e, jit))
                        fileout.close()

                        if getattr(mc.models[model_name], 'systematic_model', False):
                            continue

                        if getattr(mc.models[model_name], 'jitter_model', False):
                            continue

                        fileout = open(dir_models + dataset_name +
                                       '_' + model_name + '_full.dat', 'w')

                        if model_name + '_std' in bjd_plot[plot_x_keyword][dataset_name]:
                            fileout.write(
                                'descriptor BJD Tc_folded phase mod,+- \n')
                            for x, tfc, pha, mod, std in zip(
                                    bjd_plot[dataset_name]['x_plot'],
                                    tc_folded_plot,
                                    phase_plot,
                                    bjd_plot[plot_x_keyword][dataset_name][model_name],
                                    bjd_plot[plot_x_keyword][dataset_name][model_name + '_std']):
                                fileout.write('{0:f} {1:f} {2:f} {3:12f} {4:12f} \n'.format(
                                    x, tcf, pha, mod, std))
                            fileout.close()
                        else:
                            fileout.write(
                                'descriptor BJD Tc_folded phase mod \n')
                            for x, tcf, pha, mod in zip(bjd_plot[dataset_name]['x_plot'],
                                                        tc_folded_plot,
                                                        phase_plot,
                                                        bjd_plot[plot_x_keyword][dataset_name][model_name]):
                                fileout.write(
                                    '{0:f} {1:f} {2:f} {3:12f}\n'.format(x, tcf, pha, mod))
                            fileout.close()

                        if getattr(mc.models[model_name], 'model_class', False) in oversampled_models:
                            """
                            Additional output to deal with under-sampled lightcurves, i.e. when folding
                            the light curve from the model file is not good enough. Something similar is performed later
                            with the planetary RVs, but here we must keep into account the differences  between datasets
                            due to limb darkening, exposure times, etc.
                            """

                            variable_values = {}
                            for common_ref in mc.models[model_name].common_ref:
                                variable_values.update(
                                    mc.common_models[common_ref].convert(chain_ref))
                            variable_values.update(
                                mc.models[model_name].convert(chain_ref, dataset_name))

                            fileout = open(
                                dir_models + dataset_name + '_' + model_name + '_oversampled.dat', 'w')

                            x_range = np.arange(
                                -variable_values['P']/2., variable_values['P']/2., 0.001)
                            try:
                                delta_T = variable_values['Tc']-dataset.Tref
                            except KeyError:
                                delta_T = kepler_exo.kepler_phase2Tc_Tref(variable_values['P'],
                                                                          variable_values['mean_long'],
                                                                          variable_values['e'],
                                                                          variable_values['omega'])

                            y_plot = mc.models[model_name].compute(
                                variable_values, dataset, x_range+delta_T)

                            fileout.write('descriptor Tc_folded  mod \n')
                            for x, mod in zip(x_range, y_plot):
                                fileout.write('{0:f} {1:f} \n'.format(x, mod))
                            fileout.close()

                    fileout = open(dir_models + dataset_name +
                                   '_full.dat', 'w')
                    fileout.write('descriptor BJD mod \n')
                    for x, mod in zip(bjd_plot[dataset_name]['x_plot'],
                                      bjd_plot[plot_x_keyword][dataset_name]['complete']):
                        fileout.write('{0:f} {1:f} \n'.format(x, mod))
                    fileout.close()

                for model in planet_vars:
                    try:

                        RV_out = kepler_exo.kepler_RV_T0P(bjd_plot['full']['x_plot']-mc.Tref,
                                                          planet_vars[model]['mean_long'],
                                                          planet_vars[model]['P'],
                                                          planet_vars[model]['K'],
                                                          planet_vars[model]['e'],
                                                          planet_vars[model]['omega'])
                        fileout = open(
                            dir_models + 'RV_planet_' + model + '_kep.dat', 'w')
                        fileout.write('descriptor x_range  m_kepler \n')
                        for x, y in zip(bjd_plot['full']['x_plot'], RV_out):
                            fileout.write('{0:f} {1:f} \n'.format(x, y))
                        fileout.close()

                        x_range = np.arange(-1.50, 1.50, 0.001)
                        RV_out = kepler_exo.kepler_RV_T0P(x_range * planet_vars[model]['P'],
                                                          planet_vars[model]['mean_long'],
                                                          planet_vars[model]['P'],
                                                          planet_vars[model]['K'],
                                                          planet_vars[model]['e'],
                                                          planet_vars[model]['omega'])
                        fileout = open(
                            dir_models + 'RV_planet_' + model + '_pha.dat', 'w')
                        fileout.write('descriptor x_phase m_phase \n')
                        for x, y in zip(x_range, RV_out):
                            fileout.write('{0:f} {1:f} \n'.format(x, y))
                        fileout.close()

                        x_range = np.arange(-1.50, 1.50, 0.001)
                        if 'Tc' in planet_vars[model]:
                            Tc_range = x_range * \
                                planet_vars[model]['P'] + \
                                planet_vars[model]['Tc'] - mc.Tref
                            RV_out = kepler_exo.kepler_RV_T0P(Tc_range,
                                                              planet_vars[model]['mean_long'],
                                                              planet_vars[model]['P'],
                                                              planet_vars[model]['K'],
                                                              planet_vars[model]['e'],
                                                              planet_vars[model]['omega'])
                            fileout = open(
                                dir_models + 'RV_planet_' + model + '_Tcf.dat', 'w')
                            fileout.write('descriptor Tc_phase m_phase \n')
                            for x, y in zip(x_range, RV_out):
                                fileout.write('{0:f} {1:f} \n'.format(x, y))
                            fileout.close()

                    except:
                        pass

        print()
        print('****************************************************************************************************')
        print()
        sys.stdout.flush()

    if plot_dictionary['veuz_corner_files']:

        import csv

        print(' Writing Veusz-compatible files for personalized corner plots')

        # Transit times are too lenghty for the 'tiny' corner plot, so we apply a reduction to their value
        variable_with_offset = {}

        veusz_dir = dir_output + '/Veuz_plot/'
        if not os.path.exists(veusz_dir):
            os.makedirs(veusz_dir)

        all_variables_list = {}
        for dataset_name, dataset in mc.dataset_dict.items():
            variable_values = dataset.convert(flat_chain)

            for variable_name, variable in variable_values.items():
                all_variables_list[dataset_name +
                                   '_' + variable_name] = variable

            for model_name in dataset.models:
                variable_values = mc.models[model_name].convert(
                    flat_chain, dataset_name)
                for variable_name, variable in variable_values.items():
                    all_variables_list[dataset_name + '_' +
                                       model_name + '_' + variable_name] = variable

        for model_name, model in mc.common_models.items():
            variable_values = model.convert(flat_chain)

            for variable_name, variable in variable_values.items():

                all_variables_list[model.common_ref +
                                   '_' + variable_name] = variable
                #print(variable_name, variable)
                # Special treatment for transit time, since ti can be very long but yet very precise, making
                # the axis of corner plot quite messy
                if variable_name == 'Tc':
                    offset = np.median(variable)
                    variable_with_offset[model.common_ref +
                                         '_' + variable_name] = offset
                    all_variables_list[model.common_ref +
                                       '_' + variable_name] -= offset

                # Let's save omega in degrees, in the range 0-360
                if variable_name == 'o':
                    odeg = variable * 180 / np.pi
                    try:
                        sel = (odeg < 0.000)
                        odeg[sel] += 360.00
                    except TypeError:
                        if odeg < 0.000:
                            odeg += 360.00
                    all_variables_list[model.common_ref +
                                       '_' + variable_name + 'deg'] = odeg

        for common_ref, variable_values in planet_variables.items():
            for variable_name, variable in variable_values.items():

                # Skipping the variables that have been already included in all_variables_list
                if common_ref + '_' + variable_name in all_variables_list:
                    continue

                all_variables_list[common_ref + '_' + variable_name] = variable

                if variable_name == 'Tc':
                    offset = np.median(variable)
                    variable_with_offset[common_ref +
                                         '_' + variable_name] = offset
                    all_variables_list[common_ref +
                                       '_' + variable_name] -= offset

                if variable_name == 'o':
                    odeg = variable * 180 / np.pi
                    sel = (odeg < 0.000)
                    odeg[sel] += 360.00
                    all_variables_list[common_ref + '_' +
                                       variable_name + 'deg'] = odeg

        text_file = open(veusz_dir + "veusz_offsets.txt", "w")
        for variable_name, offset_value in variable_with_offset.items():
            text_file.write('{0:s} {1:16.9f}'.format(
                variable_name, offset_value))
        text_file.close()

        n_int = len(all_variables_list)
        output_plan = np.zeros([n_samplings, n_int], dtype=np.double)
        output_names = []
        for var_index, variable_name in enumerate(all_variables_list):
            output_plan[:, var_index] = all_variables_list[variable_name]
            output_names.extend([variable_name])

        plot_truths = np.percentile(
            output_plan[:, :], [15.865, 50, 84.135], axis=0)
        n_bins = 30 + 1

        h5f = h5py.File(veusz_dir + '_hist1d.hdf5', "w")
        data_grp = h5f.create_group("hist1d")

        data_lim = np.zeros([n_int, 2], dtype=np.double)
        data_edg = np.zeros([n_int, n_bins], dtype=np.double)
        data_skip = np.zeros(n_int, dtype=bool)

        sigma_minus = plot_truths[1, :] - plot_truths[0, :]
        sigma_plus = plot_truths[2, :] - plot_truths[1, :]
        median_vals = plot_truths[1, :]

        for ii in range(0, n_int):

            if sigma_minus[ii] == 0. and sigma_plus[ii] == 0.:
                data_skip[ii] = True
                continue

            sigma5_selection = (output_plan[:, ii] > median_vals[ii] - 5 * sigma_minus[ii]) & \
                               (output_plan[:, ii] <
                                median_vals[ii] + 5 * sigma_plus[ii])

            try:
                data_lim[ii, :] = [np.amin(output_plan[sigma5_selection, ii]), np.amax(
                    output_plan[sigma5_selection, ii])]
            except:
                continue

            if data_lim[ii, 0] == data_lim[ii, 1]:
                data_lim[ii, :] = [
                    np.amin(output_plan[:, ii]), np.amax(output_plan[:, ii])]
            if data_lim[ii, 0] == data_lim[ii, 1]:
                data_skip[ii] = True
                continue

            data_edg[ii, :] = np.linspace(
                data_lim[ii, 0], data_lim[ii, 1], n_bins)

        veusz_workaround_descriptor = 'descriptor'
        veusz_workaround_values = ''

        for ii in range(0, n_int):

            if data_skip[ii]:
                continue

            x_data = output_plan[:, ii]
            x_edges = data_edg[ii, :]

            for jj in range(0, n_int):

                if data_skip[jj]:
                    continue

                y_data = output_plan[:, jj]
                y_edges = data_edg[jj, :]

                if ii != jj:
                    hist2d = np.histogram2d(
                        x_data, y_data, bins=[x_edges, y_edges])
                    #hist1d_y = np.histogram(y_data, bins=y_edges)

                    Hflat = hist2d[0].flatten()
                    inds = np.argsort(Hflat)[::-1]
                    Hflat = Hflat[inds]
                    sm = np.cumsum(Hflat)
                    sm /= sm[-1]

                    x_edges_1d = (x_edges[1:] + x_edges[:-1]) / 2
                    y_edges_1d = (y_edges[1:] + y_edges[:-1]) / 2
                    h2d_out = np.zeros([n_bins, n_bins])
                    h2d_out[0, 1:] = x_edges_1d
                    h2d_out[1:, 0] = y_edges_1d
                    h2d_out[1:, 1:] = hist2d[0].T * 1. / np.amax(hist2d[0])

                    h2d_list = h2d_out.tolist()
                    h2d_list[0][0] = ''
                    csvfile = veusz_dir + '_hist2d___' + \
                        output_names[ii] + '___' + output_names[jj] + '.csv'
                    with open(csvfile, "w") as output:
                        writer = csv.writer(output, lineterminator='\n')
                        writer.writerows(h2d_list)

            hist1d = np.histogram(x_data, bins=x_edges)
            hist1d_norm = hist1d[0] * 1. / n_samplings
            x_edges_1d = (x_edges[1:] + x_edges[:-1]) / 2
            data_grp.create_dataset(
                output_names[ii] + '_x', data=x_edges_1d, compression="gzip")
            data_grp.create_dataset(
                output_names[ii] + '_y', data=hist1d_norm, compression="gzip")

            # data_grp.create_dataset(output_names[ii]+'_val', data=median_vals[ii])
            # data_grp.create_dataset(output_names[ii]+'_val_-', data=sigma_minus[ii])
            # data_grp.create_dataset(output_names[ii]+'_val_+', data=sigma_plus[ii])
            # data_grp.attrs[output_names[ii]+'_val'] = median_vals[ii]

            veusz_workaround_descriptor += ' ' + output_names[ii] + ',+,-'
            veusz_workaround_values += ' ' + repr(median_vals[ii]) + ' ' + repr(sigma_plus[ii]) + ' ' + repr(
                sigma_minus[ii])

        text_file = open(veusz_dir + "veusz_median_sigmas.txt", "w")
        text_file.write('%s \n' % veusz_workaround_descriptor)
        text_file.write('%s \n' % veusz_workaround_values)
        text_file.close()

        print()
        print('****************************************************************************************************')
        print()

    if plot_dictionary['P_versus_lnprob']:

        fig = plt.figure(figsize=(10, 10))

        ii = 0
        for common_ref, variable_values in planet_variables.items():
            if 'P' in variable_values:

                plt.scatter(variable_values['P'],
                            flat_lnprob, s=2, c='C'+repr(ii))
                ii += 1

        rad_filename = dir_output + 'P_versus_lnprob'

        plt.savefig(rad_filename + '.png', bbox_inches='tight', dpi=300)
        plt.close(fig)
