from __future__ import print_function

from pyorbit.classes.model_container_multinest import ModelContainerMultiNest
from pyorbit.classes.model_container_polychord import ModelContainerPolyChord
from pyorbit.classes.model_container_emcee import ModelContainerEmcee

from pyorbit.classes.input_parser import pars_input
from pyorbit.classes.io_subroutines import *
import numpy as np
import os
import matplotlib as mpl
import sys

mpl.use('Agg')

import re
import csv
import h5py
import pyorbit.classes.results_analysis as results_analysis
import pyorbit.classes.common as common
import pyorbit.classes.kepler_exo as kepler_exo
import corner
from matplotlib.ticker import AutoMinorLocator
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from pyorbit.classes.common import *

__all__ = ["pyorbit_getresults"]


def pyorbit_getresults(config_in, sampler, plot_dictionary):
    try:
        use_tex = config_in['parameters']['use_tex']
    except:
        use_tex = True

    if use_tex is False:
        print(' LaTeX disabled')

    if plot_dictionary['use_getdist']:
        from getdist import plots, MCSamples

    # plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rc('text', usetex=use_tex)

    sample_keyword = {
        'multinest': ['multinest', 'MultiNest', 'multi'],
        'polychord': ['polychord', 'PolyChord', 'polychrod', 'poly'],
        'emcee': ['emcee', 'MCMC', 'Emcee']
    }

    if sampler in sample_keyword['emcee']:

        dir_input = './' + config_in['output'] + '/emcee/'
        dir_output = './' + config_in['output'] + '/emcee_plot/'
        os.system('mkdir -p ' + dir_output)

        mc, starting_point, population, prob, state, \
            sampler_chain, sampler_lnprobability, sampler_acceptance_fraction, _, _ = \
            emcee_load_from_cpickle(dir_input)

        pars_input(config_in, mc, reload_emcee=True)

        if hasattr(mc.emcee_parameters, 'version'):
            emcee_version = mc.emcee_parameters['version'][0]
        else:
            import emcee
            emcee_version = emcee.__version__[0]

        mc.model_setup()
        """ Required to create the right objects inside each class - if defined inside """
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

    if sampler in sample_keyword['multinest']:
        plot_dictionary['lnprob_chain'] = False
        plot_dictionary['chains'] = False
        plot_dictionary['traces'] = False

        dir_input = './' + config_in['output'] + '/multinest/'
        dir_output = './' + config_in['output'] + '/multinest_plot/'
        os.system('mkdir -p ' + dir_output)

        mc = nested_sampling_load_from_cpickle(dir_input)

        mc.model_setup()
        mc.initialize_logchi2()
        results_analysis.results_resumen(mc, None, skip_theta=True)

        """ Required to create the right objects inside each class - if defined inside """
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

    if sampler in sample_keyword['polychord']:
        plot_dictionary['lnprob_chain'] = False
        plot_dictionary['chains'] = False
        plot_dictionary['traces'] = False

        dir_input = './' + config_in['output'] + '/polychord/'
        dir_output = './' + config_in['output'] + '/polychord_plot/'
        os.system('mkdir -p ' + dir_output)

        mc = nested_sampling_load_from_cpickle(dir_input)

        # pars_input(config_in, mc)

        mc.model_setup()
        mc.initialize_logchi2()
        results_analysis.results_resumen(mc, None, skip_theta=True)

        """ Required to create the right objects inside each class - if defined inside """
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

    print()
    print(' LN posterior: {0:12f}   {1:12f} {2:12f} (15-84 p) '.format(
        lnprob_med[0], lnprob_med[2], lnprob_med[1]))

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
    star_variables_MAP = results_analysis.get_stellar_parameters(mc, chain_MAP)

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
            corner_plot['labels'].append(re.sub('_', '-', var))
            corner_plot['truths'].append(chain_med[var_dict, 0])
            i_corner += 1

        corner_plot['samples'][:, -1] = flat_lnprob[:]
        corner_plot['labels'].append('ln-prob')
        corner_plot['truths'].append(lnprob_med[0])

        if plot_dictionary['use_getdist']:
            print(' Plotting full_correlation plot with GetDist')
            print()
            print(' Ignore the no burn in error warning from getdist')
            print(' since burn in has been already removed from the chains')

            plt.rc('text', usetex=False)

            samples = MCSamples(samples=corner_plot['samples'], names=corner_plot['labels'],
                                labels=corner_plot['labels'])

            g = plots.getSubplotPlotter()
            g.settings.num_plot_contours = 6
            g.triangle_plot(samples, filled=True)
            g.export(dir_output + "all_internal_variables_corner_getdist.pdf")

            print()

        else:
            # plotting mega-corner plot
            print('Plotting full_correlation plot with Corner')
            plt.rc('text', usetex=False)

            fig = corner.corner(
                corner_plot['samples'], labels=corner_plot['labels'], truths=corner_plot['truths'])
            fig.savefig(dir_output + "all_internal_variables_corner_dfm.pdf",
                        bbox_inches='tight', dpi=300)
            plt.close(fig)
            plt.rc('text', usetex=use_tex)

        print()
        print('****************************************************************************************************')
        print()
        sys.stdout.flush()

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

        plt.rc('text', usetex=False)
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

            """ Check if the semi-amplitude K is among the parameters that have been fitted. 
                If so, it computes the correpsing planetary mass with uncertainty """

            fig = corner.corner(np.asarray(corner_plot['samples']).T, labels=corner_plot['labels'],
                                truths=corner_plot['truths'])
            fig.savefig(dir_output + common_name + "_corners.pdf",
                        bbox_inches='tight', dpi=300)
            plt.close(fig)

        print()
        print('****************************************************************************************************')
        print()
        sys.stdout.flush()

    if plot_dictionary['dataset_corner']:

        print(' Dataset + models corner plots ')
        print()

        for dataset_name, dataset in mc.dataset_dict.items():

            for model_name in dataset.models:

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

                corner_plot['samples'] = []
                corner_plot['labels'] = []
                corner_plot['truths'] = []
                for var_i, var in enumerate(variable_values):
                    if np.size(variable_values[var]) <= 1:
                        continue
                    corner_plot['samples'].extend([variable_values[var]])
                    corner_plot['labels'].append(var)
                    corner_plot['truths'].append(variable_median[var])

                fig = corner.corner(np.asarray(corner_plot['samples']).T,
                                    labels=corner_plot['labels'], truths=corner_plot['truths'])
                fig.savefig(dir_output + dataset_name + '_' + model_name +
                            "_corners.pdf", bbox_inches='tight', dpi=300)
                plt.close(fig)

                print('     Dataset: ', dataset_name, '    model: ',
                      model_name, ' corner plot  done ')

        print()
        print('****************************************************************************************************')
        print()
        sys.stdout.flush()

    if plot_dictionary['write_planet_samples']:

        print(' Saving the planet variable samplings to files (with plots)')

        samples_dir = dir_output + '/planet_samples/'
        os.system('mkdir -p ' + samples_dir)

        for common_ref, variable_values in planet_variables.items():
            for variable_name, variable in variable_values.items():

                rad_filename = samples_dir + common_ref + '_' + variable_name
                fileout = open(rad_filename + '.dat', 'w')
                for val in variable:
                    fileout.write('{0:f} \n'.format(val))
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

        for theta_name, th in theta_dictionary.items():

            rad_filename = samples_dir + repr(th) + '_' + theta_name
            fileout = open(rad_filename + '.dat', 'w')
            for val in flat_chain[:, th]:
                fileout.write('{0:f} \n'.format(val))
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

                        for x, tcf, pha, y, e, sys, mod, com, obs_mod, res, jit in zip(
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
                                          '\n'.format(x, tcf, pha, y, e, sys, mod, com, obs_mod, e, res, e, jit))
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
                                fileout.write('{0:f} {1:f} {2:f} {3:f} {4:f} \n'.format(
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
                                    '{0:f} {1:f} {2:f} {3:f}\n'.format(x, tcf, pha, mod))
                            fileout.close()

                        if getattr(mc.models[model_name], 'model_class', False) == 'transit':
                            """
                            Exceptional model writing to deal with under-sampled lightcurves, i.e. when folding the 
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
                                dir_models + dataset_name + '_' + model_name + '_transit.dat', 'w')

                            x_range = np.arange(
                                -variable_values['P']/2., variable_values['P']/2., 0.001)
                            try:
                                delta_T = variable_values['Tc']-dataset.Tref
                            except KeyError:
                                delta_T = kepler_exo.kepler_phase2Tc_Tref(variable_values['P'],
                                                                          variable_values['f'],
                                                                          variable_values['e'],
                                                                          variable_values['o'])

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
                                                          planet_vars[model]['f'],
                                                          planet_vars[model]['P'],
                                                          planet_vars[model]['K'],
                                                          planet_vars[model]['e'],
                                                          planet_vars[model]['o'])
                        fileout = open(
                            dir_models + 'RV_planet_' + model + '_kep.dat', 'w')
                        fileout.write('descriptor x_range  m_kepler \n')
                        for x, y in zip(bjd_plot['full']['x_plot'], RV_out):
                            fileout.write('{0:f} {1:f} \n'.format(x, y))
                        fileout.close()

                        x_range = np.arange(-1.50, 1.50, 0.001)
                        RV_out = kepler_exo.kepler_RV_T0P(x_range * planet_vars[model]['P'],
                                                          planet_vars[model]['f'],
                                                          planet_vars[model]['P'],
                                                          planet_vars[model]['K'],
                                                          planet_vars[model]['e'],
                                                          planet_vars[model]['o'])
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
                                                              planet_vars[model]['f'],
                                                              planet_vars[model]['P'],
                                                              planet_vars[model]['K'],
                                                              planet_vars[model]['e'],
                                                              planet_vars[model]['o'])
                            fileout = open(
                                dir_models + 'RV_planet_' + model + '_Tcf.dat', 'w')
                            fileout.write('descriptor Tc_phase m_phase \n')
                            for x, y in zip(x_range, RV_out):
                                fileout.write('{0:f} {1:f} \n'.format(x, y))
                            fileout.close()

                    except:
                        print('** ERROR **')
                        pass

        print()
        print('****************************************************************************************************')
        print()
        sys.stdout.flush()

    if plot_dictionary['veuz_corner_files']:

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
                    hist2d = np.histogram2d(x_data, y_data, bins=[
                                            x_edges, y_edges], density=True)
                    hist1d_y = np.histogram(y_data, bins=y_edges, density=True)

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

    if sampler in sample_keyword['multinest'] \
       or sampler in sample_keyword['polychord']:

        fig = plt.figure(figsize=(10, 10))

        ii = 0
        for common_ref, variable_values in planet_variables.items():
            if 'P' in variable_values:

                plt.scatter(variable_values['P'],
                            flat_lnprob, s=2, c='C'+repr(ii))
                ii += 1

        rad_filename = dir_output + 'lnprob_P'

        plt.savefig(rad_filename + '.png', bbox_inches='tight', dpi=300)
        plt.close(fig)
