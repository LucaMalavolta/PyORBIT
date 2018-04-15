from classes.model_container_multinest import ModelContainerMultiNest
from classes.model_container_polychord import ModelContainerPolyChord
from classes.model_container_emcee import ModelContainerEmcee

from classes.input_parser import yaml_parser, pars_input
from classes.io_subroutines import *
import numpy as np
import os
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
import sys
mpl.use('Agg')
from matplotlib import pyplot as plt
import corner
import classes.constants as constants
import classes.kepler_exo as kepler_exo
import classes.common as common

__all__ = ["pyorbit_getresults"]


def pyorbit_getresults(config_in, sampler, plot_dictionary):

    #plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
    plt.rcParams["font.family"] = "Times New Roman"
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
        'multinest': ['multinest', 'MultiNest', 'multi'],
        'polychord':['polychord', 'PolyChord', 'polychrod', 'poly'],
        'emcee': ['emcee', 'MCMC', 'Emcee']
    }

    if sampler in sample_keyword['emcee']:

        dir_input = './' + config_in['output'] + '/emcee/'
        dir_output = './' + config_in['output'] + '/emcee_plot/'
        os.system('mkdir -p ' + dir_output)

        mc, starting_point, population, prob, state, \
        sampler_chain, sampler_lnprobability, sampler_acceptance_fraction = \
            emcee_load_from_cpickle(dir_input)

        pars_input(config_in, mc, reload_emcee=True)

        # for retrocompatibility with V5.0 of PyORBIT, the default emcee is set to 2.2.1
        if hasattr(mc.emcee_parameters, 'version'):
            emcee_version = mc.emcee_parameters['version'][0]
        else:
            if sampler_lnprobability.shape[0] > sampler_lnprobability.shape[1]:
                emcee_version = '3'
            else:
                emcee_version = '2'

        mc.model_setup()
        """ Required to create the right objects inside each class - if defined inside """
        theta_dictionary = mc.get_theta_dictionary()

        nburnin = mc.emcee_parameters['nburn']
        nthin = mc.emcee_parameters['thin']
        nsteps = sampler_chain.shape[1] * nthin

        flat_chain = emcee_flatchain(sampler_chain, nburnin, nthin)
        flat_lnprob = emcee_flatlnprob(sampler_lnprobability, nburnin, nthin, emcee_version)

        flat_BiC = -2*flat_lnprob + mc.ndim * np.log(mc.ndata)

        lnprob_med = common.compute_value_sigma(flat_lnprob)
        chain_med = common.compute_value_sigma(flat_chain)
        chain_MAP, lnprob_MAP = common.pick_MAP_parameters(flat_chain, flat_lnprob)

        print
        print 'Reference Time Tref: ', mc.Tref
        print
        print 'Dimensions = ', mc.ndim
        print 'Nwalkers = ', mc.emcee_parameters['nwalkers']
        print
        print 'Steps: ', nsteps
        print

    if sampler in sample_keyword['multinest']:

        plot_dictionary['lnprob_chain'] = False
        plot_dictionary['chains'] = False
        plot_dictionary['traces'] = False

        dir_input = './' + config_in['output'] + '/multinest/'
        dir_output = './' + config_in['output'] + '/multinest_plot/'
        os.system('mkdir -p ' + dir_output)

        mc = polychord_load_from_cpickle(dir_input)

        print mc.bounds
        #pars_input(config_in, mc)

        mc.model_setup()
        mc.initialize_logchi2()
        mc.results_resumen(None, skip_theta=True)

        """ Required to create the right objects inside each class - if defined inside """
        theta_dictionary = mc.get_theta_dictionary()
        print
        print theta_dictionary

        data_in = np.genfromtxt(dir_input + 'post_equal_weights.dat')
        flat_lnprob = data_in[:, -1]
        flat_chain = data_in[:, :-1]
        nsample = np.size(flat_lnprob)


        lnprob_med = common.compute_value_sigma(flat_lnprob)
        chain_med = common.compute_value_sigma(flat_chain)
        chain_MAP, lnprob_MAP = common.pick_MAP_parameters(flat_chain, flat_lnprob)

        print
        print 'Reference Time Tref: ', mc.Tref
        print
        print 'Dimensions = ', mc.ndim
        print
        print 'Samples: ', nsample
        print

    if sampler in sample_keyword['polychord']:

        plot_dictionary['lnprob_chain'] = False
        plot_dictionary['chains'] = False
        plot_dictionary['traces'] = False

        dir_input = './' + config_in['output'] + '/polychord/'
        dir_output = './' + config_in['output'] + '/polychord_plot/'
        os.system('mkdir -p ' + dir_output)

        mc = polychord_load_from_cpickle(dir_input)

        print mc.bounds
        #pars_input(config_in, mc)

        mc.model_setup()
        mc.initialize_logchi2()
        mc.results_resumen(None, skip_theta=True)

        """ Required to create the right objects inside each class - if defined inside """
        theta_dictionary = mc.get_theta_dictionary()
        print theta_dictionary

        data_in = np.genfromtxt(dir_input + 'pyorbit_equal_weights.txt')
        flat_lnprob = data_in[:, 1]
        flat_chain = data_in[:, 2:]
        nsample = np.size(flat_lnprob)


        lnprob_med = common.compute_value_sigma(flat_lnprob)
        chain_med = common.compute_value_sigma(flat_chain)
        chain_MAP, lnprob_MAP = common.pick_MAP_parameters(flat_chain, flat_lnprob)

        print
        print 'Reference Time Tref: ', mc.Tref
        print
        print 'Dimensions = ', mc.ndim
        print
        print 'Samples: ', nsample
        print

    print
    print ' LN posterior: %12f   %12f %12f (15-84 p) ' % (lnprob_med[0], lnprob_med[2], lnprob_med[1])

    MAP_log_priors, MAP_log_likelihood = mc.log_priors_likelihood(chain_MAP)
    BIC = -2.0 * MAP_log_likelihood + np.log(mc.ndata) * mc.ndim
    AIC = -2.0 * MAP_log_likelihood + 2.0 * mc.ndim
    AICc = AIC +  (2.0 + 2.0*mc.ndim) * mc.ndim / (mc.ndata - mc.ndim - 1.0)
    # AICc for small sample

    print
    print ' MAP log_priors     = ', MAP_log_priors
    print ' MAP log_likelihood = ', MAP_log_likelihood
    print ' MAP BIC  (using likelihood)  = ', BIC
    print ' MAP AIC  (using likelihood)  = ', AIC
    print ' MAP AICc (using likelihood) = ', AICc

    MAP_log_posterior = MAP_log_likelihood + MAP_log_priors
    BIC = -2.0 * MAP_log_posterior + np.log(mc.ndata) * mc.ndim
    AIC = -2.0 * MAP_log_posterior + 2.0 * mc.ndim
    AICc = AIC +  (2.0 + 2.0*mc.ndim) * mc.ndim / (mc.ndata - mc.ndim - 1.0)

    print
    print ' MAP BIC  (using posterior)  = ', BIC
    print ' MAP AIC  (using posterior)  = ', AIC
    print ' MAP AICc (using posterior) = ', AICc

    if mc.ndata < 40 * mc.ndim:
        print
        print ' AICc suggested over AIC because NDATA ( %12f ) < 40 * NDIM ( %12f )' % (mc.ndata, mc.ndim)
    else:
        print
        print ' AIC suggested over AICs because NDATA ( %12f ) > 40 * NDIM ( %12f )' % (mc.ndata, mc.ndim)

    print
    print '****************************************************************************************************'
    print
    print ' Print MEDIAN result '
    print

    mc.results_resumen(flat_chain)

    print
    print '****************************************************************************************************'
    print
    print ' Print MAP result (', lnprob_MAP, ')'
    print

    mc.results_resumen(chain_MAP)

    print
    print '****************************************************************************************************'
    print
    print ' Plot FLAT chain '
    print

    if plot_dictionary['lnprob_chain'] or plot_dictionary['chains']:

        #mc.results_resumen(flat_chain)

        if emcee_version == '2':
            fig = plt.figure(figsize=(12, 12))
            plt.xlabel('$\ln \mathcal{L}$')
            plt.plot(sampler_lnprobability.T, '-', alpha=0.5)
            plt.axhline(lnprob_med[0])
            plt.axvline(nburnin/nthin, c='r')
            plt.savefig(dir_output + 'LNprob_chain.png', bbox_inches='tight', dpi=300)
            plt.close(fig)
        else:
            fig = plt.figure(figsize=(12, 12))
            plt.xlabel('$\ln \mathcal{L}$')
            plt.plot(sampler_lnprobability, '-', alpha=0.5)
            plt.axhline(lnprob_med[0])
            plt.axvline(nburnin/nthin, c='r')
            plt.savefig(dir_output + 'LNprob_chain.png', bbox_inches='tight', dpi=300)
            plt.close(fig)

        print
        print '****************************************************************************************************'
        print

    if plot_dictionary['full_correlation']:

        print 'full_correlation plot'
        # plotting mega-corner plot
        plt.rc('text', usetex=False)

        corner_plot = {
            'samples': [],
            'labels': [],
            'truths': []
        }
        for var, var_dict in theta_dictionary.iteritems():
            corner_plot['samples'].extend([flat_chain[:, var_dict]])
            corner_plot['labels'].append(var)
            corner_plot['truths'].append(chain_med[var_dict, 0])

        corner_plot['samples'].extend([flat_lnprob])
        corner_plot['labels'].append('ln prob')
        corner_plot['truths'].append(lnprob_med[0])

        fig = corner.corner(np.asarray(corner_plot['samples']).T,
                            labels=corner_plot['labels'], truths=corner_plot['truths'])
        fig.savefig(dir_output + "all_internal_variables_corner.pdf", bbox_inches='tight', dpi=300)
        plt.close(fig)
        plt.rc('text', usetex=True)

        print
        print '****************************************************************************************************'
        print

    if plot_dictionary['chains']:
        print 'plotting the chains... '

        os.system('mkdir -p ' + dir_output + 'chains')
        for theta_name, ii in theta_dictionary.iteritems():
            file_name = dir_output + 'chains/' + repr(ii) + '_' + theta_name + '.png'
            fig = plt.figure(figsize=(12, 12))
            plt.plot(sampler_chain[:, :, ii].T, '-', alpha=0.5)
            plt.axvline(nburnin/nthin, c='r')
            plt.savefig(file_name, bbox_inches='tight', dpi=300)
            plt.close(fig)

        print
        print '****************************************************************************************************'
        print


    #print sampler.acor.max()
    if plot_dictionary['traces']:
        print 'plotting the Gelman-Rubin traces... '
        print
        os.system('mkdir -p ' + dir_output + 'gr_traces')
        #for theta_name, th in theta_dictionary.iteritems():

        step_sampling = np.arange(nburnin/nthin, nsteps/nthin, 1)

        for theta_name, th in theta_dictionary.iteritems():
            rhat = np.array([GelmanRubin_v2(sampler_chain[:, :steps, th]) for steps in step_sampling])
            print ' Gelman-Rubin: %5i %12f %s ' % (th, rhat[-1], theta_name)
            file_name = dir_output + 'gr_traces/v2_' + repr(th) + '_' + theta_name + '.png'
            fig = plt.figure(figsize=(12, 12))
            plt.plot(step_sampling, rhat[:], '-', color='k')
            #plt.axhline(1.01)
            plt.savefig(file_name, bbox_inches='tight', dpi=300)
            plt.close(fig)

        for theta_name, th in theta_dictionary.iteritems():
            file_name = dir_output + 'gr_traces/' + repr(th) + '_' + theta_name + '.png'
            out_absc = np.arange(0, nburnin/nthin, 1)
            out_lines = np.zeros(nburnin/nthin)

            #for ii in xrange(20, nburnin/nthin):
            #    out_lines[ii] = GelmanRubin(sampler_chain[:, :ii, th].T)
            #fig = plt.figure(figsize=(12, 12))
            #plt.plot(out_absc[20:], out_lines[20:], '-', color='k')
            #plt.axhline(1.01)

            rhat = np.array([GelmanRubin(sampler_chain[:, :steps, th].T) for steps in step_sampling])
            print ' Gelman-Rubin: %5i %12f %s ' % (th, rhat[-1], theta_name)
            fig = plt.figure(figsize=(12, 12))
            plt.plot(step_sampling, rhat[:], '-', color='k')


            plt.savefig(file_name, bbox_inches='tight', dpi=300)
            plt.close(fig)

        print
        print '****************************************************************************************************'
        print

    print

    print
    print '****************************************************************************************************'
    print
    print ' Common models analysis '
    print

    plt.rc('text', usetex=False)

    """ in this variable we store the physical variables of """
    planet_variables = {}
    planet_variables_MAP = {}

    for common_name, common_model in mc.common_models.iteritems():

        print '     Common model: ', common_name

        corner_plot = {
            'var_list': [],
            'samples': [],
            'labels': [],
            'truths': []
        }
        variable_values = common_model.convert(flat_chain)
        variable_median = common_model.convert(chain_med[:, 0])
        variable_MAP = common_model.convert(chain_MAP)

        if len(variable_median) < 1.:
            print
            continue

        n_samplings, n_pams = np.shape(flat_chain)

        print '     variable_median: ', variable_median
        print '     n_samplings, n_pams: ',  n_samplings, n_pams
        print

        """ Sometimes the user forget to include the value for the inclination among the parameters
        We do our check and then we fix it
        """
        i_is_missing = True

        """ 
        Check if the eccentricity and argument of pericenter were set as free parameters or fixed by simply 
        checking the size of their distribution
        """
        for var in variable_values.iterkeys():
            if np.size(variable_values[var]) == 1:
                    variable_values[var] = variable_values[var] * np.ones(n_samplings)
            else:
                corner_plot['var_list'].append(var)

            if var == 'i':
                i_is_missing = False

        corner_plot['samples'] = []
        corner_plot['labels'] = []
        corner_plot['truths'] = []
        for var_i, var in enumerate(corner_plot['var_list']):
            corner_plot['samples'].extend([variable_values[var]])
            corner_plot['labels'].append(var)
            corner_plot['truths'].append(variable_median[var])

        if common_model.model_class == 'planet':
            if i_is_missing:
                if 'i' in common_model.fix_list:
                    variable_values['i'] = np.random.normal(common_model.fix_list['i'][0], common_model.fix_list['i'][1], size=n_samplings)
                else:
                    variable_values['i'] = 90.00 * np.ones(n_samplings)
                    variable_median['i'] = 90.00

            Mass_E = np.empty(n_samplings)
            var_index = np.arange(0, n_samplings, dtype=int)
            star_mass_randomized = np.random.normal(mc.star_mass[0], mc.star_mass[1], size=n_samplings)

            fileout = open(dir_output + common_name + '_orbital_pams.dat', 'w')

            for P, K, e, i, star, ii in zip(
                    variable_values['P'],
                    variable_values['K'],
                    variable_values['e'],
                    variable_values['i'],
                    star_mass_randomized,
                    var_index):
                Mass_E[ii] = M_SEratio / np.sin(np.radians(i)) * \
                             kepler_exo.get_planet_mass(P, K, e, star, Minit=0.0065)

                fileout.write('{0:f} {1:f} {2:f} {3:f} {4:f} {5:f}  \n'.format(P, K, e, i, star, Mass_E[ii]))
            fileout.close()

            Mass_E_med = common.compute_value_sigma(Mass_E)
            variable_median['Me'] = Mass_E_med[0]
            variable_MAP['Me'], _ = common.pick_MAP_parameters(Mass_E, flat_lnprob)

            planet_variables[common_name] = variable_median
            planet_variables_MAP[common_name] = variable_MAP

            fig = plt.figure(figsize=(12, 12))
            plt.hist(Mass_E, bins=50)
            plt.savefig(dir_output + 'planet_mass_' + common_name + '.png', bbox_inches='tight', dpi=300)
            plt.close(fig)

            print 'Planet', common_name, '   Mass (Earths): %12f   %12f %12f (15-84 p) ' % (Mass_E_med[0], Mass_E_med[2], Mass_E_med[1])
            print
            corner_plot['samples'].extend([Mass_E])
            corner_plot['labels'].append('M [M$_\oplus$]')
            corner_plot['truths'].append(Mass_E_med[0])

        if plot_dictionary['common_corner']:

            fig = corner.corner(np.asarray(corner_plot['samples']).T, labels=corner_plot['labels'], truths=corner_plot['truths'])
            fig.savefig(dir_output + common_name + "_corners.pdf", bbox_inches='tight', dpi=300)
            plt.close(fig)

            print 'Common model: ', common_name , '  corner plot done.'
            print

    if plot_dictionary['dataset_corner']:

        print '****************************************************************************************************'
        print
        print ' Dataset + models corner plots '
        print
        for dataset_name, dataset in mc.dataset_dict.iteritems():

            for model_name in dataset.models:

                variable_values = dataset.convert(flat_chain)
                variable_median = dataset.convert(chain_med[:, 0])

                if mc.models[model_name].common_ref:
                    common_ref = mc.models[model_name].common_ref
                    variable_values.update(mc.common_models[common_ref].convert(flat_chain))
                    variable_median.update(mc.common_models[common_ref].convert(chain_med[:, 0]))

                variable_values.update(mc.models[model_name].convert(flat_chain, dataset_name))
                variable_median.update(mc.models[model_name].convert(chain_med[:, 0], dataset_name))

                corner_plot['samples'] = []
                corner_plot['labels'] = []
                corner_plot['truths'] = []
                for var_i, var in enumerate(variable_values):
                    if np.size(variable_values[var]) <= 1: continue
                    corner_plot['samples'].extend([variable_values[var]])
                    corner_plot['labels'].append(var)
                    corner_plot['truths'].append(variable_median[var])

                fig = corner.corner(np.asarray(corner_plot['samples']).T,
                                    labels=corner_plot['labels'], truths=corner_plot['truths'])
                fig.savefig(dir_output + dataset_name + '_' + model_name + "_corners.pdf", bbox_inches='tight', dpi=300)
                plt.close(fig)

                print 'Dataset: ', dataset_name , '    model: ', model_name, ' corner plot  done '

    print
    print '****************************************************************************************************'
    print
    print ' Writing all the files '

    bjd_plot = {
        'full': {
            'start': None, 'end': None, 'range': None
        }
    }

    kinds = {}
    for dataset_name, dataset in mc.dataset_dict.iteritems():
        if dataset.kind in kinds.keys():
            kinds[dataset.kind].extend([dataset_name])
        else:
            kinds[dataset.kind] = [dataset_name]

        bjd_plot[dataset_name] = {
            'start': np.amin(dataset.x0),
            'end': np.amax(dataset.x0),
            'range': np.amax(dataset.x0)-np.amin(dataset.x0),
        }

        if bjd_plot[dataset_name]['range'] < 0.1 : bjd_plot[dataset_name]['range'] = 0.1

        bjd_plot[dataset_name]['start'] -= bjd_plot[dataset_name]['range'] * 0.10
        bjd_plot[dataset_name]['end'] += bjd_plot[dataset_name]['range'] * 0.10
        bjd_plot[dataset_name]['x0_plot'] = \
            np.arange(bjd_plot[dataset_name]['start'], bjd_plot[dataset_name]['end'], 0.1)

        if bjd_plot['full']['range']:
            bjd_plot['full']['start'] = min(bjd_plot['full']['start'], np.amin(dataset.x0))
            bjd_plot['full']['end'] = max(bjd_plot['full']['end'], np.amax(dataset.x0))
            bjd_plot['full']['range'] = bjd_plot['full']['end']-bjd_plot['full']['start']
        else:
            bjd_plot['full']['start'] = np.amin(dataset.x0)
            bjd_plot['full']['end'] = np.amax(dataset.x0)
            bjd_plot['full']['range'] = bjd_plot['full']['end']-bjd_plot['full']['start']

    bjd_plot['full']['start'] -= bjd_plot['full']['range']*0.10
    bjd_plot['full']['end'] += bjd_plot['full']['range']*0.10
    bjd_plot['full']['x0_plot'] = np.arange(bjd_plot['full']['start'], bjd_plot['full']['end'],0.1)

    for dataset_name, dataset in mc.dataset_dict.iteritems():
            bjd_plot[dataset_name] = bjd_plot['full']

    bjd_plot['model_out'], bjd_plot['model_x0'] = mc.get_model(chain_med[:, 0], bjd_plot)
    bjd_plot['MAP_model_out'], bjd_plot['MAP_model_x0'] = mc.get_model(chain_MAP, bjd_plot)

    for kind_name, kind in kinds.iteritems():
        for dataset_name in kind:
            fig = plt.figure(figsize=(12, 12))
            plt.errorbar(mc.dataset_dict[dataset_name].x0,
                         mc.dataset_dict[dataset_name].y - bjd_plot['model_out'][dataset_name]['systematics'],
                         yerr=mc.dataset_dict[dataset_name]. e,
                         fmt='o', zorder=2)
            plt.plot(bjd_plot[dataset_name]['x0_plot'], bjd_plot['model_x0'][dataset_name]['complete'], zorder=2, c='b')
            plt.plot(bjd_plot[dataset_name]['x0_plot'], bjd_plot['MAP_model_x0'][dataset_name]['complete'], zorder=1, c='r')

            plt.savefig(dir_output + 'model_' + kind_name + '_' + dataset_name + '.png', bbox_inches='tight', dpi=300)
            plt.close(fig)
    print

    for prepend_keyword in ['', 'MAP_']:
        plot_out_keyword = prepend_keyword + 'model_out'
        plot_x0_keyword = prepend_keyword + 'model_x0'
        file_keyword = prepend_keyword + 'model_files'

        if prepend_keyword == '':
            planet_vars = planet_variables
        elif prepend_keyword == 'MAP_':
            planet_vars = planet_variables_MAP

        if plot_dictionary[file_keyword]:
            dir_models = dir_output + file_keyword + '/'
            os.system('mkdir -p ' + dir_models)

            for dataset_name, dataset in mc.dataset_dict.items():
                for model_name in dataset.models:
                    fileout = open(dir_models + dataset_name + '_' + model_name + '.dat', 'w')
                    common_ref = mc.models[model_name].common_ref

                    if common_ref in planet_vars:
                        phase = (dataset.x0 / planet_vars[common_ref]['P']) % 1
                    else:
                        phase = dataset.x0 * 0.00

                    fileout.write('descriptor BJD BJD0 pha val,+- sys mod full val_compare,+- res,+- \n')
                    for x, x0, pha, y, e, sys, mod, com, obs_mod, res in zip(
                        dataset.x, dataset.x0, phase, dataset.y, dataset.e,
                            bjd_plot[plot_out_keyword][dataset_name]['systematics'],
                            bjd_plot[plot_out_keyword][dataset_name][model_name],
                            bjd_plot[plot_out_keyword][dataset_name]['complete'],
                            dataset.y - bjd_plot[plot_out_keyword][dataset_name]['complete'] +
                                    bjd_plot[plot_out_keyword][dataset_name][model_name],
                            dataset.y - bjd_plot[plot_out_keyword][dataset_name]['complete']):

                        fileout.write('{0:f} {1:f} {2:f} {3:f} {4:f} {5:f} {6:1f} {7:f} {8:f} {9:f} {10:f} {11:f}'
                                      '\n'.format(x, x0, pha, y, e, sys, mod, com, obs_mod, e, res, e))
                    fileout.close()

                    fileout = open(dir_models + dataset_name + '_' + model_name + '_full.dat', 'w')

                    if model_name+'_std' in bjd_plot[plot_x0_keyword][dataset_name]:
                        fileout.write('descriptor BJD BJD0 mod,+- \n')
                        for x0, mod, std in zip(bjd_plot[dataset_name]['x0_plot'],
                                           bjd_plot[plot_x0_keyword][dataset_name][model_name],
                                           bjd_plot[plot_x0_keyword][dataset_name][model_name+'_std']):
                            fileout.write('{0:f} {1:f} {2:f} {3:f} \n'.format(x0+mc.Tref, x0, mod, std))
                        fileout.close()
                    else:
                        fileout.write('descriptor BJD BJD0 mod \n')
                        for x0, mod in zip(bjd_plot[dataset_name]['x0_plot'],
                                           bjd_plot[plot_x0_keyword][dataset_name][model_name]):
                            fileout.write('{0:f} {1:f} {2:f} \n'.format(x0+mc.Tref, x0, mod))
                        fileout.close()

                fileout = open(dir_models + dataset_name + '_full.dat', 'w')
                fileout.write('descriptor BJD BJD0 mod \n')
                for x0, mod in zip(bjd_plot[dataset_name]['x0_plot'],
                                   bjd_plot[plot_x0_keyword][dataset_name]['complete']):
                    fileout.write('{0:f} {1:f} {2:f} \n'.format(x0+mc.Tref, x0, mod))
                fileout.close()

            for model in planet_vars:

                RV_out =  kepler_exo.kepler_RV_T0P(bjd_plot['full']['x0_plot'],
                                                   planet_vars[model]['f'],
                                                   planet_vars[model]['P'],
                                                   planet_vars[model]['K'],
                                                   planet_vars[model]['e'],
                                                   planet_vars[model]['o'])
                fileout = open(dir_models + 'RV_planet_' + model + '_kep.dat', 'w')
                fileout.write('descriptor x_range x_range0 m_kepler \n')
                for x, y in zip(bjd_plot['full']['x0_plot'], RV_out):
                    fileout.write('{0:f} {1:f} {2:f} \n'.format(x+mc.Tref, x, y))
                fileout.close()

                x_range = np.arange(-0.50, 1.50, 0.001)
                RV_out = kepler_exo.kepler_RV_T0P(x_range*planet_vars[model]['P'],
                                                   planet_vars[model]['f'],
                                                   planet_vars[model]['P'],
                                                   planet_vars[model]['K'],
                                                   planet_vars[model]['e'],
                                                   planet_vars[model]['o'])
                fileout = open(dir_models + 'RV_planet_' + model + '_pha.dat', 'w')
                fileout.write('descriptor x_phase m_phase \n')
                for x, y in zip(x_range, RV_out):
                    fileout.write('{0:f} {1:f} \n'.format(x, y))
                fileout.close()





