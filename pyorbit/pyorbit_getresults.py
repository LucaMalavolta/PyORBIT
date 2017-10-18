from classes.model_container import ModelContainer
from classes.input_parser import yaml_parser, pars_input
from classes.io_subroutines import pyde_save_to_pickle, pyde_load_from_cpickle, \
    emcee_save_to_cpickle, emcee_load_from_cpickle, emcee_flatchain, emcee_flatlnprob, \
    GelmanRubin, GelmanRubin_v2
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

        mc.model_setup()
        """ Required to create the right objects inside each class - if defined inside """
        theta_dictionary = mc.get_theta_dictionary()

        nburnin = mc.emcee_parameters['nburn']
        nthin = mc.emcee_parameters['thin']
        nsteps = sampler_chain.shape[1] * nthin

        flat_chain = emcee_flatchain(sampler_chain, nburnin, nthin)
        flat_lnprob = emcee_flatlnprob(sampler_lnprobability, nburnin, nthin)

        print
        print 'Reference Time Tref: ', mc.Tref
        print
        print 'Dimensions = ', mc.ndim
        print 'Nwalkers = ', mc.emcee_parameters['nwalkers']
        print
        print 'Steps: ', nsteps

        chain_med = common.compute_value_sigma(flat_chain)
        mc.results_resumen(flat_chain)

        lnprob_med = common.compute_value_sigma(flat_lnprob)
        print ' LN probability: %12f   %12f %12f (15-84 p) ' % (lnprob_med[0], lnprob_med[2], lnprob_med[1])
        #mc.results_resumen(flat_chain)

        fig = plt.figure(figsize=(12, 12))
        plt.xlabel('$\ln \mathcal{L}$')
        plt.plot(sampler_lnprobability.T, '-', alpha=0.5)
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
                corner_plot['truths'].append(chain_med[var_dict,0])

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


        print sampler.acor.max()
        if plot_dictionary['traces']:
            print 'plotting the Gelman-Rubin traces... '
            print
            os.system('mkdir -p ' + dir_output + 'gr_traces')
            #for theta_name, th in theta_dictionary.iteritems():

            step_sampling = np.arange(nburnin/nthin, nsteps/nthin, 1)

            for theta_name, th in theta_dictionary.iteritems():
                rhat = np.array([GelmanRubin_v2(sampler_chain[:, :steps, th]) for steps in step_sampling])
                print ' LN probability: %s  %5i %12f ' % (theta_name, th, rhat[-1])
                print theta_name, th,
                file_name = dir_output + 'gr_traces/v2_' + repr(th) + '_' + theta_name + '.png'
                fig = plt.figure(figsize=(12, 12))
                plt.plot(step_sampling, rhat[:], '-', color='k')
                plt.axhline(1.01)
                plt.savefig(file_name, bbox_inches='tight', dpi=300)
                plt.close(fig)

            for theta_name, th in theta_dictionary.iteritems():
                file_name = dir_output + 'gr_traces/' + repr(th) + '_' + theta_name + '.png'
                out_absc = np.arange(0, nburnin/nthin, 1)
                out_lines = np.zeros(nburnin/nthin)

                for ii in xrange(20, nburnin/nthin):
                    out_lines[ii] = GelmanRubin(sampler_chain[:, :ii, th].T)
                fig = plt.figure(figsize=(12, 12))
                plt.plot(out_absc[20:], out_lines[20:], '-', color='k')
                plt.axhline(1.01)
                plt.savefig(file_name, bbox_inches='tight', dpi=300)
                plt.close(fig)

            print
            print '****************************************************************************************************'
            print

        print

        print
        print '****************************************************************************************************'
        print
        print ' Common models corner plots '
        print

        plt.rc('text', usetex=False)

        """ in this variable we store the physical variables of """
        planet_variables = {}

        for common_name, common_model in mc.common_models.iteritems():

            corner_plot = {
                'var_list': [],
                'samples': [],
                'labels': [],
                'truths': []
            }
            variable_values = common_model.convert(flat_chain)
            variable_median = common_model.convert(chain_med[:,0])

            n_samplings, n_pams = np.shape(flat_chain)

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
                    variable_values['i'] = 90.00 * np.ones(n_samplings)
                    variable_median['i'] = 90.00

                Mass_E = np.empty(n_samplings)
                var_index = np.arange(0, n_samplings, dtype=int)
                star_mass_randomized = np.random.normal(mc.star_mass[0], mc.star_mass[1], size=n_samplings)
                for P, K, e, i, star, ii in zip(
                        variable_values['P'],
                        variable_values['K'],
                        variable_values['e'],
                        variable_values['i'],
                        star_mass_randomized,
                        var_index):
                    Mass_E[ii] = M_SEratio * np.sin(np.radians(i)) * \
                                 kepler_exo.get_planet_mass(P, K, e, star, Minit=0.0065)

                Mass_E_med = common.compute_value_sigma(Mass_E)
                variable_median['Me'] = Mass_E_med[0]

                planet_variables[common_name] = variable_median

                fig = plt.figure(figsize=(12, 12))
                plt.hist(Mass_E, bins=50)
                plt.savefig(dir_output + 'planet_mass_' + common_name + '.png', bbox_inches='tight', dpi=300)
                plt.close(fig)

                print ' Mass (Earths): %12f   %12f %12f (15-84 p) ' % (Mass_E_med[0], Mass_E_med[2], Mass_E_med[1])
                print

                corner_plot['samples'].extend([Mass_E])
                corner_plot['labels'].append('M [M$_\oplus$]')
                corner_plot['truths'].append(Mass_E_med[0])

            fig = corner.corner(np.asarray(corner_plot['samples']).T, labels=corner_plot['labels'], truths=corner_plot['truths'])
            fig.savefig(dir_output + common_name + "_corners.pdf", bbox_inches='tight', dpi=300)
            plt.close(fig)

            print 'Common model: ', common_name , '  corner plot done.'

        print
        print '****************************************************************************************************'
        print
        print ' Dataset + models corner plots '
        print
        for dataset_name, dataset in mc.dataset_dict.iteritems():

            for model_name in dataset.models:

                common_ref = mc.models[model_name].common_ref
                variable_values = dataset.convert(flat_chain)
                variable_values.update(mc.common_models[common_ref].convert(flat_chain))
                variable_values.update(mc.models[model_name].convert(flat_chain, dataset_name))

                variable_median = dataset.convert(chain_med[:, 0])
                variable_median.update(mc.common_models[common_ref].convert(chain_med[:, 0]))
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
        print ' Dataset + models corner plots '

        bjd_plot = {
            'full': {
                'start':None, 'end': None, 'range': None
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

            bjd_plot[dataset_name]['start'] -= bjd_plot[dataset_name]['range'] * 0.05
            bjd_plot[dataset_name]['end'] += bjd_plot[dataset_name]['range'] * 0.05
            bjd_plot[dataset_name]['x0_plot'] = \
                np.arange(bjd_plot[dataset_name]['start'], bjd_plot[dataset_name]['end'], 0.1)

            if bjd_plot['full']['range']:
                bjd_plot['full']['start'] = min(bjd_plot['full']['start'], np.amin(dataset.x0))
                bjd_plot['full']['end'] = min(bjd_plot['full']['start'], np.amax(dataset.x0))
                bjd_plot['full']['range'] = bjd_plot['full']['end']-bjd_plot['full']['start']
            else:
                bjd_plot['full']['start'] = np.amin(dataset.x0)
                bjd_plot['full']['end'] = np.amax(dataset.x0)
                bjd_plot['full']['range'] = bjd_plot['full']['end']-bjd_plot['full']['start']

        bjd_plot['full']['start'] -= bjd_plot['full']['range']*0.05
        bjd_plot['full']['end'] += bjd_plot['full']['range']*0.05
        bjd_plot['full']['x0_plot'] = np.arange(bjd_plot['full']['start'],bjd_plot['full']['end'],0.1)

        for dataset_name, dataset in mc.dataset_dict.iteritems():
            if dataset.dynamical:
                bjd_plot[dataset_name] = bjd_plot['full']

        bjd_plot['model_out'], bjd_plot['model_x0'] = mc.get_model(chain_med[:, 0], bjd_plot)

        for kind_name, kind in kinds.iteritems():
            for dataset_name in kind:
                fig = plt.figure(figsize=(12, 12))
                plt.errorbar(mc.dataset_dict[dataset_name].x0,
                             mc.dataset_dict[dataset_name].y - bjd_plot['model_out'][dataset_name]['systematics'],
                             yerr=mc.dataset_dict[dataset_name]. e,
                             fmt='o', zorder=2)
                plt.plot(bjd_plot[dataset_name]['x0_plot'], bjd_plot['model_x0'][dataset_name]['complete'], zorder=1)

                plt.savefig(dir_output + 'model_' + kind_name + '_' + dataset_name + '.png', bbox_inches='tight', dpi=300)
                plt.close(fig)
        print

        if plot_dictionary['model_files']:
            os.system('mkdir -p ' + dir_output + 'model_files')

            for dataset_name, dataset in mc.dataset_dict.items():
                for model_name in dataset.models:
                    fileout = open(dir_output + 'model_files/' + dataset_name + '_' + model_name + '.dat', 'w')

                    common_ref = mc.models[model_name].common_ref

                    if common_ref in planet_variables:
                        phase = (dataset.x0 / planet_variables[common_ref]['P']) % 1
                    else:
                        phase = dataset.x0 * 0.00

                    fileout.write('descriptor BJD BJD0 pha val,+- sys mod full val_compare,+- res,+- \n')
                    for x, x0, pha, y, e, sys, mod, com, obs_mod, res in zip(
                        dataset.x, dataset.x0, phase, dataset.y, dataset.e,
                            bjd_plot['model_out'][dataset_name]['systematics'],
                            bjd_plot['model_out'][dataset_name][model_name],
                            bjd_plot['model_out'][dataset_name]['complete'],
                            dataset.y - bjd_plot['model_out'][dataset_name]['complete'] +
                                    bjd_plot['model_out'][dataset_name][model_name],
                            dataset.y - bjd_plot['model_out'][dataset_name]['complete']):

                        fileout.write('{0:f} {1:f} {2:f} {3:f} {4:f} {5:f} {6:1f} {7:f} {8:f} {9:f} {10:f} {11:f}'
                                      '\n'.format(x, x0, pha, y, e, sys, mod, com, obs_mod, e, res, e))
                    fileout.close()

                    fileout = open(dir_output + 'model_files/' + dataset_name + '_' + model_name + '_full.dat', 'w')
                    fileout.write('descriptor BJD BJD0 mod \n')
                    for x0, mod in zip(bjd_plot[dataset_name]['x0_plot'],
                                       bjd_plot['model_x0'][dataset_name]['complete']):
                        fileout.write('{0:f} {1:f} {2:f} \n'.format(x0+mc.Tref, x0, mod))
                    fileout.close()

            for model in planet_variables:
                RV_out =  kepler_exo.kepler_RV_T0P(bjd_plot['full']['x0_plot'],
                                                   planet_variables[model]['f'],
                                                   planet_variables[model]['P'],
                                                   planet_variables[model]['K'],
                                                   planet_variables[model]['e'],
                                                   planet_variables[model]['o'])
                fileout = open(dir_output + 'model_files/' + 'RV_' + model_name + '_kep.dat', 'w')
                fileout.write('descriptor x_range x_range0 m_kepler \n')
                for x, y in zip(bjd_plot['full']['x0_plot'], RV_out):
                    fileout.write('{0:f} {1:f} {2:f} \n'.format(x+mc.Tref, x, y))
                fileout.close()

                x_range = np.arange(-0.50, 1.50, 0.001)
                RV_out =  kepler_exo.kepler_RV_T0P(x_range*planet_variables[model]['P'],
                                                   planet_variables[model]['f'],
                                                   planet_variables[model]['P'],
                                                   planet_variables[model]['K'],
                                                   planet_variables[model]['e'],
                                                   planet_variables[model]['o'])
                fileout = open(dir_output + 'model_files/' + 'RV_' + model_name + '_pha.dat', 'w')
                fileout.write('descriptor x_phase0 m_phase \n')
                for x, y in zip(x_range, RV_out):
                    fileout.write('{0:f} {1:f} \n'.format(x, y))
                fileout.close()
