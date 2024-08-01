from __future__ import print_function
from pyorbit.subroutines.common import *
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo

from tqdm import tqdm

__all__ = ["print_bayesian_info", "results_summary", "results_derived", "get_planet_parameters", "get_theta_dictionary", "get_model",
           "print_theta_bounds", "print_dictionary", "get_stellar_parameters", "print_integrated_ACF"]


def print_bayesian_info(mc):

    print()
    print('====================================================================================================')
    print('     Ids, spaces (s), boundaries (b) and priors (p) of the sampler parameters     ')
    print('====================================================================================================')
    print()
    for dataset_name, dataset in mc.dataset_dict.items():

        print('----- dataset: ', dataset_name)

        print_analysis_info(dataset.sampler_parameters,
                            mc.bounds,
                            mc.spaces,
                            mc.priors,
                            dataset.prior_kind,
                            dataset.prior_pams)

        for model_name in dataset.models:
            if len(mc.models[model_name].sampler_parameters[dataset_name])==0: continue
            print('----- dataset: {0:25s} ----- model: {1:30s}'.format(dataset_name,model_name))
            print_analysis_info(mc.models[model_name].sampler_parameters[dataset_name],
                                mc.bounds,
                                mc.spaces,
                                mc.priors,
                                mc.models[model_name].prior_kind[dataset_name],
                                mc.models[model_name].prior_pams[dataset_name])

    for model_name, model in mc.common_models.items():
        print('----- common model: ', model_name)
        print_analysis_info(model.sampler_parameters,
                            mc.bounds,
                            mc.spaces,
                            mc.priors,
                            model.prior_kind,
                            model.prior_pams)


def results_summary(mc, theta,
                    skip_theta=False,
                    compute_lnprob=False,
                    chain_med=False,
                    return_samples=False,
                    is_starting_point=False,
                    is_MAP=False):
    # Function with two goals:
    # * Unfold and print out the output from theta
    # * give back a parameter name associated to each value in the result array

    if is_starting_point or is_MAP:
        fixed_warning = False
    else:
        fixed_warning = True


    print('====================================================================================================')
    if is_starting_point:
        print('     Starting point of the sample/optimization routines    ')
    else:
        print('     Statistics on the posterior of the sampler parameters     ')

    print('====================================================================================================')
    print()

    for dataset_name, dataset in mc.dataset_dict.items():

        print('----- dataset: ', dataset_name)
        print_theta_bounds(dataset.sampler_parameters, theta, mc.bounds)


        for model_name in dataset.models:
            if len(mc.models[model_name].sampler_parameters[dataset_name])==0: continue
            print('----- dataset: {0:25s} ----- model: {1:30s}'.format(dataset_name,model_name))
            print_theta_bounds(mc.models[model_name].sampler_parameters[dataset_name], theta, mc.bounds)

    for model_name, model in mc.common_models.items():
        print('----- common model: ', model_name)
        print_theta_bounds(model.sampler_parameters, theta, mc.bounds)


    print('====================================================================================================')
    if is_starting_point:
        print('     Starting point projected onto the model parameter space     ')
    else:
        print('     Statistics on the model parameters obtained from the posteriors samples     ')
    print('====================================================================================================')
    print()

    for dataset_name, dataset in mc.dataset_dict.items():
        print('----- dataset: ', dataset_name)
        parameter_values = dataset.convert(theta)
        print_dictionary(parameter_values, fixed_warning=fixed_warning)

        print()
        for model_name in dataset.models:
            parameter_values = mc.models[model_name].convert(
                theta, dataset_name)
            if len(parameter_values)==0: continue
            print('----- dataset: {0:25s} ----- model: {1:30s}'.format(dataset_name,model_name))
            print_dictionary(parameter_values, fixed_warning=fixed_warning)

    for model_name, model in mc.common_models.items():
        print('----- common model: ', model_name)
        parameter_values = model.convert(theta)
        if chain_med is not False:
            recenter_pams = {}
            parameter_values_med = model.convert(chain_med)

            for par in list(OrderedSet(model.recenter_pams) & OrderedSet(parameter_values_med)):
                recenter_pams[par] = [parameter_values_med[par],
                                      model.default_bounds[par][1] - model.default_bounds[par][0]]
            print_dictionary(parameter_values, recenter=recenter_pams, fixed_warning=fixed_warning)

        else:
            print_dictionary(parameter_values, fixed_warning=fixed_warning)

    print('====================================================================================================')
    if is_starting_point:
        print('     Derived parameters obtained from starting point     ')
    else:
        print('     Statistics on the derived parameters obtained from the posteriors samples     ')
    print('====================================================================================================')
    print()

    returned_samples = get_planet_parameters(mc, theta, verbose=True)

    if compute_lnprob:
        print()
        print('====================================================================================================')
        print('     Statistics on the log-likelihood     ')
        print('====================================================================================================')
        print()

        if len(np.shape(theta)) == 2:
            n_samples, n_values = np.shape(theta)
            logchi2_collection = np.zeros(n_samples)
            for i in range(0, n_samples):
                logchi2_collection[i] = mc(theta[i, :])
            perc0, perc1, perc2 = np.percentile(
                logchi2_collection, [15.865, 50, 84.135], axis=0)
            print(' LN probability: %12f   %12f %12f (15-84 p) ' %
                  (perc1, perc0 - perc1, perc2 - perc1))
        else:
            print(' LN probability: %12f ' % (mc(theta)))

    print()
    print('====================================================================================================')
    print('     ------------------------------------------------------------------------------------------     ')
    print('====================================================================================================')
    print()
    print()

    if return_samples:
        return returned_samples


def get_stellar_parameters(mc, theta, warnings=True, stellar_ref=None):
    try:
        n_samplings, n_pams = np.shape(theta)
    except:
        n_samplings = 1
        n_pams = np.shape(theta)

    "Stellar mass, radius and density are pre-loaded since they are are required by most of the common models"

    """
    #TODO associate to each planet the corresponding stellar parameters
    """
    try:
        stellar_model = mc.common_models[stellar_ref]
    except:
        for model_name, model_obj in mc.common_models.items():
            if getattr(model_obj,'model_class', None) == 'star_parameters':
                stellar_model = mc.common_models[model_name]

    stellar_values = stellar_model.convert(theta)

    if 'density' not in stellar_values:

        if 'radius' not in stellar_values:
            try:
                if stellar_model.prior_kind['radius'] == 'Gaussian':
                    stellar_values['radius'] = np.random.normal(stellar_model.prior_pams['radius'][0],
                                                                stellar_model.prior_pams['radius'][1],
                                                                size=n_samplings)
            except:
                if warnings:
                    print(' *** Please provide a prior on stellar mass, radius or density *** ')
                    print()

        if 'mass' not in stellar_values:
            try:
                if stellar_model.prior_kind['mass'] == 'Gaussian':
                    stellar_values['mass'] = np.random.normal(stellar_model.prior_pams['mass'][0],
                                                              stellar_model.prior_pams['mass'][1],
                                                              size=n_samplings)
            except:
                if warnings:
                    print(' *** Please provide a prior on stellar mass, radius, or density *** ')
                    print()

        if 'mass' in stellar_values.keys() and 'radius' in stellar_values.keys():
            stellar_values['density'] = stellar_values['mass'] / \
                stellar_values['radius'] ** 3
            if warnings:
                print('Note: stellar density derived from stellar mass and radius')
                print()

    else:
        if 'mass' in stellar_values and 'radius' in stellar_values:
            if warnings:
                print('Note: stellar mass, radius and density provided independently, no check for consistency is performed')
                print()
        elif 'mass' in stellar_values:
            stellar_values['radius'] = (
                stellar_values['mass'] / stellar_values['density']) ** (1. / 3.)
            if warnings:
                print('Note: stellar radius derived from stellar mass and density')
                print()
        elif 'radius' in stellar_values:
            stellar_values['mass'] = stellar_values['radius'] ** 3. * stellar_values['density']
            if warnings:
                print('Note: stellar mass derived from stellar radius and density')
                print()

        else:
            if 'mass' in stellar_model.prior_pams and 'radius' in stellar_model.prior_pams:
                if warnings:
                    print('Note: priors on stellar mass and radius provided independently from the measured density, no check for consistency is performed')
                    print()
                if stellar_model.prior_kind['mass'] == 'Gaussian':
                        stellar_values['mass'] = np.random.normal(stellar_model.prior_pams['mass'][0],
                                                                  stellar_model.prior_pams['mass'][1],
                                                                  size=n_samplings)
                if stellar_model.prior_kind['radius'] == 'Gaussian':
                    stellar_values['radius'] = np.random.normal(stellar_model.prior_pams['radius'][0],
                                                                stellar_model.prior_pams['radius'][1],
                                                                size=n_samplings)
            elif 'mass' in stellar_model.prior_pams:
                if stellar_model.prior_kind['mass'] == 'Gaussian':
                    stellar_values['mass'] = np.random.normal(stellar_model.prior_pams['mass'][0],
                                                              stellar_model.prior_pams['mass'][1],
                                                              size=n_samplings)
                    stellar_values['radius'] = (
                        stellar_values['mass'] / stellar_values['density']) ** (1. / 3.)
                if warnings:
                    print('Note: stellar radius derived from its measured density and its prior on mass')
                    print()

            elif 'radius' in stellar_model.prior_pams:

                if stellar_model.prior_kind['radius'] == 'Gaussian':
                    stellar_values['radius'] = np.random.normal(stellar_model.prior_pams['radius'][0],
                                                                stellar_model.prior_pams['radius'][1],
                                                                size=n_samplings)
                    stellar_values['mass'] = stellar_values['radius'] ** 3. * \
                        stellar_values['density']
                if warnings:
                    print('Note: stellar radius derived from its measured density and its prior on mass')
                    print()

            else:
                if warnings:
                    print(
                        ' *** Please provide a prior either on stellar Mass or stellar Radius *** ')
                    print()

    return stellar_values


def results_derived(mc, theta):
    _ = get_planet_parameters(mc, theta, verbose=True)


def get_planet_parameters(mc, theta, verbose=False):
    """
    Derived parameters from the Common Models are listed

    :param mc:
    :return:
    """

    try:
        n_samplings, n_pams = np.shape(theta)
    except:
        n_samplings = 1

    planet_parameters = {}

    for common_name, common_model in mc.common_models.items():
        parameter_values = common_model.convert(theta)
        derived_parameters = {}

        if common_model.model_class == 'planet':

            stellar_ref = getattr(common_model, 'stellar_ref', None)

            stellar_parameters = get_stellar_parameters(mc,
                                                    theta,
                                                    warnings=False,
                                                    stellar_ref=stellar_ref)

            remove_i = False
            if verbose:
                print('----- common model: ', common_name)

            """
            Check if the eccentricity and argument of pericenter were set as free parameters or fixed by simply
            checking the size of their distribution
            """

            for par in parameter_values.keys():
                if np.size(parameter_values[par]) == 1:
                    parameter_values[par] = parameter_values[par] * \
                        np.ones(n_samplings)

            if 'a_Rs' not in parameter_values.keys() and 'density' in stellar_parameters.keys():
                derived_parameters['a_Rs'] = True
                parameter_values['a_Rs'] = convert_rho_to_ars(parameter_values['P'],
                                                        stellar_parameters['density'])

            if 'i' not in parameter_values.keys():
                derived_parameters['i'] = True
                if 'i' in common_model.fix_list:

                    if verbose:
                        print('Inclination randomized to {0:3.2f} +- {1:3.2f} deg'.format(
                              common_model.fix_list['i'][0], common_model.fix_list['i'][1]))
                    parameter_values['i'] = np.random.normal(common_model.fix_list['i'][0],
                                                            common_model.fix_list['i'][1],
                                                            size=n_samplings)
                elif 'b' in parameter_values.keys() and 'a_Rs' in parameter_values.keys():
                    parameter_values['i'] = convert_b_to_i(parameter_values['b'],
                                                          parameter_values['e'],
                                                          parameter_values['omega'],
                                                          parameter_values['a_Rs'])
                else:
                    print('Inclination fixed to 90 deg!')
                    parameter_values['i'] = 90.00 * np.ones(n_samplings)
                    remove_i = True

            if 'K' in parameter_values.keys() and 'mass' in stellar_parameters.keys():
                M_Msun = kepler_exo.get_planet_mass(parameter_values['P'],
                                                                  parameter_values['K'],
                                                                  parameter_values['e'],
                                                                  stellar_parameters['mass']) \
                    / np.sin(np.radians(parameter_values['i']))

                parameter_values['M_Mj'] =M_Msun * constants.Msjup
                derived_parameters['M_Mj'] = True

                parameter_values['M_Me']  =M_Msun * constants.Msear
                derived_parameters['M_Me'] = True

            elif 'M_Me' in parameter_values.keys() and 'mass' in stellar_parameters.keys():
                derived_parameters['K'] = True
                derived_parameters['K'] = kepler_exo.kepler_K1(stellar_parameters['mass'],
                                                              parameter_values['M_Me'] /
                                                              constants.Msear,
                                                              parameter_values['P'],
                                                              parameter_values['i'],
                                                              parameter_values['e'])

                parameter_values['M_Mj'] = parameter_values['M_Me'] * \
                    (constants.Msjup/constants.Msear)
                derived_parameters['M_Mj'] = True


            if 'Tc' in parameter_values.keys():
                if 'e' in parameter_values:
                    derived_parameters['mean_long'] = True
                    parameter_values['mean_long'] = kepler_exo.kepler_Tc2phase_Tref(parameter_values['P'],
                                                                           parameter_values['Tc'] -
                                                                           mc.Tref,
                                                                           parameter_values['e'],
                                                                           parameter_values['omega'])

            elif 'mean_long' in parameter_values.keys():
                derived_parameters['Tc'] = True
                parameter_values['Tc'] = mc.Tref + kepler_exo.kepler_phase2Tc_Tref(parameter_values['P'],
                                                                                  parameter_values['mean_long'],
                                                                                  parameter_values['e'],
                                                                                  parameter_values['omega'])

            if 'R_Rs' in parameter_values.keys() and 'radius' in stellar_parameters.keys():
                parameter_values['R_Rj'] = parameter_values['R_Rs'] * \
                    constants.Rsjup * stellar_parameters['radius']
                derived_parameters['R_Rj'] = True

                parameter_values['R_Re'] = parameter_values['R_Rs'] * \
                    constants.Rsear * stellar_parameters['radius']
                derived_parameters['R_Re'] = True

            if remove_i:
                del parameter_values['i']

            try:
                k = parameter_values['R_Rs']

                parameter_values['T_41'] = parameter_values['P'] / np.pi \
                    * np.arcsin(1./parameter_values['a_Rs'] *
                                np.sqrt((1. + k)**2 - parameter_values['b']**2)
                                / np.sin(parameter_values['i']*constants.deg2rad))
                derived_parameters['T_41'] = True

                parameter_values['T_32'] = parameter_values['P'] / np.pi \
                    * np.arcsin(1./parameter_values['a_Rs'] *
                                np.sqrt((1. - k)**2 - parameter_values['b']**2)
                                / np.sin(parameter_values['i']*constants.deg2rad))
                derived_parameters['T_32'] = True

            except:
                pass

            try:
                derived_parameters['a_AU_(M)'] = True
                parameter_values['a_AU_(M)'] = convert_PMsMp_to_a(
                    parameter_values['P'],
                    stellar_parameters['mass'],
                    parameter_values['M_Me'])
            except (KeyError, ValueError):
                pass

            try:
                derived_parameters['a_AU_(rho,R)'] = True
                parameter_values['a_AU_(rho,R)'] = convert_ars_to_a(
                    parameter_values['a_Rs'],
                    stellar_parameters['radius'])
            except (KeyError, ValueError):
                pass

            try:
                derived_parameters['insol'] = True
                parameter_values['insol]'] = \
                    convert_RTaAU_to_insol(
                        stellar_parameters['radius'],
                        stellar_parameters['temperature'],
                        parameter_values['a_AU_(rho,R)']
                    )
            except (KeyError, ValueError):
                pass

            planet_parameters[common_name] = parameter_values.copy()

            for par in planet_parameters[common_name].keys():
                if par not in derived_parameters.keys():
                    del parameter_values[par]

            if verbose:
                print_dictionary(parameter_values, fixed_warning=False)

    return planet_parameters


def get_theta_dictionary(mc):
    # * give back a parameter name associated to each value in the result array

    theta_dictionary = {}
    for dataset_name, dataset in mc.dataset_dict.items():
        for par, i in dataset.sampler_parameters.items():
            try:
                theta_dictionary[dataset_name + '_' + par] = i
            except:
                theta_dictionary[repr(dataset_name) + '_' + par] = i

        for model_name in dataset.models:
            for par, i in mc.models[model_name].sampler_parameters[dataset_name].items():
                try:
                    theta_dictionary[dataset_name +
                                     '_' + model_name + '_' + par] = i
                except:
                    theta_dictionary[repr(dataset_name) +
                                     '_' + model_name + '_' + par] = i

    for model_name, model in mc.common_models.items():
        for par, i in model.sampler_parameters.items():
            theta_dictionary[model.common_ref + '_' + par] = i

    return theta_dictionary


def get_model(mc, theta, bjd_dict, **kwargs):
    model_out = {}
    model_x0 = {}

    delayed_lnlk_computation = {}

    progress_bar = kwargs.get('progress_bar', True)

    if mc.dynamical_model is not None:
        """ check if any keyword ahas get the output model from the dynamical tool
        we must do it here because all the planet are involved"""
        dynamical_output_x0 = mc.dynamical_model.compute(
            mc, theta, bjd_dict['full']['x_plot'])
        dynamical_output = mc.dynamical_model.compute(mc, theta)

    for dataset_name, dataset in mc.dataset_dict.items():

        if not getattr(dataset, 'compute_plot', True):
            continue

        x0_plot = bjd_dict[dataset_name]['x0_plot']

        n_input = np.size(x0_plot)
        model_out[dataset_name] = {}
        model_x0[dataset_name] = {}
        dataset.model_reset()

        additive_model = np.zeros(np.size(x0_plot))
        unitary_model = np.zeros(np.size(x0_plot))
        external_model = np.zeros(np.size(x0_plot))
        normalization_model = None

        parameter_values = dataset.convert(theta)
        dataset.compute(parameter_values)

        for model_name in dataset.models:
            parameter_values = {}

            for common_ref in mc.models[model_name].common_ref:
                parameter_values.update(
                    mc.common_models[common_ref].convert(theta))

            #TODO: remove try-except starting from version 11 !!
            try:
                for planet_name in mc.models[model_name].multiple_planets:
                    parameter_values.update(
                        mc.common_models[planet_name].convert_with_name(theta, planet_name))
            except:
                pass

            parameter_values.update(
                mc.models[model_name].convert(theta, dataset_name))

            if getattr(mc.models[model_name], 'jitter_model', False):
                dataset.jitter += mc.models[model_name].compute(
                    parameter_values, dataset)
                continue

            if getattr(mc.models[model_name], 'systematic_model', False):
                dataset.additive_model += mc.models[model_name].compute(
                    parameter_values, dataset)

        model_out[dataset_name]['systematics'] = dataset.additive_model.copy()
        model_out[dataset_name]['jitter'] = dataset.jitter.copy()
        model_out[dataset_name]['complete'] = np.zeros(
            dataset.n, dtype=np.double)
        model_out[dataset_name]['time_independent'] = np.zeros(
            dataset.n, dtype=np.double)

        model_x0[dataset_name]['complete'] = np.zeros(n_input, dtype=np.double)

        if 'none' in dataset.models or 'None' in dataset.models:
            continue
        if not dataset.models:
            continue

        logchi2_gp_model = None

        for model_name in dataset.models:

            parameter_values = {}
            for common_ref in mc.models[model_name].common_ref:
                parameter_values.update(
                    mc.common_models[common_ref].convert(theta))

            #TODO: remove try-except starting from version 11 !!
            try:
                for planet_name in mc.models[model_name].multiple_planets:
                    parameter_values.update(
                        mc.common_models[planet_name].convert_with_name(theta, planet_name))
            except:
                pass

            parameter_values.update(
                mc.models[model_name].convert(theta, dataset_name))

            if getattr(mc.models[model_name], 'internal_likelihood', False):
                logchi2_gp_model = model_name
                continue

            if getattr(dataset, 'dynamical', False):
                dataset.external_model = dynamical_output[dataset_name]
                external_model = dynamical_output_x0[dataset_name].copy()
                model_out[dataset_name]['dynamical'] = dynamical_output[dataset_name].copy()
                model_x0[dataset_name]['dynamical'] = dynamical_output_x0[dataset_name].copy()

            model_out[dataset_name][model_name] = mc.models[model_name].compute(
                parameter_values, dataset)

            if getattr(mc.models[model_name], 'time_independent_model', False):
                model_x0[dataset_name][model_name] = np.zeros(
                    np.size(x0_plot), dtype=np.double)
                model_out[dataset_name]['time_independent'] += mc.models[model_name].compute(
                    parameter_values, dataset)
            else:
                model_x0[dataset_name][model_name] = mc.models[model_name].compute(
                    parameter_values, dataset, x0_plot)

            if getattr(mc.models[model_name], 'systematic_model', False):
                continue

            if getattr(mc.models[model_name], 'jitter_model', False):
                continue

            if getattr(mc.models[model_name], 'unitary_model', False):
                dataset.unitary_model += model_out[dataset_name][model_name]
                unitary_model += model_x0[dataset_name][model_name]

                if dataset.normalization_model is None:
                    dataset.normalization_model = np.ones(
                        dataset.n, dtype=np.double)
                    normalization_model = np.ones(
                        np.size(x0_plot), dtype=np.double)

            elif getattr(mc.models[model_name], 'normalization_model', False):
                if dataset.normalization_model is None:
                    dataset.normalization_model = np.ones(
                        dataset.n, dtype=np.double)
                    normalization_model = np.ones(
                        np.size(x0_plot), dtype=np.double)
                dataset.normalization_model *= model_out[dataset_name][model_name]
                normalization_model *= model_x0[dataset_name][model_name]

            else:
                dataset.additive_model += model_out[dataset_name][model_name]
                additive_model += model_x0[dataset_name][model_name]

        dataset.compute_model()
        dataset.compute_residuals()

        model_x0[dataset_name]['complete'] += dataset.compute_model_from_arbitrary_datasets(additive_model,
                                                                                            unitary_model,
                                                                                            normalization_model,
                                                                                            external_model)
        model_out[dataset_name]['complete'] += dataset.model

        """ Gaussian Process check MUST be the last one or the program will fail
         that's because for the GP to work we need to know the _deterministic_ part of the model
         (i.e. the theoretical values you get when you feed your model with the parameter values) """
        if logchi2_gp_model:
            parameter_values = {}
            try:
                for common_ref in mc.models[logchi2_gp_model].common_ref:
                    parameter_values.update(
                        mc.common_models[common_ref].convert(theta))
            except:
                pass

            parameter_values.update(
                mc.models[logchi2_gp_model].convert(theta, dataset.name_ref))

            if hasattr(mc.models[logchi2_gp_model], 'delayed_lnlk_computation'):
                mc.models[logchi2_gp_model].add_internal_dataset(parameter_values, dataset)
                delayed_lnlk_computation[dataset.name_ref] = logchi2_gp_model

            else:

                model_out[dataset_name][logchi2_gp_model] = \
                    mc.models[logchi2_gp_model].sample_predict(parameter_values, dataset)
                model_out[dataset_name]['complete'] += model_out[dataset_name][logchi2_gp_model]

                """ Attempt to avoid RAM overflow
                """
                x0_len = len(x0_plot)

                plot_split_threshold = kwargs.get('plot_split_threshold', 10000)
                low_ram_plot = kwargs.get('low_ram_plot', False)

                #if (x0_len * dataset.n) > ram_occupancy:
                if (x0_len > plot_split_threshold) or low_ram_plot:

                    x0_out = np.empty(x0_len)
                    x0_var = np.empty(x0_len)

                    if dataset.n < plot_split_threshold/10.:

                        array_length = int(x0_len / dataset.n)

                        if low_ram_plot:
                            array_length = dataset.n
                        elif plot_split_threshold < array_length:
                            array_length = plot_split_threshold
                    else:
                        if low_ram_plot:
                            array_length = int(plot_split_threshold / 10.0)
                        else:
                            array_length = plot_split_threshold

                    id_start = 0
                    id_end = array_length

                    max_iterations = x0_len//array_length + 1

                    print('     Splitting the plot array to allow GP prediction of extended datasets, it may take a while...')
                    print('     # {0:d} chunks of {1:d} times each'.format(max_iterations, array_length))
                    print('     Check the documentation if the code is taking too long or if it crashes...')

                    if progress_bar:

                        for i_gp in tqdm(range(max_iterations)):

                            if (id_end+array_length >= x0_len):
                                id_end = -1

                            x0_temp = x0_plot[id_start:id_end]

                            x0_out[id_start:id_end], x0_var[id_start:id_end] = \
                                mc.models[logchi2_gp_model].sample_predict(
                                parameter_values, dataset, x0_temp, return_variance=True)

                            if id_end < 0: break

                            id_start += array_length
                            id_end += array_length
                    else:
                        for i_gp in range(max_iterations):

                            if (id_end+array_length >= x0_len):
                                id_end = -1

                            x0_temp = x0_plot[id_start:id_end]

                            x0_out[id_start:id_end], x0_var[id_start:id_end] = \
                                mc.models[logchi2_gp_model].sample_predict(
                                parameter_values, dataset, x0_temp, return_variance=True)

                            if id_end < 0: break

                            id_start += array_length
                            id_end += array_length

                else:
                    print('     computing GP prediction for the whole temporal range, it may take a while...')

                    x0_out, x0_var = \
                    mc.models[logchi2_gp_model].sample_predict(
                        parameter_values, dataset, x0_plot, return_variance=True)

                model_x0[dataset_name][logchi2_gp_model] = x0_out

                model_x0[dataset_name][logchi2_gp_model +
                                       '_std'] = np.sqrt(x0_var)
                model_x0[dataset_name]['complete'] += model_x0[dataset_name][logchi2_gp_model]

    for dataset_name, logchi2_gp_model in delayed_lnlk_computation.items():

        """ parameter_values is embedded in mc.dataset_dict[dataset_name]
            dictionary
        """
        model_out[dataset_name][logchi2_gp_model] = \
            mc.models[logchi2_gp_model].sample_predict(mc.dataset_dict[dataset_name])

        model_out[dataset_name]['complete'] += model_out[dataset_name][logchi2_gp_model]

        """ Attempt to avoid RAM overflow
        """

        x0_plot = bjd_dict[dataset_name]['x0_plot']
        x0_len = len(x0_plot)

        plot_split_threshold = kwargs.get('plot_split_threshold', 10000)
        low_ram_plot = kwargs.get('low_ram_plot', False)

        #if (6* x0_len * dataset.n * len(delayed_lnlk_computation)**2) > ram_occupancy:
        if (x0_len > plot_split_threshold) or low_ram_plot:

            x0_out = np.empty(x0_len)
            x0_var = np.empty(x0_len)


            if dataset.n < plot_split_threshold/10.:

                array_length = int(x0_len / dataset.n)

                if low_ram_plot:
                    array_length = dataset.n
                elif plot_split_threshold < array_length:
                    array_length = plot_split_threshold
            else:
                if low_ram_plot:
                    array_length = int(plot_split_threshold / 10.0)
                else:
                    array_length = plot_split_threshold

            id_start = 0
            id_end = array_length

            max_iterations = x0_len//array_length + 1

            print('     Splitting the plot array to allow GP prediction of extended datasets, it may take a while...')
            print('     # {0:d} chunks of {1:d} times each'.format(max_iterations, array_length))
            print('     Check the documentation if the code is taking too long or if it crashes...')

            if progress_bar:

                for i_gp in tqdm(range(max_iterations)):

                    #print("{0:4.1f}%".format(100/max_iterations*i_gp), end = '')

                    if (id_end+array_length >= x0_len):
                        id_end = -1

                    x0_temp = x0_plot[id_start:id_end]

                    x0_out[id_start:id_end], x0_var[id_start:id_end] = \
                        mc.models[logchi2_gp_model].sample_predict(
                        mc.dataset_dict[dataset_name], x0_temp, return_variance=True)

                    if id_end < 0: break

                    id_start += array_length
                    id_end += array_length

            else:
                for i_gp in range(max_iterations):

                    #print("{0:4.1f}%".format(100/max_iterations*i_gp), end = '')

                    if (id_end+array_length >= x0_len):
                        id_end = -1

                    x0_temp = x0_plot[id_start:id_end]

                    x0_out[id_start:id_end], x0_var[id_start:id_end] = \
                        mc.models[logchi2_gp_model].sample_predict(
                        mc.dataset_dict[dataset_name], x0_temp, return_variance=True)

                    if id_end < 0: break

                    id_start += array_length
                    id_end += array_length

            #print('   ...Done')
        else:
            print('     computing GP prediction for the whole temporal range, it may take a while...')

            x0_out, x0_var = \
            mc.models[logchi2_gp_model].sample_predict(
                mc.dataset_dict[dataset_name], x0_plot, return_variance=True)



        #model_x0[dataset_name][logchi2_gp_model], var = \
        #    mc.models[logchi2_gp_model].sample_predict(
        #        mc.dataset_dict[dataset_name], x0_plot, return_variance=True)

        model_x0[dataset_name][logchi2_gp_model] = x0_out
        model_x0[dataset_name][logchi2_gp_model + '_std'] = np.sqrt(x0_var)
        model_x0[dataset_name]['complete'] += model_x0[dataset_name][logchi2_gp_model]

    # workaround to avoid memory leaks from GP module
    # gc.collect()

    if np.shape(model_out) == 1:
        return model_out * np.ones(), model_x0

    else:
        return model_out, model_x0


def print_theta_bounds(i_dict, theta, bounds):
    format_string = '{0:12s}  {1:4d}  {2:12f} [{3:10.4f}, {4:10.4f}]'

    for par, i in i_dict.items():
        if len(np.shape(theta)) == 2:

            theta_med = compute_value_sigma(theta[:, i])
            s0, s1 = return_significant_figures(theta_med[0], theta_med[2], theta_med[1])
            format_string_long = '{0:12s}  {1:4d}  {2:12.'+repr(max(s0,s1))+'f} {3:12.'+repr(s0)+'f} {4:12.'+repr(s1)+'f}   (15-84 p)   [{5:10.4f}, {6:10.4f}]'

            print(
                format_string_long.format(par, i, theta_med[0], theta_med[2], theta_med[1], bounds[i, 0], bounds[i, 1]))
        else:
            print(format_string.format(
                par, i, theta[i], bounds[i, 0], bounds[i, 1]))
    print()


def print_analysis_info(i_dict, bounds, spaces, priors, additional_kind, additonal_pams):
    format_string_v1 = '{0:12s}  id:{1:4d}  s:{2:11s} b:[{3:12.4f}, {4:12.4f}]   p:{5:s}  '
    format_string_v2 = '{0:12s}  derived (no id, space, bound) {1:25s} p:{2:s}  '

    for par, i in i_dict.items():
        print(format_string_v1.format(
                par, i, spaces[i], bounds[i, 0], bounds[i, 1],priors[i][0]), priors[i][1])

    for par in additional_kind:
        if par not in i_dict:
            print(format_string_v2.format(par, '', additional_kind[par]), additonal_pams[par])

    print()



def return_significant_figures(perc0, perc1=None, perc2=None, are_percentiles=False):

        if perc1==None and perc2==None:
            if np.isnan(perc0):
                return 0

            abs_value = np.abs(perc0)
            if abs_value > 10.:
                return 6
            elif  abs_value > 1:
                return 6
            elif abs_value == 0:
                return 2
            else:
                return 6
                #x1 = np.log10(abs_value)
                #return int(np.ceil(abs(x1))+1)

        elif perc2==None:
            if np.isnan(perc0) or np.isnan(perc1):
                return 0
            value_err = np.abs(perc1)

            if value_err > 10. or value_err > 10.:
                return 0
            elif  value_err > 1. or value_err > 1.:
                return 1
            else:
                x0 = np.log10(value_err)
                return int(np.ceil(abs(x0))+1)
        else:
            if np.isnan(perc0) or np.isnan(perc1) or np.isnan(perc2):
                return 2, 2

            if are_percentiles:
                minus_err = np.abs(perc0 - perc1)
                plus_err = np.abs(perc2 - perc1)
            else:
                minus_err = np.abs(perc1)
                plus_err = np.abs(perc2)


            if minus_err > 10.:
                sig_minus = 0
            elif minus_err > 1.:
                sig_minus = 1
            else:
                try:
                    x0 = np.log10(minus_err)
                    sig_minus = int(np.ceil(abs(x0))+1)
                except OverflowError:
                    sig_minus = 6

            if plus_err > 10.:
                sig_plus = 0
            elif plus_err > 1.:
                sig_plus = 1
            else:
                try:
                    x0 = np.log10(plus_err)
                    sig_plus = int(np.ceil(abs(x0))+1)
                except OverflowError:
                    sig_plus = 6

            return sig_minus, sig_plus

def print_dictionary(parameter_values, recenter=[], fixed_warning=True):
    format_boundary = 0.00001

    for par_names, par_vals in parameter_values.items():
        if np.size(par_vals) > 1:
            if par_names in recenter:
                move_back = (
                    par_vals > recenter[par_names][0] + recenter[par_names][1] / 2.)
                move_forw = (
                    par_vals < recenter[par_names][0] - recenter[par_names][1] / 2.)
                par_vals_recentered = par_vals.copy()
                par_vals_recentered[move_back] -= recenter[par_names][1]
                par_vals_recentered[move_forw] += recenter[par_names][1]
                perc0, perc1, perc2 = np.percentile(
                    par_vals_recentered, [15.865, 50, 84.135], axis=0)

            else:
                perc0, perc1, perc2 = np.percentile(
                    par_vals, [15.865, 50, 84.135], axis=0)

            s0, s1 = return_significant_figures(perc0, perc1, perc2, are_percentiles=True)

            """
            if np.abs(perc1) < format_boundary or \
                    np.abs(perc0 - perc1) < format_boundary or \
                    np.abs(perc2 - perc1) < format_boundary:#0

                format_string_long_exp = '{0:12s}   {1:15e} {2:12.2e} {3:12.2e} (15-84 p)'

                print(format_string_long_exp.format(
                    par_names, perc1, perc0 - perc1, perc2 - perc1))
            else:

                format_string_long = '{0:12s}   {1:15.'+repr(max(s0,s1))+'f} {2:12.'+repr(s0)+'f} {3:12.'+repr(s1)+'f} (15-84 p)'

                print(format_string_long.format(
                    par_names, perc1, perc0 - perc1, perc2 - perc1))

            """
            format_string_long = '{0:12s}   {1:15.'+repr(max(s0,s1))+'f} {2:12.'+repr(s0)+'f} {3:12.'+repr(s1)+'f}    (15-84 p)'
            print(format_string_long.format(par_names, perc1, perc0 - perc1, perc2 - perc1))

        else:

            try:
                par_temp = par_vals[0]
            except:
                par_temp = par_vals

            s0 = return_significant_figures(par_temp)
            if fixed_warning:
                format_string = '{0:12s}   {1:15.'+repr(s0)+'f}                              [fixed]'
            else:
                format_string = '{0:12s}   {1:15.'+repr(s0)+'f} '

            print(format_string.format(par_names, par_temp))

    print()


def print_integrated_ACF(sampler_chain, theta_dict, nthin):

    from emcee.autocorr import integrated_time, function_1d, AutocorrError, auto_window
    # if not emcee.__version__[0] == '3': return

    swapped_chains = np.swapaxes(sampler_chain, 1, 0)

    print()
    print('Computing the autocorrelation time of the chains')


    try:
        tolerance = 50
        integrated_ACF = integrated_time(
            swapped_chains, tol=tolerance, quiet=False)
        print()
        print('The chains are more than 50 times longer than the ACF, the estimate can be trusted')
    except (AutocorrError):
        tolerance = 20
        print()
        print('***** WARNING ******')
        print('The integrated autocorrelation time cannot be reliably estimated')
        print('likely the chains are too short, and ACF analysis is not fully reliable')
        print('emcee.autocorr.integrated_time tolerance lowered to 20')
        print('If you still get a warning, you should drop these results entirely')
        integrated_ACF = integrated_time(
            swapped_chains, tol=tolerance, quiet=True)
    except (NameError, TypeError):
        print()
        print('Old version of emcee, this function is not implemented')
        print()
        return None, None, None, None

    """ computing the autocorrelation time every 1000 steps, skipping the first 5000 """

    n_sam = swapped_chains.shape[0]
    n_cha = swapped_chains.shape[1]
    n_dim = swapped_chains.shape[2]
    try:
        acf_len = int(np.nanmax(integrated_ACF))  # 1000//nthin
    except:
        print('Error in computing max integrated ACF, skipped ')
        print()
        return None, None, None, None

    if acf_len==0:
        print('Error in computing integrated ACF, chains too short, skipped ')
        print()
        return None, None, None, None

    c = 5

    if n_sam > acf_len*tolerance:

        acf_previous = np.ones(n_dim)
        acf_current = np.ones(n_dim)
        acf_counter = np.zeros(n_dim, dtype=np.int64)
        acf_converged_at = np.zeros(n_dim, dtype=np.int64) - 1

        i_sampler = np.arange(acf_len*3, n_sam, acf_len, dtype=np.int64)

        acf_trace = np.zeros([len(i_sampler), n_dim])
        acf_diff = np.zeros([len(i_sampler), n_dim])

        for id_sam, i_sam in enumerate(i_sampler):

            acf_previous = 1.*acf_current
            integrated_part = np.zeros(i_sam)

            for i_dim in range(0, n_dim):
                integrated_part *= 0.
                for i_cha in range(0, n_cha):
                    integrated_part += function_1d(
                        swapped_chains[:i_sam, i_cha, i_dim])

                c = 5
                integrated_part /= (1.*n_cha)
                taus = 2.0 * np.cumsum(integrated_part) - 1.0
                window = auto_window(taus, c)
                acf_current[i_dim] = taus[window]

            acf_trace[id_sam, :] = acf_current * 1.
            acf_diff[id_sam, :] = np.abs(acf_current-acf_previous)/acf_current

            converged = (i_sam > integrated_ACF*tolerance) \
                & (acf_diff[id_sam, :] < 0.01) \
                & (acf_converged_at < 0.)

            acf_counter[converged] += 1

            converged_stable = (converged) & (acf_counter == 5)

            acf_converged_at[converged_stable] = i_sam * nthin

            #if np.sum((acf_converged_at > 0), dtype=np.int64) == n_dim:
            #    break


        acf_converged_flag = (acf_converged_at>0)
        how_many_ACT = (n_sam - acf_converged_at/nthin)/integrated_ACF
        how_many_ACT[(acf_converged_at < 0)] = -1

        still_required_050 = (50-how_many_ACT)*(nthin*integrated_ACF)
        still_required_050_flag = (still_required_050 > 0)
        still_required_050[~still_required_050_flag] = 0
        still_required_100 = (100-how_many_ACT)*(nthin*integrated_ACF)
        still_required_100_flag = (still_required_100 > 0)
        still_required_100[~still_required_100_flag] = 0

        print()
        print('Reference thinning used in the analysis:', nthin)
        print(
            'Step length used in the analysis: {0:d}*nthin = {1:d}'.format(acf_len, acf_len*nthin))
        print()
        print('Convergence criteria: less than 1% variation in ACF after {0:d} times the integrated ACF'.format(
            tolerance))
        print('At least 50*ACF after convergence, 100*ACF would be ideal')
        print('Negative values: not converged yet')
        print()
        print('   sampler parameter              |    ACF   | ACF*nthin | converged at | nsteps/ACF | to 100*ACF')
        print('                                  |          |           |              |            | ')

        for key_name, key_val in theta_dict.items():
            print('   {0:30s} | {1:7.2f}  | {2:8.0f}  |  {3:7.0f}     |  {4:6.0f}    | {5:7.0f} '.format(key_name,
                                                                                                   integrated_ACF[key_val],
                                                                                                    integrated_ACF[key_val] *
                                                                                                    nthin,
                                                                                                    acf_converged_at[key_val],
                                                                                                    how_many_ACT[key_val],
                                                                                                    still_required_100[key_val]))

        print()

        if np.sum((acf_converged_at > 0), dtype=np.int64) == n_dim:
            if np.sum(how_many_ACT > 100, dtype=np.int64) == n_dim:
                print('All the chains are longer than 100*ACF ')
            elif (np.sum(how_many_ACT > 50, dtype=np.int64) == n_dim):
                print("""All the chains are longer than 50*ACF, but some are shorter than 100*ACF
PyORBIT should keep running for about {0:.0f} more steps to reach 100*ACF""".format(np.average(still_required_100[still_required_100_flag])))
            else:
                print("""All the chains have converged, but PyORBIT should keep running for about:
{0:.0f} more steps to reach 50*ACF,
{1:.0f} more steps to reach 100*ACF""".format(np.average(still_required_050[still_required_050_flag]), np.average(still_required_100[still_required_100_flag])))

            print('Suggested value for burnin: {0:.0f}'.format(np.average(acf_converged_at)))

        else:
            print(' {0:5.0f} chains have not converged yet, keep going '.format(
                np.sum((acf_converged_at < 0), dtype=np.int64)))

        print()

        return i_sampler, acf_trace, acf_diff, acf_converged_at

    else:
        print("Chains too short to apply convergence criteria")
        print(
            "They should be at least {0:d}*nthin = {1:d}".format(50*acf_len, 50*acf_len*nthin))
        print()

        return None, None, None, None
