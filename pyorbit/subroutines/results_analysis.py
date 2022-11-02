from __future__ import print_function
from pyorbit.subroutines.common import *
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo

from tqdm import tqdm

__all__ = ["results_resumen", "results_derived", "get_planet_variables", "get_theta_dictionary", "get_model",
           "print_theta_bounds", "print_dictionary", "get_stellar_parameters", "print_integrated_ACF"]


def results_resumen(mc, theta,
                    skip_theta=False,
                    compute_lnprob=False,
                    chain_med=False,
                    return_samples=False,
                    is_starting_point=False):
    # Function with two goals:
    # * Unfold and print out the output from theta
    # * give back a parameter name associated to each value in the result array

    print()
    print('====================================================================================================')
    if skip_theta:
        print('     Boundaries of the sampler variables     ')
    elif is_starting_point:
        print('     Starting point of the sample/optimization routines    ')
    else:
        print('     Statistics on the posterior of the sampler variables     ')

    print('====================================================================================================')
    print()
    for dataset_name, dataset in mc.dataset_dict.items():
        print('----- dataset: ', dataset_name)
        print_theta_bounds(dataset.variable_sampler,
                           theta, mc.bounds, skip_theta)

        for model_name in dataset.models:
            print('---------- ', dataset_name,
                  '     ----- model: ', model_name)
            print_theta_bounds(
                mc.models[model_name].variable_sampler[dataset_name], theta, mc.bounds, skip_theta)

    for model_name, model in mc.common_models.items():
        print('----- common model: ', model_name)
        print_theta_bounds(model.variable_sampler,
                           theta, mc.bounds, skip_theta)

    if skip_theta:
        return

    print('====================================================================================================')
    if is_starting_point:
        print('     Starting point projected onto the physical space     ')
    else:
        print('     Statistics on the physical parameters obtained from the posteriors samples     ')
    print('====================================================================================================')
    print()

    for dataset_name, dataset in mc.dataset_dict.items():
        print('----- dataset: ', dataset_name)
        variable_values = dataset.convert(theta)
        print_dictionary(variable_values)

        print()
        for model_name in dataset.models:
            print('---------- ', dataset_name,
                  '     ----- model: ', model_name)
            variable_values = mc.models[model_name].convert(
                theta, dataset_name)
            print_dictionary(variable_values)

    for model_name, model in mc.common_models.items():
        print('----- common model: ', model_name)
        variable_values = model.convert(theta)
        if chain_med is not False:
            recenter_pams = {}
            variable_values_med = model.convert(chain_med)

            # for var in list(set(mc.recenter_pams_dataset) & set(mc.variable_sampler[dataset_name])):
            for var in list(set(model.recenter_pams) & set(variable_values_med)):
                recenter_pams[var] = [variable_values_med[var],
                                      model.default_bounds[var][1] - model.default_bounds[var][0]]
            print_dictionary(variable_values, recenter=recenter_pams)

        else:
            print_dictionary(variable_values)

    print('====================================================================================================')
    if is_starting_point:
        print('     Derived parameters obtained from starting point     ')
    else:
        print('     Statistics on the derived parameters obtained from the posteriors samples     ')
    print('====================================================================================================')
    print()

    returned_samples = get_planet_variables(mc, theta, verbose=True)

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
        print(' Trying to identify the stellar parameters')
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
                print('Note: stellar density derived from its mass and radius')
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
                print('Note: stellar radius derived from its mass and density')
                print()
        elif 'radius' in stellar_values:
            stellar_values['mass'] = stellar_values['radius'] ** 3. * stellar_values['density']
            if warnings:
                print('Note: stellar mass derived from its radius and density')
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
    _ = get_planet_variables(mc, theta, verbose=True)


def get_planet_variables(mc, theta, verbose=False):
    """
    Derived parameters from the Common Models are listed

    :param mc:
    :return:
    """

    try:
        n_samplings, n_pams = np.shape(theta)
    except:
        n_samplings = 1

    planet_variables = {}

    for common_name, common_model in mc.common_models.items():
        variable_values = common_model.convert(theta)
        derived_variables = {}

        if common_model.model_class == 'planet':

            stellar_ref = getattr(common_model, 'stellar_ref', None)

            stellar_values = get_stellar_parameters(mc,
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

            for var in variable_values.keys():
                if np.size(variable_values[var]) == 1:
                    variable_values[var] = variable_values[var] * \
                        np.ones(n_samplings)

            if 'a_Rs' not in variable_values.keys() and 'density' in stellar_values.keys():
                derived_variables['a_Rs'] = True
                variable_values['a_Rs'] = convert_rho_to_a(variable_values['P'],
                                                        stellar_values['density'])

            if 'i' not in variable_values.keys():
                derived_variables['i'] = True
                if 'i' in common_model.fix_list:

                    if verbose:
                        print('Inclination randomized to {0:3.2f} +- {1:3.2f} deg'.format(
                              common_model.fix_list['i'][0], common_model.fix_list['i'][1]))
                    variable_values['i'] = np.random.normal(common_model.fix_list['i'][0],
                                                            common_model.fix_list['i'][1],
                                                            size=n_samplings)
                elif 'b' in variable_values.keys() and 'a_Rs' in variable_values.keys():
                    variable_values['i'] = convert_b_to_i(variable_values['b'],
                                                          variable_values['e'],
                                                          variable_values['omega'],
                                                          variable_values['a_Rs'])
                else:
                    print('Inclination fixed to 90 deg!')
                    variable_values['i'] = 90.00 * np.ones(n_samplings)
                    remove_i = True

            if 'K' in variable_values.keys() and 'mass' in stellar_values.keys():
                M_Msun = kepler_exo.get_planet_mass(variable_values['P'],
                                                                  variable_values['K'],
                                                                  variable_values['e'],
                                                                  stellar_values['mass']) \
                    / np.sin(np.radians(variable_values['i']))

                variable_values['M_Mj'] =M_Msun * constants.Msjup
                derived_variables['M_Mj'] = True

                variable_values['M_Me']  =M_Msun * constants.Msear
                derived_variables['M_Me'] = True

            elif 'M_Me' in variable_values.keys() and 'mass' in stellar_values.keys():
                derived_variables['K'] = True
                derived_variables['K'] = kepler_exo.kepler_K1(stellar_values['mass'],
                                                              variable_values['M_Me'] /
                                                              constants.Msear,
                                                              variable_values['P'],
                                                              variable_values['i'],
                                                              variable_values['e'])

                variable_values['M_Mj'] = variable_values['M_Me'] * \
                    (constants.Msjup/constants.Msear)
                derived_variables['M_Mj'] = True


            if 'Tc' in variable_values.keys():
                if 'e' in variable_values:
                    derived_variables['mean_long'] = True
                    variable_values['mean_long'] = kepler_exo.kepler_Tc2phase_Tref(variable_values['P'],
                                                                           variable_values['Tc'] -
                                                                           mc.Tref,
                                                                           variable_values['e'],
                                                                           variable_values['omega'])

            elif 'mean_long' in variable_values.keys():
                derived_variables['Tc'] = True
                variable_values['Tc'] = mc.Tref + kepler_exo.kepler_phase2Tc_Tref(variable_values['P'],
                                                                                  variable_values['mean_long'],
                                                                                  variable_values['e'],
                                                                                  variable_values['omega'])

            if 'R_Rs' in variable_values.keys() and 'radius' in stellar_values.keys():
                variable_values['R_Rj'] = variable_values['R_Rs'] * \
                    constants.Rsjup * stellar_values['radius']
                derived_variables['R_Rj'] = True

                variable_values['R_Re'] = variable_values['R_Rs'] * \
                    constants.Rsear * stellar_values['radius']
                derived_variables['R_Re'] = True

            if remove_i:
                del variable_values['i']

            try:
                k = variable_values['R_Rs']

                variable_values['T_41'] = variable_values['P'] / np.pi \
                    * np.arcsin(1./variable_values['a_Rs'] *
                                np.sqrt((1. + k)**2 - variable_values['b']**2)
                                / np.sin(variable_values['i']*constants.deg2rad))
                derived_variables['T_41'] = True

                variable_values['T_32'] = variable_values['P'] / np.pi \
                    * np.arcsin(1./variable_values['a_Rs'] *
                                np.sqrt((1. - k)**2 - variable_values['b']**2)
                                / np.sin(variable_values['i']*constants.deg2rad))
                derived_variables['T_32'] = True

            except:
                pass

            try:
                derived_variables['a_AU_(M)'] = True
                variable_values['a_AU_(M)'] = convert_PMsMp_to_a(
                    variable_values['P'],
                    stellar_values['mass'],
                    variable_values['M_Me'])
            except (KeyError, ValueError):
                pass

            try:
                derived_variables['a_AU_(rho,R)'] = True
                variable_values['a_AU_(rho,R)'] = convert_ars_to_a(
                    variable_values['a_Rs'],
                    stellar_values['radius'])
            except (KeyError, ValueError):
                pass

            try:
                derived_variables['insol'] = True
                variable_values['insol]'] = \
                    convert_RTaAU_to_insol(
                        stellar_values['radius'],
                        stellar_values['temperature'],
                        variable_values['a_AU_(rho,R)']
                    )
            except (KeyError, ValueError):
                pass

            planet_variables[common_name] = variable_values.copy()

            for var in planet_variables[common_name].keys():
                if var not in derived_variables.keys():
                    del variable_values[var]

            if verbose:
                print_dictionary(variable_values)

    return planet_variables


def get_theta_dictionary(mc):
    # * give back a parameter name associated to each value in the result array

    theta_dictionary = {}
    for dataset_name, dataset in mc.dataset_dict.items():
        for var, i in dataset.variable_sampler.items():
            try:
                theta_dictionary[dataset_name + '_' + var] = i
            except:
                theta_dictionary[repr(dataset_name) + '_' + var] = i

        for model_name in dataset.models:
            for var, i in mc.models[model_name].variable_sampler[dataset_name].items():
                try:
                    theta_dictionary[dataset_name +
                                     '_' + model_name + '_' + var] = i
                except:
                    theta_dictionary[repr(dataset_name) +
                                     '_' + model_name + '_' + var] = i

    for model_name, model in mc.common_models.items():
        for var, i in model.variable_sampler.items():
            theta_dictionary[model.common_ref + '_' + var] = i

    return theta_dictionary


def get_model(mc, theta, bjd_dict):
    model_out = {}
    model_x0 = {}

    delayed_lnlk_computation = {}

    if mc.dynamical_model is not None:
        """ check if any keyword ahas get the output model from the dynamical tool
        we must do it here because all the planet are involved"""
        dynamical_output_x0 = mc.dynamical_model.compute(
            mc, theta, bjd_dict['full']['x_plot'])
        dynamical_output = mc.dynamical_model.compute(mc, theta)

    for dataset_name, dataset in mc.dataset_dict.items():

        x0_plot = bjd_dict[dataset_name]['x0_plot']

        n_input = np.size(x0_plot)
        model_out[dataset_name] = {}
        model_x0[dataset_name] = {}
        dataset.model_reset()

        additive_model = np.zeros(np.size(x0_plot))
        unitary_model = np.zeros(np.size(x0_plot))
        external_model = np.zeros(np.size(x0_plot))
        normalization_model = None

        variable_values = dataset.convert(theta)
        dataset.compute(variable_values)

        for model_name in dataset.models:
            variable_values = {}

            for common_ref in mc.models[model_name].common_ref:
                variable_values.update(
                    mc.common_models[common_ref].convert(theta))

            # try:
            #    for common_ref in mc.models[model_name].common_ref:
            #        variable_values.update(mc.common_models[common_ref].convert(theta))
            # except:
            #    continue
            variable_values.update(
                mc.models[model_name].convert(theta, dataset_name))

            if getattr(mc.models[model_name], 'jitter_model', False):
                dataset.jitter += mc.models[model_name].compute(
                    variable_values, dataset)
                continue

            if getattr(mc.models[model_name], 'systematic_model', False):
                dataset.additive_model += mc.models[model_name].compute(
                    variable_values, dataset)

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

            variable_values = {}
            for common_ref in mc.models[model_name].common_ref:
                variable_values.update(
                    mc.common_models[common_ref].convert(theta))
            variable_values.update(
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
                variable_values, dataset)

            if getattr(mc.models[model_name], 'time_independent_model', False):
                model_x0[dataset_name][model_name] = np.zeros(
                    np.size(x0_plot), dtype=np.double)
                model_out[dataset_name]['time_independent'] += mc.models[model_name].compute(
                    variable_values, dataset)
            else:
                model_x0[dataset_name][model_name] = mc.models[model_name].compute(
                    variable_values, dataset, x0_plot)

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
            variable_values = {}
            try:
                for common_ref in mc.models[logchi2_gp_model].common_ref:
                    variable_values.update(
                        mc.common_models[common_ref].convert(theta))
            except:
                pass

            variable_values.update(
                mc.models[logchi2_gp_model].convert(theta, dataset.name_ref))

            if hasattr(mc.models[logchi2_gp_model], 'delayed_lnlk_computation'):
                mc.models[logchi2_gp_model].add_internal_dataset(variable_values, dataset)
                delayed_lnlk_computation[dataset.name_ref] = logchi2_gp_model

            else:

                model_out[dataset_name][logchi2_gp_model] = \
                    mc.models[logchi2_gp_model].sample_predict(variable_values, dataset)
                    #mc.models[logchi2_gp_model].sample_conditional(variable_values, dataset)
                model_out[dataset_name]['complete'] += model_out[dataset_name][logchi2_gp_model]

                """ Attempt to avoid RAM overflow
                    A float64 takes 64 bits -> 8 kilobytes
                    here we set it to take 0.5 Gb
                """
                ram_occupancy= (0.5 *  1024**3) / 8
                x0_len = len(x0_plot)
                if (x0_len * dataset.n) > ram_occupancy:

                    print('     Splitting the plot array to allow GP prediction of extended datasets, it may take a while...')

                    x0_out = np.empty(x0_len)
                    x0_var = np.empty(x0_len)

                    array_length = int(x0_len / dataset.n)

                    id_start = 0
                    id_end = array_length

                    max_iterations = x0_len//array_length + 1
                    for i_gp in tqdm(range(max_iterations+1)):

                        #print("{0:4.1f}%".format(100/max_iterations*i_gp), end = '')

                        if (id_end+array_length >= x0_len):
                            id_end = -1

                        x0_temp = x0_plot[id_start:id_end]

                        x0_out[id_start:id_end], x0_var[id_start:id_end] = \
                            mc.models[logchi2_gp_model].sample_predict(
                            variable_values, dataset, x0_temp, return_variance=True)

                        if id_end < 0: break

                        id_start += array_length
                        id_end += array_length

                    #print('   ...Done')
                else:
                    print('     computing GP prediction for the whole temporal range, it may take a while...')

                    x0_out, x0_var = \
                    mc.models[logchi2_gp_model].sample_predict(
                        variable_values, dataset, x0_plot, return_variance=True)

                model_x0[dataset_name][logchi2_gp_model] = x0_out

                model_x0[dataset_name][logchi2_gp_model +
                                       '_std'] = np.sqrt(x0_var)
                model_x0[dataset_name]['complete'] += model_x0[dataset_name][logchi2_gp_model]

    for dataset_name, logchi2_gp_model in delayed_lnlk_computation.items():

        """ variable_values is embedded in mc.dataset_dict[dataset_name]
            dictionary
        """
        model_out[dataset_name][logchi2_gp_model] = \
            mc.models[logchi2_gp_model].sample_predict(mc.dataset_dict[dataset_name])

        model_out[dataset_name]['complete'] += model_out[dataset_name][logchi2_gp_model]

        """ Attempt to avoid RAM overflow
            A float64 takes 64 bits -> 8 kilobytes
            here we set it to take 0.5 Gb
        """
        ram_occupancy= (0.5 *  1024**3) / 8

        x0_plot = bjd_dict[dataset_name]['x0_plot']
        x0_len = len(x0_plot)
        #print(ram_occupancy,x0_len, dataset.n )
        #print((6* x0_len * dataset.n * len(delayed_lnlk_computation)**2))
        if (6* x0_len * dataset.n * len(delayed_lnlk_computation)**2) > ram_occupancy:

            print('     Splitting the plot array to allow GP prediction of extended datasets, it may take a while...')

            x0_out = np.empty(x0_len)
            x0_var = np.empty(x0_len)

            array_length = int(x0_len / dataset.n)
            #print(array_length)

            id_start = 0
            id_end = array_length

            max_iterations = x0_len//array_length + 1
            #print(max_iterations)
            #print()
            for i_gp in tqdm(range(max_iterations+1)):

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


def print_theta_bounds(i_dict, theta, bounds, skip_theta=False):
    format_string = '{0:12s}  {1:4d}  {2:12f} ([{3:10f}, {4:10f}])'
    format_string_long = '{0:12s}  {1:4d}  {2:12f}   {3:12f}  {4:12f} (15-84 p) ([{5:9f}, {6:9f}])'
    format_string_notheta = '{0:12s}  {1:4d}  ([{2:10f}, {3:10f}])'

    for var, i in i_dict.items():

        if skip_theta:
            print(format_string_notheta.format(
                var, i, bounds[i, 0], bounds[i, 1]))
        elif len(np.shape(theta)) == 2:

            theta_med = compute_value_sigma(theta[:, i])
            print(
                format_string_long.format(var, i, theta_med[0], theta_med[2], theta_med[1], bounds[i, 0], bounds[i, 1]))
        else:
            print(format_string.format(
                var, i, theta[i], bounds[i, 0], bounds[i, 1]))
    print()


def print_dictionary(variable_values, recenter=[]):
    format_string = '{0:12s}   {1:15f} '
    format_string_long = '{0:12s}   {1:15f}   {2:15f}  {3:15f} (15-84 p)'
    format_string_exp = '{0:12s}   {1:15e} '
    format_string_long_exp = '{0:12s}   {1:15e}   {2:15e}  {3:15e} (15-84 p)'
    format_boundary = 0.000010

    for var_names, var_vals in variable_values.items():
        if np.size(var_vals) > 1:
            if var_names in recenter:
                move_back = (
                    var_vals > recenter[var_names][0] + recenter[var_names][1] / 2.)
                move_forw = (
                    var_vals < recenter[var_names][0] - recenter[var_names][1] / 2.)
                var_vals_recentered = var_vals.copy()
                var_vals_recentered[move_back] -= recenter[var_names][1]
                var_vals_recentered[move_forw] += recenter[var_names][1]
                perc0, perc1, perc2 = np.percentile(
                    var_vals_recentered, [15.865, 50, 84.135], axis=0)

            else:
                perc0, perc1, perc2 = np.percentile(
                    var_vals, [15.865, 50, 84.135], axis=0)

            if np.abs(perc1) < format_boundary or \
                    np.abs(perc0 - perc1) < format_boundary or \
                    np.abs(perc2 - perc1) < format_boundary:
                print(format_string_long_exp.format(
                    var_names, perc1, perc0 - perc1, perc2 - perc1))
            else:
                print(format_string_long.format(
                    var_names, perc1, perc0 - perc1, perc2 - perc1))

        else:
            try:
                if np.abs(var_vals[0]) < format_boundary:
                    print(format_string_exp.format(var_names, var_vals[0]))
                else:
                    print(format_string.format(var_names, var_vals[0]))

            except:
                if np.abs(var_vals) < format_boundary:
                    print(format_string_exp.format(var_names, var_vals))
                else:
                    print(format_string.format(var_names, var_vals))

    print()


def print_integrated_ACF(sampler_chain, theta_dict, nthin):

    from emcee.autocorr import integrated_time, function_1d, AutocorrError, auto_window
    # if not emcee.__version__[0] == '3': return

    swapped_chains = np.swapaxes(sampler_chain, 1, 0)

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
        return

    """ computing the autocorrelation time every 1000 steps, skipping the first 5000 """

    n_sam = swapped_chains.shape[0]
    n_cha = swapped_chains.shape[1]
    n_dim = swapped_chains.shape[2]
    try:
        acf_len = int(np.nanmax(integrated_ACF))  # 1000//nthin
    except:
        print('Error in computing max integrated ACF, skipped ')
        print()
        return

    if acf_len==0:
        print('Error in computing integrated ACF, chains too short, skipped ')
        print()
        return

    c = 5

    if n_sam > acf_len*tolerance:

        acf_previous = np.ones(n_dim)
        acf_current = np.ones(n_dim)
        acf_converged_at = np.zeros(n_dim, dtype=np.int16) - 1

        for i_sam in range(acf_len*3, n_sam, acf_len):

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

            acf_current[acf_current < 0.1] = 0.1
            sel = (i_sam > integrated_ACF*tolerance) & (np.abs(acf_current -
                                                               acf_previous)/acf_current < 0.01) & (acf_converged_at < 0.)
            acf_converged_at[sel] = i_sam * nthin

            if np.sum((acf_converged_at > 0), dtype=np.int16) == n_dim:
                break

        how_many_ACT = (n_sam - acf_converged_at/nthin)/integrated_ACF
        how_many_ACT[(acf_converged_at < 0)] = -1

        still_required_050 = (50-how_many_ACT)*(nthin*integrated_ACF)
        still_required_100 = (100-how_many_ACT)*(nthin*integrated_ACF)

        print()
        print('Computing the autocorrelation time of the chains')
        print('Reference thinning used in the analysis:', nthin)
        print(
            'Step length used in the analysis: {0:d}*nthin = {1:d}'.format(acf_len, acf_len*nthin))
        print()
        print('Convergence criteria: less than 1% variation in ACF after {0:d} times the integrated ACF'.format(
            tolerance))
        print('At least 50*ACF after convergence, 100*ACF would be ideal')
        print('Negative values: not converged yet')
        print()
        print('          sample variable      |    ACF   | ACF*nthin | converged at | nteps/ACF ')
        print('                               |          |           |              |           ')

        for key_name, key_val in theta_dict.items():
            print('          {0:20s} | {1:7.3f}  | {2:8.1f}  |  {3:7.0f}     |  {4:7.1f}   '.format(key_name,
                                                                                                    integrated_ACF[key_val],
                                                                                                    integrated_ACF[key_val] *
                                                                                                    nthin,
                                                                                                    acf_converged_at[key_val],
                                                                                                    how_many_ACT[key_val]))

        print()

        if np.sum((acf_converged_at > 0), dtype=np.int16) == n_dim:
            if np.sum(how_many_ACT > 100, dtype=np.int16) == n_dim:
                print('All the chains are longer than 100*ACF ')
            elif (np.sum(how_many_ACT > 50, dtype=np.int16) == n_dim):
                print("""All the chains are longer than 50*ACF, but some are shorter than 100*ACF
PyORBIT should keep running for at least {0:9.0f} more steps to reach 100*ACF""".format(np.amax(still_required_100)))
            else:
                print("""All the chains have converged, but PyORBIT should keep running for at least:
{0:9.0f} more steps to reach 50*ACF
{1:9.0f} more steps to reach 100*ACF""".format(np.amax(still_required_050), np.amax(still_required_100)))

            print('Suggested value for burnin: ', np.amax(acf_converged_at))

        else:
            print(' {0:5.0f} chains have not converged yet, keep going '.format(
                np.sum((acf_converged_at < 0), dtype=np.int16)))

        print()

    else:
        print("Chains too shoort to apply convergence criteria")
        print(
            "They should be at least {0:d}*nthin = {1:d}".format(50*acf_len, 50*acf_len*nthin))
        print()
