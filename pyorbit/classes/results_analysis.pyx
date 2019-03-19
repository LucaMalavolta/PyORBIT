from __future__ import print_function
from common import *

__all__ = ["results_resumen", "results_derived", "get_planet_variables", "get_theta_dictionary", "get_model",
           "print_theta_bounds", "print_dictionary"]


def results_resumen(mc, theta, skip_theta=False, compute_lnprob=False, chain_med=False):
    # Function with two goals:
    # * Unfold and print out the output from theta
    # * give back a parameter name associated to each value in the result array

    print()
    print('====================================================================================================')
    if skip_theta:
        print('     Boundaries of the sampler variables     ')
    else:
        print('     Statistics on the posterior of the sampler variables     ')

    print('====================================================================================================')
    print()
    for dataset_name, dataset in mc.dataset_dict.iteritems():
        print('----- dataset: ', dataset_name)
        print_theta_bounds(dataset.variable_sampler, theta, mc.bounds, skip_theta)

        for model_name in dataset.models:
            print('---------- ', dataset_name, '     ----- model: ', model_name)
            print_theta_bounds(mc.models[model_name].variable_sampler[dataset_name], theta, mc.bounds, skip_theta)

    for model in mc.common_models.itervalues():
        print('----- common model: ', model.common_ref)
        print_theta_bounds(model.variable_sampler, theta, mc.bounds, skip_theta)

    if skip_theta:
        return

    print('====================================================================================================')
    print('     Statistics on the physical parameters obtained from the posteriors samples     ')
    print('====================================================================================================')
    print()

    for dataset_name, dataset in mc.dataset_dict.iteritems():
        print('----- dataset: ', dataset_name)
        variable_values = dataset.convert(theta)
        print_dictionary(variable_values)

        print()
        for model_name in dataset.models:
            print('---------- ', dataset_name, '     ----- model: ', model_name)
            variable_values = mc.models[model_name].convert(theta, dataset_name)
            print_dictionary(variable_values)

    for model in mc.common_models.itervalues():
        print('----- common model: ', model.common_ref)
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
    print('     Statistics on the derived parameters obtained from the posteriors samples     ')
    print('====================================================================================================')
    print()

    _ = get_planet_variables(mc, theta, verbose=True)

    if compute_lnprob:
        print()
        print('====================================================================================================')
        print('     Statistics on the log-likelihood     ')
        print('====================================================================================================')
        print()

        if len(np.shape(theta)) == 2:
            n_samples, n_values = np.shape(theta)
            logchi2_collection = np.zeros(n_samples)
            for i in xrange(0, n_samples):
                logchi2_collection[i] = mc(theta[i, :])
            perc0, perc1, perc2 = np.percentile(logchi2_collection, [15.865, 50, 84.135], axis=0)
            print(' LN probability: %12f   %12f %12f (15-84 p) ' % (perc1, perc0 - perc1, perc2 - perc1))
        else:
            print(' LN probability: %12f ' % (mc(theta)))

    print()
    print('====================================================================================================')
    print('     ------------------------------------------------------------------------------------------     ')
    print('====================================================================================================')
    print()
    print()


def get_stellar_parameters(mc, theta):
    try:
        n_samplings, n_pams = np.shape(theta)
    except:
        n_samplings = 1
        n_pams = np.shape(theta)

    "Stellar mass, radius and density are pre-oaded since they are are required by most of the common models"
    stellar_model = mc.common_models['star_parameters']
    stellar_values = stellar_model.convert(theta)
    stellar_values['provided'] = True

    if 'rho' not in stellar_values:

        if 'radius' not in stellar_values:
            try:
                if stellar_model.prior_kind['radius'] == 'Gaussian':
                    stellar_values['radius'] = np.random.normal(stellar_model.prior_pams['radius'][0],
                                                                stellar_model.prior_pams['radius'][1],
                                                                size=n_samplings)
            except:
                print(' *** Please provide a prior on stellar Radius *** ')
                stellar_values['provided'] = False

        if 'mass' not in stellar_values:
            try:
                if stellar_model.prior_kind['mass'] == 'Gaussian':
                    stellar_values['mass'] = np.random.normal(stellar_model.prior_pams['mass'][0],
                                                              stellar_model.prior_pams['mass'][1],
                                                              size=n_samplings)
            except:
                print(' *** Please provide a prior on stellar Mass *** ')
                stellar_values['provided'] = False

        stellar_values['rho'] = stellar_values['mass'] / stellar_values['radius'] ** 3

    else:
        if 'mass' in stellar_values:
            stellar_values['radius'] = (stellar_values['mass'] / stellar_values['rho']) ** (1. / 3.)
        elif 'radius' in stellar_values:
            stellar_values['mass'] = stellar_values['radius'] ** 3. * stellar_values['rho']
        else:
            if 'mass' in stellar_model.prior_pams:
                if stellar_model.prior_kind['mass'] == 'Gaussian':
                    stellar_values['mass'] = np.random.normal(stellar_model.prior_pams['mass'][0],
                                                              stellar_model.prior_pams['mass'][1],
                                                              size=n_samplings)
                    stellar_values['radius'] = (stellar_values['mass'] / stellar_values['rho']) ** (1. / 3.)

            elif 'radius' in stellar_model.prior_pams:

                if stellar_model.prior_kind['radius'] == 'Gaussian':
                    stellar_values['radius'] = np.random.normal(stellar_model.prior_pams['mass'][0],
                                                                stellar_model.prior_pams['mass'][1],
                                                                size=n_samplings)
                    stellar_values['mass'] = stellar_values['radius'] ** 3. * stellar_values['rho']
            else:
                print(' *** Please provide a prior either on stellar Mass or stellar Radius *** ')
                stellar_values['provided'] = False


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
    stellar_values = get_stellar_parameters(mc, theta)

    planet_variables = {}

    for common_name, common_model in mc.common_models.iteritems():
        variable_values = common_model.convert(theta)
        derived_variables = {}

        if common_model.model_class == 'planet':

            remove_i = False
            if verbose:
                print('----- common model: ', common_name)

            """
            Check if the eccentricity and argument of pericenter were set as free parameters or fixed by simply
            checking the size of their distribution
            """

            var_index = np.arange(0, n_samplings, dtype=int)
            for var in variable_values.iterkeys():
                if np.size(variable_values[var]) == 1:
                    variable_values[var] = variable_values[var] * np.ones(n_samplings)

            if 'a' not in variable_values.keys():
                derived_variables['a'] = True
                variable_values['a'] = convert_rho_to_a(variable_values['P'],
                                                        stellar_values['rho'])

            if 'i' not in variable_values.keys():
                derived_variables['i'] = True
                if 'i' in common_model.fix_list:
                    if verbose:
                        print('Inclination fixed to ', common_model.fix_list['i'][0])
                    variable_values['i'] = np.random.normal(common_model.fix_list['i'][0],
                                                            common_model.fix_list['i'][1],
                                                            size=n_samplings)
                elif 'b' in variable_values.keys():
                    variable_values['i'] = convert_b_to_i(variable_values['b'],
                                                          variable_values['e'],
                                                          variable_values['o'],
                                                          variable_values['a'])
                else:
                    variable_values['i'] = 90.00 * np.ones(n_samplings)
                    remove_i = True

            if 'K' in variable_values.keys():
                variable_values['M'] = np.empty(n_samplings)
                derived_variables['M'] = True

                for P, K, e, i, star, ii in zip(
                        variable_values['P'],
                        variable_values['K'],
                        variable_values['e'],
                        variable_values['i'],
                        stellar_values['mass'],
                        var_index):
                    variable_values['M'][ii] = kepler_exo.get_planet_mass(P, K, e, star, Minit=0.0065) \
                                               / np.sin(np.radians(i))

            elif 'M' in variable_values.keys():
                variable_values['K'] = np.empty(n_samplings)
                derived_variables['K'] = True
                for star, M, P, i, e, ii in zip(stellar_values['mass'],
                                                variable_values['M'],
                                                variable_values['P'],
                                                variable_values['i'],
                                                variable_values['e']):
                    derived_variables['K'][ii] = kepler_exo.kepler_K1(star, M / constants.Msear, P, i, e)

            if 'Tc' in variable_values.keys():
                variable_values['f'] = np.empty(n_samplings)
                derived_variables['f'] = True

                for P, Tc, e, o, ii in zip(
                        variable_values['P'],
                        variable_values['Tc'],
                        variable_values['e'],
                        variable_values['o'],
                        var_index):
                    variable_values['f'][ii] = kepler_exo.kepler_Tc2phase_Tref(P, Tc - mc.Tref, e, o)

            elif 'f' in variable_values.keys():
                variable_values['Tc'] = np.empty(n_samplings)
                derived_variables['Tc'] = True

                for P, f, e, o, ii in zip(
                        variable_values['P'],
                        variable_values['f'],
                        variable_values['e'],
                        variable_values['o'],
                        var_index):
                    variable_values['Tc'][ii] = mc.Tref + kepler_exo.kepler_phase2Tc_Tref(P, f, e, o)

            if 'R' in variable_values.keys() and stellar_values['provided']:
                variable_values['R_Rj'] = variable_values['R'] * constants.Rsjup * stellar_values['radius']
                derived_variables['R_Rj'] = True

                variable_values['R_Re'] = variable_values['R'] * constants.Rsear * stellar_values['radius']
                derived_variables['R_Re'] = True

            if 'M' in variable_values.keys() and stellar_values['provided']:
                variable_values['M_Mj'] = variable_values['M'] * constants.Msjup * stellar_values['mass']
                derived_variables['M_Mj'] = True

                variable_values['M_Me'] = variable_values['M'] * constants.Msear * stellar_values['mass']
                derived_variables['M_Me'] = True

            if remove_i:
                del variable_values['i']

            if 'b' in variable_values.keys() and stellar_values['provided']:

                k = variable_values['R']

                variable_values['T_41'] =  variable_values['P'] / np.pi \
                                          * np.arcsin(1./variable_values['a'] *
                                                      np.sqrt((1. + k)**2 - variable_values['b']**2)
                                                      / np.sin(variable_values['i']*constants.deg2rad))
                derived_variables['T_41'] = True

                variable_values['T_32'] =  variable_values['P'] / np.pi \
                                          * np.arcsin(1./variable_values['a'] *
                                                      np.sqrt((1. - k)**2 - variable_values['b']**2)
                                                      / np.sin(variable_values['i']*constants.deg2rad))
                derived_variables['T_32'] = True

            planet_variables[common_name] = variable_values.copy()

            for var in variable_values.keys():
                if var not in derived_variables.keys():
                    del variable_values[var]

            if verbose:
                print_dictionary(variable_values)

    return planet_variables


def get_theta_dictionary(mc):
    # * give back a parameter name associated to each value in the result array

    theta_dictionary = {}
    for dataset_name, dataset in mc.dataset_dict.iteritems():
        for var, i in dataset.variable_sampler.iteritems():
            try:
                theta_dictionary[dataset_name + '_' + var] = i
            except:
                theta_dictionary[repr(dataset_name) + '_' + var] = i

        for model_name in dataset.models:
            for var, i in mc.models[model_name].variable_sampler[dataset_name].iteritems():
                try:
                    theta_dictionary[dataset_name + '_' + model_name + '_' + var] = i
                except:
                    theta_dictionary[repr(dataset_name) + '_' + model_name + '_' + var] = i

    for model in mc.common_models.itervalues():
        for var, i in model.variable_sampler.iteritems():
            theta_dictionary[model.common_ref + '_' + var] = i

    return theta_dictionary


def get_model(mc, theta, bjd_dict):
    model_out = {}
    model_x0 = {}

    delayed_lnlk_computation = {}

    if mc.dynamical_model is not None:
        """ check if any keyword ahas get the output model from the dynamical tool
        we must do it here because all the planet are involved"""
        dynamical_output_x0 = mc.dynamical_model.compute(mc, theta, bjd_dict['full']['x0_plot'])
        dynamical_output = mc.dynamical_model.compute(mc, theta)

    for dataset_name, dataset in mc.dataset_dict.iteritems():

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
                variable_values.update(mc.common_models[common_ref].convert(theta))

            # try:
            #    for common_ref in mc.models[model_name].common_ref:
            #        variable_values.update(mc.common_models[common_ref].convert(theta))
            # except:
            #    continue
            variable_values.update(mc.models[model_name].convert(theta, dataset_name))

            if getattr(mc.models[model_name], 'model_class', None) is 'common_jitter':
                dataset.jitter = mc.models[model_name].compute(variable_values, dataset)
                continue

            if getattr(mc.models[model_name], 'systematic_model', False):
                dataset.additive_model += mc.models[model_name].compute(variable_values, dataset)

        model_out[dataset_name]['systematics'] = dataset.additive_model.copy()
        model_out[dataset_name]['jitter'] = dataset.jitter.copy()
        model_out[dataset_name]['complete'] = np.zeros(dataset.n, dtype=np.double)  # dataset.additive_model.copy()

        model_x0[dataset_name]['complete'] = np.zeros(n_input, dtype=np.double)

        if 'none' in dataset.models or 'None' in dataset.models:
            continue
        if not dataset.models:
            continue

        logchi2_gp_model = None

        for model_name in dataset.models:

            variable_values = {}
            for common_ref in mc.models[model_name].common_ref:
                variable_values.update(mc.common_models[common_ref].convert(theta))
            variable_values.update(mc.models[model_name].convert(theta, dataset_name))

            if getattr(mc.models[model_name], 'internal_likelihood', False):
                logchi2_gp_model = model_name
                continue

            if getattr(mc.models[model_name], 'systematic_model', False):
                continue

            if getattr(dataset, 'dynamical', False):
                dataset.external_model = dynamical_output[dataset_name]
                external_model = dynamical_output_x0[dataset_name].copy()
                model_out[dataset_name]['dynamical'] = dynamical_output[dataset_name].copy()
                model_x0[dataset_name]['dynamical'] = dynamical_output_x0[dataset_name].copy()

            model_out[dataset_name][model_name] = mc.models[model_name].compute(variable_values, dataset)

            if getattr(mc.models[model_name], 'time_independent_model', False):
                model_x0[dataset_name][model_name] = np.zeros(np.size(x0_plot), dtype=np.double)
            else:
                model_x0[dataset_name][model_name] = mc.models[model_name].compute(variable_values, dataset, x0_plot)

            if getattr(mc.models[model_name], 'unitary_model', False):
                dataset.unitary_model += model_out[dataset_name][model_name]
                unitary_model += model_x0[dataset_name][model_name]

                if dataset.normalization_model is None:
                    dataset.normalization_model = np.ones(dataset.n, dtype=np.double)
                    normalization_model = np.ones(np.size(x0_plot), dtype=np.double)

            elif getattr(mc.models[model_name], 'normalization_model', False):
                if dataset.normalization_model is None:
                    dataset.normalization_model = np.ones(dataset.n, dtype=np.double)
                    normalization_model = np.ones(np.size(x0_plot), dtype=np.double)
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
                    variable_values.update(mc.common_models[common_ref].convert(theta))
            except:
                pass

            variable_values.update(mc.models[logchi2_gp_model].convert(theta, dataset.name_ref))

            if hasattr(mc.models[logchi2_gp_model], 'delayed_lnlk_computation'):
                mc.models[logchi2_gp_model].add_internal_dataset(variable_values, dataset,
                                                                 reset_status=delayed_lnlk_computation)
                delayed_lnlk_computation[dataset.name_ref] = logchi2_gp_model

            else:
                model_out[dataset_name][logchi2_gp_model] = \
                    mc.models[logchi2_gp_model].sample_conditional(variable_values, dataset)
                model_out[dataset_name]['complete'] += model_out[dataset_name][logchi2_gp_model]

                model_x0[dataset_name][logchi2_gp_model], var = \
                    mc.models[logchi2_gp_model].sample_predict(variable_values, dataset, x0_plot)

                model_x0[dataset_name][logchi2_gp_model + '_std'] = np.sqrt(var)
                model_x0[dataset_name]['complete'] += model_x0[dataset_name][logchi2_gp_model]

    for dataset_name, logchi2_gp_model in delayed_lnlk_computation.iteritems():
        model_out[dataset_name][logchi2_gp_model] = \
            mc.models[logchi2_gp_model].sample_conditional(mc.dataset_dict[dataset_name])

        model_out[dataset_name]['complete'] += model_out[dataset_name][logchi2_gp_model]

        model_x0[dataset_name][logchi2_gp_model], var = \
            mc.models[logchi2_gp_model].sample_predict(mc.dataset_dict[dataset_name], x0_plot)

        model_x0[dataset_name][logchi2_gp_model + '_std'] = np.sqrt(var)
        model_x0[dataset_name]['complete'] += model_x0[dataset_name][logchi2_gp_model]

    # workaround to avoid memory leaks from GP module
    # gc.collect()

    if np.shape(model_out) == 1:
        return model_out * np.ones(), model_x0

    else:
        return model_out, model_x0


def print_theta_bounds(i_dict, theta, bounds, skip_theta=False):
    format_string = '{0:10s}  {1:4d}  {2:12f} ([{3:10f}, {4:10f}])'
    format_string_long = '{0:10s}  {1:4d}  {2:12f}   {3:12f}  {4:12f} (15-84 p) ([{5:9f}, {6:9f}])'
    format_string_notheta = '{0:10s}  {1:4d}  ([{2:10f}, {3:10f}])'

    for var, i in i_dict.iteritems():

        if skip_theta:
            print(format_string_notheta.format(var, i, bounds[i, 0], bounds[i, 1]))
        elif len(np.shape(theta)) == 2:

            theta_med = compute_value_sigma(theta[:, i])
            print(
                format_string_long.format(var, i, theta_med[0], theta_med[2], theta_med[1], bounds[i, 0], bounds[i, 1]))
        else:
            print(format_string.format(var, i, theta[i], bounds[i, 0], bounds[i, 1]))
    print()


def print_dictionary(variable_values, recenter=[]):
    format_string = '{0:10s}   {1:15f} '
    format_string_long = '{0:10s}   {1:15f}   {2:15f}  {3:15f} (15-84 p)'

    for var_names, var_vals in variable_values.iteritems():
        if np.size(var_vals) > 1:
            if var_names in recenter:
                move_back = (var_vals > recenter[var_names][0] + recenter[var_names][1] / 2.)
                move_forw = (var_vals < recenter[var_names][0] - recenter[var_names][1] / 2.)
                var_vals_recentered = var_vals.copy()
                var_vals_recentered[move_back] -= recenter[var_names][1]
                var_vals_recentered[move_forw] += recenter[var_names][1]
                perc0, perc1, perc2 = np.percentile(var_vals_recentered, [15.865, 50, 84.135], axis=0)

            else:
                perc0, perc1, perc2 = np.percentile(var_vals, [15.865, 50, 84.135], axis=0)

            print(format_string_long.format(var_names, perc1, perc0 - perc1, perc2 - perc1))
        else:
            try:
                print(format_string.format(var_names, var_vals[0]))
            except:
                print(format_string.format(var_names, var_vals))

    print()
