from __future__ import print_function

from pyorbit.model_definitions import *
# Special import for Dataset, it had to be escluded from
# model_definitons to avoid circular import
from pyorbit.common.dataset import Dataset

from pyorbit.subroutines.common import np, get_2darray_from_val

import sys
import yaml
import copy
from scipy.stats import gaussian_kde, multivariate_normal

from pyorbit.subroutines.common import np, get_2darray_from_val

__all__ = ["pars_input", "yaml_parser", "yaml_fix_nested"]


def yaml_parser(file_conf):
    stream = open(file_conf, 'r')

    try:
        config_in = yaml.load(stream, Loader=yaml.FullLoader)
    except AttributeError:
        config_in = yaml.load(stream)
        print(' Consider updating YAML')
    except OSError:
        print(' OSError while reading the configuration file')
        quit()

    if 'output' not in config_in:

        for extension in accepted_extensions:
            if file_conf.find(extension) > 0:
                output_name = file_conf.replace(extension, "")
                continue

        config_in['output'] = output_name

    return config_in

def yaml_fix_nested(config_in):

    print()
    print('Internal reformatting for Nested Sampling compatibility of priors')
    common_conf = config_in['common']
    if 'planets' in common_conf:
        for planet_name in common_conf['planets']:

            if common_conf['planets'][planet_name].get('orbit', 'keplerian') == 'circular':
                continue

            prev_parametrization = common_conf['planets'][planet_name].get('parametrization', 'Eastman2013')
            if prev_parametrization[:5] == 'Stand':
                continue

            change_parametrization = False
            if common_conf['planets'][planet_name].get('priors', False):
                if 'e' in common_conf['planets'][planet_name]['priors'] or \
                    'omega' in common_conf['planets'][planet_name]['priors']:
                    change_parametrization = True

            if change_parametrization:
                if prev_parametrization[-5:] == 'Tcent':
                    common_conf['planets'][planet_name]['parametrization'] = 'Standard_Tcent'
                else:
                    common_conf['planets'][planet_name]['parametrization'] = 'Standard'
                print('    planet {0:s} - parametrization changed from {1:s} to {2:s}'.format(planet_name,
                 prev_parametrization, common_conf['planets'][planet_name]['parametrization']))


    if 'star' not in common_conf:
        print()
        return config_in

    star_conf = config_in['common']['star']
    for starpam_name in star_conf:
        """ We assume that a prior on ld_c1 will be necessarily specified in a limb_darkening object"""
        try:
            if 'ld_c1' in star_conf[starpam_name]['priors']:
                star_conf[starpam_name]['parametrization'] = 'Standard'
                print('    LD {0:s} - LD parametrization changed to Standard'.format(starpam_name))
        except:
            pass

    print()

    return config_in


def pars_input(config_in, mc, input_datasets=None, reload_emcee=False, reload_zeus=False, reload_affine=False, shutdown_jitter=False):

    mc.output_name = config_in['output']

    conf_inputs = config_in['inputs']
    conf_models = config_in['models']
    conf_common = config_in['common']
    conf_parameters = config_in['parameters']
    conf_solver = config_in['solver']

    if conf_models is None:
        conf_models = {'dummy_model': True}


    """ Beginning of snippet dedicated to the reloading of parameters that are not involved in the fit procedure"""
    if reload_emcee or reload_zeus or reload_affine:
        if hasattr(mc, 'emcee_parameters'):
            conf = None
            for conf_name in ['zeus', 'mcmc', 'affine', 'emcee']:
                if conf_name in conf_solver: conf = conf_solver[conf_name]

            if conf:
                if 'multirun' in conf:
                    mc.emcee_parameters['multirun'] = np.asarray(
                        conf['multirun'], dtype=np.int64)

                if 'multirun_iter' in conf:
                    mc.emcee_parameters['multirun_iter'] = np.asarray(
                        conf['multirun_iter'], dtype=np.int64)

                if 'nsave' in conf:
                    mc.emcee_parameters['nsave'] = np.asarray(
                        conf['nsave'], dtype=np.double)

                if 'nsteps' in conf:
                    mc.emcee_parameters['nsteps'] = np.asarray(
                        conf['nsteps'], dtype=np.int64)

                if 'nburn' in conf:
                    mc.emcee_parameters['nburn'] = np.asarray(
                        conf['nburn'], dtype=np.int64)

                if 'thin' in conf:
                    mc.emcee_parameters['thin'] = np.asarray(
                        conf['thin'], dtype=np.int64)


        if hasattr(mc, 'zeus_parameters'):
            conf = None
            for conf_name in ['emcee', 'mcmc', 'affine', 'zeus']:
                if conf_name in conf_solver: conf = conf_solver[conf_name]

            if conf:
                if 'nsave' in conf:
                    mc.zeus_parameters['nsave'] = np.asarray(
                        conf['nsave'], dtype=np.double)

                if 'nsteps' in conf:
                    mc.zeus_parameters['nsteps'] = np.asarray(
                        conf['nsteps'], dtype=np.int64)

                if 'nburn' in conf:
                    mc.zeus_parameters['nburn'] = np.asarray(
                        conf['nburn'], dtype=np.int64)

                if 'thin' in conf:
                    mc.zeus_parameters['thin'] = np.asarray(
                        conf['thin'], dtype=np.int64)


        # Check if inclination has been updated
        for model_name, model_conf in conf_common.items():
            if not isinstance(model_name, str):
                model_name = repr(model_name)

            if model_name == 'planets':

                for planet_name, planet_conf in model_conf.items():

                    if 'star_parameters' in planet_conf:
                        mc.common_models[planet_name].stellar_ref = planet_conf['star_parameters']
                    else:
                        for common_name, common_obj in mc.common_models.items():
                            if getattr(common_obj,'model_class', None) == 'star_parameters':
                                mc.common_models[planet_name].stellar_ref = common_name

                    if 'fixed' in planet_conf:
                        fixed_conf = planet_conf['fixed']
                        for var in fixed_conf:
                            mc.common_models[planet_name].fix_list[var] = get_2darray_from_val(
                                fixed_conf[var])

            #TODO: make more general case - star_parameters with different label
            if model_conf.get('kind',model_name) == 'star_parameters':
                bounds_space_priors_starts_fixed(
                    mc, mc.common_models[model_name], model_conf)

            if model_name == 'star':
                for submodel_name, submodel_conf in model_conf.items():
                    if submodel_conf.get('kind', submodel_name) == 'star_parameters':
                        bounds_space_priors_starts_fixed(
                            mc, mc.common_models[submodel_name], submodel_conf)
        return
    """ End of snippet dedicated to the reloading of parameters that are not involved in the fit procedure"""

    for dataset_name, dataset_conf in conf_inputs.items():

        if not isinstance(dataset_name, str):
            dataset_name = repr(dataset_name)

        """ The keyword in dataset_dict and the name assigned internally to the dataset must be the same
            or everything will fall apart """
        print('Opening: ', dataset_name)

        mc.dataset_dict[dataset_name] = Dataset(dataset_name,
                                                dataset_conf['kind'],
                                                np.atleast_1d(dataset_conf['models']).tolist())

        try:
            data_input = input_datasets[dataset_name]
        except:
            try:
                data_input = mc.dataset_dict[dataset_name].convert_dataset_from_file(
                    dataset_conf['file'])
            except:
                print('Either a file or an input dataset must be provided')
                quit()

        mc.dataset_dict[dataset_name].define_dataset_base(
            data_input, False, shutdown_jitter)

        if mc.Tref:
            mc.dataset_dict[dataset_name].common_Tref(mc.Tref)
        else:
            mc.Tref = mc.dataset_dict[dataset_name].Tref

        if 'boundaries' in dataset_conf:
            bound_conf = dataset_conf['boundaries']
            for var in bound_conf:
                mc.dataset_dict[dataset_name].bounds[var] = np.asarray(
                    bound_conf[var], dtype=np.double)

        if 'spaces' in dataset_conf:
            space_conf = dataset_conf['spaces']
            for var in space_conf:
                mc.dataset_dict[dataset_name].spaces[var] = space_conf[var]

        if 'starts' in dataset_conf:
            mc.starting_point_flag = True
            starts_conf = dataset_conf['starts']
            for var in starts_conf:
                mc.dataset_dict[dataset_name].starts[var] = np.asarray(
                    starts_conf[var], dtype=np.double)

        mc.dataset_dict[dataset_name].update_bounds_spaces_priors_starts()

        """ Added in PyORBIT 9: ancillary datasets for additional datasets that
            are required for the modelling (e.g., instrumental parameters) but
            do not required any modelling themselves
            Differently from the standard pyorbit dataset, the first row must
            include the dataset names
         """
        if 'ancillary_name' in dataset_conf:
            mc.dataset_dict[dataset_name].ancillary = input_datasets[dataset_conf['ancillary_name']]

        for ancill in ['ancillary_file', 'ancillary', 'ancillary_data', 'ancillary_dataset']:
            if ancill in dataset_conf:
                mc.dataset_dict[dataset_name].convert_ancillary_from_file(dataset_conf[ancill])


    for model_name, model_conf in conf_common.items():

        """ Accept names starting with numbers (but why you should do that?"""
        if not isinstance(model_name, str):
            model_name = repr(model_name)

        if model_name == 'planets':

            for planet_name, planet_conf in model_conf.items():

                print('Adding common model of planet: ', planet_name)

                if not isinstance(planet_name, str):
                    planet_name = repr(planet_name)

                mc.common_models[planet_name] = define_common_type_to_class['planets'](
                    planet_name)

                bounds_space_priors_starts_fixed(
                    mc, mc.common_models[planet_name], planet_conf)

                try:
                    if planet_conf['orbit'] in mc.common_models[planet_name].orbit_list:

                        #mc.planet_dict[planet_name] = planet_conf['orbit']
                        mc.common_models[planet_name].orbit = planet_conf['orbit']

                        print('Using orbital model: ',
                              mc.common_models[planet_name].orbit)

                        if planet_conf['orbit'] == 'circular':
                            mc.common_models[planet_name].fix_list['e'] = np.asarray(
                                [0.000, 0.0000], dtype=np.double)
                            mc.common_models[planet_name].fix_list['omega'] = np.asarray(
                                [90.0, 0.0000], dtype=np.double)
                        if planet_conf['orbit'] == 'dynamical':
                            mc.dynamical_dict[planet_name] = True
                    else:
                        #mc.planet_dict[planet_name] = mc.common_models[planet_name].orbit
                        print('Orbital model not recognized, switching to: ',
                              mc.common_models[planet_name].orbit)
                except:
                    #mc.planet_dict[planet_name] = mc.common_models[planet_name].orbit
                    print('Using default orbital model: ',
                          mc.common_models[planet_name].orbit)

                try:
                    if planet_conf['parametrization'] in mc.common_models[planet_name].parametrization_list:
                        mc.common_models[planet_name].parametrization = planet_conf['parametrization']
                        print('Using orbital parametrization: ',
                              mc.common_models[planet_name].parametrization)
                    else:
                        print('Orbital parametrization not recognized, switching to: ',
                              mc.common_models[planet_name].parametrization)

                    if mc.common_models[planet_name].parametrization[-5:] == 'Tcent' or \
                            mc.common_models[planet_name].parametrization[-5:] == 'Tc':
                        print('Using Central Time of Transit instead of phase')
                        mc.common_models[planet_name].use_time_of_transit = True

                except:
                    print('Using default orbital parametrization: ',
                          mc.common_models[planet_name].parametrization)

                try:
                    mc.common_models[planet_name].use_inclination = planet_conf['use_inclination']
                    print('Inclination will be included as free parameter: ',
                          planet_conf['use_inclination'])
                except:
                    for key in ['boundaries', 'spaces', 'priors', 'starts', 'fixed']:
                        if key in planet_conf and 'i' in planet_conf[key]:
                            mc.common_models[planet_name].use_inclination = True
                            print(
                                'Inclination will be included as free parameter: ', True)
                    # False by default, unless the user has specified some of its properties

                try:
                    mc.common_models[planet_name].use_semimajor_axis = planet_conf['use_semimajor_axis']
                    print('Semi-major axis will be included as free parameter: {}'.format(
                        planet_conf['use_inclination']))
                except:
                    for key in ['boundaries', 'spaces', 'priors', 'starts', 'fixed']:
                        if key in planet_conf and 'a' in planet_conf[key]:
                            mc.common_models[planet_name].use_semimajor_axis = True
                            print(
                                'Semi-major axis will be included as free parameter: ', True)
                    # False by default, unless the user has specified some of its properties

                try:
                    mc.common_models[planet_name].use_time_of_transit = planet_conf['use_time_of_transit']
                    print('Using Central Time of Transit instead of phase: {}'.format(
                        planet_conf['use_time_of_transit']))
                except:
                    # False by default
                    pass

                try:
                    mc.common_models[planet_name].use_mass_for_planets = planet_conf['use_mass_for_planets']
                    print('Using planetary mass instead of RV semiamplitude: '.format(
                        planet_conf['use_mass_for_planets']))
                except:
                    # False by default
                    pass

                print()

        elif model_name == 'star':

            for conf_name, star_conf in model_conf.items():
                if 'type' in star_conf:
                    model_type = star_conf['type']
                elif 'kind' in star_conf:
                    model_type = star_conf['kind']
                elif 'model' in star_conf:
                    model_type = star_conf['model']
                else:
                    model_type = conf_name

                """
                Two ways to include the limb darkening:
                    1) just specify its properties, if all data has been gathered with the same filter
                    2) include different models, for dataset obtained with different filters
                """

                dict_copy = model_conf[conf_name].copy()

                for key in ['type', 'kind','model', 'priors', 'spaces', 'boundaries', 'starts', 'fixed', 'parametrization']:
                    if key in dict_copy:
                        del dict_copy[key]

                if len(dict_copy) == 0:
                    dict_copy = {conf_name: model_conf[conf_name].copy()}

                for key_name, key_vals in dict_copy.items():
                    mc.common_models[key_name] = define_common_type_to_class[model_type](
                        key_name)
                    bounds_space_priors_starts_fixed(
                        mc, mc.common_models[key_name], key_vals)

                    """ Automatic detection of common models without dataset-specific parameters"""
                    for dataset_name, dataset in mc.dataset_dict.items():
                        if key_name in dataset.models and not (key_name in conf_models):
                            conf_models[key_name] = {'common': key_name}

                    try:
                        if key_vals['parametrization'] in mc.common_models[key_name].parametrization_list:
                            mc.common_models[key_name].parametrization = key_vals['parametrization']
                            print('Using limb darkening coefficient parametrization: ',
                                  mc.common_models[key_name].parametrization)
                    except:
                        continue

        else:

            """
                If the object type/kind has been specified, than the section name is treated just as a label.
                Otherwise, it is assumed that the section name correspond to the object name
            """
            if 'type' in model_conf:
                model_type = model_conf['type']
            elif 'kind' in model_conf:
                model_type = model_conf['kind']
            elif 'model' in model_conf:
                model_type = model_conf['model']
            else:
                model_type = model_name

            mc.common_models[model_name] = define_common_type_to_class[model_type](
                model_name)
            bounds_space_priors_starts_fixed(
                mc, mc.common_models[model_name], model_conf)

            """ Automatic detection of common models without dataset-specific parameters
                If the required parameters are included in the "model-specific section",
                the relative dictionary is created and the keywords are copied there
            """

            for dataset_name, dataset in mc.dataset_dict.items():
                if model_name in dataset.models and (model_name in conf_models):
                    """ Using the same name for a common model and a data-specific model causes
                    the automatic assignment of the former to the latter
                    """
                    conf_models[model_name].update(model_conf)
                    conf_models[model_name]['common'] = model_name
                if model_name in dataset.models and not (model_name in conf_models):
                    """ Common model will be used as data-specific model for this dataset"""

                    try:
                        conf_models[model_name] = model_conf
                        conf_models[model_name]['common'] = model_name
                    except:
                        conf_models[model_name] = {'common': model_name}


        # Add the stellar parameter class anyway
        #if 'star_parameters' not in mc.common_models:
        #    mc.common_models['star_parameters'] = define_common_type_to_class['star_parameters'](
        #        'star_parameters')


    """ Check if there is any planet that requires dynamical computations"""
    if mc.dynamical_dict:
        mc.dynamical_model = DynamicalIntegrator()

    for model_name, model_conf in conf_models.items():

        if model_name == 'dummy_model':
            continue

        if not isinstance(model_name, str):
            model_name = repr(model_name)

        """ Check if the keplerian approximation must be used for this dataset even if the planet has the dynamical flag"""

        keplerian_approximation = False
        if 'keplerian_approximation' in model_conf:
            keplerian_approximation = model_conf['keplerian_approximation']
            print('Using Keplerian approximation')

        if 'type' in model_conf:
            model_type = model_conf['type']
        elif 'kind' in model_conf:
            model_type = model_conf['kind']
        elif 'model' in model_conf:
            model_type = model_conf['model']
        else:
            model_type = model_name

        """ some models requires one or more planets, with some specific properties"""
        if model_type in model_requires_planets or model_type in single_planet_model:

            """ radial_velocities and transits are just wrappers for the planets to be actually included in the model, so we
                substitute it with the individual planets in the list"""

            try:
                model_name_expanded = []
                model_name_original = []
                planet_list = []
                for key, val in model_conf['models'].items():
                    model_name_original.append(key)
                    model_name_expanded.append(key)
                    planet_list.append(val)

            except:
                try:
                    planet_list = np.atleast_1d(model_conf['planets']).tolist()
                except:
                    planet_list = np.atleast_1d(model_conf['common']).tolist()
                model_name_original = [model_name for pl_name in planet_list]
                model_name_expanded = [model_name +
                                       '_' + pl_name for pl_name in planet_list]

            """ Let's avoid some dumb user using the planet names to name the models"""

            if model_type in model_requires_planets:
                """ For each dataset we check if the current model is included in the list of models.
                    We then remove the generic model name  and include all the planet-specific  model names
                """
                for dataset_name, dataset in mc.dataset_dict.items():
                    if model_name in dataset.models:
                        dataset.models.remove(model_name)
                        dataset.models.extend(model_name_expanded)

                        if len(list(set(planet_list) & set(mc.dynamical_dict))) and not keplerian_approximation:
                            dataset.dynamical = True

            for model_name_exp, planet_name in zip(model_name_expanded, planet_list):

                try:
                    mc.models[model_name_exp] = \
                        define_type_to_class[model_type][mc.common_models[planet_name].orbit](
                            model_name_exp, planet_name)

                    if keplerian_approximation:
                        """ override default model from input file, to apply keplerian approximation """
                        mc.models[model_name_exp] = \
                            define_type_to_class[model_type]['keplerian'](
                                model_name_exp, planet_name)
                        mc.common_models[planet_name].use_mass_for_planets = True
                        print('Using planetary mass instead of RV semiamplitude: ',
                              planet_conf['use_mass_for_planets'])

                except:
                    mc.models[model_name_exp] = \
                        define_type_to_class[model_type](
                            model_name_exp, planet_name)

                    mc.models[model_name_exp].model_conf = model_conf.copy()

                if model_type in transit_time_model:

                    for dataset_name, dataset in mc.dataset_dict.items():
                        if planet_name in mc.dynamical_dict and \
                                model_name_exp in dataset.models and \
                                not keplerian_approximation:
                            dataset.planet_name = planet_name
                            dataset.dynamical = True
                            mc.dynamical_t0_dict[planet_name] = dataset_name

                """ This snippet will work only for transit class"""
                if mc.models[model_name_exp].model_class in model_requires_limbdarkening:

                    try:
                        common_name = mc.models[model_name_exp].model_conf['limb_darkening']
                    except:
                        common_name = 'limb_darkening'

                    print('  LC model: {0:s} is using {1:s} LD parameters'.format(
                        model_name_exp, common_name))

                    mc.models[model_name_exp].model_conf['limb_darkening_model'] = \
                        mc.common_models[common_name].ld_type

                    mc.models[model_name_exp].model_conf['limb_darkening_ncoeff'] = \
                        mc.common_models[common_name].ld_ncoeff

                    mc.models[model_name_exp].common_ref.append(common_name)

                try:
                    common_name = mc.models[model_name_exp].model_conf['star_parameters']
                except:
                    common_name = 'star_parameters'

                mc.models[model_name_exp].common_ref.append(common_name)
                mc.models[model_name_exp].stellar_ref = common_name

                for dataset_name, dataset in mc.dataset_dict.items():
                    if model_name_exp in dataset.models:
                        try:
                            if dataset_name not in model_conf:
                                model_conf[dataset_name] = {}
                                mc.models[model_name_exp].model_conf[dataset_name] = {}

                            bounds_space_priors_starts_fixed(mc, mc.models[model_name_exp], model_conf[dataset_name],
                                                            dataset=dataset_name, backup_conf=model_conf)
                        except TypeError:
                            continue

        else:

            if model_conf.get('common', False):
                common_ref = model_conf['common']
            elif hasattr(define_type_to_class[model_type], 'default_common'):

                """ New: some data-specific models may noy need a common model, e.g., local polynomial trends,
                    however the code requires that such common model is provided in order to have a fall-back for the
                    default priors, boundaries, spaces, and fixed parameters. Such common model is saved with the same
                    name of the data-specifc model.
                    Default common model only works it has been defined in the model class as a Class attribute
                """
                common_type = define_type_to_class[model_type].default_common
                common_ref = model_name

                if model_name not in mc.common_models:
                    mc.common_models[common_ref] = define_common_type_to_class[common_type](common_ref)
            else:
                common_ref = None

            mc.models[model_name] = define_type_to_class[model_type](model_name, common_ref)

            """ Check if we need to add the limb darkening parameters to the model"""
            if mc.models[model_name].model_class in model_requires_limbdarkening:

                try:
                    common_name = mc.models[model_name].model_conf['limb_darkening']
                except:
                    common_name = 'limb_darkening'

                print('  model: {0:s} is using {1:s} LD parameters'.format(
                    model_name, common_name))

                model_conf['limb_darkening_model'] = \
                    mc.common_models[common_name].ld_type

                model_conf['limb_darkening_ncoeff'] = \
                    mc.common_models[common_name].ld_ncoeff

                mc.models[model_name].common_ref.append(common_name)

            """ Adding the stlerra parameters common model by default """
            ##TODO: is it really needed?
            #try:
            #    common_name = mc.models[model_name].model_conf['star_parameters']
            #except:
            #    common_name = 'star_parameters'
            #mc.models[model_name].common_ref.append(common_name)
            #mc.models[model_name].stellar_ref = common_name







            """ A model can be exclusively unitary, additive, or normalization.
                How the individual model are combined to provided the final model is embedded in the Dataset class
            """
            try:
                if model_conf['normalization_model'] is True:
                    mc.models[model_name].normalization_model = True
                    mc.models[model_name].unitary_model = False
                    print('Model type: normalization')
            except KeyError:
                pass

            try:
                if model_conf['multiplicative_model'] is True:
                    mc.models[model_name].multiplicative_model = True
                    mc.models[model_name].normalization_model = False
                    mc.models[model_name].unitary_model = False
                    print('Model type: normalization')
            except KeyError:
                pass


            try:
                if model_conf['unitary_model'] is True:
                    mc.models[model_name].unitary_model = True
                    mc.models[model_name].normalization_model = False
                    print('Model type: unitary')
            except KeyError:
                pass

            try:
                if model_conf['additive_model'] is True:
                    mc.models[model_name].unitary_model = False
                    mc.models[model_name].normalization_model = False
                    print('Model type: additive')
            except KeyError:
                pass

            """ Using default boundaries if one dataset is missing"""
            try:
                for dataset_name, dataset in mc.dataset_dict.items():

                    if model_name in dataset.models:

                        if dataset_name not in model_conf:
                            model_conf[dataset_name] = {}

                        bounds_space_priors_starts_fixed(mc,
                                                         mc.models[model_name],
                                                         model_conf[dataset_name],
                                                         dataset=dataset_name,
                                                         backup_conf=model_conf)

            except:
                pass

            try:
                mc.models[model_name].model_conf.update(model_conf)
            except:
                mc.models[model_name].model_conf = model_conf.copy()

    if 'ordered_planets' in conf_parameters:
        mc.ordered_planets = conf_parameters['ordered_planets']

    if 'Tref' in conf_parameters:
        mc.Tref = np.asarray(conf_parameters['Tref'])
        for dataset_name in mc.dataset_dict:
            mc.dataset_dict[dataset_name].common_Tref(mc.Tref)

    if 'star_mass' in conf_parameters:
        mc.star_mass = np.asarray(
            conf_parameters['star_mass'][:], dtype=np.double)
    if 'star_radius' in conf_parameters:
        mc.star_radius = np.asarray(
            conf_parameters['star_radius'][:], dtype=np.double)

    if 'dynamical_integrator' in conf_parameters:
        try:
            mc.dynamical_model.dynamical_integrator = conf_parameters['dynamical_integrator']
        except:
            pass

    if 'dynamical_integrator' in conf_solver:
        try:
            mc.dynamical_model.dynamical_integrator = conf_solver['dynamical_integrator']
        except:
            pass

    shutdown_jitter = False
    if 'shutdown_jitter' in conf_parameters:
        shutdown_jitter = conf_parameters['shutdown_jitter']

    if 'shutdown_jitter' in conf_solver:
        shutdown_jitter = conf_solver['shutdown_jitter']

    if 'pyde' in conf_solver and hasattr(mc, 'pyde_parameters'):
        conf = conf_solver['pyde']

        if 'ngen' in conf:
            mc.pyde_parameters['ngen'] = np.asarray(
                conf['ngen'], dtype=np.double)

        if 'npop_mult' in conf:
            mc.pyde_parameters['npop_mult'] = np.asarray(
                conf['npop_mult'], dtype=np.int64)

        if 'shutdown_jitter' in conf:
            mc.pyde_parameters['shutdown_jitter'] = np.asarray(
                conf['shutdown_jitter'], dtype=bool)
        elif shutdown_jitter:
            mc.pyde_parameters['shutdown_jitter'] = shutdown_jitter

        if 'include_priors' in conf:
            mc.include_priors = np.asarray(conf['include_priors'], dtype=bool)

        if 'use_threading_pool' in conf:
            mc.pyde_parameters['use_threading_pool'] = np.asarray(conf['use_threading_pool'], dtype=bool)


    if hasattr(mc, 'emcee_parameters'):
        conf = None
        for conf_name in ['zeus', 'mcmc', 'affine', 'emcee']:
            if conf_name in conf_solver: conf = conf_solver[conf_name]
        if conf:
            if 'multirun' in conf:
                mc.emcee_parameters['multirun'] = np.asarray(
                    conf['multirun'], dtype=np.int64)

            if 'multirun_iter' in conf:
                mc.emcee_parameters['multirun_iter'] = np.asarray(
                    conf['multirun_iter'], dtype=np.int64)

            if 'nsave' in conf:
                mc.emcee_parameters['nsave'] = np.asarray(
                    conf['nsave'], dtype=np.double)

            if 'nsteps' in conf:
                mc.emcee_parameters['nsteps'] = np.asarray(
                    conf['nsteps'], dtype=np.int64)

            if 'nburn' in conf:
                mc.emcee_parameters['nburn'] = np.asarray(
                    conf['nburn'], dtype=np.int64)

            if 'npop_mult' in conf:
                mc.emcee_parameters['npop_mult'] = np.asarray(
                    conf['npop_mult'], dtype=np.int64)

            if 'thin' in conf:
                mc.emcee_parameters['thin'] = np.asarray(
                    conf['thin'], dtype=np.int64)

            if 'shutdown_jitter' in conf:
                mc.emcee_parameters['shutdown_jitter'] = np.asarray(
                    conf['shutdown_jitter'], dtype=bool)
            elif shutdown_jitter:
                mc.emcee_parameters['shutdown_jitter'] = shutdown_jitter

            if 'include_priors' in conf:
                mc.include_priors = np.asarray(conf['include_priors'], dtype=bool)

            if 'use_threading_pool' in conf:
                mc.emcee_parameters['use_threading_pool'] = np.asarray(conf['use_threading_pool'], dtype=bool)

    if hasattr(mc, 'zeus_parameters'):
        conf = None
        for conf_name in ['emcee', 'mcmc', 'affine', 'zeus']:
            if conf_name in conf_solver: conf = conf_solver[conf_name]
        if conf:
            if 'nsteps' in conf:
                mc.zeus_parameters['nsteps'] = np.asarray(
                    conf['nsteps'], dtype=np.int64)

            if 'nburn' in conf:
                mc.zeus_parameters['nburn'] = np.asarray(
                    conf['nburn'], dtype=np.int64)

            if 'npop_mult' in conf:
                mc.zeus_parameters['npop_mult'] = np.asarray(
                    conf['npop_mult'], dtype=np.int64)

            if 'thin' in conf:
                mc.zeus_parameters['thin'] = np.asarray(
                    conf['thin'], dtype=np.int64)

            if 'shutdown_jitter' in conf:
                mc.zeus_parameters['shutdown_jitter'] = np.asarray(
                    conf['shutdown_jitter'], dtype=bool)
            elif shutdown_jitter:
                mc.zeus_parameters['shutdown_jitter'] = shutdown_jitter

            if 'include_priors' in conf:
                mc.include_priors = np.asarray(conf['include_priors'], dtype=bool)

            if 'use_threading_pool' in conf:
                mc.zeus_parameters['use_threading_pool'] = np.asarray(conf['use_threading_pool'], dtype=bool)


    if 'nested_sampling' in conf_solver and hasattr(mc, 'nested_sampling_parameters'):
        conf = conf_solver['nested_sampling']

        for key_name, key_value in conf.items():
            mc.nested_sampling_parameters[key_name] = key_value

        if 'include_priors' in conf:
            mc.include_priors = np.asarray(conf['include_priors'], dtype=bool)

        if 'ordered_planets' in conf:
            mc.ordered_planets = conf['ordered_planets']

    if 'optimize' in conf_solver and hasattr(mc, 'optimize_parameters'):
        conf = conf_solver['optimize']

        for key_name, key_value in conf.items():
            mc.optimize_parameters[key_name] = key_value

    if 'recenter_bounds' in conf_solver:
        """
        required to avoid a small bug in the code
        if the dispersion of PyDE walkers around the median value is too broad,
        then emcee walkers will start outside the bounds, causing an error
        """
        mc.recenter_bounds_flag = conf_solver['recenter_bounds']

    if 'include_priors' in conf_solver:
        mc.include_priors = np.asarray(
            conf_solver['include_priors'], dtype=bool)

def bounds_space_priors_starts_fixed(mc,
                                     model_obj,
                                     input_conf,
                                     dataset=None,
                                     add_var_name='',
                                     backup_conf=None):
    # type: (object, object, object, object, object, object) -> object

    conf = input_conf.copy()
    key_list = ['boundaries', 'spaces', 'priors', 'starts', 'fixed']


    """ Copy the missing keys from the backup configuration
    Example: we have N+1 datasets, namely LCdata_transit00 ... LCdata_transitN,
    all including a local_polynomial_trend model.
    In the "models" section of the configuration file we can specify:

    models:
      ...
      local_polynomial_trend:
        common: polynomial_trend
        normalization_model: True
        order: 2
        boundaries:
          poly_c0: [-2.0, 2.0]
          poly_c1: [-1.0, 1.0]
          poly_c2: [-1.0, 1.0]
        LCdata_transit00:
          boundaries:
            poly_c0: [-4.0, 4.0]
            poly_c1: [-3.0, 3.0]
            poly_c2: [-3.0, 3.0]
    parameters:
      ...

    Datasets LCdata_transit00 will have boundaries poly_c0: [-4.0, 4.0],
    poly_c1: [-3.0, 3.0], etc
    All the other datasets will have bondaries poly_c0: [-2.0, 2.0],
    poly_c1: [-1.0, 1.0], etc

    This formatting applies to all the items in "key_list"

    Code is super-utly, and probably a real coder will faint by looking at it :-(
    """

    for key_name in key_list:
        if backup_conf is None:
            continue
        if key_name not in backup_conf:
            continue
        if key_name not in conf:
            conf[key_name] = backup_conf[key_name]
        else:
            for var_name in backup_conf[key_name]:
                if var_name not in  conf[key_name]:
                    conf[key_name][var_name] = backup_conf[key_name][var_name]

    if dataset is None:
        if 'boundaries' in conf:
            bound_conf = conf['boundaries']
            for var in bound_conf:
                model_obj.bounds[add_var_name +
                                 var] = np.asarray(bound_conf[var], dtype=np.double)

        if 'spaces' in conf:
            space_conf = conf['spaces']
            for var in space_conf:
                model_obj.spaces[add_var_name+var] = space_conf[var]

        if 'priors' in conf:
            prior_conf = conf['priors']
            for var in prior_conf:

                if var == 'multivariate':
                    model_obj.multivariate_priors = True

                    model_obj.multivariate_vars = prior_conf[var]['variables']
                    data_file = np.genfromtxt(prior_conf[var]['file'])

                    ll = []
                    for ii in range(len(model_obj.multivariate_vars)):
                        ll.append(data_file[:,ii])
                    cov_data = np.stack(ll, axis=0)
                    model_obj.multivariate_cov = np.cov(cov_data)
                    model_obj.multivariate_med = np.median(data_file, axis=0)
                    model_obj.multivariate_func = multivariate_normal(model_obj.multivariate_med,
                                                                      model_obj.multivariate_cov)
                    ll = None
                    cov_data = None
                    data_file = None

                else:
                    prior_pams = np.atleast_1d(prior_conf[var])
                    model_obj.prior_kind[add_var_name+var] = prior_pams[0]

                    if prior_pams[0] == 'File':
                        data_file = np.genfromtxt(prior_pams[1])
                        model_obj.prior_pams[add_var_name + var] = \
                            gaussian_kde(data_file)
                        data_file = None
                    elif np.size(prior_pams) > 1:
                        model_obj.prior_pams[add_var_name + var] = \
                            np.asarray(prior_pams[1:], dtype=np.double)
                    else:
                        model_obj.prior_pams[add_var_name + var] = \
                            np.asarray([0.00], dtype=np.double)

        if 'starts' in conf:
            mc.starting_point_flag = True
            starts_conf = conf['starts']
            for var in starts_conf:
                model_obj.starts[add_var_name + var] = \
                    np.asarray(starts_conf[var], dtype=np.double)

        if 'fixed' in conf:
            fixed_conf = conf['fixed']
            for var in fixed_conf:
                model_obj.fix_list[add_var_name + var] = \
                    get_2darray_from_val(fixed_conf[var])

    else:
        model_obj.bounds[dataset] = {}
        model_obj.spaces[dataset] = {}
        model_obj.starts[dataset] = {}
        model_obj.fix_list[dataset] = {}
        model_obj.prior_kind[dataset] = {}
        model_obj.prior_pams[dataset] = {}

        if 'boundaries' in conf:
            bound_conf = conf['boundaries']
            for var in bound_conf:
                model_obj.bounds[dataset][add_var_name + var] = \
                    np.asarray(bound_conf[var], dtype=np.double)

        if 'spaces' in conf:
            space_conf = conf['spaces']
            for var in space_conf:
                model_obj.spaces[dataset][add_var_name + var] = \
                    space_conf[var]

        if 'priors' in conf:
            prior_conf = conf['priors']
            for var in prior_conf:



                if var == 'multivariate':
                    model_obj.multivariate_priors[dataset] = True

                    model_obj.multivariate_vars[dataset] = prior_conf[var]['variables']
                    data_file = np.genfromtxt(prior_conf[var]['file'])

                    ll = []
                    for ii in range(len(model_obj.multivariate_vars[dataset])):
                        ll.append(data_file[:,ii])
                    cov_data = np.stack(ll, axis=0)
                    model_obj.multivariate_cov[dataset] = np.cov(cov_data)
                    model_obj.multivariate_med[dataset] = np.median(data_file, axis=0)
                    model_obj.multivariate_func[dataset] = multivariate_normal(model_obj.multivariate_med[dataset],
                                                                      model_obj.multivariate_cov[dataset])
                    ll = None
                    cov_data = None
                    data_file = None





                prior_pams = np.atleast_1d(prior_conf[var])
                model_obj.prior_kind[dataset][add_var_name + var] = \
                    prior_conf[var][0]

                if prior_conf[var][0] == 'File':
                    data_file = np.genfromtxt(prior_pams[1])
                    model_obj.prior_pams[dataset][add_var_name + var] = \
                        gaussian_kde(data_file)
                    data_file = None
                elif np.size(prior_conf[var]) > 1:
                    model_obj.prior_pams[dataset][add_var_name + var] = \
                        np.asarray(prior_pams[1:], dtype=np.double)
                else:
                    model_obj.prior_pams[dataset][add_var_name + var] = \
                        np.asarray([0.00], dtype=np.double)


        if 'starts' in conf:
            mc.starting_point_flag = True
            starts_conf = conf['starts']
            for var in starts_conf:
                model_obj.starts[dataset][add_var_name + var] = \
                    np.asarray(starts_conf[var], dtype=np.double)

        if 'fixed' in conf:
            fixed_conf = conf['fixed']
            for var in fixed_conf:
                model_obj.fix_list[dataset][add_var_name + var] = \
                    get_2darray_from_val(fixed_conf[var])


    return
