from __future__ import print_function
from pyorbit.models.dataset import Dataset
import sys
import yaml
import copy
from scipy.stats import gaussian_kde

from pyorbit.classes.common import np, get_2darray_from_val

from pyorbit.models.planets import CommonPlanets
from pyorbit.models.activity import CommonActivity
from pyorbit.models.radial_velocities import \
    RVkeplerian, RVdynamical, \
    TransitTimeKeplerian, TransitTimeDynamical, DynamicalIntegrator

from pyorbit.models.batman_transit import Batman_Transit
from pyorbit.models.batman_transit_with_ttv import Batman_Transit_With_TTV

from pyorbit.models.gp_semiperiodic_activity import \
    GaussianProcess_QuasiPeriodicActivity
from pyorbit.models.gp_semiperiodic_activity_common import \
    GaussianProcess_QuasiPeriodicActivity_Common
from pyorbit.models.gp_semiperiodic_activity_shared import \
    GaussianProcess_QuasiPeriodicActivity_Shared

from pyorbit.models.gp_semiperiodic_activity_alternative import \
    GaussianProcess_QuasiPeriodicActivity_Alternative

from pyorbit.models.gp_semiperiodic_activity_derivative import \
    GaussianProcess_QuasiPeriodicActivity_Derivative

from pyorbit.models.celerite_rotation import Celerite_Rotation
from pyorbit.models.celerite_rotation_legacy import Celerite_Rotation_Legacy

from pyorbit.models.correlations import LocalCorrelation
from pyorbit.models.correlated_jitter import LocalCorrelatedJitter
from pyorbit.models.polynomial_trend import CommonPolynomialTrend, \
    PolynomialTrend, LocalPolynomialTrend
from pyorbit.models.common_offset import CommonOffset, Offset
from pyorbit.models.common_jitter import CommonJitter, Jitter
from pyorbit.models.sinusoid_common_period import SinusoidCommonPeriod

from pyorbit.models.batman_limb_darkening import Batman_LimbDarkening_Linear, \
    Batman_LimbDarkening_Quadratic, \
    Batman_LimbDarkening_SquareRoot, Batman_LimbDarkening_Logarithmic, \
    Batman_LimbDarkening_Exponential, Batman_LimbDarkening_Power2, \
    Batman_LimbDarkening_NonLinear

from pyorbit.models.dilution_factor import CommonDilutionFactor, DilutionFactor
from pyorbit.models.normalization_factor import CommonNormalizationFactor, \
    NormalizationFactor
from pyorbit.models.star_parameters import CommonStarParameters

__all__ = ["pars_input", "yaml_parser"]

"""
 model_requires_planets: all those models that requires AT LEAST one of the planets in the system must be listed here
    this is the case for dataset that contains the signature of multiple planets, e.g., RVs or transit light curve 
 single_planet_model: the model is associated to a specific planet, e.g., time of transits 
"""

model_requires_planets = ['radial_velocities',
                          'rv_planets', 'batman_transit', 'batman_transit_with_ttv']
single_planet_model = ['Tc_planets', 'transit_times']
transit_time_model = ['Tc_planets', 'transit_times']

define_common_type_to_class = {
    'planets': CommonPlanets,
    'activity': CommonActivity,
    'polynomial_trend': CommonPolynomialTrend,
    'common_offset': CommonOffset,
    'common_jitter': CommonJitter,
    'batman_ld_linear': Batman_LimbDarkening_Linear,
    'batman_ld_quadratic': Batman_LimbDarkening_Quadratic,
    'batman_ld_square-root': Batman_LimbDarkening_SquareRoot,
    'batman_ld_logarithmic': Batman_LimbDarkening_Logarithmic,
    'batman_ld_exponential': Batman_LimbDarkening_Exponential,
    'batman_ld_power2': Batman_LimbDarkening_Power2,
    'batman_ld_nonlinear': Batman_LimbDarkening_NonLinear,
    'dilution_factor': CommonDilutionFactor,
    'normalization_factor': CommonNormalizationFactor,
    'star_parameters': CommonStarParameters
}

define_type_to_class = {
    'radial_velocities': {'circular': RVkeplerian,
                          'keplerian': RVkeplerian,
                          'dynamical': RVdynamical},
    'rv_planets': {'circular': RVkeplerian,
                   'keplerian': RVkeplerian,
                   'dynamical': RVdynamical},
    'Tc_planets': {'circular': TransitTimeKeplerian,
                   'keplerian': TransitTimeKeplerian,
                   'dynamical': TransitTimeDynamical},
    'transit_times': {'circular': TransitTimeKeplerian,
                      'keplerian': TransitTimeKeplerian,
                      'dynamical': TransitTimeDynamical},
    'batman_transit': Batman_Transit,
    'batman_transit_with_ttv': Batman_Transit_With_TTV,
    'gp_quasiperiodic': GaussianProcess_QuasiPeriodicActivity,
    'gp_quasiperiodic_common': GaussianProcess_QuasiPeriodicActivity_Common,
    'gp_quasiperiodic_shared': GaussianProcess_QuasiPeriodicActivity_Shared,
    'gp_quasiperiodic_alternative': GaussianProcess_QuasiPeriodicActivity_Alternative,
    'gp_quasiperiodic_derivative': GaussianProcess_QuasiPeriodicActivity_Derivative,
    'celerite_rotation': Celerite_Rotation,
    'celerite_rotation_legacy': Celerite_Rotation_Legacy,
    'local_correlation': LocalCorrelation,
    'polynomial_trend': PolynomialTrend,
    'local_polynomial_trend': LocalPolynomialTrend,
    'common_offset': Offset,
    'common_jitter': Jitter,
    'sinusoid_common_period': SinusoidCommonPeriod,
    'dilution_factor': DilutionFactor,
    'normalization_factor': NormalizationFactor,
    'local_correlated_jitter': LocalCorrelatedJitter,
}

accepted_extensions = ['.yaml', '.yml', '.conf', '.config', '.input', ]

star_properties_list = ['limb_darkening', 'dilution_factor']


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


def pars_input(config_in, mc, input_datasets=None, reload_emcee=False, shutdown_jitter=False):

    mc.output_name = config_in['output']

    conf_inputs = config_in['inputs']
    conf_models = config_in['models']
    conf_common = config_in['common']
    conf_parameters = config_in['parameters']
    conf_solver = config_in['solver']

    if conf_models is None:
        conf_models = {'dummy_model': True}

    if reload_emcee:
        if 'emcee' in conf_solver:
            conf = conf_solver['emcee']

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

        # Check if inclination has been updated
        for model_name, model_conf in conf_common.items():
            if not isinstance(model_name, str):
                model_name = repr(model_name)

            if model_name == 'planets':

                for planet_name, planet_conf in model_conf.items():

                    if 'fixed' in planet_conf:
                        fixed_conf = planet_conf['fixed']
                        for var in fixed_conf:
                            mc.common_models[planet_name].fix_list[var] = get_2darray_from_val(
                                fixed_conf[var])

            if model_name == 'star_parameters':
                bounds_space_priors_starts_fixed(
                    mc, mc.common_models[model_name], model_conf)

            if model_name == 'star':
                try:
                    star_conf = model_conf['star_parameters']
                    bounds_space_priors_starts_fixed(
                        mc, mc.common_models['star_parameters'], star_conf)
                except:
                    print()
                    print(" Error in reading the priors from stellar parameters ")

        return

    for dataset_name, dataset_conf in conf_inputs.items():

        if not isinstance(dataset_name, str):
            dataset_name = repr(dataset_name)

        """ The keyword in dataset_dict and the name assigned internally to the databes must be the same
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
                            mc.common_models[planet_name].fix_list['o'] = np.asarray(
                                [np.pi/2., 0.0000], dtype=np.double)
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
                else:
                    model_type = conf_name

                """
                Two ways to include the limb darkening:
                    1) just specify its properties, if all data has been gathered with the same filter
                    2) include different models, for dataset obtained with different filters
                """

                dict_copy = model_conf[conf_name].copy()

                for key in ['type', 'kind', 'priors', 'spaces', 'boundaries', 'starts', 'fixed', 'parametrization']:
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

        if 'star_parameters' not in mc.common_models:
            mc.common_models['star_parameters'] = define_common_type_to_class['star_parameters'](
                'star_parameters')

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
        else:
            model_type = model_name

        if model_type in model_requires_planets or model_type in single_planet_model:

            """ radial_velocities and transits are just wrappers for the planets to be actually included in the model, so we
                substitute it with the individual planets in the list"""

            # try:
            #    planet_list = np.atleast_1d(model_conf['planets']).tolist()
            # except:
            #    planet_list = np.atleast_1d(model_conf['common']).tolist()

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
                if mc.models[model_name_exp].model_class == 'transit':

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
                
                mc.models[model_name_exp].common_ref.append('star_parameters')

                for dataset_name, dataset in mc.dataset_dict.items():
                    if model_name_exp in dataset.models:
                        try:
                            if dataset_name not in model_conf:
                                model_conf[dataset_name] = {}
                                mc.models[model_name_exp].model_conf[dataset_name] = {}

                            bounds_space_priors_starts_fixed(mc, mc.models[model_name_exp], model_conf[dataset_name],
                                                            dataset_1=dataset_name)
                        except TypeError:
                            continue

        elif model_type == 'local_correlation' or model_type == 'correlation':
            mc.models[model_name] = \
                define_type_to_class[model_type](model_name, None)
            mc.models[model_name].model_conf = model_conf.copy()
            bounds_space_priors_starts_fixed(
                mc, mc.models[model_name], model_conf, dataset_1=model_conf['reference'])

        else:

            try:
                mc.models[model_name] = \
                    define_type_to_class[model_type](
                        model_name, model_conf['common'])
            except KeyError:
                mc.models[model_name] = \
                    define_type_to_class[model_type](model_name, None)

            try:
                if model_conf['normalization_model'] is True:
                    mc.models[model_name].normalization_model = True
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

                        bounds_space_priors_starts_fixed(mc, mc.models[model_name], model_conf[dataset_name],
                                                         dataset_1=dataset_name)
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

    if 'emcee' in conf_solver and hasattr(mc, 'emcee_parameters'):
        conf = conf_solver['emcee']

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

def bounds_space_priors_starts_fixed(mc, model_obj, conf, dataset_1=None, dataset_2=None, add_var_name=''):
    # type: (object, object, object, object, object, object) -> object

    if dataset_1 is None:
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
                prior_pams = np.atleast_1d(prior_conf[var])
                model_obj.prior_kind[add_var_name+var] = prior_pams[0]

                if prior_pams[0] == 'File':
                    data_file = np.genfromtxt(prior_pams[1])
                    model_obj.prior_pams[add_var_name + var] = \
                        gaussian_kde(data_file)
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

    elif dataset_2 is None:
        model_obj.bounds[dataset_1] = {}
        model_obj.spaces[dataset_1] = {}
        model_obj.starts[dataset_1] = {}
        model_obj.fix_list[dataset_1] = {}
        model_obj.prior_kind[dataset_1] = {}
        model_obj.prior_pams[dataset_1] = {}

        if 'boundaries' in conf:
            bound_conf = conf['boundaries']
            for var in bound_conf:
                model_obj.bounds[dataset_1][add_var_name + var] = \
                    np.asarray(bound_conf[var], dtype=np.double)

        if 'spaces' in conf:
            space_conf = conf['spaces']
            for var in space_conf:
                model_obj.spaces[dataset_1][add_var_name + var] = \
                    space_conf[var]

        if 'priors' in conf:
            prior_conf = conf['priors']
            for var in prior_conf:
                model_obj.prior_kind[dataset_1][add_var_name + var] = \
                    prior_conf[var][0]

                if prior_conf[var][0] == 'File':
                    data_file = np.genfromtxt(prior_conf[var][1])
                    model_obj.prior_pams[dataset_1][add_var_name + var] = \
                        gaussian_kde(data_file)
                elif np.size(prior_conf[var]) > 1:
                    model_obj.prior_pams[dataset_1][add_var_name + var] = \
                        np.asarray(prior_conf[var][1:], dtype=np.double)
                else:
                    model_obj.prior_pams[dataset_1][add_var_name + var] = \
                        np.asarray([0.00], dtype=np.double)
                
                #model_obj.prior_pams[dataset_1][add_var_name +
                #                                var] = np.asarray(prior_conf[var][1:], dtype=np.double)

        if 'starts' in conf:
            mc.starting_point_flag = True
            starts_conf = conf['starts']
            for var in starts_conf:
                model_obj.starts[dataset_1][add_var_name + var] = \
                    np.asarray(starts_conf[var], dtype=np.double)

        if 'fixed' in conf:
            fixed_conf = conf['fixed']
            for var in fixed_conf:
                model_obj.fix_list[dataset_1][add_var_name + var] = \
                    get_2darray_from_val(fixed_conf[var])

    else:

        if 'boundaries' in conf:
            bound_conf = conf['boundaries']
            for var in bound_conf:
                model_obj.bounds[dataset_1][dataset_2][add_var_name + var] = \
                    np.asarray(bound_conf[var], dtype=np.double)

        if 'spaces' in conf:
            space_conf = conf['spaces']
            for var in space_conf:
                model_obj.spaces[dataset_1][dataset_2][add_var_name + var] = \
                    space_conf[var]

        if 'priors' in conf:
            prior_conf = conf['priors']
            for var in prior_conf:
                model_obj.prior_kind[dataset_1][dataset_2][add_var_name + var] = \
                    prior_conf[var][0]

                if prior_conf[var][0] == 'File':
                    data_file = np.genfromtxt(prior_conf[var][1])
                    model_obj.prior_pams[dataset_1][dataset_2][add_var_name + var] = \
                        gaussian_kde(data_file)
                elif np.size(prior_conf[var]) > 1:
                    model_obj.prior_pams[dataset_1][dataset_2][add_var_name + var] = \
                        np.asarray(prior_conf[var][1:], dtype=np.double)
                else:
                    model_obj.prior_pams[dataset_1][dataset_2][add_var_name + var] = \
                        np.asarray([0.00], dtype=np.double)

                #model_obj.prior_pams[dataset_1][dataset_2][add_var_name + var] = \
                #    np.asarray(prior_conf[var][1:], dtype=np.double)

        if 'starts' in conf:
            mc.starting_point_flag = True
            starts_conf = conf['starts']
            for var in starts_conf:
                model_obj.starts[dataset_1][dataset_2][add_var_name + var] = \
                    np.asarray(starts_conf[var], dtype=np.double)

        if 'fixed' in conf:
            fixed_conf = conf['fixed']
            for var in fixed_conf:
                model_obj.fix_list[dataset_1][dataset_2][add_var_name + var] = \
                    get_2darray_from_val(fixed_conf[var])
    return
