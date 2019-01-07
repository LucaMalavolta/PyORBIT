from common import *
from ..models.dataset import *
from ..models.planets import CommonPlanets
from ..models.activity import CommonActivity
from ..models.radial_velocities import RVkeplerian, RVdynamical, TransitTimeKeplerian, TransitTimeDynamical, DynamicalIntegrator, RVkeplerianMass
from ..models.batman_transit import Batman_Transit
from ..models.gp_semiperiodic_activity import GaussianProcess_QuasiPeriodicActivity
from ..models.gp_semiperiodic_activity_common import GaussianProcess_QuasiPeriodicActivity_Common
from ..models.gp_semiperiodic_activity_shared import GaussianProcess_QuasiPeriodicActivity_Shared

if 'celerite' in sys.modules:
  from ..models.celerite_semiperiodic_activity import Celerite_QuasiPeriodicActivity

from ..models.correlations import Correlation_SingleDataset
from ..models.polynomial_trend import CommonPolynomialTrend, PolynomialTrend
from ..models.common_offset import CommonOffset, Offset
from ..models.common_jitter import CommonJitter, Jitter
from ..models.sinusoid_common_period import SinusoidCommonPeriod
__all__ = ["pars_input", "yaml_parser"]

define_common_type_to_class = {
    'planets': CommonPlanets,
    'activity': CommonActivity,
    'polynomial_trend': CommonPolynomialTrend,
    'common_offset': CommonOffset,
    'common_jitter': CommonJitter
}

define_type_to_class = {
    'radial_velocities': {'circular': RVkeplerian,
                          'keplerian': RVkeplerian,
                          'dynamical': RVdynamical,
                          'keplerian_mass': RVkeplerianMass},
    'rv_planets': {'circular': RVkeplerian,
                   'keplerian': RVkeplerian,
                   'dynamical': RVdynamical,
                   'keplerian_mass': RVkeplerianMass},
    'transit_time': {'circular': TransitTimeKeplerian,
                     'keplerian': TransitTimeKeplerian,
                     'dynamical': TransitTimeDynamical},
    'batman_transit': Batman_Transit,
    'gp_quasiperiodic': GaussianProcess_QuasiPeriodicActivity,
    'gp_quasiperiodic_common': GaussianProcess_QuasiPeriodicActivity_Common,
    'gp_quasiperiodic_shared': GaussianProcess_QuasiPeriodicActivity_Shared,
    'celerite_quasiperiodic': Celerite_QuasiPeriodicActivity,
    'correlation_singledataset': Correlation_SingleDataset,
    'polynomial_trend': PolynomialTrend,
    'common_offset': Offset,
    'common_jitter': Jitter,
    'sinusoid_common_period': SinusoidCommonPeriod
}

accepted_extensions = ['.yaml', '.yml', '.conf', '.config', '.input', ]


def yaml_parser(file_conf):
    stream = file(file_conf, 'r')
    config_in = yaml.load(stream)

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

    if reload_emcee:
        if 'star_mass' in conf_parameters:
            mc.star_mass = np.asarray(conf_parameters['star_mass'][:], dtype=np.double)
        if 'star_radius' in config_in:
            mc.star_radius = np.asarray(conf_parameters['star_radius'][:], dtype=np.double)

        if 'emcee' in conf_solver:
            conf = conf_solver['emcee']

            if 'multirun' in conf:
                mc.emcee_parameters['multirun'] = np.asarray(conf['multirun'], dtype=np.int64)

            if 'multirun_iter' in conf:
                mc.emcee_parameters['multirun_iter'] = np.asarray(conf['multirun_iter'], dtype=np.int64)

            if 'nsave' in conf:
                mc.emcee_parameters['nsave'] = np.asarray(conf['nsave'], dtype=np.double)

            if 'nsteps' in conf:
                mc.emcee_parameters['nsteps'] = np.asarray(conf['nsteps'], dtype=np.int64)

            if 'nburn' in conf:
                mc.emcee_parameters['nburn'] = np.asarray(conf['nburn'], dtype=np.int64)

            if 'thin' in conf:
                mc.emcee_parameters['thin'] = np.asarray(conf['thin'], dtype=np.int64)

        # Check if inclination has been updated
        for model_name, model_conf in conf_common.iteritems():
            if not isinstance(model_name, str):
                model_name = repr(model_name)

            if model_name == 'planets':

                for planet_name, planet_conf in model_conf.iteritems():

                    if 'fixed' in planet_conf:
                        fixed_conf = planet_conf['fixed']
                        for var in fixed_conf:
                            mc.common_models[planet_name].fix_list[var] = get_2darray_from_val(fixed_conf[var])


        return

    for dataset_name, dataset_conf in conf_inputs.iteritems():

        if not isinstance(dataset_name, str):
            dataset_name = repr(dataset_name)

        """ The keyword in dataset_dict and the name assigned internally to the databes must be the same
            or everything will fall apart """
        mc.dataset_dict[dataset_name] = Dataset(dataset_name,
                                                dataset_conf['kind'],
                                                np.atleast_1d(dataset_conf['models']).tolist())

        try:
            data_input = input_datasets[dataset_name]
        except:
            try:
                data_input = mc.dataset_dict[dataset_name].convert_dataset_from_file(dataset_conf['file'])
            except:
                print 'Either a file or an input dataset must be provided'
                quit()

        mc.dataset_dict[dataset_name].define_dataset_base(data_input, False, shutdown_jitter)

        if mc.Tref:
            mc.dataset_dict[dataset_name].common_Tref(mc.Tref)
        else:
            mc.Tref = mc.dataset_dict[dataset_name].Tref

        if 'boundaries' in dataset_conf:
            bound_conf = dataset_conf['boundaries']
            for var in bound_conf:
                mc.dataset_dict[dataset_name].bounds[var] = np.asarray(bound_conf[var], dtype=np.double)

        if 'spaces' in dataset_conf:
            space_conf = dataset_conf['spaces']
            for var in space_conf:
                mc.dataset_dict[dataset_name].spaces[var] = space_conf[var]

        if 'starts' in dataset_conf:
            mc.starting_point_flag = True
            starts_conf = dataset_conf['starts']
            for var in starts_conf:
                mc.dataset_dict[dataset_name].starts[var] = np.asarray(starts_conf[var], dtype=np.double)

        mc.dataset_dict[dataset_name].update_bounds_spaces_priors_starts()

    for model_name, model_conf in conf_common.iteritems():

        if not isinstance(model_name, str):
            model_name = repr(model_name)

        if model_name == 'planets':

            for planet_name, planet_conf in model_conf.iteritems():

                if not isinstance(planet_name, str):
                    planet_name = repr(planet_name)

                mc.common_models[planet_name] = define_common_type_to_class['planets'](planet_name)

                bounds_space_priors_starts_fixed(mc, mc.common_models[planet_name], planet_conf)

                if 'orbit' in planet_conf:
                    mc.planet_dict[planet_name] = planet_conf['orbit']

                    if planet_conf['orbit'] == 'circular':
                        mc.common_models[planet_name].fix_list['e'] = np.asarray([0.000, 0.0000], dtype=np.double)
                        mc.common_models[planet_name].fix_list['o'] = np.asarray([np.pi/2., 0.0000], dtype=np.double)
                    if planet_conf['orbit'] == 'dynamical':
                        mc.dynamical_dict[planet_name] = True
                else:
                    mc.planet_dict[planet_name] = 'keplerian'

                if 'parametrization' in planet_conf:
                    if planet_conf['parametrization'] in mc.common_models[planet_name].parametrization_list:
                        mc.common_models[planet_name].parametrization = planet_conf['parametrization']
                    else:
                        print 'Orbital parametrization not recognized, switching to ', \
                            mc.common_models[planet_name].parametrization, ' for planet ', planet_name
                else:
                    print 'Using standard orbital parametrization ', mc.common_models[planet_name].parametrization, \
                        ' for planet ', planet_name

        else:
            if 'type' in model_conf:
                model_type = model_conf['type']
            elif 'kind' in model_conf:
                model_type = model_conf['kind']
            else:
                model_type = model_name

            mc.common_models[model_name] = define_common_type_to_class[model_type](model_name)
            bounds_space_priors_starts_fixed(mc, mc.common_models[model_name], model_conf)

            """ Automatic detection of common models without dataset-specific parameters"""
            for dataset in mc.dataset_dict.itervalues():
                if model_name in dataset.models and not (model_name in conf_models):
                    conf_models[model_name] = {'common': model_name}


    """ Check if there is any planet that requires dynamical computations"""
    if mc.dynamical_dict:
        mc.dynamical_model = DynamicalIntegrator()

    for model_name, model_conf in conf_models.iteritems():

        if not isinstance(model_name, str):
            model_name = repr(model_name)

        """ Check if the keplerian approximation must be used for this dataset even if the planet has the dynamical flag"""

        if 'keplerian_approximation' in model_conf:
            keplerian_approximation = model_conf['keplerian_approximation']
            print 'Using Keplerian approximation'
        else:
            keplerian_approximation = False

        if 'type' in model_conf:
            model_type = model_conf['type']
        elif 'kind' in model_conf:
            model_type = model_conf['kind']
        else:
            model_type = model_name

        if model_type == 'radial_velocities' or model_type == 'rv_planets' or model_type == 'batman_transit':

            """ radial_velocities and transits are just wrappers for the planets to be actually included in the model, so we
                substitue it with the individual planets in the list"""

            try:
                planet_list = np.atleast_1d(model_conf['planets']).tolist()
            except:
                planet_list = np.atleast_1d(model_conf['common']).tolist()

            model_name_expanded = [model_name + '_' + pl_name for pl_name in planet_list]
            """ Let's avoid some dumb user using the planet names to name the models"""

            for dataset in mc.dataset_dict.itervalues():
                if model_name in dataset.models:
                    dataset.models.remove(model_name)
                    dataset.models.extend(model_name_expanded)

                    if len(list(set(planet_list) & set(mc.dynamical_dict))) and not keplerian_approximation:
                        dataset.dynamical = True

            for model_name_exp, planet_name in zip(model_name_expanded, planet_list):

                try:
                    """ This snippet will work only for RV class"""
                    if keplerian_approximation:
                        """ override default model from input file, to apply keplerian approximation """
                        mc.models[model_name_exp] = \
                            define_type_to_class[model_type]['keplerian_mass'](model_name_exp, planet_name)
                    else:
                        mc.models[model_name_exp] = \
                            define_type_to_class[model_type][mc.planet_dict[planet_name]](model_name_exp, planet_name)

                except:
                    mc.models[model_name_exp] = \
                            define_type_to_class[model_type](model_name_exp, planet_name)

                    mc.models[model_name_exp].model_conf = model_conf.copy()

                for dataset_name in list(set(model_name_exp) & set(mc.dataset_dict)):
                    bounds_space_priors_starts_fixed(mc, mc.models[model_name_exp], model_conf[dataset_name],
                                                   dataset_1=dataset_name)

        elif model_type == 'transit_time':
            """ Only one planet for each file with transit times... mixing them would cause HELL"""

            try:
                planet_name = np.atleast_1d(model_conf['planet']).tolist()[0]
            except:
                planet_name = np.atleast_1d(model_conf['common']).tolist()[0]

            if keplerian_approximation:
                """ override default model from input file, to apply keplerian approximation """
                mc.models[model_name] = \
                    define_type_to_class[model_type]['keplerian'](model_name_exp, planet_name)
            else:
                mc.models[model_name] = \
                    define_type_to_class[model_type][mc.planet_dict[planet_name]](model_name, planet_name)

            #mc.models[model_name] = \
            #        define_type_to_class[model_type][mc.planet_dict[planet_name]](model_name, planet_name)

            """  CHECK THIS!!!!
            bounds_space_priors_starts_fixed(mc, mc.models[model_name], model_conf)
            """
            for dataset_name, dataset in mc.dataset_dict.iteritems():
                if planet_name in mc.dynamical_dict and model_name in dataset.models and not keplerian_approximation:
                    dataset.planet_name = planet_name
                    dataset.dynamical = True
                    mc.dynamical_t0_dict[planet_name] = dataset_name

        elif model_type == 'correlation_singledataset':
            mc.models[model_name] = \
                    define_type_to_class[model_type](model_name, None)
            mc.models[model_name].model_conf = model_conf.copy()
            bounds_space_priors_starts_fixed(mc, mc.models[model_name], model_conf, dataset_1=model_conf['reference'])

        else:

            try:
                mc.models[model_name] = \
                    define_type_to_class[model_type](model_name, model_conf['common'])
            except:
                mc.models[model_name] = \
                    define_type_to_class[model_type](model_name, None)

            mc.models[model_name].model_conf = model_conf.copy()

            #for dataset_name in list(set(model_conf) & set(mc.dataset_dict)):
            #    bounds_space_priors_starts_fixed(mc, mc.models[model_name], model_conf[dataset_name], dataset_1=dataset_name)

            """ Using default boundaries if one dataset is missing"""
            if not mc.models[model_name].list_pams_dataset:
                continue
            for dataset_name, dataset in mc.dataset_dict.iteritems():
                if model_name in dataset.models:

                    if dataset_name not in model_conf:
                        model_conf[dataset_name] = {}
                        print
                        print '*********************************** WARNING *********************************** '
                        print 'Using default boundaries for dataset-specific parameters '
                        print 'for model: ', model_name, ' dataset: ', dataset_name
                        print

                    bounds_space_priors_starts_fixed(mc, mc.models[model_name], model_conf[dataset_name],
                                                   dataset_1=dataset_name)

                #mc.models[model_name].setup_dataset(mc.dataset_dict[dataset_name])

    if 'Tref' in conf_parameters:
        mc.Tref = np.asarray(conf_parameters['Tref'])
        for dataset_name in mc.dataset_dict:
            mc.dataset_dict[dataset_name].common_Tref(mc.Tref)

    if 'star_mass' in conf_parameters:
        mc.star_mass = np.asarray(conf_parameters['star_mass'][:], dtype=np.double)
    if 'star_radius' in config_in:
        mc.star_radius = np.asarray(conf_parameters['star_radius'][:], dtype=np.double)

    if 'dynamical_integrator' in conf_solver:
        mc.dynamical_model.dynamical_integrator = conf_solver['dynamical_integrator']

    if 'pyde' in conf_solver and hasattr(mc, 'pyde_parameters'):
        conf = conf_solver['pyde']

        if 'ngen' in conf:
            mc.pyde_parameters['ngen'] = np.asarray(conf['ngen'], dtype=np.double)

        if 'npop_mult' in conf:
            mc.pyde_parameters['npop_mult'] = np.asarray(conf['npop_mult'], dtype=np.int64)

        if 'shutdown_jitter' in conf:
            mc.pyde_parameters['shutdown_jitter'] = np.asarray(conf['shutdown_jitter'], dtype=bool)

        if 'include_priors' in conf:
            mc.include_priors = np.asarray(conf['include_priors'], dtype=bool)

    if 'emcee' in conf_solver and hasattr(mc, 'emcee_parameters'):
        conf = conf_solver['emcee']

        if 'multirun' in conf:
            mc.emcee_parameters['multirun'] = np.asarray(conf['multirun'], dtype=np.int64)

        if 'multirun_iter' in conf:
            mc.emcee_parameters['multirun_iter'] = np.asarray(conf['multirun_iter'], dtype=np.int64)

        if 'nsave' in conf:
            mc.emcee_parameters['nsave'] = np.asarray(conf['nsave'], dtype=np.double)

        if 'nsteps' in conf:
            mc.emcee_parameters['nsteps'] = np.asarray(conf['nsteps'], dtype=np.int64)

        if 'nburn' in conf:
            mc.emcee_parameters['nburn'] = np.asarray(conf['nburn'], dtype=np.int64)

        if 'npop_mult' in conf:
            mc.emcee_parameters['npop_mult'] = np.asarray(conf['npop_mult'], dtype=np.int64)

        if 'thin' in conf:
            mc.emcee_parameters['thin'] = np.asarray(conf['thin'], dtype=np.int64)

        if 'shutdown_jitter' in conf:
            mc.emcee_parameters['shutdown_jitter'] = np.asarray(conf['shutdown_jitter'], dtype=bool)

        if 'include_priors' in conf:
            mc.include_priors = np.asarray(conf['include_priors'], dtype=bool)

    if 'nested_sampling' in conf_solver and hasattr(mc, 'nested_sampling_parameters'):
        conf = conf_solver['nested_sampling']

        for key_name, key_value in conf.items():
            mc.nested_sampling_parameters[key_name] = key_value

        if 'include_priors' in conf:
            mc.include_priors = np.asarray(conf['include_priors'], dtype=bool)

    if 'recenter_bounds' in conf_solver:
        """
        required to avoid a small bug in the code
        if the dispersion of PyDE walkers around the median value is too broad,
        then emcee walkers will start outside the bounds, causing an error
        """
        mc.recenter_bounds_flag = conf_solver['recenter_bounds']

    if 'include_priors' in conf_solver:
        mc.include_priors = np.asarray(conf_solver['include_priors'], dtype=bool)

    if 'use_threading_pool' in conf_solver:
        mc.use_threading_pool = conf_solver['use_threading_pool']


def bounds_space_priors_starts_fixed(mc, model_obj, conf, dataset_1=None, dataset_2=None, add_var_name =''):
    # type: (object, object, object, object, object, object) -> object

    if dataset_1 is None:
        if 'boundaries' in conf:
            bound_conf = conf['boundaries']
            for var in bound_conf:
                model_obj.bounds[add_var_name+var] = np.asarray(bound_conf[var], dtype=np.double)

        if 'spaces' in conf:
            space_conf = conf['spaces']
            for var in space_conf:
                model_obj.spaces[add_var_name+var] =space_conf[var]

        if 'priors' in conf:
            prior_conf = conf['priors']
            for var in prior_conf:
                prior_pams = np.atleast_1d(prior_conf[var])
                model_obj.prior_kind[add_var_name+var] = prior_pams[0]
                if np.size(prior_pams) > 1:
                    model_obj.prior_pams[add_var_name+var] = np.asarray(prior_pams[1:], dtype=np.double)
                else:
                    model_obj.prior_pams[add_var_name+var] = np.asarray([0.00], dtype=np.double)

        if 'starts' in conf:
            mc.starting_point_flag = True
            starts_conf = conf['starts']
            for var in starts_conf:
                model_obj.starts[add_var_name+var] = np.asarray(starts_conf[var], dtype=np.double)

        if 'fixed' in conf:
            fixed_conf = conf['fixed']
            for var in fixed_conf:
                model_obj.fix_list[add_var_name+var] = get_2darray_from_val(fixed_conf[var])

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
                model_obj.bounds[dataset_1][add_var_name+var] = np.asarray(bound_conf[var], dtype=np.double)

        if 'spaces' in conf:
            space_conf = conf['spaces']
            for var in space_conf:
                model_obj.spaces[dataset_1][add_var_name+var] = space_conf[var]

        if 'priors' in conf:
            prior_conf = conf['priors']
            for var in prior_conf:
                model_obj.prior_kind[dataset_1][add_var_name+var] = prior_conf[var][0]
                model_obj.prior_pams[dataset_1][add_var_name+var] = np.asarray(prior_conf[var][1:], dtype=np.double)

        if 'starts' in conf:
            mc.starting_point_flag = True
            starts_conf = conf['starts']
            for var in starts_conf:
                print dataset_1, add_var_name+var, starts_conf[var]
                model_obj.starts[dataset_1][add_var_name+var] = np.asarray(starts_conf[var], dtype=np.double)

        if 'fixed' in conf:
            fixed_conf = conf['fixed']
            for var in fixed_conf:
                model_obj.fix_list[dataset_1][add_var_name+var] = get_2darray_from_val(fixed_conf[var])

    else:

        if 'boundaries' in conf:
            bound_conf = conf['boundaries']
            for var in bound_conf:
                model_obj.bounds[dataset_1][dataset_2][add_var_name + var] = \
                    np.asarray(bound_conf[var], dtype=np.double)

        if 'spaces' in conf:
            space_conf = conf['spaces']
            for var in space_conf:
                model_obj.spaces[dataset_1][dataset_2][add_var_name + var] = space_conf[var]

        if 'priors' in conf:
            prior_conf = conf['priors']
            for var in prior_conf:
                model_obj.prior_kind[dataset_1][dataset_2][add_var_name + var] = prior_conf[var][0]
                model_obj.prior_pams[dataset_1][dataset_2][add_var_name + var] = \
                    np.asarray(prior_conf[var][1:], dtype=np.double)

        if 'starts' in conf:
            mc.starting_point_flag = True
            starts_conf = conf['starts']
            for var in starts_conf:
                model_obj.starts[dataset_1][dataset_2][add_var_name + var] = \
                    np.asarray(starts_conf[var], dtype=np.double)

        if 'fixed' in conf:
            fixed_conf = conf['fixed']
            for var in fixed_conf:
                model_obj.fix_list[dataset_1][dataset_2][add_var_name + var] = get_2darray_from_val(fixed_conf[var])
    return
