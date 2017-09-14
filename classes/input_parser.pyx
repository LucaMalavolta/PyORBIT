from common import *
from dataset import *
from abstract_common import CommonPlanets, CommonActivity
from radial_velocities import RVkeplerian, RVdynamical, TransitTimeKeplerian, TransitTimeDynamical, DynamicalIntegrator
from gp_semiperiodic_activity import GaussianProcess_QuasiPeriodicActivity
#TD from gaussian import GaussianProcess_QuasiPeriodicActivity
#TD from curvature import CurvatureCommonVariables
#TD from correlations import CorrelationsCommonVariables
#TD from sinusoids import SinusoidsCommonVariables

define_common_type_to_class = {
    'planets': CommonPlanets,
    'activity': CommonActivity
    #TD 'gp_quasiperiodic': GaussianProcess_QuasiPeriodicActivity,
    #TD 'curvature': CurvatureCommonVariables,
    #TD 'correlation': CorrelationsCommonVariables
}

define_type_to_class = {
    'radial_velocities': {'keplerian': RVkeplerian,
                          'dynamical': RVdynamical},
    'transit_time': {'keplerian': TransitTimeKeplerian,
                     'dynamical': TransitTimeDynamical},
    'gp_quasiperiodic': GaussianProcess_QuasiPeriodicActivity
    #TD 'gp_quasiperiodic': GaussianProcess_QuasiPeriodicActivity,
    #TD 'curvature': CurvatureCommonVariables,
    #TD 'correlation': CorrelationsCommonVariables
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


def pars_input(config_in, mc):

    mc.output_name = config_in['output']

    conf_inputs = config_in['inputs']
    conf_models = config_in['models']
    conf_common = config_in['common']
    conf_parameters = config_in['parameters']
    conf_solver = config_in['solver']

    for counter in conf_inputs:
        dataset_name = conf_inputs[counter]['file']

        """ The keyword in dataset_dict and the name assigned internally to the databes must be the same
            or everything will fall apart """
        mc.dataset_dict[dataset_name] = Dataset(conf_inputs[counter]['kind'], dataset_name, conf_inputs[counter]['models'])

        mc.dataset_index[counter] = dataset_name

        if counter == 0:
            mc.Tref = mc.dataset_dict[dataset_name].Tref
        else:
            mc.dataset_dict[dataset_name].common_Tref(mc.Tref)

        if 'name' in conf_inputs[counter]:
            mc.dataset_dict[dataset_name].name = conf_inputs[counter]['name']
        else:
            mc.dataset_dict[dataset_name].name = repr(counter)

        if 'boundaries' in conf_inputs[counter]:
            bound_conf = conf_inputs[counter]['boundaries']
            for var in bound_conf:
                mc.dataset_dict[dataset_name].bounds[var] = np.asarray(bound_conf[var], dtype=np.double)

        if 'starts' in conf_inputs[counter]:
            mc.starting_point_flag = True
            starts_conf = conf_inputs[counter]['starts']
            for var in starts_conf:
                mc.dataset_dict[dataset_name].starts[var] = np.asarray(starts_conf[var], dtype=np.double)

        mc.dataset_dict[dataset_name].update_priors_starts_bounds()

    for model_name, model_conf in conf_common.iteritems():

        if model_name == 'planets':

            for planet_name, planet_conf in model_conf.iteritems():

                mc.common_models[planet_name] = define_common_type_to_class['planets'](planet_name)
                boundaries_fixed_priors_stars(mc, mc.common_models[planet_name], planet_conf)

                if 'orbit' in planet_conf:
                    mc.planet_dict[planet_name] = planet_conf['orbit']
                else:
                    mc.planet_dict[planet_name] = 'keplerian'

        else:
            if 'type' in model_conf:
                model_type = model_conf['type']
            elif 'kind' in model_conf:
                model_type = model_conf['kind']
            else:
                model_type = model_name

            mc.common_models[model_name] = define_common_type_to_class[model_type](model_name)
            boundaries_fixed_priors_stars(mc, mc.common_models[model_name], model_conf)


    """ Check if there is any planet that requires dynamical computations"""
    for planet_name, orbit in mc.planet_dict.iteritems():
        if orbit == 'dynamical':
            mc.dynamical_dict[planet_name] = True

        if mc.dynamical_dict:
            mc.dynamical_model = DynamicalIntegrator()

    for model_name, model_conf in conf_models.iteritems():

        if 'type' in model_conf:
            model_type = model_conf['type']
        elif 'kind' in model_conf:
            model_type = model_conf['kind']
        else:
            model_type = model_name

        if model_type == 'radial_velocities':

            """ radial_velocities is just a wrapper for the planets to be actually included in the model, so we
                substitue it with the individual planets in the list"""

            model_name_expanded = [model_name + '_' + pl_name for pl_name in model_conf['planets']]
            """ Let's avoid some dumb user using the planet names to name the models"""

            for dataset in mc.dataset_dict.itervalues():
                if model_name in dataset.models:
                    dataset.models.remove(model_name)
                    dataset.models.extend(model_name_expanded)

                    if len(list(set(model['planets']) & set(mc.dynamical_dict))):
                        dataset.dynamical = True

            for model_name_exp, planet_name in zip(model_name_expanded, model_conf['planets']):

                mc.models[model_name_exp] = \
                    define_type_to_class[model_type][mc.planet_dict[planet_name]](model_name_exp, planet_name)
                boundaries_fixed_priors_stars(mc, mc.models[model_name_exp], model_conf[planet_name])

        elif model_type == 'transit_time':
            """ Only one planet for each file with transit times... mixing them would cause HELL"""

            planet_name = model_conf['planet']
            mc.models[model_name] = \
                    define_type_to_class[model_type][mc.planet_dict[planet_name]](model_name, planet_name)
            boundaries_fixed_priors_stars(mc, mc.models[model_name], model_conf)

            for dataset_name, dataset in mc.dataset_dict.iteritems():
                if planet_name in mc.dynamical_dict and planet_name in dataset.models:
                    dataset.planet_name = planet_name
                    dataset.dynamical = True
                    mc.t0_dict[planet_name] = dataset_name

        else:

            mc.models[model_name] = \
                    define_type_to_class[model_type](model_name, model_conf['common'])
            boundaries_fixed_priors_stars(mc, mc.models[model_name], model_conf)

    if 'Tref' in conf_parameters:
        mc.Tref = np.asarray(conf_parameters['Tref'])
        for dataset_name in mc.dataset_dict:
            mc.dataset_dict[dataset_name].common_Tref(mc.Tref)

    if 'star_mass' in conf_parameters:
        mc.star_mass = np.asarray(conf_parameters['star_mass'][:], dtype=np.double)
    if 'star_radius' in config_in:
        mc.star_radius = np.asarray(conf_parameters['star_radius'][:], dtype=np.double)

    if 'dynamical_integrator' in conf_solver:
        mc.dynamical_model.dynamical_integrator = conf_parameters['dynamical_integrator']

    if 'pyde' in conf_solver:
        conf = conf_solver['pyde']

        if 'ngen' in conf:
            mc.pyde_parameters['ngen'] = np.asarray(conf['ngen'], dtype=np.double)

        if 'npop_mult' in conf:
            mc.pyde_parameters['npop_mult'] = np.asarray(conf['npop_mult'], dtype=np.int64)

    if 'emcee' in conf_solver:
        conf = conf_solver['emcee']

        if 'multirun' in conf:
            mc.emcee_parameters['multirun'] = np.asarray(conf['multirun'], dtype=np.int64)

        if 'MultiRun_iter' in conf:
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

    if 'polychord' in conf_solver:
        conf = conf_solver['polychord']

        if 'nlive' in conf:
            mc.polychord_parameters['nlive'] = np.asarray(conf['nlive'], dtype=np.int64)

        if 'nlive_mult' in conf:
            mc.polychord_parameters['nlive_mult'] = np.asarray(conf['nlive_mult'], dtype=np.int64)

        if 'num_repeats_mult' in conf:
            mc.polychord_parameters['num_repeats_mult'] = np.asarray(conf['num_repeats_mult'], dtype=np.int64)

        if 'feedback' in conf:
            mc.polychord_parameters['feedback'] = np.asarray(conf['feedback'], dtype=np.int64)

        if 'precision_criterion' in conf:
            mc.polychord_parameters['precision_criterion'] = np.asarray(conf['precision_criterion'], dtype=np.double)

        if 'max_ndead' in conf:
            mc.polychord_parameters['max_ndead'] = np.asarray(conf['max_ndead'], dtype=np.int64)

        if 'boost_posterior' in conf:
            mc.polychord_parameters['boost_posterior'] = np.asarray(conf['boost_posterior'], dtype=np.double)

        if 'read_resume' in conf:
            mc.polychord_parameters['read_resume'] = np.asarray(conf['read_resume'], dtype=bool)

        if 'base_dir' in conf:
            mc.polychord_parameters['base_dir'] = np.asarray(conf['base_dir'], dtype=str)

        #if 'file_root' in conf:
        #    mc.polychord_parameters['file_root'] = np.asarray(conf['file_root'], dtype=str)

        if 'shutdown_jitter' in conf:
            mc.polychord_parameters['shutdown_jitter'] = np.asarray(conf['shutdown_jitter'], dtype=bool)

    if 'recenter_bounds' in conf_solver:
        """ 
        required to avoid a small bug in the code
        if the dispersion of PyDE walkers around the median value is too broad,
        then emcee walkers will start outside the bounds, causing an error
        """
        mc.recenter_bounds_flag = conf_solver['recenter_bounds']





def dump():

    for model in conf_models:
        if 'type' in conf_models[model]:
            mc.models[model] = define_type_to_class[conf_models[model]['type']](model)
        else:
            mc.models[model] = define_type_to_class[model](model)

        if mc.models[model].model_class is 'planets':
            conf = conf_models[model]
            for counter in conf:
                planet_name = 'Planet_' + repr(counter)
                planet_conf = conf[counter]

                print mc.models[model]
                mc.models[model].add_planet(planet_name)

                if 'orbit' in planet_conf:
                    # By default orbits are keplerians
                    if planet_conf['orbit'] == 'circular':
                        mc.models['planets'].switch_to_circular(planet_name)
                    if planet_conf['orbit'] == 'dynamical':
                        mc.models['planets'].switch_to_dynamical(planet_name)

                if 'transit' in planet_conf:
                    if planet_conf['transit']:
                        mc.models['planets'].switch_on_transit(planet_name)

                if 'inclination' in planet_conf:
                    mc.models['planets'].inclination[planet_name] = planet_conf['inclination']

                if 'radius' in planet_conf:
                    mc.models['planets'].radius[planet_name] = planet_conf['radius']

                boundaries_fixed_priors_stars(mc, model, planet_conf, planet_name)

        if mc.models[model].model_class is 'correlation':

            conf = conf_models[model]
            correlation_common = False

            """ When including the specific values for each dataset association, the existence of common variables must have
                been already checked, just to avoid problems to those distracted users that include the Common block after
                the dataset-specific ones
            """
            for counter_ref in conf:
                if counter_ref is 'common':
                    correlation_common = True

            for counter_ref in conf:
                if counter_ref is 'common' or counter_ref is 'type':
                    continue

                dataname_ref = mc.dataset_index[counter_ref]
                mc.models[model].add_dataset(dataname_ref)

                for counter_asc in conf[counter_ref]:
                    dataname_asc = mc.dataset_index[counter_asc]
                    mc.models[model].add_associated_dataset(mc, dataname_ref, dataname_asc)
                    free_zeropoint = False

                    """ Apply common settings (if present) 
                        before overriding them with the specific values (if provided)
                    """
                    if correlation_common:
                        common_conf = conf[counter_ref]['common']

                        """ The xero point of the correlation plot is already included as a free parameter
                         as the offset of the associated dataset, so it is disabled by default. However there may be 
                         situations in which it is still needed to have it as a free parameter, so this option is given"""

                        if 'free_zeropoint' not in common_conf or common_conf['Free_ZeroPoint'] is False:
                            mc.models[model].fix_list[dataname_ref][dataname_asc]['correlation_0'] = 0.0000
                            free_zeropoint = True

                        """ By default the origin of the x axis of the independent parameter (usually an activity indicator)
                        is set to the median value of the parameter. The user has te option to specify a value"""
                        if 'abscissa_zero' in common_conf:
                            mc.models[model].x_zero[dataname_ref][dataname_asc] = common_conf['abscissa_zero']

                        if 'order' in common_conf:
                            mc.models[model].order[dataname_ref][dataname_asc] = \
                                np.asarray(common_conf['order'], dtype=np.int64)

                        boundaries_fixed_priors_stars(mc, model, common_conf,
                                                      dataname_ref, dataname_asc, 'correlation_')

                    if free_zeropoint is False and \
                        ('free_zeropoint' not in conf[counter_ref][counter_asc] or
                         conf[counter_ref][counter_asc]['free_zeropoint'] is False):
                        mc.models[model].fix_list[dataname_ref][dataname_asc]['correlation_0'] = 0.0000

                    if 'abscissa_zero' in conf[counter_ref][counter_asc]:
                        mc.models[model].x_zero[dataname_ref][dataname_asc] = common_conf['abscissa_zero']

                    if 'order' in conf[counter_ref][counter_asc]:
                        mc.models[model].order[dataname_ref][dataname_asc] = \
                            np.asarray(conf[counter_ref][counter_asc]['order'], dtype=np.int64)

                    boundaries_fixed_priors_stars(mc, model, conf[counter_ref][counter_asc],
                                                  dataname_ref, dataname_asc, 'correlation_')

        if mc.models[model].model_class is 'gaussian_process':
            conf = conf_models[model]

            for name_ref in conf:
                if name_ref == 'type':
                    continue

                print name_ref
                if name_ref == 'common':
                    dataset_name = mc.models[model].common_ref
                else:
                    dataset_name = mc.dataset_index[name_ref]

                mc.models[model].add_dataset(dataset_name)

                boundaries_fixed_priors_stars(mc, model, conf[name_ref], dataset_name)

        if mc.models[model].model_class is 'curvature':

            conf = conf_models[model]

            if 'order' in conf:
                mc.models[model].order = np.asarray(conf['order'], dtype=np.int64)

            boundaries_fixed_priors_stars(mc, model, conf)

                #if 'Sinusoids' in config_in:
    #    conf = config_in['Sinusoids']
    #    mc.scv.Prot_bounds = np.asarray(conf['Prot'], dtype=np.double)
    #    if 'Priors' in conf:
    #        mc.scv.prior_kind[var] = conf['Priors'][var][0]
    #        mc.scv.prior_pams[var] = np.asarray(conf['Priors'][var][1:], dtype=np.double)
    #
    #    for counter in conf['Seasons']:
    #        mc.scv.add_season_range(np.asarray(conf['Seasons'][counter][:2], dtype=np.double),
    #                                conf['Seasons'][counter][2:])

            # if 'Phase_dataset' in conf:
            #    # Additional activity indicators associated to RVs must the same order and number of sinusoids,
            #    #
            #    mc.scv.phase = np.asarray(conf['Phase_dataset'], dtype=np.int64)

    if 'Tref' in conf_parameters:
        mc.Tref = np.asarray(conf_parameters['Tref'])
        for dataset_name in mc.dataset_dict:
            mc.dataset_dict[dataset_name].common_Tref(mc.Tref)

    if 'star_mass' in conf_parameters:
        mc.star_mass = np.asarray(conf_parameters['star_mass'][:], dtype=np.double)
    if 'star_radius' in config_in:
        mc.star_radius = np.asarray(conf_parameters['star_radius'][:], dtype=np.double)

    if 'dynamical_integrator' in conf_solver:
        mc.dynamical_model.dynamical_integrator = conf_parameters['dynamical_integrator']

    if 'pyde' in config_in:
        conf = config_in['pyde']

        if 'ngen' in conf:
            mc.pyde_parameters['ngen'] = np.asarray(conf['ngen'], dtype=np.double)

        if 'npop_mult' in conf:
            mc.pyde_parameters['npop_mult'] = np.asarray(conf['npop_mult'], dtype=np.int64)

    if 'emcee' in conf_solver:
        conf = conf_solver['emcee']

        if 'multirun' in conf:
            mc.emcee_parameters['multirun'] = np.asarray(conf['multirun'], dtype=np.int64)

        if 'MultiRun_iter' in conf:
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

    if 'polychord' in conf_solver:
        conf = conf_solver['polychord']

        if 'nlive' in conf:
            mc.polychord_parameters['nlive'] = np.asarray(conf['nlive'], dtype=np.int64)

        if 'nlive_mult' in conf:
            mc.polychord_parameters['nlive_mult'] = np.asarray(conf['nlive_mult'], dtype=np.int64)

        if 'num_repeats_mult' in conf:
            mc.polychord_parameters['num_repeats_mult'] = np.asarray(conf['num_repeats_mult'], dtype=np.int64)

        if 'feedback' in conf:
            mc.polychord_parameters['feedback'] = np.asarray(conf['feedback'], dtype=np.int64)

        if 'precision_criterion' in conf:
            mc.polychord_parameters['precision_criterion'] = np.asarray(conf['precision_criterion'], dtype=np.double)

        if 'max_ndead' in conf:
            mc.polychord_parameters['max_ndead'] = np.asarray(conf['max_ndead'], dtype=np.int64)

        if 'boost_posterior' in conf:
            mc.polychord_parameters['boost_posterior'] = np.asarray(conf['boost_posterior'], dtype=np.double)

        if 'read_resume' in conf:
            mc.polychord_parameters['read_resume'] = np.asarray(conf['read_resume'], dtype=bool)

        if 'base_dir' in conf:
            mc.polychord_parameters['base_dir'] = np.asarray(conf['base_dir'], dtype=str)

        #if 'file_root' in conf:
        #    mc.polychord_parameters['file_root'] = np.asarray(conf['file_root'], dtype=str)

        if 'shutdown_jitter' in conf:
            mc.polychord_parameters['shutdown_jitter'] = np.asarray(conf['shutdown_jitter'], dtype=bool)

    if 'recenter_bounds' in conf_solver:
        """ 
        required to avoid a small bug in the code
        if the dispersion of PyDE walkers around the median value is too broad,
        then emcee walkers will start outside the bounds, causing an error
        """
        mc.recenter_bounds_flag = conf_solver['recenter_bounds']


def boundaries_fixed_priors_stars(mc, model_obj, conf, dataset_1=None, dataset_2=None, add_var_name =''):
    # type: (object, object, object, object, object, object) -> object

    if dataset_1 is None:
        if 'boundaries' in conf:
            bound_conf = conf['boundaries']
            for var in bound_conf:
                model_obj.bounds[add_var_name+var] = np.asarray(bound_conf[var], dtype=np.double)

        if 'fixed' in conf:
            fixed_conf = conf['fixed']
            for var in fixed_conf:
                model_obj.fix_list[add_var_name+var] = np.asarray(fixed_conf[var], dtype=np.double)

        if 'priors' in conf:
            prior_conf = conf['priors']
            for var in prior_conf:
                model_obj.prior_kind[add_var_name+var] = prior_conf[var][0]
                model_obj.prior_pams[add_var_name+var] = np.asarray(prior_conf[var][1:], dtype=np.double)

        if 'starts' in conf:
            mc.starting_point_flag = True
            starts_conf = conf['starts']
            for var in starts_conf:
                model_obj.starts[add_var_name+var] = np.asarray(starts_conf[var], dtype=np.double)

    elif dataset_2 is None:
        if 'boundaries' in conf:
            bound_conf = conf['boundaries']
            for var in bound_conf:
                model_obj.bounds[dataset_1][add_var_name+var] = np.asarray(bound_conf[var], dtype=np.double)

        if 'fixed' in conf:
            fixed_conf = conf['fixed']
            for var in fixed_conf:
                model_obj.fix_list[dataset_1][add_var_name+var] = np.asarray(fixed_conf[var], dtype=np.double)

        if 'priors' in conf:
            prior_conf = conf['priors']
            for var in prior_conf:
                model_obj.prior_kind[dataset_1][add_var_name+var] = prior_conf[var][0]
                model_obj.prior_pams[dataset_1][add_var_name+var] = np.asarray(prior_conf[var][1:], dtype=np.double)

        if 'starts' in conf:
            mc.starting_point_flag = True
            starts_conf = conf['starts']
            for var in starts_conf:
                model_obj.starts[dataset_1][add_var_name+var] = np.asarray(starts_conf[var], dtype=np.double)

    else:

        if 'boundaries' in conf:
            bound_conf = conf['boundaries']
            for var in bound_conf:
                model_obj.bounds[dataset_1][dataset_2][add_var_name + var] = \
                    np.asarray(bound_conf[var], dtype=np.double)

        if 'fixed' in conf:
            fixed_conf = conf['fixed']
            for var in fixed_conf:
                model_obj.fix_list[dataset_1][dataset_2][add_var_name + var] = \
                    np.asarray(fixed_conf[var], dtype=np.double)

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

    return