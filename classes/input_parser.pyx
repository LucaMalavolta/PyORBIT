from common import *
from dataset import *


def yaml_parser(file_conf, mc):
    stream = file(file_conf, 'r')
    config_in = yaml.load(stream)

    conf = config_in['Inputs']
    for counter in conf:
        print conf[counter]['Kind'], conf[counter]['File'], conf[counter]['Models']
        dataset_name = conf[counter]['File']

        """ The keywird in dataset_dict and the name assigned internally to the databes must be the same
            or everything will fall apart """
        if 'Tcent' in conf[counter]['Kind']:
            planet_name = 'Planet_' + repr(conf[counter]['Planet'])
            mc.dataset_dict[dataset_name] = \
                TransitCentralTimes(counter, conf[counter]['Kind'], dataset_name, conf[counter]['Models'])
            mc.dataset_index[counter] = dataset_name
            mc.dataset_dict[dataset_name].set_planet(planet_name)
            mc.t0_list[planet_name] = mc.dataset_dict[dataset_name]
        else:
            mc.dataset_dict[dataset_name] = \
                Dataset(counter, conf[counter]['Kind'], dataset_name, conf[counter]['Models'])

        if counter == 0:
            mc.Tref = mc.dataset_dict[dataset_name].Tref
        else:
            mc.dataset_dict[dataset_name].common_Tref(mc.Tref)

        if 'Name' in conf[counter]:
            mc.dataset_dict[dataset_name].name = conf[counter]['Name']
        else:
            mc.dataset_dict[dataset_name].name = 'Unspecified'

        if 'Boundaries' in conf[counter]:
            bound_conf = conf[counter]['Boundaries']
            for var in bound_conf:
                mc.dataset_dict[dataset_name].bounds[var] = np.asarray(bound_conf[var], dtype=np.double)

        if 'Starts' in conf[counter]:
            mc.starting_point_flag = True
            starts_conf = conf[counter]['Starts']
            for var in starts_conf:
                mc.dataset_dict[dataset_name].starts[var] = np.asarray(starts_conf[var], dtype=np.double)

    mc.planet_name = config_in['Output']

    if 'Planets' in config_in:
        conf = config_in['Planets']
        for counter in conf:
            planet_name = 'Planet_' + repr(counter)
            planet_conf = conf[counter]
            mc.pcv.add_planet(planet_name)

            if 'Boundaries' in planet_conf:
                bound_conf = planet_conf['Boundaries']
                for var in bound_conf:
                    mc.pcv.bounds[planet_name][var] = np.asarray(bound_conf[var], dtype=np.double)

            if 'Fixed' in planet_conf:
                fixed_conf = planet_conf['Fixed']
                for var in fixed_conf:
                    mc.pcv.fix_list[planet_name][var] = np.asarray(fixed_conf[var], dtype=np.double)

            if 'Priors' in planet_conf:
                prior_conf = planet_conf['Priors']
                for var in prior_conf:
                    mc.pcv.prior_kind[planet_name][var] = prior_conf[var][0]
                    mc.pcv.prior_pams[planet_name][var] = np.asarray(prior_conf[var][1:], dtype=np.double)
                print mc.pcv.prior_kind
                print mc.pcv.prior_pams

            if 'Starts' in planet_conf:
                mc.starting_point_flag = True
                starts_conf = planet_conf['Starts']
                for var in starts_conf:
                    mc.pcv.starts[planet_name][var] = np.asarray(starts_conf[var], dtype=np.double)

            if 'Orbit' in planet_conf:
                # By default orbits are keplerians
                if planet_conf['Orbit'] == 'circular':
                    mc.pcv.switch_to_circular(planet_name)
                if planet_conf['Orbit'] == 'dynamical':
                    mc.pcv.switch_to_dynamical(planet_name)

            if 'Transit' in planet_conf:
                if planet_conf['Transit']:
                    mc.pcv.switch_on_transit(planet_name)

            if 'Inclination' in planet_conf:
                mc.pcv.inclination[planet_name] = planet_conf['Inclination']
            if 'Radius' in planet_conf:
                mc.pcv.radius[planet_name] = planet_conf['Radius']

    if 'Correlations' in config_in:
        conf = config_in['Correlations']
        for counter in conf:
            dataset_name = mc.dataset_index[counter]
            mc.cov.add_dataset(dataset_name)

            #### AAAAAA Fix this
            if 'Association' in conf[name_ref]:
                mc.dataset_list[name_ref].associate_to_other_dataset(mc.dataset_list[conf[name_ref]['Association']])

            if 'Order' in conf[name_ref]:
                mc.cov.order[mc.dataset_list[name_ref].name_ref] = np.asarray(conf[name_ref]['Order'], dtype=np.int64)

            if 'Boundaries' in conf[name_ref]:
                bound_conf = conf[name_ref]['Boundaries']
                for var in bound_conf:
                    mc.cov.bounds[name_ref]['correlation_' + var] = np.asarray(bound_conf[var], dtype=np.double)

            if 'Fixed' in conf[name_ref]:
                fixed_conf = conf[name_ref]['Fixed']
                for var in fixed_conf:
                    mc.cov.fix_list[name_ref]['correlation_' + var] = np.asarray(fixed_conf[var], dtype=np.double)

            if 'Priors' in conf[name_ref]:
                prior_conf = conf[name_ref]['Priors']
                for var in prior_conf:
                    mc.cov.prior_kind[name_ref]['correlation_' + var] = prior_conf[var][0]
                    mc.cov.prior_pams[name_ref]['correlation_' + var] = np.asarray(prior_conf[var][1:], dtype=np.double)

            if 'Starts' in conf[name_ref]:
                mc.starting_point_flag = True
                starts_conf = conf[name_ref]['Starts']
                for var in starts_conf:
                    mc.cov.starts[name_ref]['correlation_' + var] = np.asarray(starts_conf[var], dtype=np.double)

    if 'Sinusoids' in config_in:
        conf = config_in['Sinusoids']
        mc.scv.Prot_bounds = np.asarray(conf['Prot'], dtype=np.double)
        if 'Priors' in conf:
            mc.scv.prior_kind[var] = conf['Priors'][var][0]
            mc.scv.prior_pams[var] = np.asarray(conf['Priors'][var][1:], dtype=np.double)

        for counter in conf['Seasons']:
            mc.scv.add_season_range(np.asarray(conf['Seasons'][counter][:2], dtype=np.double),
                                    conf['Seasons'][counter][2:])

            # if 'Phase_dataset' in conf:
            #    # Additional activity indicators associated to RVs must the same order and number of sinusoids,
            #    #
            #    mc.scv.phase = np.asarray(conf['Phase_dataset'], dtype=np.int64)

    if 'Gaussian' in config_in:
        conf = config_in['Gaussian']
        for name_ref in conf:
            if name_ref == 'Common':
                dataset_name = 'Common'
            else:
                dataset_name = mc.dataset_index[name_ref]
            mc.gcv.add_dataset(dataset_name)

            if 'Boundaries' in conf[name_ref]:
                bound_conf = conf[name_ref]['Boundaries']
                for var in bound_conf:
                    mc.gcv.bounds[dataset_name][var] = np.asarray(bound_conf[var], dtype=np.double)

            if 'Fixed' in conf[name_ref]:
                fixed_conf = conf[name_ref]['Fixed']
                for var in fixed_conf:
                    mc.gcv.fix_list[dataset_name][var] = np.asarray(fixed_conf[var], dtype=np.double)

            if 'Priors' in conf[name_ref]:
                prior_conf = conf[name_ref]['Priors']
                for var in prior_conf:
                    mc.gcv.prior_kind[dataset_name][var] = prior_conf[var][0]
                    mc.gcv.prior_pams[dataset_name][var] = np.asarray(prior_conf[var][1:], dtype=np.double)

            if 'Starts' in conf[name_ref]:
                mc.starting_point_flag = True
                starts_conf = conf[name_ref]['Starts']
                for var in starts_conf:
                    mc.gcv.starts[dataset_name][var] = np.asarray(starts_conf[var], dtype=np.double)

    if 'Curvature' in config_in:
        conf = config_in['Curvature']

        if 'Order' in conf:
            mc.ccv.order = np.asarray(conf['Order'], dtype=np.int64)

        if 'Boundaries' in conf:
            bound_conf = conf['Boundaries']
            for var in bound_conf:
                mc.ccv.bounds[var] = np.asarray(bound_conf[var], dtype=np.double)

        if 'Fixed' in conf:
            fixed_conf = conf['Fixed']
            for var in fixed_conf:
                mc.ccv.fix_list[var] = np.asarray(fixed_conf[var], dtype=np.double)

        if 'Priors' in conf:
            prior_conf = conf['Priors']
            for var in prior_conf:
                mc.ccv.prior_kind[var] = prior_conf[var][0]
                mc.ccv.prior_pams[var] = np.asarray(prior_conf[var][1:], dtype=np.double)

        if 'Starts' in conf:
            mc.starting_point_flag = True
            starts_conf = conf['Starts']
            for var in starts_conf:
                mc.ccv.starts[var] = np.asarray(starts_conf[var], dtype=np.double)

    if 'Tref' in config_in:
        mc.Tref = np.asarray(config_in['Tref'])
        for dataset_name in mc.dataset_dict:
            mc.dataset_dict[dataset_name].common_Tref(mc.Tref)

    if 'pyDE' in config_in:
        conf = config_in['pyDE']

        if 'Ngen' in conf:
            mc.pyde_parameters['ngen'] = np.asarray(conf['Ngen'], dtype=np.double)

        if 'Npop_mult' in conf:
            mc.pyde_parameters['npop_mult'] = np.asarray(conf['Npop_mult'], dtype=np.int64)

    if 'emcee' in config_in:
        conf = config_in['emcee']

        if 'MultiRun' in conf:
            mc.emcee_parameters['MultiRun'] = np.asarray(conf['MultiRun'], dtype=np.int64)

        if 'MultiRun_iter' in conf:
            mc.emcee_parameters['MultiRun_iter'] = np.asarray(conf['MultiRun_iter'], dtype=np.int64)

        if 'Nsave' in conf:
            mc.emcee_parameters['nsave'] = np.asarray(conf['Nsave'], dtype=np.double)

        if 'Nsteps' in conf:
            mc.emcee_parameters['nsteps'] = np.asarray(conf['Nsteps'], dtype=np.int64)

        if 'Nburn' in conf:
            mc.emcee_parameters['nburn'] = np.asarray(conf['Nburn'], dtype=np.int64)

        if 'Npop_mult' in conf:
            mc.emcee_parameters['npop_mult'] = np.asarray(conf['Npop_mult'], dtype=np.int64)

        if 'Thin' in conf:
            mc.emcee_parameters['thin'] = np.asarray(conf['Thin'], dtype=np.int64)

    if 'PolyChord' in config_in:
        conf = config_in['PolyChord']

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

    if 'Recenter_Bounds' in config_in:
        # required to avoid a small bug in the code
        # if the dispersion of PyDE walkers around the median value is too broad,
        # then emcee walkers will start outside the bounds, causing an error
        mc.recenter_bounds_flag = config_in['Recenter_Bounds']

    if 'Star_Mass' in config_in:
        mc.star_mass = np.asarray(config_in['Star_Mass'][:], dtype=np.double)
    if 'Star_Radius' in config_in:
        mc.star_radius = np.asarray(config_in['Star_Radius'][:], dtype=np.double)

    if 'Dynamical_Integrator' in config_in:
        mc.dynamical_model.dynamical_integrator = config_in['Dynamical_Integrator']

