from common import *
from dataset import *

def yaml_parser(file_conf, mc):
    stream = file(file_conf, 'r')
    config_in = yaml.load(stream)

    conf = config_in['Inputs']
    for counter in conf:
        print conf[counter]['Kind'], conf[counter]['File'], conf[counter]['Models']
        mc.dataset_list.append(Dataset(counter, conf[counter]['Kind'], conf[counter]['File'], conf[counter]['Models']))

        if counter == 0:
            mc.Tref = mc.dataset_list[0].Tref
        else:
            mc.dataset_list[counter].common_Tref(mc.Tref)

        if 'Boundaries' in conf[counter]:
            bound_conf = conf[counter]['Boundaries']
            for var in bound_conf:
                mc.dataset_list[counter].bounds[var] = np.asarray(bound_conf[var], dtype=np.double)

        if 'Starts' in conf[counter]:
            mc.starting_point_flag = True
            starts_conf = conf[counter]['Starts']
            for var in bound_conf:
                mc.dataset_list[counter].starts[var] = np.asarray(starts_conf[var], dtype=np.double)

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

            if 'Tcent' in planet_conf:
                mc.dataset_list.append(TransitCentralTimes(planet_name, planet_conf['Tcent']))
                mc.dataset_list[-1].common_Tref(mc.Tref)
                mc.t0_list[planet_name] = mc.dataset_list[-1]

            if 'Inclination' in planet_conf:
                mc.pcv.inclination[planet_name] = planet_conf['Inclination']
            if 'Radius' in planet_conf:
                mc.pcv.radius[planet_name] = planet_conf['Radius']

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
                dataset_name = mc.dataset_list[name_ref].name_ref
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
        for dataset in mc.dataset_list:
            dataset.common_Tref(mc.Tref)

    if 'MultiRun' in config_in['emcee']:
        mc.MultiRun = np.asarray(config_in['emcee']['MultiRun'], dtype=np.int64)

    if 'Ngen' in config_in['emcee']:
        mc.ngen = np.asarray(config_in['emcee']['Ngen'], dtype=np.double)

    if 'Nsave' in config_in['emcee']:
        mc.nsave = np.asarray(config_in['emcee']['Nsave'], dtype=np.double)

    if 'Nsteps' in config_in['emcee']:
        mc.nsteps = np.asarray(config_in['emcee']['Nsteps'], dtype=np.int64)

    if 'Nburn' in config_in['emcee']:
        mc.nburn = np.asarray(config_in['emcee']['Nburn'], dtype=np.int64)

    if 'Npop_mult' in config_in['emcee']:
        mc.npop_mult = np.asarray(config_in['emcee']['Npop_mult'], dtype=np.int64)

    if 'Thin' in config_in['emcee']:
        mc.thin = np.asarray(config_in['emcee']['Thin'], dtype=np.int64)

    if 'Recenter_Bounds' in config_in['emcee']:
        # required to avoid a small bug in the code
        # if the dispersion of PyDE walkers around the median value is too broad,
        # then emcee walkers will start outside the bounds, causing an error
        mc.recenter_bounds_flag = config_in['emcee']['Recenter_Bounds']

    if 'Star_Mass' in config_in:
        mc.star_mass = np.asarray(config_in['Star_Mass'][:], dtype=np.double)
    if 'Star_Radius' in config_in:
        mc.star_radius = np.asarray(config_in['Star_Radius'][:], dtype=np.double)

    if 'Dynamical_Integrator' in config_in:
        mc.pcv.dynamical_integrator = config_in['Dynamical_Integrator']

    mc.model_setup()
