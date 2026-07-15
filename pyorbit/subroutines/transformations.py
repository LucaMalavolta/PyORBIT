from pyorbit.keywords_definitions import *

__all__ = ['_prepare_planet_parametrization',
           '_prepare_planet_scaled_semimajor_axis',
           '_prepare_planet_mass',
           '_prepare_stellar_mass',
           '_prepare_planet_inclination',
           '_prepare_planet_time_inferior_conjunction']


def _prepare_planet_parametrization(main_object, mc, **kwargs):

    """ Check the orbital parametrization of the planet and update
        the list of parameters accordingly """

    if mc.common_models[main_object.planet_ref].parametrization[:8] == 'Ford2006' \
        and mc.common_models[main_object.planet_ref].orbit != 'circular':
        main_object.list_pams_common.discard('e')
        main_object.list_pams_common.discard('omega')

        main_object.list_pams_common.update(['e_coso'])
        main_object.list_pams_common.update(['e_sino'])

    elif mc.common_models[main_object.planet_ref].parametrization[:8] != 'Standard' \
        and mc.common_models[main_object.planet_ref].orbit != 'circular':
            # 'Eastman2013' is the standard choice
        main_object.list_pams_common.discard('e')
        main_object.list_pams_common.discard('omega')

        main_object.list_pams_common.update(['sre_coso'])
        main_object.list_pams_common.update(['sre_sino'])


def _prepare_planet_scaled_semimajor_axis(main_object, mc, **kwargs):

    """ Check if the mass of the star is computed using a multivariate prior.
    This is a star property and not a planet one, but we need it defined here and now
    because we need to know how to compute the stellar density if
    it replaces the semimajor axis as a free parameter in the transit fit
    """

    try:
        multivariate_pams = mc.common_models[main_object.stellar_ref].multivariate_pams
    except AttributeError:
        multivariate_pams = []

    """ Default parametrization uses the stellar density and the impact
        parameter, it is possible to switch back to scaled semi-major axis and
        inclination respectively by activating the proper flag """

    if mc.common_models[main_object.planet_ref].use_scaled_semimajor_axis:
        """ a is the semi-major axis (in units of stellar radii) """
        main_object.list_pams_common.update(['a_Rs'])
        main_object.compute_scaled_semimajor_axis = False

    else:
        if 'mass' in multivariate_pams and 'radius' in multivariate_pams:
            main_object.list_pams_common.update(['mass'])
            main_object.list_pams_common.update(['radius'])
            main_object.multivariate_mass_radius = True
        elif mc.common_models[main_object.stellar_ref].compute_density:
            main_object.list_pams_common.update(['mass'])
            main_object.list_pams_common.update(['radius'])
            main_object.multivariate_mass_radius = True
        else:
            """ this is the density of the star (in solar units) """
            main_object.list_pams_common.update(['density'])
            main_object.multivariate_mass_radius = False


def _prepare_planet_semimajor_axis(main_object, mc, **kwargs):

    """ Check if the mass of the star is computed using a multivariate prior.
    This is a star property and not a planet one, but we need it defined here and now
    because we need to know how to compute the stellar density if
    it replaces the semimajor axis as a free parameter in the transit fit
    """

    """ Default parametrization uses the stellar density and the impact
        parameter, it is possible to switch back to scaled semi-major axis and
        inclination respectively by activating the proper flag """

    if mc.common_models[main_object.planet_ref].use_scaled_semimajor_axis:
        """ a is the semi-major axis (in units of stellar radii) """
        main_object.list_pams_common.update(['a_Rs'])
        main_object.compute_scaled_semimajor_axis = False

    if mc.common_models[main_object.planet_ref].use_semimajor_axis:
        """ a is the semi-major axis (in units of stellar radii) """
        main_object.list_pams_common.update(['a_AU'])
        main_object.compute_semimajor_axis = False


def _prepare_planet_mass(main_object, mc, **kwargs):

    if not mc.common_models[main_object.planet_ref].use_mass \
        and not mc.common_models[main_object.planet_ref].use_scaled_mass \
        and not mc.common_models[main_object.planet_ref].use_stellar_scaled_mass:

        print("UNRECOVERABLE ERROR model {0:s} :".format(main_object.model_name))
        print('    Dynamical modelling requires the mass or the scaled mass of the planet as free parameters')
        print('    for efficient exploration of parameter space')
        quit()

    if mc.common_models[main_object.planet_ref].use_mass:
        main_object.list_pams_common.update(['M_Me'])

    if mc.common_models[main_object.planet_ref].use_scaled_mass:
        main_object.list_pams_common.update(['Me_Ms'])

    if mc.common_models[main_object.planet_ref].use_stellar_scaled_mass:
        main_object.list_pams_common.update(['M_Ms'])


def _prepare_stellar_mass(main_object, mc, **kwargs):

        try:
            multivariate_pams = mc.common_models[main_object.stellar_ref].multivariate_pams
        except AttributeError:
            multivariate_pams = []

        if 'mass' in multivariate_pams and 'radius' in multivariate_pams:
            main_object.list_pams_common.update(['mass'])
            main_object.list_pams_common.update(['radius'])
        elif mc.common_models[main_object.stellar_ref].compute_density:
            main_object.list_pams_common.update(['mass'])
            main_object.list_pams_common.update(['radius'])
        elif mc.common_models[main_object.stellar_ref].compute_density:
            main_object.list_pams_common.update(['mass'])
            main_object.list_pams_common.update(['radius'])
        elif mc.common_models[main_object.stellar_ref].compute_mass:
            main_object.list_pams_common.update(['density'])
            main_object.list_pams_common.update(['radius'])
        elif mc.common_models[main_object.stellar_ref].compute_radius:
            main_object.list_pams_common.update(['density'])
            main_object.list_pams_common.update(['mass'])


def _prepare_planet_inclination(main_object, mc, **kwargs):

    if mc.common_models[main_object.planet_ref].use_inclination:
        """ i is the orbital inclination (in degrees) """
        main_object.list_pams_common.update(['i'])
        main_object.compute_inclination = False
    else:
        """ b is the impact parameter """
        main_object.list_pams_common.update(['b'])
        main_object.compute_inclination = True

def _prepare_planet_time_inferior_conjunction(main_object, mc, **kwargs):

        if mc.common_models[main_object.planet_ref].use_time_inferior_conjunction:
            main_object.list_pams_common.update(['Tc'])
        else:
            main_object.list_pams_common.update(['mean_long'])
            main_object.compute_time_inferior_conjunction = True
            # mean longitude = argument of pericenter + mean anomaly at Tref


