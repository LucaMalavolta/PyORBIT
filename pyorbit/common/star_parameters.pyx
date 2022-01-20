from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *


class CommonStarParameters(AbstractCommon):
    ''' This class must be used by each planet in the system
            model_name is the way the planet is identified
    '''
    model_class = 'star_parameters'

    list_pams = {
        'radius',  # radius of the star, in Solar radii
        'mass',  # mass of the star, in Solar masses
        'density',  # density of the star, in Solar density units
        'o_star', # Sky-projected angle between stellar rotation axis and normal of orbit plane [deg]
        'i_star', # Inclination of the star
        'v_sini' # projected rotational velocity of the star
        'temperature' #effective temperature of the star, in K
    }

    default_bounds = {
        'radius': [0.0000, 2.0000],
        'mass': [0.0000, 2.0000],
        'density': [0.0000, 5.0000],
        'o_star': [0.0000, 180.0000],
        'i_star': [0.0000, 180.0000],
        'v_sini': [0.0000, 200.0000],
        'temperature': [2000., 11000.]
    }

    """ Must be the same parameters as in list_pams, because priors are applied only to _physical_ parameters """
    default_priors = {
        'radius': ['Uniform', []],
        'mass': ['Uniform', []],
        'density': ['Uniform', []],
        'o_star': ['Uniform', []],
        'i_star': ['Uniform', []],
        'v_sini': ['Uniform', []],
        'temperature': ['Uniform', []]
    }

    default_spaces = {
        'radius': 'Linear',
        'mass': 'Linear',
        'density': 'Linear',
        'o_star': 'Linear',
        'i_star': 'Linear',
        'v_sini': 'Linear',
        'temperature': 'Linear'
    }

    default_fixed = {
        'radius': 1.0000,
        'mass': 1.0000,
        'density': 1.0000,
        'o_star': 0.0000,
        'i_star': 90.0000,
        'v_sini': 1.6000,
        'temperature': 5777
    }

    recenter_pams = {}

    #def __init__(self, *args, **kwargs):
    #    super(CommonStarParameters, self).__init__(*args, **kwargs)

    def define_special_variable_properties(self, ndim, output_lists, var):

        if not(var == "mass" or var == "radius") or \
            not ('mass' in self.multivariate_vars
                and 'radius' in self.multivariate_vars
                and self.multivariate_priors) \
            or 'mass' in self.fix_list \
            or 'radius' in self.fix_list \
            or 'mass' in self.variable_sampler \
            or 'radius' in self.variable_sampler:
            return ndim, output_lists, False

        self.transformation['mass'] = get_var_val
        self.variable_index['mass'] = ndim
        self.transformation['radius'] = get_var_val
        self.variable_index['radius'] = ndim + 1

        self.transformation['density'] = get_2var_rho
        self.variable_index['density'] = [ndim, ndim + 1]
        variable_list = ['mass', 'radius']

        for var in variable_list:

            if var not in self.bounds:
                self.bounds[var] = self.default_bounds[var]

            if var not in self.spaces:
                self.spaces[var] = self.default_spaces[var]

            output_lists['bounds'].append(self.bounds[var])

            if var not in self.prior_pams:
                self.prior_kind[var] = self.default_priors[var][0]
                self.prior_pams[var] = self.default_priors[var][1]

            nested_coeff = nested_sampling_prior_prepare(self.prior_kind[var],
                                                          output_lists['bounds'][-1],
                                                          self.prior_pams[var],
                                                          self.spaces[var])

            output_lists['spaces'].append(self.spaces[var])
            output_lists['priors'].append([self.prior_kind[var], self.prior_pams[var], nested_coeff])

            self.variable_sampler[var] = ndim
            ndim += 1

        for var in ['density']:
            if var not in self.bounds:
                self.bounds[var] = self.default_bounds[var]

            if var not in self.prior_pams:

                if var in self.bounds:
                    self.prior_pams[var] = self.bounds[var]
                else:
                    self.prior_pams[var] = self.default_bounds[var]

                self.prior_kind[var] = 'Uniform'

        return ndim, output_lists, True
