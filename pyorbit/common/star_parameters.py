from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *


class CommonStarParameters(AbstractCommon):
    ''' This class must be used by each planet in the system
            model_name is the way the planet is identified
    '''
    model_class = 'star_parameters'

    parameters_dictionary = {
        'radius': # Radius of the star, in Solar radii
            {
                'bounds': [0., 2.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 1.0,
                'unit': 'solar_unit',
            },
        'mass': # Mass of the star, in Solar masses
            {
                'bounds': [0., 2.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 1.0,
                'unit': 'solar_unit',
            },
        'density': # Density of the star, in Solar density units
            {
                'bounds': [0., 5.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 1.0,
                'unit': 'solar_unit',
            },
        'i_star': # Inclination of the star
            {
                'bounds': [0., 180.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 90,
                'unit': 'degree',
            },
        'v_sini': # Projected rotational velocity of the star
            {
                'bounds': [0., 200.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 1.6,
                'unit': 'km/s',
            },
        'rotation_period': # Rotation period of the star
            {
                'bounds': [1., 1000.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 27,
                'unit': 'days',
            },
        'temperature': # Effective temperature of the photosphere
            {
                'bounds': [2000., 11000.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 5777,
                'unit': 'kelvin',
            },
        'line_contrast':
            {
                'bounds': [0., 100.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 50,
                'unit': 'percentual',
            },
        'line_fwhm':
            {
                'bounds': [0., 12.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 6,
                'unit': 'km/s',
            },
        'rv_center':
            {
                'bounds': [-3e2, 3e2],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0,
                'unit': 'km/s',
            },
    }

    recenter_pams = set()

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
