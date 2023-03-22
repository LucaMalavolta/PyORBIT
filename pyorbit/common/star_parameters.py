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
        'veq_star': #equatorial velocity of the star
            {
                'bounds': [0.00, 70.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 1.6,
                'unit': 'km/s',
            },
        'alpha_rotation':
            {
                'bounds': [0.00, 1.00],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 1.6,
                'unit': 'unit',
            },

    }
    recenter_pams = set()


    def __init__(self, *args, **kwargs):
        super(CommonStarParameters, self).__init__(*args, **kwargs)

        self.use_equatorial_velocity = False
        self.use_stellar_rotation = False
        self.use_stellar_inclination = False
        self.use_differential_rotation = False

    def initialize_model(self, mc, **kwargs):

        self.use_equatorial_velocity = kwargs.get('use_equatorial_velocity', self.use_equatorial_velocity)

        """ check if the stellar_inclination has to be used as parameter """
        self.use_stellar_inclination = kwargs.get('use_stellar_inclination', self.use_stellar_inclination)

        """ check if the stellar_rotation has to be used as parameter """
        self.use_stellar_rotation = kwargs.get('use_stellar_rotation', self.use_stellar_rotation)
        self.use_differential_rotation = kwargs.get('use_differential_rotation', self.use_differential_rotation)

        self.use_differential_rotation = kwargs.get('use_differential_rotation', self.use_differential_rotation)
        if self.use_differential_rotation:
            print('Differential rotation in model switched on: ')
            print(' - Equatorial velocity and stellar inclination as independent parameters ')
            print(' - Check on rotational period, unless otherwise specified ')

            self.use_equatorial_velocity = True
            self.use_stellar_inclination = True

            self.use_stellar_rotation = kwargs.get('use_stellar_rotation', True)

    def define_special_parameter_properties(self, ndim, output_lists, pam):

        skip_first_parametrization = False
        skip_second_parametrization = False


        if self.use_equatorial_velocity and self.use_stellar_rotation:
            if not(pam == "veq_star" or pam == 'i_star' or pam == 'rotation_period'):
                skip_first_parametrization = True

            for var_check in ['v_sini', 'veq_star', 'i_star', 'radius', 'rotation_period']:
                if var_check in self.sampler_parameters:
                    skip_first_parametrization = True

            if 'v_sini' in self.fix_list or 'radius' in self.fix_list:
                skip_first_parametrization = True

        elif self.use_equatorial_velocity:
            if not(pam == "veq_star" or pam == 'i_star'):
                skip_first_parametrization = True

            for var_check in ['v_sini', 'veq_star', 'i_star']:
                if var_check in self.sampler_parameters:
                    skip_first_parametrization = True

            if 'v_sini' in self.fix_list:
                skip_first_parametrization = True


        if not(pam == "mass" or pam == "radius") or \
            not ('mass' in self.multivariate_pams
                and 'radius' in self.multivariate_pams
                and self.multivariate_priors) \
            or 'mass' in self.fix_list \
            or 'radius' in self.fix_list \
            or 'mass' in self.sampler_parameters \
            or 'radius' in self.sampler_parameters:
            skip_second_parametrization = True


        if skip_first_parametrization and skip_second_parametrization:
            return ndim, output_lists, False

        if not skip_first_parametrization:
            self.transformation['veq_star'] = get_var_val
            self.parameter_index['veq_star'] = ndim
            self.transformation['i_star'] = get_var_val
            self.parameter_index['i_star'] = ndim + 1

            self.transformation['v_sini'] = get_2var_vsini
            self.parameter_index['v_sini'] = [ndim, ndim + 1]

            if self.use_stellar_rotation:
                self.transformation['rotation_period'] = get_var_val
                self.parameter_index['rotation_period'] = ndim + 2

                self.transformation['radius'] = get_2var_veq_rot_radius
                self.parameter_index['radius'] = [ndim, ndim + 2]

                parameter_list = ['veq_star', 'i_star', 'rotation_period']
                derived_list = ['v_sini', 'radius']
            else:
                parameter_list = ['veq_star', 'i_star']
                derived_list = ['v_sini']

        if not skip_second_parametrization:
            self.transformation['mass'] = get_var_val
            self.parameter_index['mass'] = ndim
            self.transformation['radius'] = get_var_val
            self.parameter_index['radius'] = ndim + 1

            self.transformation['density'] = get_2var_rho
            self.parameter_index['density'] = [ndim, ndim + 1]
            parameter_list = ['mass', 'radius']
            derived_list = ['density']

        for pam in parameter_list:

            if pam not in self.bounds:
                self.bounds[pam] = self.default_bounds[pam]

            if pam not in self.spaces:
                self.spaces[pam] = self.default_spaces[pam]

            output_lists['bounds'].append(self.bounds[pam])

            if pam not in self.prior_pams:
                self.prior_kind[pam] = self.default_priors[pam][0]
                self.prior_pams[pam] = self.default_priors[pam][1]

            nested_coeff = nested_sampling_prior_prepare(self.prior_kind[pam],
                                                          output_lists['bounds'][-1],
                                                          self.prior_pams[pam],
                                                          self.spaces[pam])

            output_lists['spaces'].append(self.spaces[pam])
            output_lists['priors'].append([self.prior_kind[pam], self.prior_pams[pam], nested_coeff])

            self.sampler_parameters[pam] = ndim
            ndim += 1

        for pam in derived_list:
            if pam not in self.bounds:
                self.bounds[pam] = self.default_bounds[pam]

            if pam not in self.prior_pams:

                if pam in self.bounds:
                    self.prior_pams[pam] = self.bounds[pam]
                else:
                    self.prior_pams[pam] = self.default_bounds[pam]

                self.prior_kind[pam] = 'Uniform'

        return ndim, output_lists, True
