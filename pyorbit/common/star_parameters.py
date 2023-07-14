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
                'fixed' : 0.6,
                'unit': 'unit',
            },
        'convective_c1':
            {
                'bounds': [0.00, 5.00],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0,
                'unit': 'unit',
            },
        'convective_c2':
            {
                'bounds': [-5.00, 0.00],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0,
                'unit': 'unit',
            },
        'convective_c3':
            {
                'bounds': [-5.00, 5.00],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0,
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
        self.use_stellar_radius = False
        self.convective_order = 0
        self.preserve_density = True

    def initialize_model(self, mc, **kwargs):

        """ check if the stellar_rotation has to be used as parameter
            when the stellar rotation period is given as a prior, then the equatorial velcoity
            is computed using the rotational period and the radius of the star.
            The stellar inclination must be used as a free parameter
            and the veq_sini as a prior to be checked a posteriori, as the determination of the
            stellar inclination from the veq_sini could bias its distribution
        """
        self.use_stellar_rotation = kwargs.get('use_stellar_rotation', self.use_stellar_rotation)
        if self.use_stellar_rotation:
            self.use_equatorial_velocity = False
            self.use_stellar_inclination = True
            self.use_stellar_radius = True


        """ The user can force the use of the equatorial velocity, the stellar inclination, 
            and the stellar radius, by activating the corresponding flags in the model"""
        self.use_equatorial_velocity = kwargs.get('use_equatorial_velocity', self.use_equatorial_velocity)
        self.use_stellar_inclination = kwargs.get('use_stellar_inclination', self.use_stellar_inclination)
        self.use_stellar_radius = kwargs.get('use_stellar_radius', self.use_stellar_radius)

        """ check if the differential rotation should be included in the model"""
        self.use_differential_rotation = kwargs.get('use_differential_rotation', self.use_differential_rotation)

        """ the user can decide how to deal with the mass-radius-density correlation
            density and (sometimes) radius can be involved in transit fit, while there is no way to measure the 
            mass from the star from radial velocities or photometry, so the mass should have lower priority 
            as a free parameter.
        """
        self.preserve_density = kwargs.get('preserve_density', self.preserve_density)
        self.preserve_density = ~ kwargs.get('use_mass_radius', ~ self.preserve_density)

        self.convective_order = kwargs.get('convective_order', self.convective_order)

        if self.use_equatorial_velocity and self.use_stellar_inclination and self.use_stellar_rotation and self.use_stellar_radius:
            print('Possible source of unexpected behaviour, I will quit')
            print('These parameters are correlated and should not be all free simultaneously:')
            print('- stellar rotation period ')
            print('- stellar radius ')
            print('- stellar inclination ')
            print('- equatorial velocity')
            print()
            quit()


    def define_derived_parameters(self):

        # prima controlla se la variabile è una variabile del sample o derivata 
        # Prima controlla le variabili derivate 
        # se tutti i parametri necessari sono stati inseriti allora procedi con le variabili derivate 
        # il controllo va fatto ogni volt aperchè non sappiamo quali sono le variabili da derivare 


        if not self.use_equatorial_velocity and \
            'rotation_period' in self.sampler_parameters and \
            'radius' in self.sampler_parameters and \
            'veq_star' not in self.parameter_index:

            pam00_index = self.parameter_index['rotation_period']
            pam01_index = self.parameter_index['radius']

            self.transformation['veq_star'] = get_2var_prot_rstar_veq
            self.parameter_index['veq_star'] = [pam00_index, pam01_index]

            parameter_list = ['v_sini', 'rotation_period', 'radius']
            derived_list.a = ['veq_star', 'i_star']



        if self.use_stellar_rotation and pam=='rotation_period':
            self.transformation['rotation_period'] = get_var_val
            self.parameter_index['rotation_period'] = ndim




        if self.use_stellar_inclination and pam=='i_star':
            self.transformation['i_star'] = get_var_val
            self.parameter_index['i_star'] = ndim



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



        if  'rotation_period' in self.sampler_parameters and \
            'radius' in self.sampler_parameters and \
            'veq_star' not in self.parameter_index:

            pam00_index = self.parameter_index['rotation_period']
            pam01_index = self.parameter_index['radius']

            self.transformation['veq_star'] = get_2var_prot_rstar_veq
            self.parameter_index['veq_star'] = [pam00_index, pam01_index]

            parameter_list = ['v_sini', 'rotation_period', 'radius']
            derived_list.a = ['veq_star', 'i_star']

            get_2var_prot_rstar_veq


        parameter_list = ['rotation_period']
        derived_list = []









    def define_special_parameter_properties(self, ndim, output_lists, pam):

        skip_first_parametrization = True
        skip_second_parametrization = True
        skip_third_parametrization = True
        skip_fourth_parametrization = True

        # prima controlla se la variabile è una variabile del sample o derivata 
        # Prima controlla le variabili derivate 
        # se tutti i parametri necessari sono stati inseriti allora procedi con le variabili derivate 
        # il controllo va fatto ogni volt aperchè non sappiamo quali sono le variabili da derivare 


        if not self.use_equatorial_velocity and \
            'rotation_period' in self.sampler_parameters and \
            'radius' in self.sampler_parameters and \
            'veq_star' not in self.parameter_index:

            pam00_index = self.parameter_index['rotation_period']
            pam01_index = self.parameter_index['radius']

            self.transformation['veq_star'] = get_2var_prot_rstar_veq
            self.parameter_index['veq_star'] = [pam00_index, pam01_index]

            parameter_list = ['v_sini', 'rotation_period', 'radius']
            derived_list.a = ['veq_star', 'i_star']



        if self.use_stellar_rotation and pam=='rotation_period':
            self.transformation['rotation_period'] = get_var_val
            self.parameter_index['rotation_period'] = ndim




        if self.use_stellar_inclination and pam=='i_star':
            self.transformation['i_star'] = get_var_val
            self.parameter_index['i_star'] = ndim


        if 

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

        else:
    


        if  'rotation_period' in self.sampler_parameters and \
            'radius' in self.sampler_parameters and \
            'veq_star' not in self.parameter_index:

            pam00_index = self.parameter_index['rotation_period']
            pam01_index = self.parameter_index['radius']

            self.transformation['veq_star'] = get_2var_prot_rstar_veq
            self.parameter_index['veq_star'] = [pam00_index, pam01_index]

            parameter_list = ['v_sini', 'rotation_period', 'radius']
            derived_list.a = ['veq_star', 'i_star']

            get_2var_prot_rstar_veq


        parameter_list = ['rotation_period']
            derived_list = []





        self.use_equatorial_velocity = False
            self.use_stellar_radius = True





        #print(self.parameter_index)

        if self.use_equatorial_velocity and self.use_stellar_inclination:
            if self.use_stellar_rotation:
                if pam == "veq_star" or pam == 'i_star' or pam == 'rotation_period':
                    skip_first_parametrization = False

                for var_check in ['v_sini', 'veq_star', 'i_star', 'radius', 'rotation_period']:
                    if var_check in self.sampler_parameters:
                        skip_first_parametrization = True

                if 'v_sini' in self.fix_list or 'radius' in self.fix_list:
                    skip_first_parametrization = True

            elif self.use_stellar_radius:
                if pam == "veq_star" or pam == 'i_star' or pam == 'radius':
                    skip_first_parametrization = False

                for var_check in ['v_sini', 'veq_star', 'i_star', 'radius', 'rotation_period']:
                    if var_check in self.sampler_parameters:
                        skip_first_parametrization = True

                if 'v_sini' in self.fix_list or 'rotation_period' in self.fix_list:
                    skip_first_parametrization = True

            else:
                if pam == "veq_star" or pam == 'i_star':
                    skip_first_parametrization = False

                for var_check in ['v_sini', 'veq_star', 'i_star']:
                    if var_check in self.sampler_parameters:
                        skip_first_parametrization = True

                if 'v_sini' in self.fix_list:
                    skip_first_parametrization = True

        elif self.use_stellar_rotation or self.use_stellar_radius:
            if pam == "v_sini" or pam == 'rotation_period' or pam=='radius':
                    skip_fourth_parametrization = False

            for var_check in ['v_sini', 'i_star', 'radius', 'veq_star','rotation_period']:
                    if var_check in self.sampler_parameters:
                        skip_fourth_parametrization = True

            if 'i_star' in self.fix_list or 'veq_star' in self.fix_list:
                skip_fourth_parametrization = True


        if pam == "mass" or pam == "radius":
            if self.preserve_density:
                skip_second_parametrization = True
                skip_third_parametrization = False
            else:
                skip_second_parametrization = False

        if ('mass' in self.multivariate_pams
                and 'radius' in self.multivariate_pams
                and self.multivariate_priors):
            skip_second_parametrization = False
            skip_third_parametrization = True

        if ('mass' in self.fix_list \
            or 'radius' in self.fix_list \
            or 'mass' in self.sampler_parameters \
            or 'radius' in self.sampler_parameters):
            skip_second_parametrization = True
            skip_third_parametrization = True

        if ('density' in self.fix_list
            or 'density' in self.sampler_parameters):
            skip_third_parametrization = True

        if skip_first_parametrization and skip_second_parametrization and skip_third_parametrization and skip_fourth_parametrization:
            return ndim, output_lists, False
        
        if (skip_first_parametrization +
              skip_second_parametrization +  skip_third_parametrization +
               skip_fourth_parametrization) == 2:
            print('You got a random error that happens only once every five times',)
            print('The code will quit, you will need to relaunching it again',)
            quit()

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

            elif self.use_stellar_radius:
                self.transformation['radius'] = get_var_val
                self.parameter_index['radius'] = ndim + 2

                self.transformation['rotation_period'] = get_2var_veq_radius_rot
                self.parameter_index['rotation_period'] = [ndim, ndim + 2]

                parameter_list = ['veq_star', 'i_star', 'radius']
                derived_list = ['v_sini', 'rotation_period']

            else:
                parameter_list = ['veq_star', 'i_star']
                derived_list = ['v_sini']

        if not skip_fourth_parametrization:
            self.transformation['v_sini'] = get_var_val
            self.parameter_index['v_sini'] = ndim
            self.transformation['rotation_period'] = get_var_val
            self.parameter_index['rotation_period'] = ndim + 1
            self.transformation['radius'] = get_var_val
            self.parameter_index['radius'] = ndim + 2

            self.transformation['i_star'] = get_3var_vsini_prot_rstar_istar
            self.parameter_index['i_star'] = [ndim, ndim+1, ndim+2]

            self.transformation['veq_star'] = get_2var_prot_rstar_veq
            self.parameter_index['veq_star'] = [ndim + 1, ndim+2]

            parameter_list = ['v_sini', 'rotation_period', 'radius']
            derived_list = ['veq_star', 'i_star']



        if not skip_second_parametrization:
            self.transformation['mass'] = get_var_val
            self.parameter_index['mass'] = ndim
            self.transformation['radius'] = get_var_val
            self.parameter_index['radius'] = ndim + 1

            self.transformation['density'] = get_2var_rho
            self.parameter_index['density'] = [ndim, ndim + 1]
            parameter_list = ['mass', 'radius']
            derived_list = ['density']

        if not skip_third_parametrization:
            self.transformation['density'] = get_var_val
            self.parameter_index['density'] = ndim
            self.transformation['radius'] = get_var_val
            self.parameter_index['radius'] = ndim + 1

            self.transformation['mass'] = get_2var_rho
            self.parameter_index['mass'] = [ndim, ndim + 1]
            parameter_list = ['density', 'radius']
            derived_list = ['mass']

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
