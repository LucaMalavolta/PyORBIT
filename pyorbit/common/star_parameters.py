from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *
from pyorbit.keywords_definitions import *

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
        'cosi_star': # Inclination of the star
            {
                'bounds': [0., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 1.,
                'unit': 'unitary',
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
        'activity_decay': # Rotation period of the star
            {
                'bounds': [10., 10000.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 1000,
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
    recenter_pams = OrderedSet()


    def __init__(self, *args, **kwargs):
        super(CommonStarParameters, self).__init__(*args, **kwargs)

        self.use_equatorial_velocity = False
        self.use_stellar_rotation_period = False
        self.use_stellar_inclination = False
        self.use_cosine_stellar_inclination = False
        self.use_differential_rotation = False
        self.use_stellar_radius = False
        self.use_projected_velocity = True
        self.convective_order = 0
        self.compute_mass = True
        self.compute_radius = False
        self.compute_density = False

    def initialize_model(self, mc, **kwargs):

        """ check if the stellar_rotation_period has to be used as parameter
            when the stellar rotation period is given as a prior, then the equatorial velcoity
            is computed using the rotational period and the radius of the star.
            The stellar inclination must be used as a free parameter
            and the veq_sini as a prior to be checked a posteriori, as the determination of the
            stellar inclination from the veq_sini could bias its distribution
            From several experiments, I determined that PyDE/MCMC convergence is
            much more difficult if v_eq and i_star are left as free parameters and
            the output is compared with too many priors
        """

        for keyword in keywords_stellar_rotation:
            self.use_stellar_rotation_period = kwargs.get(keyword, self.use_stellar_rotation_period)
        if self.use_stellar_rotation_period:
            self.use_equatorial_velocity = False
            self.use_stellar_inclination = True
            self.use_stellar_radius = True
            self.use_projected_velocity = False

        """ check if the differemntial rotation should be included in the model"""
        for keyword in keywords_differential_rotation:
            self.use_differential_rotation = kwargs.get(keyword, self.use_differential_rotation)
        if self.use_differential_rotation and not self.use_stellar_rotation_period:
            self.use_equatorial_velocity = True
            self.use_stellar_inclination = True
            self.use_projected_velocity = False


        """ The user can force the use of the equatorial velocity, the stellar inclination,
            and the stellar radius, by activating the corresponding flags in the model"""
        self.use_equatorial_velocity = kwargs.get('use_equatorial_velocity', self.use_equatorial_velocity)
        self.use_stellar_inclination = kwargs.get('use_stellar_inclination', self.use_stellar_inclination)
        self.use_cosine_stellar_inclination = kwargs.get('use_cosine_stellar_inclination', False) # Directly addressed here for back-compatibility
        if self.use_cosine_stellar_inclination:
            self.use_stellar_inclination = True
        self.use_stellar_radius = kwargs.get('use_stellar_radius', self.use_stellar_radius)
        self.use_projected_velocity = kwargs.get('use_projected_velocity', self.use_projected_velocity)

        """ the user can decide how to deal with the mass-radius-density correlation
            Density and (sometimes) radius can be involved in transit fit, while there is no way to measure the
            mass from the star from radial velocities or photometry, so the mass should have lower priority
            as a free parameter.
        """
        self.compute_mass = kwargs.get('compute_mass', self.compute_mass)
        self.compute_radius = kwargs.get('compute_radius', self.compute_radius)
        self.compute_density = kwargs.get('compute_density', self.compute_density)

        try:
            multivariate_pams = self.multivariate_pams
            if len(multivariate_pams) > 0:
                self.compute_density = True
        except AttributeError:
            pass

        if self.compute_density:
            self.compute_radius = False
            self.compute_mass = False
        if self.compute_radius:
            self.compute_density = False
            self.compute_mass = False
        if self.compute_mass:
            self.compute_density = False
            self.compute_radius = False

        if not (self.compute_mass or self.compute_radius or self.compute_density):
            self.compute_mass = True

        self.convective_order = kwargs.get('convective_order', self.convective_order)

        if self.use_equatorial_velocity and self.use_stellar_inclination and self.use_stellar_rotation_period and self.use_stellar_radius:
            print('Possible source of unexpected behaviour, I will quit')
            print('These parameters are correlated and should not be all free simultaneously:')
            print('- stellar rotation period ')
            print('- stellar radius ')
            print('- stellar inclination ')
            print('- equatorial velocity')
            print()
            quit()


    def define_derived_parameters(self):

        derived_list = []

        if self.use_cosine_stellar_inclination and \
            'cosi_star' in self.sampler_parameters and \
            'i_star' not in self.parameter_index:

            pam00_index = self.sampler_parameters['cosi_star']

            self.transformation['i_star'] = get_var_arccosine
            self.parameter_index['i_star'] = pam00_index

            derived_list.append('i_star')

        if not self.use_equatorial_velocity and \
            'rotation_period' in self.sampler_parameters and \
            'radius' in self.sampler_parameters and \
            'veq_star' not in self.parameter_index:

            pam00_index = self.sampler_parameters['rotation_period']
            pam01_index = self.sampler_parameters['radius']

            self.transformation['veq_star'] = get_2var_prot_rstar_veq
            self.parameter_index['veq_star'] = [pam00_index, pam01_index]

            derived_list.append('veq_star')

        if not self.use_equatorial_velocity and \
            'rotation_period' in self.sampler_parameters and \
            'radius' in self.sampler_parameters and \
            'i_star' in self.sampler_parameters and  \
            'v_sini' not in self.sampler_parameters:

            pam00_index = self.sampler_parameters['rotation_period']
            pam01_index = self.sampler_parameters['radius']
            pam02_index = self.sampler_parameters['i_star']

            self.transformation['v_sini'] = get_3var_prot_rstar_istar_veq
            self.parameter_index['v_sini'] = [pam00_index, pam01_index, pam02_index]

            derived_list.append('v_sini')

        if not self.use_equatorial_velocity and \
            'rotation_period' in self.sampler_parameters and \
            'radius' in self.sampler_parameters and \
            'cosi_star' in self.sampler_parameters and  \
            'v_sini' not in self.sampler_parameters:

            pam00_index = self.sampler_parameters['rotation_period']
            pam01_index = self.sampler_parameters['radius']
            pam02_index = self.sampler_parameters['cosi_star']

            self.transformation['v_sini'] = get_3var_prot_rstar_cosistar_veq
            self.parameter_index['v_sini'] = [pam00_index, pam01_index, pam02_index]

            derived_list.append('v_sini')

        if 'veq_star' in self.sampler_parameters and \
            'i_star' in self.sampler_parameters and  \
            'v_sini' not in self.parameter_index:

            pam00_index = self.sampler_parameters['veq_star']
            pam01_index = self.sampler_parameters['i_star']

            self.transformation['v_sini'] = get_2var_veq_istar_vsini
            self.parameter_index['v_sini'] = [pam00_index, pam01_index]

            derived_list.append('v_sini')

        if 'veq_star' in self.sampler_parameters and \
            'cosi_star' in self.sampler_parameters and  \
            'v_sini' not in self.parameter_index:

            pam00_index = self.sampler_parameters['veq_star']
            pam01_index = self.sampler_parameters['cosi_star']

            self.transformation['v_sini'] = get_2var_veq_cosi_vsini
            self.parameter_index['v_sini'] = [pam00_index, pam01_index]

            derived_list.append('v_sini')

        if 'veq_star' in self.sampler_parameters and \
            'radius' in self.sampler_parameters and  \
            'rotation_period' not in self.parameter_index:

            pam00_index = self.sampler_parameters['veq_star']
            pam01_index = self.sampler_parameters['radius']

            self.transformation['rotation_period'] = get_2var_veq_radius_rot
            self.parameter_index['rotation_period'] = [pam00_index, pam01_index]

            derived_list.append('rotation_period')

        if 'density' in self.sampler_parameters and  \
            'radius' in self.sampler_parameters and \
            'mass' not in self.parameter_index:

            pam00_index = self.sampler_parameters['density']
            pam01_index = self.sampler_parameters['radius']

            self.transformation['mass'] = get_2var_mass
            self.parameter_index['mass'] = [pam00_index, pam01_index]

            derived_list.append('mass')

        if 'mass' in self.sampler_parameters and  \
            'radius' in self.sampler_parameters and \
            'density' not in self.parameter_index:

            pam00_index = self.sampler_parameters['mass']
            pam01_index = self.sampler_parameters['radius']

            self.transformation['density'] = get_2var_rho
            self.parameter_index['density'] = [pam00_index, pam01_index]

            derived_list.append('density')

        if 'mass' in self.sampler_parameters and  \
            'density' in self.sampler_parameters and \
            'radius' not in self.parameter_index:

            pam00_index = self.sampler_parameters['mass']
            pam01_index = self.sampler_parameters['density']

            self.transformation['radius'] = get_2var_radius
            self.parameter_index['radius'] = [pam00_index, pam01_index]

            derived_list.append('radius')


        for pam in derived_list:
            if pam not in self.bounds:
                self.bounds[pam] = self.default_bounds[pam]

            if pam not in self.prior_pams:

                if pam in self.bounds:
                    self.prior_pams[pam] = self.bounds[pam]
                else:
                    self.prior_pams[pam] = self.default_bounds[pam]

                self.prior_kind[pam] = 'Uniform'

        return

