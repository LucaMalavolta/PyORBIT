from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *

from pyorbit.subroutines.common import np, convert_rho_to_ars, convert_b_to_i
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.keywords_definitions import *


class CommonPlanets(AbstractCommon):
    """
    Inherited class from AbstractCommon

    For computational reason it is better to fit :math:`\sqrt{e}\sin{\omega}` and :math:`\sqrt{e}\cos{\omega}`.

    Attributes:
        :model_class (string): identify the kind of class
        :list_pams: all the possible parameters that can be assigned to a planet are listed here
        :default_bounds: these default boundaries are used when the user does not define them in the yaml file
        :recenter_pams: circular parameters that may need a recentering around the most likely value after the global
            optimization run
        :period_average: variable used only by TRADES
    """

    model_class = 'planet'

    """ choice to parametrize the eccentricity and argument of pericenter:
        Standard: $e$ and $\omega$
        Ford2006: $e \cos{\omega }$ and $e \sin{\omega}$
        Eastman2013: $\sqrt{e} \cos{\omega }$ and $\sqrt{e} \sin{\omega}$
    """
    parametrization_list = ['Ford2006', 'Eastman2013', 'Standard',
                            'Ford2006_Tcent', 'Eastman2013_Tcent', 'Standard_Tcent',
                            'Ford2006_Tc', 'Eastman2013_Tc', 'Standard_Tc']
    orbit_list = ['circular', 'keplerian', 'dynamical']


    parameters_dictionary = {
        'P': # Orbital period of the planet
            {
                'bounds': [0.4, 100000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
                'unit': 'days',
            },
        'K': # RV semi-amplitude, in m/s
            {
                'bounds': [0.001, 2000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
            },
        'Tc':
            {
                'bounds': [0.0, 1000.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
            },
        'mean_long':
            {
                'bounds': [0.0, 360.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0000,
            },
        'e_coso':
            {
                'bounds': [-1.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0000,
            },
        'e_sino':
            {
                'bounds': [-1.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0000,
            },
        'sre_coso':
            {
                'bounds': [-1.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0000,
            },
        'sre_sino':
            {
                'bounds': [-1.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0000,
            },
        'e':
            {
                'bounds': [0.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0000,
            },
        'omega':
            {
                'bounds': [0.0, 360.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 90.,
            },
        'M_Me':
            {
                'bounds': [0.05, 1000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
            },
        'i':
            {
                'bounds':  [0.0, 180],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 90.,
            },
        'Omega':
            {
                'bounds':  [0.0, 180],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 90.,
            },
        'R_Rs':
            {
                'bounds': [0.00001, 0.5],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.05,
            },
        'a_Rs':
            {
                'bounds': [0.00001, 500.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 1.0,
            },
        'b':
            {
                'bounds': [0.0, 2.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0,
            },
        'lambda':
            {
                'bounds': [-180.0000, 180.0000],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0,
            },
        'phase_amp':
            {
                'bounds':  [0.00, 0.50],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0,
            },
        'delta_occ':
            {
                'bounds': [0.00, 0.50],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0,
            },
        'phase_off':
            {
                'bounds': [-180.0, 180.0],  #FIX
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0,
            },
        'albedo':
            {
                'bounds': [0., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0,
            },
        'redist':
            {
                'bounds': [0., 1.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0,
            },
        'insol':
            {
                'bounds': [0, 1000000000],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : 1.00000,
            },
    }

    recenter_pams = {'mean_long', 'omega', 'Omega'}

    def __init__(self, *args, **kwargs):
        super(CommonPlanets, self).__init__(*args, **kwargs)

        self.orbit = 'keplerian'
        self.parametrization = 'Eastman2013'

        self.use_inclination = False
        self.use_semimajor_axis = False
        self.use_time_inferior_conjunction = False
        self.use_mass = False
        self.use_shared_ttvs = False

        self.omega_star = True

        self.period_average = None
        # Variable used only by TRADES

        self.compute_inclination = True
        self.compute_time_inferior_conjunction = True
        self.compute_mass = False
        self.compute_mean_longitude = False
        self.compute_semimajor_axis = True

    def initialize_model(self, mc, **kwargs):

        self.orbit = kwargs.get('orbit', self.orbit)
        if self.orbit in self.orbit_list:
            print('Using orbital model: ', self.orbit)

            if self.orbit == 'circular':
                self.fix_list['e'] = np.asarray([0.000, 0.0000], dtype=np.double)
                self.fix_list['omega'] = np.asarray([90.0, 0.0000], dtype=np.double)

        else:
            print('ERROR in configuration file - orbital model: not supported')
            quit()


        self.parametrization = kwargs.get('parametrization', self.parametrization)
        if self.parametrization in self.parametrization_list:
            print('Using orbital parametrization: ', self.parametrization)

            if self.parametrization[-5:] == 'Tcent' or self.parametrization[-5:] == 'Tc':
                self.use_time_inferior_conjunction = True
        else:
            print('ERROR in configuration file - orbital model: not supported')
            quit()

        self.use_semimajor_axis = kwargs.get('use_semimajor_axis', self.use_semimajor_axis)
        if self.use_semimajor_axis:
            self.compute_semimajor_axis = False
            print('Semi-major axis replacing stellar density as a free parameter: ', True)

        self.use_inclination = kwargs.get('use_inclination', self.use_inclination)
        if self.use_inclination:
            self.compute_inclination = False
            print('Inclination replacing impact parameter as a free parameter: ', True)

        self.use_mass = kwargs.get('use_mass', self.use_mass)
        if self.use_mass :
            print('Planetary mass replacing RV semi-amplitude as a free parameter: ', True)

        self.use_time_inferior_conjunction = kwargs.get('use_time_inferior_conjunction', self.use_time_inferior_conjunction)
        if self.use_time_inferior_conjunction:
            self.compute_time_inferior_conjunction = False
            self.compute_mean_longitude = True
            print('Time of inferior conjunction replacing mean longitude as a free parameter: ', True)
        else:
            self.compute_time_inferior_conjunction = True
            self.compute_mean_longitude = False
        

        for use_shared_ttvs in keywords_shared_ttv:
            self.use_shared_ttvs = kwargs.get(use_shared_ttvs, self.use_shared_ttvs)
            if self.use_shared_ttvs:
                break

        for tc_list in keywords_tc_list:
            self.tc_list = kwargs.get(tc_list, None)
            if self.tc_list:
                print('List of times of inferior conjuctions: ', self.tc_list)
                break

        for tc_flag in keywords_tc_flag:
            self.tc_flag = kwargs.get(tc_flag, None)
            if self.tc_flag:
                print('Dataset flag of times of inferior conjuctions: ', self.tc_flag)
                break

        print()


    def define_derived_parameters(self):

        derived_list = []

        if 'e_coso' in self.sampler_parameters and  \
            'e_sino' in self.sampler_parameters:

            pam00_index = self.sampler_parameters['e_coso']
            pam01_index = self.sampler_parameters['e_sino']

            try:
                del self.parameter_index['e_coso']
                del self.parameter_index['e_sino']
            except:
                pass

            if 'e' not in self.parameter_index:
                self.transformation['e'] = get_2var_e
                self.parameter_index['e'] = [pam00_index, pam01_index]
                derived_list.append('e')

            if 'omega' not in self.parameter_index:
                self.transformation['omega'] = get_2var_o
                self.parameter_index['omega'] = [pam00_index, pam01_index]
                derived_list.append('omega')

        if 'sre_coso' in self.sampler_parameters and  \
            'sre_sino' in self.sampler_parameters:

            pam00_index = self.sampler_parameters['sre_coso']
            pam01_index = self.sampler_parameters['sre_sino']

            try:
                del self.parameter_index['sre_coso']
                del self.parameter_index['sre_sino']
            except:
                pass

            if 'e' not in self.parameter_index:
                self.transformation['e'] = get_2var_sre
                self.parameter_index['e'] = [pam00_index, pam01_index]
                derived_list.append('e')

            if 'omega' not in self.parameter_index:
                self.transformation['omega'] = get_2var_o
                self.parameter_index['omega'] = [pam00_index, pam01_index]
                derived_list.append('omega')


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


    def define_starting_point_from_derived(self, starting_point, var_sampler):
        """
        Eccentricity and argument of pericenter require a special treatment

        since they can be provided as fixed individual values or may need to be combined
        in :math:`\sqrt{e}\sin{\omega}` and :math:`\sqrt{e}\cos{\omega}` if are both free variables

        Args:
            :starting_point:
            :var_sampler:
        Returns:
            :bool:

        """
        if var_sampler == 'sre_coso' or var_sampler=='sre_sino':

            if 'e' in self.starts and 'omega' in self.starts:

                starting_point[self.sampler_parameters['sre_coso']] = \
                    np.sqrt(self.starts['e']) * np.cos(self.starts['omega'])
                starting_point[self.sampler_parameters['sre_sino']] = \
                    np.sqrt(self.starts['e']) * np.sin(self.starts['omega'])

            elif 'sre_coso' in self.starts and 'sre_sino' in self.starts:
                starting_point[self.sampler_parameters['sre_coso']] = self.starts['sre_coso']
                starting_point[self.sampler_parameters['sre_coso']] = self.starts['sre_sino']

            return True

        if var_sampler == 'e_coso' or var_sampler=='e_sino':

            if 'e' in self.starts and 'omega' in self.starts:
                starting_point[self.sampler_parameters['e_coso']] = \
                    self.starts['e'] * np.cos(self.starts['omega'])
                starting_point[self.sampler_parameters['e_sino']] = \
                    self.starts['e'] * np.sin(self.starts['omega'])

            elif 'e_coso' in self.starts and 'e_sino' in self.starts:
                starting_point[self.sampler_parameters['e_coso']] = self.starts['e_coso']
                starting_point[self.sampler_parameters['e_coso']] = self.starts['e_sino']

            return True

        return False


    def update_parameter_values_for_dynamical(self, parameter_values, Tref, prepend=''):

        if self.compute_inclination:
            if self.compute_semimajor_axis:
                parameter_values[prepend+'a_Rs'] = convert_rho_to_ars(parameter_values[prepend+'P'], parameter_values['density'])
            parameter_values[prepend+'i'] = convert_b_to_i(
                parameter_values[prepend+'b'], parameter_values[prepend+'e'], parameter_values[prepend+'omega'], parameter_values[prepend+'a_Rs'])

        if self.compute_time_inferior_conjunction:
            parameter_values[prepend+'Tc']= kepler_exo.kepler_phase2Tc_Tref(
                parameter_values[prepend+'P'],
                parameter_values[prepend+'mean_long'],
                parameter_values[prepend+'e'],
                parameter_values[prepend+'omega']) + Tref

        if self.compute_mean_longitude:
            parameter_values[prepend+'mean_long'] = kepler_exo.kepler_Tc2phase_Tref(
                parameter_values[prepend+'P'],
                parameter_values[prepend+'Tc'] - Tref,
                parameter_values[prepend+'e'],
                parameter_values[prepend+'omega'])