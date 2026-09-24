from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *

from pyorbit.subroutines.common import np, convert_rho_to_ars, convert_b_to_i
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.keywords_definitions import *
from packaging.version import Version


class CommonExomoons(AbstractCommon):

    """
    Inherited class from AbstractCommon
    """
    model_class = 'exomoon'

    """ choice to parametrize the eccentricity and argument of pericenter:
        Standard: $e$ and $\omega$
        Ford2006: $e \cos{\omega }$ and $e \sin{\omega}$
        Eastman2013: $\sqrt{e} \cos{\omega }$ and $\sqrt{e} \sin{\omega}$
    """
    parametrization_list = ['Ford2006', 'Eastman2013', 'Standard',
                            'Ford2006_Tcent', 'Eastman2013_Tcent', 'Standard_Tcent',
                            'Ford2006_Tc', 'Eastman2013_Tc', 'Standard_Tc']
    orbit_list = ['circular', 'keplerian']

    parameters_dictionary = {
        'em_P': # Orbital period of the planet
            {
                'bounds': [0.4, 100000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
                'unit': 'days',
            },
        'em_Tc':
            {
                'bounds': [0.0, 1000.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
            },
        'em_tau':
            {
                'bounds': [0.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
            },
        'em_mean_long':
            {
                'bounds': [0.0, 360.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0000,
            },
        'em_coso':
            {
                'bounds': [-1.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0000,
            },
        'em_sino':
            {
                'bounds': [-1.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0000,
            },
        'em_sre_coso':
            {
                'bounds': [-1.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0000,
            },
        'em_sre_sino':
            {
                'bounds': [-1.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0000,
            },
        'em_e':
            {
                'bounds': [0.0, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.0000,
            },
        'em_omega':
            {
                'bounds': [0.0, 360.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 90.,
            },
        'em_M_Me':
            {
                'bounds': [0.05, 1000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
            },
        'em_M_Ms':
            {
                'bounds': [0.000000001, 0.5],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
            },
        'em_Me_Ms':
            {
                'bounds': [0.1, 2000.],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
            },
        'em_i':
            {
                'bounds':  [0.0, 180],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 90.,
            },
        'em_Omega':
            {
                'bounds':  [0.0, 360.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 180.,
            },
        'em_R_Rs':
            {
                'bounds': [0.00001, 0.5],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.05,
            },
        'em_a_Rs':
            {
                'bounds': [0.00001, 500.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 1.0,
            },
    }

    recenter_pams = {'em_mean_long', 'em_omega', 'em_Omega'}

    def __init__(self, *args, **kwargs):
        super(CommonPlanets, self).__init__(*args, **kwargs)

        self.orbit = 'keplerian'
        self.parametrization = 'Standard'

        self.use_inclination = True
        self.use_time_inferior_conjunction = True
        self.use_mass = True
        self.use_longitude_of_nodes = False

        self.compute_inclination = False
        self.compute_time_inferior_conjunction = False
        self.compute_mass = False
        self.compute_mean_longitude = True
        self.compute_planet_mass = False

    def initialize_model(self, mc, **kwargs):

        self.Tref = mc.Tref

        self.use_circular_orbit = kwargs.get('use_circular_orbit', False)
        self.orbit = kwargs.get('orbit', self.orbit)

        if self.orbit in self.orbit_list:

            if self.orbit == 'circular' or self.use_circular_orbit:
                self.fix_list['em_e'] = np.asarray([0.000, 0.0000], dtype=np.double)
                self.fix_list['em_omega'] = np.asarray([90.0, 0.0000], dtype=np.double)

        else:
            print("UNRECOVERABLE ERROR model {0:s} :".format(self.common_ref))
            print('    {0:s} orbital model not supported, check configuration file'.format(self.orbit))
            quit()


        self.parametrization = kwargs.get('parametrization', self.parametrization)
        if self.parametrization in self.parametrization_list:
            if self.parametrization[-5:] == 'Tcent' or self.parametrization[-5:] == 'Tc':
                self.use_time_inferior_conjunction = True
        else:
            print("UNRECOVERABLE ERROR model {0:s} :".format(self.common_ref))
            print('    {0:s} parametrization not supported, check configuration file'.format(self.parametrization))
            quit()


        self.use_inclination = kwargs.get('use_inclination', self.use_inclination)
        if self.use_inclination:
            self.compute_inclination = False

        self.use_mass = kwargs.get('use_mass', self.use_mass)


        self.use_time_inferior_conjunction = kwargs.get('use_time_inferior_conjunction', self.use_time_inferior_conjunction)
        if self.use_time_inferior_conjunction:
            self.compute_time_inferior_conjunction = False
            self.compute_mean_longitude = True
        else:
            self.compute_time_inferior_conjunction = True
            self.compute_mean_longitude = False

        self.use_longitude_of_nodes = kwargs.get('use_longitude_of_nodes', self.use_longitude_of_nodes)



        self.default_Omega = 180.00


    def print_info(self):
        print("*** exomoon {0:s} global parameters:".format(self.common_ref))

        if self.orbit == 'circular' or self.use_circular_orbit:
            print('    orbital model: circular')
        else:
            print("    orbital model: ", self.orbit)
        print("    orbital parametrization: ", self.parametrization)
        if self.use_time_inferior_conjunction:
            print("    time of inferior conjunction replacing mean longitude as a free parameter: ", self.use_time_inferior_conjunction)
        else:
            print('    mean longitude as free parameter (non-transiting planets or dynamical model): ', True)

        print("    default longitude of ascending node (when fixed): ", self.default_Omega)

        if self.use_inclination:
            print("    using inclination - no alternatives -: ", self.use_inclination)
        if self.use_mass :
            print("    using exomoon mass - no alternatives: ", self.use_mass)

        if self.compute_semimajor_axis_from_mass:
            print("    planetary semi-major axis computed from planet/star masses: ", self.compute_semimajor_axis_from_mass)
            print("       WARNING: this flag should be true only if you are fitting for planetary masses")
            print("                and stellar density is not a free/fixed parameter (e.g. as in transit fitting)")

        if self.use_longitude_of_nodes:
            print("    longitude of ascending node as a free parameter: ", self.use_longitude_of_nodes)

        print()

    def define_derived_parameters(self):

        derived_list = []

        if 'em_e_coso' in self.sampler_parameters and  \
            'em_e_sino' in self.sampler_parameters:

            pam00_index = self.sampler_parameters['em_e_coso']
            pam01_index = self.sampler_parameters['em_e_sino']

            try:
                del self.parameter_index['em_e_coso']
                del self.parameter_index['em_e_sino']
            except:
                pass

            if 'em_e' not in self.parameter_index:
                self.transformation['em_e'] = get_2var_e
                self.parameter_index['em_e'] = [pam00_index, pam01_index]
                derived_list.append('em_e')

            if 'em_omega' not in self.parameter_index:
                self.transformation['em_omega'] = get_2var_o
                self.parameter_index['em_omega'] = [pam00_index, pam01_index]
                derived_list.append('em_omega')

        if 'em_sre_coso' in self.sampler_parameters and  \
            'em_sre_sino' in self.sampler_parameters:

            pam00_index = self.sampler_parameters['em_sre_coso']
            pam01_index = self.sampler_parameters['em_sre_sino']

            try:
                del self.parameter_index['em_sre_coso']
                del self.parameter_index['em_sre_sino']
            except:
                pass

            if 'em_e' not in self.parameter_index:
                self.transformation['em_e'] = get_2var_sre
                self.parameter_index['em_e'] = [pam00_index, pam01_index]
                derived_list.append('em_e')

            if 'em_omega' not in self.parameter_index:
                self.transformation['em_omega'] = get_2var_o
                self.parameter_index['em_omega'] = [pam00_index, pam01_index]
                derived_list.append('em_omega')

            #if 'mean_long' in self.sampler_parameters:
            #    pam00_index = self.sampler_parameters['P']
            #    pam01_index = self.sampler_parameters['mean_long']

            #    if 'Tc_Tref' not in self.parameter_index:
            #        self.transformation['Tc_Tref'] = kepler_exo.kepler_compute_deltaTc_from_meanlong
            #        self.parameter_index['Tc_Tref'] = [pam00_index, pam01_index]
            #        derived_list.append('Tc_Tref')

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
        if var_sampler == 'em_sre_coso' or var_sampler=='em_sre_sino':

            if 'em_e' in self.starts and 'em_omega' in self.starts:

                starting_point[self.sampler_parameters['em_sre_coso']] = \
                    np.sqrt(self.starts['em_e']) * np.cos(self.starts['em_omega'])
                starting_point[self.sampler_parameters['em_sre_sino']] = \
                    np.sqrt(self.starts['em_e']) * np.sin(self.starts['em_omega'])

            elif 'em_sre_coso' in self.starts and 'em_sre_sino' in self.starts:
                starting_point[self.sampler_parameters['em_sre_coso']] = self.starts['em_sre_coso']
                starting_point[self.sampler_parameters['em_sre_coso']] = self.starts['em_sre_sino']

            return True

        if var_sampler == 'em_e_coso' or var_sampler=='em_e_sino':

            if 'em_e' in self.starts and 'em_omega' in self.starts:
                starting_point[self.sampler_parameters['em_e_coso']] = \
                    self.starts['em_e'] * np.cos(self.starts['em_omega'])
                starting_point[self.sampler_parameters['em_e_sino']] = \
                    self.starts['em_e'] * np.sin(self.starts['em_omega'])

            elif 'em_e_coso' in self.starts and 'em_e_sino' in self.starts:
                starting_point[self.sampler_parameters['em_e_coso']] = self.starts['em_e_coso']
                starting_point[self.sampler_parameters['em_e_coso']] = self.starts['em_e_sino']

            return True

        return False


    def update_parameter_values_for_dynamical(self, parameter_values, input_prepend=''):

        if input_prepend == '':
            prepend = ''
        else:
            prepend = input_prepend + '__'

        if self.compute_time_inferior_conjunction:
            parameter_values[prepend+'em_Tc']= kepler_exo.kepler_compute_deltaTc_from_meanlong(
                parameter_values[prepend+'em_P'],
                parameter_values[prepend+'em_mean_long'],
                parameter_values[prepend+'em_e'],
                parameter_values[prepend+'em_omega'],
                parameter_values[prepend+'em_Omega']) + self.Tref

        if self.compute_mean_longitude:
            parameter_values[prepend+'em_mean_long'] = kepler_exo.kepler_compute_meanlong_from_deltaTc(
                parameter_values[prepend+'em_P'],
                parameter_values[prepend+'em_Tc'] - self.Tref,
                parameter_values[prepend+'em_e'],
                parameter_values[prepend+'em_omega'],
                parameter_values[prepend+'em_Omega'])

            parameter_values[prepend+'em_Tperi'] = kepler_exo.kepler_compute_deltaTperi_from_deltaTc(
                    parameter_values[prepend+'em_P'],
                    parameter_values[prepend+'em_Tc'] - self.Tref,
                    parameter_values[prepend+'em_e'],
                    parameter_values[prepend+'em_omega'])
        else:
            parameter_values[prepend+'em_Tperi'] = kepler_exo.kepler_compute_deltaTperi_from_meanlong(
                    parameter_values[prepend+'em_P'],
                    parameter_values[prepend+'em_mean_long'],
                    parameter_values[prepend+'em_e'],
                    parameter_values[prepend+'em_omega'],
                    parameter_values[prepend+'em_Omega'])