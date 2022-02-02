from pyorbit.subroutines.common import *
from pyorbit.common.abstract_common import *


class CommonPlanets(AbstractCommon):
    """
    Inherited class from AbstractCommon

    For computational reason it is better to fit :math:`\sqrt{e}\sin{\omega}` and :math:`\sqrt{e}\cos{\omega}`.
    :func:`define_special_variable_properties` and :func:`define_special_starting_point` must be redefined

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

    list_pams = {
        'P',  # Period, log-uniform prior
        'K',  # RV semi-amplitude, log-uniform prior
        'Tc', # central time of transit
        'mean_long',  # mean longitude = argument of pericenter + mean anomaly at Tref
        'e',  # eccentricity, uniform prior - to be fixed
        'omega',  # argument of pericenter (in degrees)
        'M_Me',  # Mass in Earth masses
        'i',  # orbital inclination (in degrees)
        'Omega', # longitude of ascending node (usually 180 degrees when unknown)
        'R_Rs',  # planet radius (in units of stellar radii)
        'a_Rs',  # semi-major axis (in units of stellar radii)
        'b',   # impact parameter
        'lambda', # Sky-projected angle between stellar rotation axis and normal of orbit plane [deg]
        'phase_amp', #Amplitude of the phase light curve
        'delta_occ', #depth of the occultation, as measured at the maximum value of the reflected light curve
        'phase_off', #offset of the bright spot
        'albedo',
        'redist',
        'insol'
    }

    default_bounds = {
        'P': [0.4, 100000.0],
        'K': [0.5, 2000.0],
        'Tc': [0.0, 1000.0],
        'mean_long': [0.0, 360.],
        'e_coso': [-1.0, 1.0],
        'e_sino': [-1.0, 1.0],
        'sre_coso': [-1.0, 1.0],
        'sre_sino': [-1.0, 1.0],
        'e': [0.0, 1.0],
        'omega': [0.0, 360.],
        # Used by TTVfast/TRADES
        'M_Me': [0.5, 1000.0],  # Fix the unit
        'i': [0.0, 180.0],
        'Omega': [0.0, 180],
        # Used by BATMAN
        'R_Rs': [0.00001, 0.5],  # Fix the unit
        'a_Rs': [0.00001, 500.],  # Fix the unit
        'b': [0.0, 2.0],
        'lambda': [-90.0000, 90.0000],
        'phase_amp': [0.00, 0.50],
        'delta_occ': [0.00, 0.50],
        'phase_off': [-np.pi, np.pi],
        'albedo': [0., 1.],
        'redist': [0., 1.],
        'insol': [0, 1000000000]
    }

    """ Must be the same parameters as in list_pams, because priors are
    applied only to _physical_ parameters """
    default_priors = {
        'P': ['Uniform', []],
        'K': ['Uniform', []],
        'Tc': ['Uniform', []],
        'mean_long': ['Uniform', []],
        #'e': ['BetaDistribution', [0.71, 2.57]],
        'e': ['Uniform', []],
        'e_coso': ['Uniform', []],
        'e_sino': ['Uniform', []],
        'sre_coso': ['Uniform', []],
        'sre_sino': ['Uniform', []],
        'omega': ['Uniform', []],
        'M_Me': ['Uniform', []],  # Fix the unit
        'i': ['Uniform', []],
        'Omega': ['Uniform', []],
        'R_Rs': ['Uniform', []],  # Fix the unit
        'a_Rs': ['Uniform', []],  # Fix the unit
        'b': ['Uniform', []],  # Fix the unit
        'lambda': ['Uniform', []],  # Fix the unit
        'phase_amp': ['Uniform', []],  # Fix the unit
        'delta_occ': ['Uniform', []],  # Fix the unit
        'phase_off': ['Uniform', []],  # Fix the unit
        'albedo': ['Uniform', []],  # Fix the unit
        'redist': ['Uniform', []],  # Fix the unit
        'insol': ['Uniform', []],  # Fix the unit
        }

    default_spaces = {
        'P': 'Log_Base2',
        'K': 'Log_Base2',
        'Tc': 'Linear',
        'mean_long': 'Linear',
        'e_coso': 'Linear',
        'e_sino': 'Linear',
        'sre_coso': 'Linear',
        'sre_sino': 'Linear',
        'e': 'Linear',
        'omega': 'Linear',
        'M_Me': 'Log_Base2',
        'i': 'Linear',
        'Omega': 'Linear',
        'R_Rs': 'Linear',
        'a_Rs': 'Linear',
        'b': 'Linear',
        'lambda': 'Linear',
        'phase_amp': 'Linear',
        'delta_occ': 'Linear',
        'phase_off': 'Linear',
        'albedo': 'Linear',
        'redist': 'Linear',
        'insol': 'Log_Base2',
    }

    default_fixed = {
        'e_coso': 0.0000,
        'e_sino': 0.0000,
        'sre_coso': 0.0000,
        'sre_sino': 0.0000,
        'e': 0.0000,
        'omega': 90.,
        'i': 90.000000,
        'Omega': 90.,
        'R_Rs': 0.05,
        'a_Rs': 1.0,
        'b': 1.0,
        'lambda': 0.0000,
        'phase_amp': 0.000,
        'delta_occ': 0.000,
        'phase_off': 0.000,
        'albedo': 0.000,
        'redist': 0.000,
        'insol': 1000000.000,
    }

    recenter_pams = {'mean_long', 'omega', 'Omega'}

    def __init__(self, *args, **kwargs):
        super(CommonPlanets, self).__init__(*args, **kwargs)

        self.orbit = 'keplerian'
        self.parametrization = 'Eastman2013'

        self.use_inclination = False
        self.use_semimajor_axis = False
        self.use_time_of_transit = False
        self.use_mass_for_planets = False

        self.omega_star = True

        self.period_average = None
        # Variable used only by TRADES

    def define_special_variable_properties(self, ndim, output_lists, var):
        """ Boundaries definition for eccentricity :math:`e` and argument of pericenter :math:`\omega`

        The internal variable to be fitted are :math:`\sqrt{e}\sin{\omega}` and :math:`\sqrt{e}\cos{\omega}`.
        With this parametrization it is not possible to naturally put a boundary to :math:`e` without affecting the
        :math:`\omega`.
        Additionally the subroutine will check if either :math:`e` or :math:`\omega` have been provided as fixed values.
        If true, the parametrization will consist of :math:`e` or :math:`\omega`  instead of
        :math:`\sqrt{e}\sin{\omega}` and :math:`\sqrt{e}\cos{\omega}`

        Args:
            :ndim: number of parameters already processed by other models
            :var: input variable, either :math:`e` or :math:`\omega`
        Returns:
            :ndim: updated dimensionality of the problem
            :bounds_list: additional boundaries to be added to the original list
        """

        if var == 'P':
            if var in self.fix_list:
                self.period_average = self.fix_list['P'][0]
            else:
                if var in self.bounds:
                    self.period_average = np.average(self.bounds[var])
                else:
                    self.period_average = np.average(self.default_bounds[var])
            return ndim, output_lists, False

        if not(var == "e" or var == "omega"):
            return ndim, output_lists, False

        if 'e' in self.fix_list or \
           'omega' in self.fix_list:
            return ndim, output_lists, False

        for var_check in ['e', 'omega', 'e_coso', 'e_sino', 'sre_coso', 'sre_sino']:
            if var_check in self.variable_sampler:
                return ndim, output_lists, False

        if self.parametrization[:8] == 'Standard':
            self.transformation['e'] = get_var_val
            self.variable_index['e'] = ndim
            self.transformation['omega'] = get_var_val
            self.variable_index['omega'] = ndim + 1
            variable_list = ['e', 'omega']

        else:
            if self.parametrization[:8] == 'Ford2006':
                self.transformation['e'] = get_2var_e
                self.variable_index['e'] = [ndim, ndim + 1]
                variable_list = ['e_coso', 'e_sino']
            else:
                # 'Eastman2013' is the standard choice
                self.transformation['e'] = get_2var_sre
                self.variable_index['e'] = [ndim, ndim + 1]
                variable_list = ['sre_coso', 'sre_sino']

            self.transformation['omega'] = get_2var_o
            self.variable_index['omega'] = [ndim, ndim + 1]

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

        for var in ['e', 'omega']:
            if var not in self.bounds:
                self.bounds[var] = self.default_bounds[var]

            if var not in self.prior_pams:

                if var in self.bounds:
                    self.prior_pams[var] = self.bounds[var]
                else:
                    self.prior_pams[var] = self.default_bounds[var]

                self.prior_kind[var] = 'Uniform'

        return ndim, output_lists, True

    def define_special_starting_point(self, starting_point, var_sampler):
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

                starting_point[self.variable_sampler['sre_coso']] = \
                    np.sqrt(self.starts['e']) * np.cos(self.starts['omega'])
                starting_point[self.variable_sampler['sre_sino']] = \
                    np.sqrt(self.starts['e']) * np.sin(self.starts['omega'])

            elif 'sre_coso' in self.starts and 'sre_sino' in self.starts:
                starting_point[self.variable_sampler['sre_coso']] = self.starts['sre_coso']
                starting_point[self.variable_sampler['sre_coso']] = self.starts['sre_sino']

            return True

        if var_sampler == 'e_coso' or var_sampler=='e_sino':

            if 'e' in self.starts and 'omega' in self.starts:
                starting_point[self.variable_sampler['e_coso']] = \
                    self.starts['e'] * np.cos(self.starts['omega'])
                starting_point[self.variable_sampler['e_sino']] = \
                    self.starts['e'] * np.sin(self.starts['omega'])

            elif 'e_coso' in self.starts and 'e_sino' in self.starts:
                starting_point[self.variable_sampler['e_coso']] = self.starts['e_coso']
                starting_point[self.variable_sampler['e_coso']] = self.starts['e_sino']

            return True

        return False

    def special_fix_population(self, population):

        n_pop = np.size(population, axis=0)
        if 'e_sino' in self.variable_sampler and \
           'e_coso' in self.variable_sampler:
                e_sino_list = self.variable_sampler['e_sino']
                e_coso_list = self.variable_sampler['e_coso']
                e_pops = np.sqrt(population[:, e_sino_list] ** 2 + population[:, e_coso_list] ** 2)
                o_pops = np.arctan2(population[:, e_sino_list], population[:, e_coso_list], dtype=np.double)
                # e_mean = (self[planet_name].bounds['e'][0] +
                # self[planet_name].bounds['e'][1]) / 2.
                for ii in range(0, n_pop):
                    if not self.bounds['e'][0] + 0.02 <= e_pops[ii] < \
                                    self.bounds['e'][1] - 0.02:
                        e_random = np.random.uniform(self.bounds['e'][0],
                                                     self.bounds['e'][1])
                        population[ii, e_sino_list] = e_random * np.sin(o_pops[ii])
                        population[ii, e_coso_list] = e_random * np.cos(o_pops[ii])

        if 'sre_sino' in self.variable_sampler and \
           'sre_coso' in self.variable_sampler:
                sre_sino_list = self.variable_sampler['sre_sino']
                sre_coso_list = self.variable_sampler['sre_coso']
                e_pops = population[:, sre_sino_list] ** 2 + population[:, sre_coso_list] ** 2
                o_pops = np.arctan2(population[:, sre_sino_list], population[:, sre_coso_list], dtype=np.double)
                # e_mean = (self[planet_name].bounds['e'][0] +
                # self[planet_name].bounds['e'][1]) / 2.
                for ii in range(0, n_pop):
                    if not self.bounds['e'][0] + 0.02 <= e_pops[ii] < \
                                    self.bounds['e'][1] - 0.02:
                        e_random = np.random.uniform(self.bounds['e'][0],
                                                     self.bounds['e'][1])
                        population[ii, sre_sino_list] = np.sqrt(e_random) * np.sin(o_pops[ii])
                        population[ii, sre_coso_list] = np.sqrt(e_random) * np.cos(o_pops[ii])