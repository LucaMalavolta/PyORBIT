from abstract_common import *


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
    parametrization = 'Eastman2013'
    parametrization_list = ['Ford2006', 'Eastman2013', 'Standard']

    list_pams = {
        'P',  # Period, log-uniform prior
        'K',  # RV semi-amplitude, log-uniform prior
        'f',  # mean longitude = argument of pericenter + mean anomaly at Tref
        'e',  # eccentricity, uniform prior - to be fixed
        'o',  # argument of pericenter (in radians)
        'M',  # Mass in Earth masses
        'i',  # orbital inclination (in degrees)
        'lN',  # longitude of ascending node (usually 180 degrees when unknown)
        'R',  # planet radius (in units of stellar radii)
        'a'  # semi-major axis (in units of stellar radii)
    }

    default_bounds = {
        'P': [0.4, 100000.0],
        'K': [0.5, 2000.0],
        'f': [0.0, 2 * np.pi],
        'e_coso': [-1.0, 1.0],
        'e_sino': [-1.0, 1.0],
        'sre_coso': [-1.0, 1.0],
        'sre_sino': [-1.0, 1.0],
        'e': [0.0, 1.0],
        'o': [0.0, 2 * np.pi],
        # Used by TTVfast/TRADES
        'M': [0.5, 1000.0],  # Fix the unit
        'i': [0.0, 180.0],
        'lN': [0.0, 2 * np.pi],
        # Used by BATMAN
        'R': [0.00001, 0.5],  # Fix the unit
        'a': [0.00001, 50.]  # Fix the unit
    }

    """ Must be the same parameters as in list_pams, because priors are applied only to _physical_ parameters """
    default_priors = {
        'P': ['Uniform', []],
        'K': ['Uniform', []],
        'f': ['Uniform', []],
        'e': ['BetaDistribution', [0.71, 2.57]],
        'e_coso': ['Uniform', []],
        'e_sino': ['Uniform', []],
        'sre_coso': ['Uniform', []],
        'sre_sino': ['Uniform', []],
        'o': ['Uniform', []],
        'M': ['Uniform', []],  # Fix the unit
        'i': ['Uniform', []],
        'lN': ['Uniform', []],
        'R': ['Uniform', []],  # Fix the unit
        'a': ['Uniform', []]  # Fix the unit
    }

    default_spaces = {
        'P': 'Logarithmic',
        'K': 'Logarithmic',
        'f': 'Linear',
        'e_coso': 'Linear',
        'e_sino': 'Linear',
        'sre_coso': 'Linear',
        'sre_sino': 'Linear',
        'e': 'Linear',
        'o': 'Linear',
        'M': 'Linear',
        'i': 'Linear',
        'lN': 'Linear',
        'R': 'Linear',
        'a': 'Linear'
    }

    omega_star = True

    recenter_pams = {'f', 'o', 'lN'}

    # Variable used only by TRADES
    period_average = None

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

        if not(var == "e" or var == "o"):
            return ndim, output_lists, False

        if 'e' in self.fix_list or \
           'o' in self.fix_list:
            return ndim, output_lists, False

        for var_check in ['e', 'o', 'e_coso', 'e_sino', 'sre_coso', 'sre_sino']:
            if var_check in self.variable_sampler:
                return ndim, output_lists, False

        if self.parametrization == 'Standard:':
            self.transformation['e'] = get_var_val
            self.variable_index['e'] = ndim
            self.transformation['o'] = get_var_val
            self.variable_index['o'] = ndim + 1
            variable_list = ['e', 'o']

        else:
            if self.parametrization == 'Ford2006':
                self.transformation['e'] = get_2var_e
                self.variable_index['e'] = [ndim, ndim + 1]
                variable_list = ['e_coso', 'e_sino']
            else:
                # 'Eastman2013' is the standard choice
                self.transformation['e'] = get_2var_sre
                self.variable_index['e'] = [ndim, ndim + 1]
                variable_list = ['sre_coso', 'sre_sino']

            self.transformation['o'] = get_2var_o
            self.variable_index['o'] = [ndim, ndim + 1]

        for var in variable_list:

            if var not in self.bounds:
                self.bounds[var] = self.default_bounds[var]

            self.spaces[var] = self.default_spaces[var]

            output_lists['bounds'].append(self.bounds[var])

            if var not in self.prior_pams:
                self.prior_kind[var] = self.default_priors[var][0]
                self.prior_pams[var] = self.default_priors[var][1]

            output_lists['spaces'].append(self.spaces[var])
            output_lists['priors'].append([self.prior_kind[var], self.prior_pams[var]])
            output_lists['nested'].append(nested_sampling_prior_transformation(
                    self.prior_kind[var],
                    output_lists['bounds'][-1],
                    self.prior_pams[var],
                )
            )

            self.variable_sampler[var] = ndim
            ndim += 1

        for var in ['e', 'o']:
            if var not in self.prior_pams:

                if var in self.bounds:
                    self.prior_pams[var] = self.bounds[var]
                else:
                    self.prior_pams[var] = self.default_bounds[var]

                self.prior_kind[var] = 'Uniform'

        return ndim, output_lists, True

    def define_special_starting_point(self, starting_point, var):
        """
        Eccentricity and argument of pericenter require a special treatment

        since they can be provided as fixed individual values or may need to be combined
        in :math:`\sqrt{e}\sin{\omega}` and :math:`\sqrt{e}\cos{\omega}` if are both free variables

        Args:
            :starting_point:
            :var:
        Returns:
            :bool:

        """

        if not (var == "e" or var == "o"):
            return False

        if 'sre_coso' in self.variable_sampler and \
                        'sre_sino' in self.variable_sampler:

            if 'e' in self.starts and 'o' in self.starts:
                starting_point[self.variable_sampler['sre_coso']] = \
                    np.sqrt(self.starts['e']) * np.cos(self.starts['o'])
                starting_point[self.variable_sampler['sre_sino']] = \
                    np.sqrt(self.starts['e']) * np.sin(self.starts['o'])

            elif 'sre_coso' in self.starts and 'sre_sino' in self.starts:
                starting_point[self.variable_sampler['sre_coso']] = self.starts['sre_coso']
                starting_point[self.variable_sampler['sre_coso']] = self.starts['sre_sino']

        if 'e_coso' in self.variable_sampler and \
                        'e_sino' in self.variable_sampler:

            if 'e' in self.starts and 'o' in self.starts:
                starting_point[self.variable_sampler['e_coso']] = \
                    self.starts['e'] * np.cos(self.starts['o'])
                starting_point[self.variable_sampler['e_sino']] = \
                    self.starts['e'] * np.sin(self.starts['o'])

            elif 'e_coso' in self.starts and 'e_sino' in self.starts:
                starting_point[self.variable_sampler['e_coso']] = self.starts['e_coso']
                starting_point[self.variable_sampler['e_coso']] = self.starts['e_sino']


        return True

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
                for ii in xrange(0, n_pop):
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
                for ii in xrange(0, n_pop):
                    if not self.bounds['e'][0] + 0.02 <= e_pops[ii] < \
                                    self.bounds['e'][1] - 0.02:
                        e_random = np.random.uniform(self.bounds['e'][0],
                                                     self.bounds['e'][1])
                        population[ii, sre_sino_list] = np.sqrt(e_random) * np.sin(o_pops[ii])
                        population[ii, sre_coso_list] = np.sqrt(e_random) * np.cos(o_pops[ii])