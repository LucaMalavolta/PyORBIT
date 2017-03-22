# A class for common variables (number of planets, number of sinusoids...)
# Info to unpack the variables inside emcee must be included here
# Physical effects must be included here
from common import *


class PlanetsCommonVariables:
    def __init__(self):
        self.n_planets = 0
        self.planet_name = []

        self.bounds = {}
        self.starts = {}

        self.variables = {}
        self.var_list = {}
        self.fix_list = {}

        self.prior_kind = {}
        self.prior_pams = {}

        self.inclination = {}
        self.radius = {}

        self.circular = {}
        self.dynamical = {}
        self.transit = {}
        self.dynamical_integrator = 'TRADES'

        self.list_pams = {}

        self.fixed = []
        self.nfix = 0

        self.dynamical_set = {}

        self.list_pams_kepler = {
            'P': 'LU',  # Period, log-uniform prior
            'K': 'LU',  # RV semi-amplitude, log-uniform prior
            'f': 'U',  # RV vurve phase, log-uniform prior
            'e': 'U',  # eccentricity, uniform prior - to be fixed
            'o': 'U'}  # argument of pericenter

        ''' Orbital parameters to be used in the dynamical fit '''
        self.list_pams_dynamical = {
            'P': 'LU',  # Period in days
            'M': 'LU',  # Mass in Earth masses
            'i': 'U',  # inclination in degrees
            'f': 'U',  # phase - as defined by Malavolta+2016
            'lN': 'U',  # longitude of ascending node
            'e': 'U',  # eccentricity, uniform prior - to be fixed
            'o': 'U'}  # argument of pericenter

        """
        From the batman guide:
            params = batman.TransitParams()  # object to store transit parameters
            params.t0 = 0.  # time of inferior conjunction
            params.per = 1.  # orbital period
            params.rp = 0.1  # planet radius (in units of stellar radii)
            params.a = 15.  # semi-major axis (in units of stellar radii)
            params.inc = 87.  # orbital inclination (in degrees)
            params.ecc = 0.  # eccentricity
            params.w = 90.  # longitude of periastron (in degrees)
            params.limb_dark = "nonlinear"  # limb darkening model
            params.u = [0.5, 0.1, 0.1, -0.1]  # limb darkening coefficients

            If some of these parameters are included among the planet parameters
            they will be taken from the PlanetCommonVariables class instead that from here
            Variables need to be converted from the PyORBIT standard to BATMAN standard
        """
        self.list_pams_transit = {
            'P': 'LU', # Period, log-uniform prior
            'f': 'U',  # RV curve phase, log-uniform prior
            'e': 'U',  # eccentricity, uniform prior - to be fixed
            'o': 'U',  # argument of pericenter (in radians)
            'R': 'U',  # planet radius (in units of stellar radii)
            'a': 'U',  # semi-major axis (in units of stellar radii)
            'i': 'U'   # orbital inclination (in degrees)
        }

        """These default boundaries are used when the user does not define them in the yaml file"""
        self.default_bounds = {
            'P': [0.4, 100000.0],
            'K': [0.5, 2000.0],
            'f': [0.0, 2 * np.pi],
            'ecoso': [-1.0, 1.0],
            'esino': [-1.0, 1.0],
            'e': [0.0, 1.0],
            'o': [0.0, 2 * np.pi],
            # Used by TTVfast/TRADES
            'M': [0.5, 10000],  # Fix the unit
            'i': [0.0, 180.0],
            'lN': [0.0, 2 * np.pi],
            # Used by BATMAN
            'R': [0.00001, 0.5],  # Fix the unit
            'a': [0.00001, 50.]  # Fix the unit
        }

        """
        Addition for the transit part
        """
        self.limbdark = "nonlinear"
        self.list_pams_limbdark = {
            'uniform': {'U0': 'U'},
            'linear': {'U0': 'U'},
            'quadratic': {'U0': 'U', 'U1': 'U'},
            'square-root': {'U0': 'U', 'U1': 'U'},
            'logarithmic': {'U0': 'U', 'U1': 'U'},
            'exponential': {'U0': 'U', 'U1': 'U'},
            'nonlinear': {'U0': 'U', 'U1': 'U', 'U2': 'U', 'U3': 'U'}
        }
        # ld_options = ["uniform", "linear", "quadratic", "nonlinear"]
        # params.limb_dark = "nonlinear"  # limb darkening model
        # params.u = [0.5, 0.1, 0.1, -0.1]  # limb darkening coefficients
        self.limbdark_bounds = {
            'uniform': {'U0': [-1.00,1.00]},
            'linear': {'U0': [-1.00,1.00]},
            'quadratic': {'U0': [-1.00,1.00], 'U1': [-1.00,1.00]},
            'square-root': {'U0': [-1.00,1.00], 'U1': [-1.00,1.00]},
            'logarithmic': {'U0': [-1.00,1.00], 'U1': [-1.00,1.00]},
            'exponential': {'U0': [-1.00,1.00], 'U1': [-1.00,1.00]},
            'nonlinear': {'U0': [-1.00,1.00], 'U1': [-1.00,1.00], 'U2': [-1.00,1.00], 'U3': [-1.00,1.00]}
        }

    def add_planet(self, name_ref):
        self.n_planets += 1
        self.planet_name.append(name_ref)

        self.list_pams[name_ref] = self.list_pams_kepler
        self.starts[name_ref] = {}
        self.fix_list[name_ref] = {}
        self.var_list[name_ref] = {}
        self.variables[name_ref] = {}

        self.bounds[name_ref] = {'e': [0.00, 1.00]}

        self.prior_kind[name_ref] = {}
        self.prior_pams[name_ref] = {}

        self.inclination[name_ref] = [90.000, 0.000]

    def switch_to_circular(self, name_ref):
        self.fix_list[name_ref]['e'] = 0.00000
        self.fix_list[name_ref]['o'] = np.pi / 2.
        self.circular[name_ref] = True

    def switch_to_dynamical(self, name_ref):
        self.list_pams[name_ref] = self.list_pams_dynamical
        self.dynamical[name_ref] = True

    # AAAAAAHHHHHH
    def switch_on_transit(self, name_ref):
        for var in self.list_pams_transit:
            self.list_pams[var] = self.list_pams_transit[var]
        for var in self.list_pams_limbdark[self.limbdark]:
            self.list_pams[var] = self.list_pams_limbdark[self.limbdark][var]
            self.default_bounds[var] = self.limbdark_bounds[self.limbdark][var]

        self.transit[name_ref] = True

    def define_bounds(self, mc):
        """ Bounds are defined in this class, where all the Planet-related variables are stored"""

        for pl_name in self.planet_name:
            mc.variable_list[pl_name] = {}
            ndim_buffer = mc.ndim

            for var in self.list_pams[pl_name]:
                '''We check for each parameter (except eccentricity and omega) if the variable is a
                    fixed value or a free variable, and move the parameter into the requested space
                    Notice that 'e' and 'w' are not yet included in list_pams[pl_name] at this stage
                '''
                if "e" == var or var == "o": continue
                if var in self.fix_list[pl_name]:
                    self.variables[pl_name][var] = get_fix_val
                    self.var_list[pl_name][var] = self.nfix
                    self.fixed.append(self.fix_list[pl_name][var])
                    self.nfix += 1
                else:
                    '''If no bounds have been specified in the input file, we use the default ones
                        Bounds must be provided in any case to avoid a failure of PyDE '''
                    if var in self.bounds[pl_name]:
                        bounds_tmp = self.bounds[pl_name][var]
                    else:
                        bounds_tmp = self.default_bounds[var]

                    if self.list_pams[pl_name][var] == 'U':
                        self.variables[pl_name][var] = get_var_val
                        mc.bounds_list.append(bounds_tmp)
                    elif self.list_pams[pl_name][var] == 'LU':
                        self.variables[pl_name][var] = get_var_exp
                        mc.bounds_list.append(np.log2(bounds_tmp))

                    self.var_list[pl_name][var] = mc.ndim
                    mc.variable_list[pl_name][var] = mc.ndim
                    mc.ndim += 1

            '''eccentricity and argument of pericenter require a special treatment
             since they can be provided as fixed individual values or may need to be combined
             in ecosw and esinw if are both free variables'''
            if 'e' in self.fix_list[pl_name] and 'o' in self.fix_list[pl_name]:
                self.variables[pl_name]['e'] = get_fix_val
                self.var_list[pl_name]['e'] = self.nfix
                self.fixed.append(self.fix_list[pl_name]['e'])
                self.nfix += 1
                self.variables[pl_name]['o'] = get_fix_val
                self.var_list[pl_name]['o'] = self.nfix
                self.fixed.append(self.fix_list[pl_name]['o'])
                self.nfix += 1
            elif 'e' in self.fix_list[pl_name]:
                self.variables[pl_name]['e'] = get_fix_val
                self.var_list[pl_name]['e'] = self.nfix
                self.fixed.append(self.fix_list[pl_name]['e'])
                self.nfix += 1
                self.variables[pl_name]['o'] = get_var_val
                self.var_list[pl_name]['o'] = mc.ndim
                mc.variable_list[pl_name]['o'] = mc.ndim
                if 'o' in self.bounds[pl_name]:
                    bounds_tmp = self.bounds[pl_name]['o']
                else:
                    bounds_tmp = self.default_bounds['o']
                mc.bounds_list.append(bounds_tmp)
                mc.ndim += 1
            elif 'o' in self.fix_list[pl_name]:
                self.variables[pl_name]['o'] = get_fix_val
                self.var_list[pl_name]['o'] = self.nfix
                self.fixed.append(self.fix_list[pl_name]['o'])
                self.nfix += 1
                self.variables[pl_name]['e'] = get_var_val
                self.var_list[pl_name]['e'] = mc.ndim
                mc.variable_list[pl_name]['e'] = mc.ndim
                if 'e' in self.bounds[pl_name]:
                    bounds_tmp = self.bounds[pl_name]['e']
                else:
                    bounds_tmp = self.default_bounds['e']
                mc.bounds_list.append(bounds_tmp)
                mc.ndim += 1
            else:
                self.variables[pl_name]['e'] = get_2var_e
                self.var_list[pl_name]['e'] = [mc.ndim, mc.ndim + 1]
                self.variables[pl_name]['o'] = get_2var_o
                self.var_list[pl_name]['o'] = [mc.ndim, mc.ndim + 1]
                mc.variable_list[pl_name]['ecoso'] = mc.ndim
                mc.variable_list[pl_name]['esino'] = mc.ndim + 1
                mc.bounds_list.append(self.default_bounds['ecoso'])
                mc.bounds_list.append(self.default_bounds['esino'])
                mc.ndim += 2

            mc.variable_list[pl_name]['kepler_pams'] = np.arange(ndim_buffer, mc.ndim, 1)

    def starting_point(self, mc):

        """Default values are already set in the array"""
        for pl_name in self.planet_name:
            if pl_name not in self.starts: continue

            for var in self.list_pams[pl_name]:
                if "e" == var or var == "o": continue
                if var in self.starts[pl_name]:
                    if self.list_pams[pl_name][var] == 'U':
                        start_converted = self.starts[pl_name][var]
                    elif self.list_pams[pl_name][var] == 'LU':
                        start_converted = np.log2(self.starts[pl_name][var])

                    mc.starting_point[mc.variable_list[pl_name][var]] = start_converted

            '''eccentricity and argument of pericenter require a special treatment
             since they can be provided as fixed individual values or may need to be combined
             in ecosw and esinw if are both free variables'''
            if 'e' in self.fix_list[pl_name] and 'o' in self.fix_list[pl_name]:
                pass
            elif 'e' in self.fix_list[pl_name]:
                if 'o' in self.starts[pl_name]:
                    mc.starting_point[mc.variable_list[pl_name]['o']] = self.starts[pl_name]['o']

            elif 'o' in self.fix_list[pl_name]:
                if 'e' in self.starts[pl_name]:
                    mc.starting_point[mc.variable_list[pl_name]['e']] = self.starts[pl_name]['e']
            else:
                if 'e' in self.starts[pl_name] and 'o' in self.starts[pl_name]:
                    mc.starting_point[mc.variable_list[pl_name]['ecoso']] = \
                        np.sqrt(self.starts[pl_name]['e']) * np.cos(self.starts[pl_name]['o'])
                    mc.starting_point[mc.variable_list[pl_name]['esino']] = \
                        np.sqrt(self.starts[pl_name]['e']) * np.sin(self.starts[pl_name]['o'])
                elif 'ecoso' in self.starts[pl_name] and 'esino' in self.starts[pl_name]:
                    mc.starting_point[mc.variable_list[pl_name]['ecoso']] = self.starts[pl_name]['ecoso']
                    mc.starting_point[mc.variable_list[pl_name]['esino']] = self.starts[pl_name]['esino']

    def return_priors(self, pl_name, theta):
        prior_out = 0.00
        kep_pams = self.convert(pl_name, theta)
        for key in self.prior_pams[pl_name]:
            prior_out += giveback_priors(self.prior_kind[pl_name][key], self.prior_pams[pl_name][key], kep_pams[key])
        return prior_out

    def convert(self, pl_name, theta):
        dict_out = {}
        for key in self.list_pams[pl_name]:
            dict_out[key] = (self.variables[pl_name][key](theta, self.fixed, self.var_list[pl_name][key]))
        return dict_out

    def initialize(self, mc):
        for pl_name in self.planet_name:
            for key in mc.variable_list[pl_name]:
                if key != 'kepler_pams':
                    mc.pam_names[mc.variable_list[pl_name][key]] = key

    def print_vars(self, mc, theta):
        for pl_name in self.planet_name:
            out_list = self.convert(pl_name, theta)

            print pl_name, ' vars: ', np.asarray(theta[mc.variable_list[pl_name]['kepler_pams']])
            print pl_name, ' pams: ', out_list
        print

    def kepler_from_dynamical(self, mc, theta, planet_name):
        dict_pams = self.convert(planet_name, theta)
        return {
            'P': dict_pams['P'],
            'K': kp.kepler_K1(
                mc.star_mass[0], dict_pams['M']/mc.M_SEratio, dict_pams['P'], dict_pams['i'], dict_pams['e']),
            'f': dict_pams['f'],
            'e': dict_pams['e'],
            'o': dict_pams['o']}
