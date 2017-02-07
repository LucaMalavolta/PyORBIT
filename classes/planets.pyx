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
        self.dynamical_integrator = 'TRADES'

        self.list_pams = {}

        self.fixed = []
        self.nfix = 0

        self.dynamical_set = {}

        self.list_pams_default = {
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
            'lN': [0.0, 2 * np.pi]}

    def add_planet(self, name_ref):
        self.n_planets += 1
        self.planet_name.append(name_ref)

        self.list_pams[name_ref] = self.list_pams_default
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

    def define_bounds(self, mc):

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

        # return self.variables[pl_name]['P'](theta, fixed, ), self.variables[pl_name]['K'](theta, fixed, 2), \
        #       self.variables[pl_name]['f'](theta, fixed, 2), self.variables[pl_name]['e'](theta, fixed, 2), \
        #       self.variables[pl_name]['o'](theta, fixed, 2)

    def compute(self, theta, dataset, pl_name):
        dict_pams = self.convert(pl_name, theta)
        return kp.kepler_RV_T0P(dataset.x0, dict_pams['f'], dict_pams['P'], dict_pams['K'], dict_pams['e'],
                                dict_pams['o'])

    def print_vars(self, mc, theta):
        for pl_name in self.planet_name:
            out_list = self.convert(pl_name, theta)

            for key in mc.variable_list[pl_name]:
                if key != 'kepler_pams':
                    mc.pam_names[mc.variable_list[pl_name][key]] = key
            print pl_name, ' vars: ', np.asarray(theta[mc.variable_list[pl_name]['kepler_pams']])
            print pl_name, ' pams: ', out_list

    # Compute the expected RV from a single planet, assuming non-interacting orbits in multi-planet systems
    @staticmethod
    def model_kepler(orbit_pams, x0):
        e = orbit_pams['e']
        if e > 1.0:
            e = 1.0
        rv_out = kp.kepler_RV_T0P(x0,
                                  orbit_pams['f'],
                                  orbit_pams['P'],
                                  orbit_pams['K'],
                                  e,
                                  orbit_pams['o'])
        return rv_out

    def prepare_dynamical(self, mc):
        """
        Prepare the input parameters according to chosen dynamical integrator
        :param mc:
        :return:
        """

        if self.dynamical_integrator == 'TRADES':
            self.prepare_dynamical_trades(mc)
        if self.dynamical_integrator == 'ttvfast':
            self.prepare_dynamical_ttvfast(mc)
        return

    def compute_dynamical(self, mc,  *args, **kwargs):
        """
        Run the appropriate subroutine according to chosen dynamical integrator
        :param mc:
        :return:
        """
        if self.dynamical_integrator == 'TRADES':
            output = self.compute_dynamical_trades(mc, *args, **kwargs)
        if self.dynamical_integrator == 'ttvfast':
            output = self.compute_dynamical_ttvfast(mc, *args, **kwargs)
        return output

    def prepare_dynamical_trades(self, mc):
        """
        :param mc:
        :return:
        """
        dataset_rv = 0
        int_buffer = dict(rv_times=[], t0_times=[], rv_ref=[], t0_ref=[], key_ref={})

        """ Putting all the RV epochs in the same array, flagging in the temporary buffer
            the stored values according to their dataset of origin
        """
        for dataset in mc.dataset_list:
            if dataset.kind == 'RV':
                int_buffer['rv_times'].extend(dataset.x.tolist())
                int_buffer['rv_ref'].extend(dataset.x * 0 + dataset_rv)
                int_buffer['key_ref'][dataset.name_ref] = dataset_rv
                dataset_rv += 1
            elif dataset.kind == 'Tcent':
                int_buffer['t0_times'].extend(dataset.x.tolist())

        """ Creating the flag array after all the RV epochs have been mixed
        """
        self.dynamical_set['data'] = {'selection': {}}
        for dataset in mc.dataset_list:
            if dataset.kind == 'RV':
                self.dynamical_set['data']['selection'][dataset.name_ref] = \
                    (np.asarray(int_buffer['rv_ref']) == int_buffer['key_ref'][dataset.name_ref])

        self.dynamical_set['data']['rv_times'] = np.float64(int_buffer['rv_times'])

        """ Substituting empty arrays with None to allow minimum determination
            We give for granted that at least one of the two vector is not null, otherwise
            it doesn't make any sense to run the program at all
        """
        if np.size(int_buffer['rv_times']) == 0:
            int_buffer['rv_times'] = int_buffer['t0_times']

        if np.size(int_buffer['t0_times']) == 0:
            int_buffer['t0_times'] = int_buffer['rv_times']

        rv_minmax = [np.amin(int_buffer['rv_times']), np.amax(int_buffer['rv_times'])]
        t0_minmax = [np.amin(int_buffer['t0_times']), np.amax(int_buffer['t0_times'])]

        self.dynamical_set['trades'] = {
            'ti_beg': np.min([rv_minmax[0], t0_minmax[0]]) - 10,
            'ti_end': np.max([rv_minmax[1], t0_minmax[1]]) + 10}
        self.dynamical_set['trades']['ti_int'] = \
            self.dynamical_set['trades']['ti_end'] - self.dynamical_set['trades']['ti_beg']
        self.dynamical_set['trades']['ti_ref'] = np.float64(mc.Tref)
        self.dynamical_set['trades']['i_step'] = np.float64(1.e-3)

        self.dynamical_set['data']['plan_ref'] = {}

        """First iteration to identify the number of transits
        stored for each planet, including the planets in the dynamical
        simulation but without observed transit
        """
        t0_ntot = [0]
        t0_flag = [False]
        n_body = 1 # The star is included among the bodies

        period_list = []
        planet_list = []
        for pl_name in self.dynamical:
            period_list.extend([(self.bounds[pl_name]['P'][1]+self.bounds[pl_name]['P'][0])/2.])
            planet_list.extend([pl_name])
        sort_planets = np.argsort(period_list)
        
        #for pl_name in self.dynamical_list:
        for pl_name in np.asarray(planet_list)[sort_planets]:
            tmp_t0_ntot = [0]
            tmp_t0_flag = [True]
            if pl_name in mc.t0_list:
                tmp_t0_ntot = [mc.t0_list[pl_name].n]
                tmp_t0_flag = [True]
            t0_ntot.extend(tmp_t0_ntot)
            t0_flag.extend(tmp_t0_flag)
            self.dynamical_set['data']['plan_ref'][pl_name] = np.asarray(n_body).astype(int)
            n_body += 1


        """ TRADES requires at least one planet to be observed transiting, since it was create primarily for
         TTV analysis. We use a workaround in case there are no T_cent observed by creating and passing to
          TRADES a fake dataset that  is not include in the chi^2 computation
        """
        if np.max(t0_ntot) == 0:
            self.dynamical_set['fake_t0s'] = True
            t0_flag[1] = True
            t0_ntot[1] = 3
        else:
            self.dynamical_set['fake_t0s'] = False

        """ Identification the maximum number of transit """
        n_max_t0 = np.max(t0_ntot)

        self.dynamical_set['data']['t0_tot'] = np.asarray(t0_ntot).astype(int)
        self.dynamical_set['data']['t0_flg'] = np.asarray(t0_flag)

        self.dynamical_set['data']['t0_num'] = np.zeros([n_max_t0, n_body]).astype(int)
        self.dynamical_set['data']['t0_obs'] = np.zeros([n_max_t0, n_body], dtype=np.float64)
        self.dynamical_set['data']['t0_err'] = np.zeros([n_max_t0, n_body], dtype=np.float64)

        for pl_name in self.dynamical:
            plan_n = self.dynamical_set['data']['plan_ref'][pl_name]
            if pl_name in mc.t0_list:
                self.dynamical_set['data']['t0_num'][0:t0_ntot[plan_n], plan_n] = mc.t0_list[pl_name].n_transit[:].astype(int)
                self.dynamical_set['data']['t0_obs'][0:t0_ntot[plan_n], plan_n] = mc.t0_list[pl_name].x
                self.dynamical_set['data']['t0_err'][0:t0_ntot[plan_n], plan_n] = mc.t0_list[pl_name].e

        if self.dynamical_set['fake_t0s']:
                self.dynamical_set['data']['t0_num'][0:3, 1] = np.arange(0, 3, 1, dtype=int)
                self.dynamical_set['data']['t0_obs'][0:3, 1] = np.arange(-1, 2, 1, dtype=np.float64)*10.0 + mc.Tref
                self.dynamical_set['data']['t0_err'][0:3, 1] = 0.1

        self.dynamical_set['trades']['n_body'] = n_body

        pytrades.args_init(self.dynamical_set['trades']['ti_beg'],
                           self.dynamical_set['trades']['ti_ref'],
                           self.dynamical_set['trades']['ti_int'],
                           self.dynamical_set['trades']['n_body'],
                           self.dynamical_set['data']['t0_tot'],
                           self.dynamical_set['data']['t0_num'],
                           self.dynamical_set['data']['t0_obs'],
                           self.dynamical_set['data']['t0_err'])

        print 'TRADES parameters', self.dynamical_set['trades']
        print

        return

    def compute_dynamical_trades(self, mc, theta, full_orbit=None):
        """ This function compute the expected TTV and RVs for dynamically interacting planets.
            The user can specify which planets are subject to interactions, e.g. long-period planets can be approximated
            with a Keplerian function"""

        """ setting up the dictionaries with the orbital parameters required by TRADES,
        """

        self.dynamical_set['pams'] = {
            'M': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64),
            'R': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64),
            'P': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64),
            'e': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64),
            'o': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64),
            'i': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64),
            'lN': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64),
            'mA': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64)
        }

        """ Adding star parameters"""
        self.dynamical_set['pams']['M'][0] = mc.star_mass[0]
        self.dynamical_set['pams']['R'][0] = mc.star_radius[0]
        
        for pl_name in self.dynamical:
            
            n_plan = self.dynamical_set['data']['plan_ref'][pl_name]
            dict_pams = self.convert(pl_name, theta)
            
            self.dynamical_set['pams']['M'][n_plan] = dict_pams['M'] / mc.M_SEratio
            self.dynamical_set['pams']['R'][n_plan] = mc.pcv.radius[pl_name][0] / mc.R_SEratio
            self.dynamical_set['pams']['P'][n_plan] = dict_pams['P']
            self.dynamical_set['pams']['e'][n_plan] = dict_pams['e']
            self.dynamical_set['pams']['o'][n_plan] = dict_pams['o'] * (180. / np.pi)
            self.dynamical_set['pams']['i'][n_plan] = dict_pams['i']
            self.dynamical_set['pams']['lN'][n_plan] = dict_pams['lN'] * (180. / np.pi)
            self.dynamical_set['pams']['mA'][n_plan] = (dict_pams['f']-dict_pams['o']) * (180. / np.pi)

        #sample_plan[:, convert_out['Tcent']] = mc.Tref + kp.kepler_Tcent_T0P(
        #    sample_plan[:, convert_out['P']], sample_plan[:, convert_out['f']],
        #    sample_plan[:, convert_out['e']], sample_plan[:, convert_out['o']])

        if full_orbit is None:
            rv_sim, t0_sim = pytrades.kelements_to_data(
                self.dynamical_set['trades']['ti_beg'],
                self.dynamical_set['trades']['ti_ref'],
                self.dynamical_set['trades']['i_step'],
                self.dynamical_set['trades']['ti_int'],
                self.dynamical_set['pams']['M'],
                self.dynamical_set['pams']['R'],
                self.dynamical_set['pams']['P'],
                self.dynamical_set['pams']['e'],
                self.dynamical_set['pams']['o'],
                self.dynamical_set['pams']['mA'],
                self.dynamical_set['pams']['i'],
                self.dynamical_set['pams']['lN'],
                self.dynamical_set['data']['rv_times'],
                self.dynamical_set['data']['t0_flg'],
                self.dynamical_set['data']['t0_tot'],
                self.dynamical_set['data']['t0_num'])
        else:
            rv_sim, t0_sim = pytrades.kelements_to_data(
                self.dynamical_set['trades']['ti_beg'],
                self.dynamical_set['trades']['ti_ref'],
                self.dynamical_set['trades']['i_step'],
                self.dynamical_set['trades']['ti_int'],
                self.dynamical_set['pams']['M'],
                self.dynamical_set['pams']['R'],
                self.dynamical_set['pams']['P'],
                self.dynamical_set['pams']['e'],
                self.dynamical_set['pams']['o'],
                self.dynamical_set['pams']['mA'],
                self.dynamical_set['pams']['i'],
                self.dynamical_set['pams']['lN'],
                full_orbit,
                self.dynamical_set['data']['t0_flg'],
                self.dynamical_set['data']['t0_tot'],
                self.dynamical_set['data']['t0_num'])

        output = {}
        t0_sim -= mc.Tref
        for dataset in mc.dataset_list:
            if dataset.kind == 'RV' and full_orbit is None:
                output[dataset.name_ref] = rv_sim[self.dynamical_set['data']['selection'][dataset.name_ref]]
                #print 'RV out', output[dataset.name_ref]
            elif dataset.kind == 'Tcent':
                n_plan = self.dynamical_set['data']['plan_ref'][dataset.planet_name]
                output[dataset.name_ref] = t0_sim[:self.dynamical_set['data']['t0_tot'][n_plan], n_plan]
        if full_orbit is not None:
            output['full_orbit'] = rv_sim

        return output

    def prepare_dynamical_ttvfast(self, mc):

        self.dynamical_set['rv_times'] = []
        self.dynamical_set['data_selection'] = {}

        dataset_rv = 0
        int_buffer = dict(rv_times=[], t0_times=[], rv_ref=[], t0_ref=[], key_ref={})

        for dataset in mc.dataset_list:
            if dataset.kind == 'RV':
                int_buffer['rv_times'].extend(dataset.x0.tolist())
                int_buffer['rv_ref'].extend(dataset.x0 * 0 + dataset_rv)
                int_buffer['key_ref'][dataset.name_ref] = dataset_rv
                dataset_rv += 1
            elif dataset.kind == 'Tcent':
                int_buffer['t0_times'].extend(dataset.x0.tolist())

        if np.size(int_buffer['t0_times']) == 0:
            int_buffer['t0_times'] = int_buffer['rv_times']

        for dataset in mc.dataset_list:
            if dataset.kind == 'RV':
                self.dynamical_set['data_selection'][dataset.name_ref] = \
                    (np.asarray(int_buffer['rv_ref']) == int_buffer['key_ref'][dataset.name_ref])

        self.dynamical_set['rv_times'] = int_buffer['rv_times']
        self.dynamical_set['len_rv'] = np.size(self.dynamical_set['rv_times'])
        if self.dynamical_set['len_rv'] == 0:
            self.dynamical_set['rv_times'] = None
            int_buffer['rv_times'] = int_buffer['t0_times']

        rv_minmax = [np.amin(int_buffer['rv_times']), np.amax(int_buffer['rv_times'])]
        t0_minmax = [np.amin(int_buffer['t0_times']), np.amax(int_buffer['t0_times'])]

        self.dynamical_set['ttvfast'] = {
            't_beg': np.min([rv_minmax[0], t0_minmax[0]]) - 10,
            't_end': np.max([rv_minmax[1], t0_minmax[1]]) + 10}

        if self.dynamical_set['ttvfast']['t_beg'] > -10. :
            """It means that both rv_minmax[0] and t0_minmax[0] are greater than zero, i.e. that both
            RV and TTV epochs start after Tref """
            self.dynamical_set['ttvfast']['t_beg'] = 0.0000

        ####printself.dynamical_set['ttvfast']['t_beg']

    def compute_dynamical_ttvfast(self, mc, theta, full_orbit=None):
        """ This function compute the expected TTV and RVs for dynamically interacting planets.
            The user can specify which planets are subject to interactions, e.g. long-period planets can be approximated
            with a Keplerian function"""

        input_flag = 0
        n_plan = 0

        P_min = None

        plan_ref = {}
        params = [mc.G_ttvfast, mc.star_mass[0]]

        for pl_name in self.dynamical:
            plan_ref[pl_name] = n_plan

            dict_pams = self.convert(pl_name, theta)

            mA = (dict_pams['f']-dict_pams['o']) * (180. / np.pi) + \
                 self.dynamical_set['ttvfast']['t_beg']/dict_pams['P'] * 360.0000000000

            params.extend([
                dict_pams['M'] / mc.M_SEratio,  # mass in Solar unit
                dict_pams['P'],
                dict_pams['e'],
                dict_pams['i'],
                dict_pams['lN'] * (180. / np.pi),
                dict_pams['o'] * (180. / np.pi),
                mA])

            n_plan += 1
            if P_min is None:
                P_min = dict_pams['P']
            else:
                P_min = np.min(np.asarray([P_min, dict_pams['P']]))

        t_step = P_min/20.
        if full_orbit is None:
            pos, rv = ttvfast._ttvfast._ttvfast(params,
                                                t_step,
                                                self.dynamical_set['ttvfast']['t_beg'],
                                                self.dynamical_set['ttvfast']['t_end'],
                                                n_plan, input_flag,
                                                self.dynamical_set['len_rv'],
                                                self.dynamical_set['rv_times'])
        else:
            pos, rv = ttvfast._ttvfast._ttvfast(params,
                                                t_step,
                                                self.dynamical_set['ttvfast']['t_beg'],
                                                self.dynamical_set['ttvfast']['t_end'],
                                                n_plan, input_flag,
                                                len(full_orbit), full_orbit.tolist())
        positions = np.asarray(pos)
        rv_meas = np.asarray(rv) * mc.AUday2ms
        output = {}
        #print 'RV  ------->  ', rv[:10]
        #print self.dynamical_set['len_rv'], rv_meas[:10]
        #print 't_beg', self.dynamical_set['ttvfast']['t_beg']
        #print 't_end', self.dynamical_set['ttvfast']['t_end']
        #print 'Full orbit flag: ', full_orbit
        #print 'Full orbit flag: ', full_orbit
        #print positions[:10,:]

        for dataset in mc.dataset_list:
            if dataset.kind == 'RV' and full_orbit is None:
                output[dataset.name_ref] = rv_meas[self.dynamical_set['data_selection'][dataset.name_ref]]

            elif dataset.kind == 'Tcent':
                t0_sel = (positions[0][:] == plan_ref[dataset.planet_name])
                t0_mod = positions[2][t0_sel]  # T0 of the transit
                nn_mod = np.asarray(positions[1][t0_sel], dtype=np.int16)  # Number of transit as stored by TTVfast

                """TTVfast transit numbering precedes the one of our dataset by construction, so
                        we need to add an offset to the transit reference number in order to have the
                        right association. We do that by identifying the first TTVfast T0 corresponding to
                        the 0th transit in our dataset"""
                min_idx = np.argmin(np.abs(t0_mod - dataset.x0[0]))
                ref_idf = nn_mod[min_idx]
                #print ref_idf,  np.sum(t0_sel), len(t0_mod), t0_mod

                if np.sum(t0_sel) > np.max(dataset.n_transit) + ref_idf:
                    output[dataset.name_ref] = t0_mod[dataset.n_transit + ref_idf]
                    """With this approach, missing T0s in the input dataset are automatically skipped """
                else:
                    """The output vector contains less T0s than supposed: collision between planets?? """
                    output[dataset.name_ref] = dataset.n_transit * 0.0000

        if full_orbit is not None:
            output['full_orbit'] = rv_meas
        #print rv_meas[:3], self.dynamical_set['rv_times'][:3]
        return output

    def kepler_from_dynamical(self, mc, theta, planet_name):
        dict_pams = self.convert(planet_name, theta)
        return {
            'P': dict_pams['P'],
            'K': kp.kepler_K1(
                mc.star_mass[0], dict_pams['M']/mc.M_SEratio, dict_pams['P'], dict_pams['i'], dict_pams['e']),
            'f': dict_pams['f'],
            'e': dict_pams['e'],
            'o': dict_pams['o']}
