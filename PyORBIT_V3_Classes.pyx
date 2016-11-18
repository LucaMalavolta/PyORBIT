import numpy as np
import kepler_exo as kp
import yaml
import george
import spotpy
import ttvfast


def get_var_log(var, fix, i):
    return np.log2(var[i], dtype=np.double)


def get_var_exp(var, fix, i):
    return np.exp2(var[i], dtype=np.double)


def get_var_val(var, fix, i):
    return var[i]


def get_fix_val(var, fix, i):
    return fix[i]


def get_2var_e(var, fix, i):
    ecoso = var[i[0]]
    esino = var[i[1]]
    return np.square(ecoso, dtype=np.double) + np.square(esino, dtype=np.double)


def get_2var_o(var, fix, i):
    ecoso = var[i[0]]
    esino = var[i[1]]
    return np.arctan2(esino, ecoso, dtype=np.double)


def giveback_priors(kind, pams, val):
    if kind == 'Gaussian':
        return -(val - pams[0]) ** 2 / (2 * pams[1] ** 2)
    if kind == 'Uniform':
        return 0.0


# A class for common variables (number of planets, number of sinusoids...)
# Info to unpack the variables inside emcee must be included here
# Physical effects must be included here
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

        self.circular = {}
        self.dynamical = {}

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

        """These default bounds are used when the user does not define them in the yaml file"""
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
                        self.starts[pl_name]['e'] * np.cos(self.starts[pl_name]['o'] * (np.pi / 180.))
                    mc.starting_point[mc.variable_list[pl_name]['esino']] = \
                        self.starts[pl_name]['e'] * np.sin(self.starts[pl_name]['o'] * (np.pi / 180.))
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

        self.dynamical_set['rv_times'] = []
        self.dynamical_set['data_selection'] = {}

        dataset_rv = 0
        int_buffer = dict(t0_times=[], rv_ref=[], t0_ref=[], key_ref={})

        for dataset in mc.dataset_list:
            if dataset.kind == 'RV':
                self.dynamical_set['rv_times'].extend(dataset.x0.tolist())
                int_buffer['rv_ref'].extend(dataset.x0 * 0 + dataset_rv)
                int_buffer['key_ref'][dataset.name_ref] = dataset_rv
                dataset_rv += 1
            elif dataset.kind == 'Tcent':
                int_buffer['t0_times'].extend(dataset.x0.tolist())

        if np.size(int_buffer['t0_times']) == 0:
            int_buffer['t0_times'] = self.dynamical_set['rv_times']

        for dataset in mc.dataset_list:
            if dataset.kind == 'RV':
                self.dynamical_set['data_selection'][dataset.name_ref] = \
                    (np.asarray(int_buffer['rv_ref']) == int_buffer['key_ref'][dataset.name_ref])

        self.dynamical_set['len_rv'] = np.size(self.dynamical_set['rv_times'])
        if self.dynamical_set['len_rv'] == 0:
            self.dynamical_set['rv_times'] = None

        '''We assumed that each planet in the dynamical integrations has well defined period boundaries
            AND that its period has not been fixed (or the whole dynamical analysis would be meaningless'''

        rv_minmax = [np.amin(self.dynamical_set['rv_times']), np.amax(self.dynamical_set['rv_times'])]
        t0_minmax = [np.amin(int_buffer['t0_times']), np.amax(int_buffer['t0_times'])]

        self.dynamical_set['ttvfast'] = {
            't_beg': np.min([rv_minmax[0], t0_minmax[0]]) - 10,
            't_end': np.max([rv_minmax[1], t0_minmax[1]]) + 10}

    def compute_dynamical(self, mc, theta, full_orbit=None):
        """ This function compute the expected TTV and RVs for dynamically interacting planets.
            The user can specify which planets are subject to interactions, e.g. long-period planets can be approximated
            with a Keplerian function"""

        input_flag = 0
        n_plan = 0

        P_min = None

        plan_ref = {}
        params = [mc.G_ttvfast, mc.star_mass_val]

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
                mc.star_mass_val, dict_pams['M']/mc.M_SEratio, dict_pams['P'], dict_pams['i'], dict_pams['e']),
            'f': dict_pams['f'],
            'e': dict_pams['e'],
            'o': dict_pams['o']}


class SinusoidsCommonVariables:
    ''' This is among the first classes I've been writing in Python, and differently from other
    classes in this suite it has not updated (only bug-fixed) while my Python skilled were improving
    So you may find many sub-optimal approaches
    '''

    def __init__(self):

        self.n_datasets = 0

        self.offset_coherence = True  # Same phase offset across seasons
        self.offset_common_id = -1
        self.offset_reference_name = ''
        self.use_offset = {}

        self.phase_coherence = False
        self.phase_sincro = False

        self.n_pha = 0

        # No 'Season' selection in the configuration file, so we make believe the code that all the data
        # is contained within one season
        self.season_sel = False
        self.n_seasons = 1
        self.season_name = ['Season_0']
        self.season_list = [0., 5000000.0]
        self.season_range = np.asarray(self.season_list, dtype=np.double)

        # self.Prot_bounded = False
        # self.pha_bounded = False

        self.Prot_bounds = [0., 30.0]
        self.pof_bounds = [0., 1.0]
        self.pha_bounds = [0., 1.0]

        self.phase_list = []

        self.prior_kind = {}
        self.prior_pams = {}

    def add_season_range(self, range_input, phase_input):
        # Reset default values at first call
        if not self.season_sel:
            self.n_seasons = 0
            self.season_list = []
            self.season_name = []
            self.season_sel = True

        self.season_name.append('Season_' + repr(self.n_seasons))
        self.season_list.append(range_input)
        self.season_range = np.asarray(self.season_list, dtype=np.double)

        self.phase_list.append(phase_input)
        self.phase = np.asarray(self.phase_list, dtype=np.int64)

        self.n_pha = np.amax(self.phase, axis=1)  # maximum value for each period
        self.n_pha_max = np.amax(self.n_pha)  # maximum value for each period

        self.n_seasons += 1
        return

    # def add_phase_offset(self, dataset, season_name):
    #    # Check if the dataset needs an offset with respect to the
    #    # reference dataset (simply the first one)
    #    add_bound = False
    #    try:
    #        value = self.use_offset[dataset.kind][season_name + '_off']
    #    except KeyError:
    #        # Key is not present
    #        if self.offset_skip_first:
    #            self.use_offset[dataset.kind][season_name] = False
    #            self.offset_skip_first = False
    #        else:
    #            self.use_offset[dataset.kind][season_name] = True
    #            add_bound = True
    #    return add_bound

    def setup_model_sinusoids(self, dataset):
        dataset.n_amp = self.phase[:, dataset.ind]
        dataset.p_mask = np.zeros([self.n_seasons, dataset.n], dtype=bool)
        dataset.n_seasons = 0
        dataset.season_flag = np.zeros(self.n_seasons, dtype=bool)

        # If this is the first dataset, it is assumed as the reference one for the offset
        # otherwise, the offset is applied
        try:
            value = self.use_offset[dataset.kind]
        except KeyError:
            if self.offset_reference_name == '':
                self.offset_reference_name = dataset.kind
                self.use_offset[dataset.kind] = False
            else:
                self.use_offset[dataset.kind] = not self.phase_sincro  # True

        for ii in xrange(0, self.n_seasons):
            p_sel = (self.season_range[ii, 0] < dataset.x) & (dataset.x < self.season_range[ii, 1])
            if np.sum(p_sel) == 0:
                'No data points within the specified range'
            else:
                dataset.season_flag[ii] = True
                dataset.p_mask[ii, :] = p_sel[:]
                dataset.n_seasons += 1

        if dataset.kind == 'RV': dataset.sinamp_bounds = np.asarray([0., 100.], dtype=np.double)
        if dataset.kind == 'Phot': dataset.sinamp_bounds = np.asarray([0., 0.5], dtype=np.double)
        if dataset.kind == 'FWHM': dataset.sinamp_bounds = np.asarray([0., 2000.], dtype=np.double)
        if dataset.kind == 'BIS': dataset.sinamp_bounds = np.asarray([0., 60.], dtype=np.double)
        if dataset.kind == 'Act': dataset.sinamp_bounds = np.asarray([0., 10.], dtype=np.double)
        return

    def compute(self, mc, theta, dataset):
        # MC = Model_Container object
        # Prot and pha_in could be brought out from the llop, but I would not work
        # for plaet-only analysis
        Prot = theta[mc.variable_list['Common']['Prot']]
        pha_in = np.zeros([self.n_seasons, self.n_pha_max], dtype=np.double)
        off_in = np.zeros([self.n_seasons], dtype=np.double)
        amp_in = np.zeros([self.n_seasons, self.n_pha_max], dtype=np.double)
        for jj in range(0, self.n_seasons):
            pha_in[jj, :self.n_pha[jj]] = theta[mc.variable_list[self.season_name[jj] + '_pha']]
            if dataset.season_flag[jj]:
                if self.use_offset[dataset.kind]:
                    off_in[:] = theta[mc.variable_list[dataset.kind][self.season_name[jj] + '_off']]
                amp_in[jj, :dataset.n_amp[jj]] = \
                    theta[mc.variable_list[dataset.name_ref][self.season_name[jj] + '_amp']]
        return self.model_sinusoids(dataset, Prot, amp_in, pha_in, off_in)

    def print_vars(self, mc, theta):
        # Prot and pha_in could be brought out from the llop, but I would not work
        # for planet-only analysis
        mc.pam_names[mc.variable_list['Prot']] = 'Prot'
        print 'Prot ', theta[mc.variable_list['Prot']]

        for jj in range(0, self.n_seasons):
            id_var = mc.variable_list[self.season_name[jj] + '_pha']
            if np.size(id_var) == 0:
                continue
            if np.size(id_var) == 1:
                mc.pam_names[id_var] = self.season_name[jj] + '_pha'
            else:
                for ii in id_var:
                    mc.pam_names[ii] = self.season_name[jj] + '_' + repr(ii - id_var[0]) + '_pha'

            print self.season_name[jj], '_pha', theta[id_var]

        for dataset in mc.dataset_list:
            for jj in range(0, self.n_seasons):
                if dataset.season_flag[jj]:
                    if self.use_offset[dataset.kind]:
                        id_var = mc.variable_list[dataset.kind][self.season_name[jj] + '_off']
                        if np.size(id_var) == 0:
                            continue
                        if np.size(id_var) == 1:
                            mc.pam_names[id_var] = dataset.kind + '_' + self.season_name[jj] + '_off'
                        else:
                            for ii in id_var:
                                mc.pam_names[ii] = dataset.kind + '_' + self.season_name[jj] + \
                                                   '_' + repr(ii - id_var[0]) + '_off'

                        print dataset.name_ref, dataset.kind, self.season_name[jj], '_off', \
                            theta[mc.variable_list[dataset.kind][self.season_name[jj] + '_off']]

                    id_var = mc.variable_list[dataset.name_ref][self.season_name[jj] + '_amp']
                    if np.size(id_var) == 0:
                        continue
                    if np.size(id_var) == 1:
                        mc.pam_names[id_var] = dataset.name_ref[:-4] + '_' + self.season_name[jj] + '_amp'
                    else:
                        for ii in id_var:
                            mc.pam_names[ii] = dataset.name_ref[:-4] + '_' + self.season_name[jj] + \
                                               '_' + repr(ii - id_var[0]) + '_amp'

                    print dataset.name_ref, self.season_name[jj], '_amp', \
                        theta[mc.variable_list[dataset.name_ref][self.season_name[jj] + '_amp']]

    def define_bounds(self, mc):
        mc.bounds_list.append(self.Prot_bounds[:])
        mc.variable_list['Prot'] = mc.ndim
        mc.ndim += 1

        for jj in range(0, self.n_seasons):
            for kk in range(0, self.n_pha[jj]):
                mc.bounds_list.append(self.pha_bounds[:])
            mc.variable_list[self.season_name[jj] + '_pha'] = np.arange(mc.ndim,
                                                                        mc.ndim + self.n_pha[jj], 1)
            mc.ndim += self.n_pha[jj]

        for dataset in mc.dataset_list:

            # two nested case:
            # 1) dataset.kind has an offset or not (if it is the reference offset)
            # 2) the offset is the same for every season or not
            # WARNING: the offset is defined for each DATASET.KIND and not for each DATASET.NAME_REF
            # since the offset is a physical effect (not an instrumental one)

            for jj in range(0, self.n_seasons):

                if dataset.season_flag[jj]:

                    if self.use_offset[dataset.kind]:
                        if (not self.offset_coherence) or \
                                (self.offset_coherence and self.offset_common_id < 0):
                            mc.bounds_list.append(self.pof_bounds[:])
                            mc.variable_list[dataset.kind][self.season_name[jj] + '_off'] = mc.ndim
                            self.offset_common_id = mc.ndim
                            mc.ndim += 1
                        else:
                            mc.variable_list[dataset.kind][self.season_name[jj] + '_off'] = \
                                self.offset_common_id

                    for kk in xrange(0, dataset.n_amp[jj]):
                        mc.bounds_list.append(dataset.sinamp_bounds)

                    mc.variable_list[dataset.name_ref][self.season_name[jj] + '_amp'] = \
                        np.arange(mc.ndim, mc.ndim + dataset.n_amp[jj], 1)
                    mc.ndim += dataset.n_amp[jj]

    def starting_point(self, mc):
        pass
        ### MUST BE FIXED ###

    def model_sinusoids(self, dataset, p_rot, amp, pha, off):
        # np.size(SinAmp)==np.sum(n_amp) ???? What if an activity dataset is missing?
        # cycle for np.sum(rva_mask[:, jj]) <= 0 ??

        # Prot= Rotational periodo of the star
        # n_amp = number of sinusoids in the fit
        # amp = amplitudes
        # sin = phase of each sinusoid
        # off = overall offset of the sinusoids
        # sel = if the model has to be restricted to a specific temporal range

        model = np.zeros(dataset.n, dtype=np.double)
        xph = (dataset.x0 / p_rot) % 1
        har = np.arange(1, np.size(amp, axis=1) + 1, 1., dtype=np.double)

        for ii in xrange(0, self.n_seasons):
            for jj in xrange(0, np.size(pha[ii, :])):
                model += dataset.p_mask[ii, :] * amp[ii, jj] * np.sin(
                    (har[jj] * xph + pha[ii, jj] + off[ii]) * 2. * np.pi)
        return model


class GaussianProcessCommonVariables:
    def __init__(self):
        self.fix_list = {}
        self.bounds = {}
        self.starts = {}

        self.var_list = {}
        self.variables = {}

        self.prior_kind = {}
        self.prior_pams = {}

        self.fixed = []
        self.nfix = 0

        ''' Three parameters out of four are the same for all the datasets, since they are related to
        the properties of the physical process rather than the observed effects on a dataset
         From Grunblatt+2015, Affer+2016
         - theta: is usually related to the rotation period of the star( or one of its harmonics);
         - lambda: is the correlation decay timescale, and it can be related to the lifetime of the active regions.
         - omega: is the length scale of the periodic component, and can be linked to the size evolution of the active regions;
         - h: represents the amplitude of the correlations '''

        self.list_pams_common = {
            'Prot': 'U',
            'Pdec': 'U',
            'Oamp': 'LU'}

        self.list_pams_dataset = {'Hamp': 'U'}

        return

    def convert_val2gp(self, input_pam):
        """
        :param input_pam: dictonary with the 'physically meaningful' parameters of the GP kernel
        :return: dictonary with the parameters to be fed to 'george'
        WARNING: this subroutine is HIGHLY specific of your choice of the kernel! I reccomend to
        create a new Class with different transformations if you are planning of using a different
        kernel combination
        """
        return {
            'ln_P': np.log(input_pam['Prot']),
            'metric': input_pam['Pdec'] ** 2,
            'gamma': 2. / (input_pam['Oamp'] ** 2),
            'amp2': np.power(input_pam['Hamp'], 2)}

    def convert_gp2val(self, input_pam):
        """
        :param input_pam: dictonary with the parameters to be fed to 'george'
        :return: dictonary with the 'physically meaningful' parameters of the GP kernel
        WARNING: this subroutine is HIGHLY specific of your choice of the kernel! I reccomend to
        create a new Class with different transformations if you are planning of using a different
        kernel combination
        """
        return {
            'Prot': np.exp(input_pam['ln_P']),
            'Pdec': np.sqrt(input_pam['metric']),
            'Oamp': np.sqrt(2. / input_pam['gamma']),
            'Hamp': np.sqrt(input_pam['amp2'])}

    def add_dataset(self, name_ref='Common'):
        # planet_name : name of the dataset, otherwise 'common'

        self.bounds[name_ref] = {}
        self.variables[name_ref] = {}
        self.fix_list[name_ref] = {}
        self.var_list[name_ref] = {}

        self.prior_kind[name_ref] = {}
        self.prior_pams[name_ref] = {}

        if name_ref == 'Common':
            for name in self.list_pams_common:
                self.bounds['Common'][name] = [0.01, 100]
        else:
            for name in self.list_pams_dataset:
                self.bounds[name_ref][name] = [0.00001, 100]

    def define_bounds(self, mc):
        for var in self.list_pams_common:
            if var in self.fix_list['Common']:
                self.variables['Common'][var] = get_fix_val
                self.var_list['Common'][var] = self.nfix
                self.fixed.append(self.fix_list['Common'][var])
                self.nfix += 1
            else:
                if self.list_pams_common[var] == 'U':
                    self.variables['Common'][var] = get_var_val
                    mc.bounds_list.append(self.bounds['Common'][var])
                if self.list_pams_common[var] == 'LU':
                    self.variables['Common'][var] = get_var_exp
                    mc.bounds_list.append(np.log2(self.bounds['Common'][var]))

                self.var_list['Common'][var] = mc.ndim
                mc.variable_list['Common'][var] = mc.ndim
                mc.ndim += 1

        for dataset in mc.dataset_list:
            if 'gaussian' in dataset.models:
                for var in self.list_pams_dataset:
                    if var in self.fix_list[dataset.name_ref]:
                        self.variables[dataset.name_ref][var] = get_fix_val
                        self.var_list[dataset.name_ref][var] = self.nfix
                        self.fixed.append(self.fix_list[dataset.name_ref][var])
                        self.nfix += 1
                    else:
                        if self.list_pams_dataset[var] == 'U':
                            self.variables[dataset.name_ref][var] = get_var_val
                            mc.bounds_list.append(self.bounds[dataset.name_ref][var])
                        if self.list_pams_dataset[var] == 'LU':
                            self.variables[dataset.name_ref][var] = get_var_exp
                            mc.bounds_list.append(np.log2(self.bounds[dataset.name_ref][var]))
                        self.var_list[dataset.name_ref][var] = mc.ndim
                        mc.variable_list[dataset.name_ref][var] = mc.ndim
                        mc.ndim += 1

    def starting_point(self, mc):
        if 'Common' in self.starts:
            for var in self.starts['Common']:
                if self.list_pams_common[var] == 'U':
                    start_converted = self.starts['Common'][var]
                if self.list_pams_common[var] == 'LU':
                    start_converted = np.log2(self.starts['Common'][var])
                mc.starting_point[mc.variable_list['Common'][var]] = start_converted

        for dataset in mc.dataset_list:
            if 'gaussian' in dataset.models and dataset.name_ref in self.starts:
                    for var in self.starts[dataset.name_ref]:
                        if self.list_pams_common[var] == 'U':
                            start_converted = self.starts[dataset.name_ref][var]
                        if self.list_pams_common[var] == 'LU':
                            start_converted = np.log2(self.starts[dataset.name_ref][var])
                        mc.starting_point[mc.variable_list[dataset.name_ref][var]] = start_converted

    def convert(self, theta, d_name=None):
        dict_out = {}
        for key in self.list_pams_common:
            dict_out[key] = self.variables['Common'][key](theta, self.fixed, self.var_list['Common'][key])
        # If we need the parameters for the prior, we are not providing any name for the dataset
        if d_name is not None:
            for key in self.list_pams_dataset:
                dict_out[key] = self.variables[d_name][key](theta, self.fixed, self.var_list[d_name][key])
        return dict_out

    def return_priors(self, theta, d_name=None):
        prior_out = 0.00
        kep_pams = self.convert(theta, d_name)
        if d_name is None:
            for key in self.prior_pams['Common']:
                prior_out += giveback_priors(self.prior_kind['Common'][key], self.prior_pams['Common'][key],
                                             kep_pams[key])
        else:
            for key in self.prior_pams[d_name]:
                prior_out += giveback_priors(self.prior_kind[d_name][key], self.prior_pams[d_name][key], kep_pams[key])
        return prior_out

    def lnlk_compute(self, theta, dataset):
        # 2 steps:
        #   1) theta parameters must be converted in physical units (e.g. from logarithmic to linear space)
        #   2) physical values must be converted to {\tt george} input parameters
        gp_pams = self.convert_val2gp(self.convert(theta, dataset.name_ref))
        # gp_pams['ln_P] = ln_theta = ln_Period -> ExpSine2Kernel(gamma, ln_period)
        # gp_pams['metric'] = metric = r^2 = lambda**2  -> ExpSquaredKernel(metric=r^2)
        # gp_pams['gamma'] = Gamma =  1/ (2 omega**2) -> ExpSine2Kernel(gamma, ln_period)
        # gp_pams['amp2] = h^2 -> h^2 * ExpSquaredKernel * ExpSine2Kernel
        kernel = gp_pams['amp2'] * george.kernels.ExpSquaredKernel(metric=gp_pams['metric']) * \
                 george.kernels.ExpSine2Kernel(gamma=gp_pams['gamma'], ln_period=gp_pams['ln_P'])

        gp = george.GP(kernel)
        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        gp.compute(dataset.x0, env)

        return gp.lnlikelihood(dataset.y - dataset.model, quiet=True)

    def sample_compute(self, theta, dataset):

        gp_pams = self.convert_val2gp(self.convert(theta, dataset.name_ref))

        kernel = gp_pams['amp2'] * george.kernels.ExpSquaredKernel(metric=gp_pams['metric']) * \
                 george.kernels.ExpSine2Kernel(gamma=gp_pams['gamma'], ln_period=gp_pams['ln_P'])

        gp = george.GP(kernel)
        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        gp.compute(dataset.x0, env)

        return gp.sample_conditional(dataset.y - dataset.model, dataset.x)

    def print_vars(self, mc, theta):

        for name in self.list_pams_common:
            if name in mc.variable_list['Common']:
                mc.pam_names[mc.variable_list['Common'][name]] = name
                var = self.variables['Common'][name](theta, self.fixed, self.var_list['Common'][name])
                print 'GaussianProcess ', name, var, self.var_list['Common'][name], '(', theta[
                    self.var_list['Common'][name]], ')'

        for dataset in mc.dataset_list:
            for name in self.list_pams_dataset:
                if name in mc.variable_list[dataset.name_ref]:
                    mc.pam_names[mc.variable_list[dataset.name_ref][name]] = name
                    var = self.variables[dataset.name_ref][name](theta, self.fixed,
                                                                 self.var_list[dataset.name_ref][name])
                    print 'GaussianProcess ', dataset.name_ref, name, var, '(', theta[
                        self.var_list[dataset.name_ref][name]], ')'


class Dataset:
    def __init__(self, ind, kind, input_file, models):
        self.ind = ind
        self.kind = kind
        # 'RV', 'PHOT', 'ACT'...
        self.models = models
        self.name_ref = input_file

        self.list_sys = ['jitter', 'offset', 'linear']
        self.bounds = {}
        self.starts = {}
        self.n_sys = {}

        """There is a mix of arrays and dictonaries here because this was one of the first
            classes I created. I will probably consider switch everything to dictionary
            without affecting the functionality in a later time"""

        print 'Opening: ', input_file
        self.data = np.loadtxt(input_file)

        n_cols = np.size(self.data, axis=1)

        self.x = np.asarray(self.data[:, 0], dtype=np.double)
        self.y = np.asarray(self.data[:, 1], dtype=np.double)
        self.e = np.asarray(self.data[:, 2], dtype=np.double)

        self.n = np.size(self.x)

        self.sys = {}
        for var in self.list_sys:
            self.sys[var] = np.zeros(self.n, dtype=np.double) - 1
        # it was zero for jitter and offset, and -1 for linear for unknown reasons

        if n_cols > 3:
            self.sys['jitter'] = np.asarray(self.data[:, 3], dtype=np.double)
        if n_cols > 4:
            self.sys['offset'] = np.asarray(self.data[:, 4], dtype=np.double)
        if n_cols > 5:
            self.sys['linear'] = np.asarray(self.data[:, 5], dtype=np.double)

        # use different offsets for the data
        # off must start from zero
        # -1 values for self.j and self.l mean that these will not be used

        # self.a = np.asarray(self.data[:, 5], dtype=np.double)
        # Flag for activity fit: we can choose to fit the same RV amplitude
        # of the activity signal for all the datasets, use dfferent values
        # for each dataset or exclude some of the datasets

        # Model for RV systematics
        for var in self.list_sys:
            self.n_sys[var] = np.max(self.sys[var].astype(np.int64)) + 1

        print 'N = ', self.n
        for var in self.list_sys:
            print 'N '+var+' = ', self.n_sys[var]
        # print 'N activ. = ', self.n_a
        print

        self.Tref = np.mean(self.x, dtype=np.double)
        self.x0 = self.x - self.Tref

        self.mask = {}
        for var in self.list_sys:
            self.mask[var] = np.zeros([self.n, self.n_sys[var]], dtype=bool)
            for ii in xrange(0, self.n_sys[var]):
                self.mask[var][(abs(self.sys[var] - ii) < 0.1), ii] = True

        self.model = np.zeros(self.n, dtype=np.double)
        self.jitter = np.zeros(self.n, dtype=np.double)

        """Default boundaries are defined according to the characteristic of the dataset"""
        self.default_bounds = {'offset': [np.min(self.y), np.max(self.y)],
                               'jitter': [0., 50 * np.max(self.e)],
                               'linear': [-1., 1.]}

    def common_Tref(self, Tref_in):
        self.Tref = Tref_in
        self.x0 = self.x - self.Tref
        return

    def model_reset(self):
        self.model[:] = 0.0
        self.jitter[:] = 0.0
        return

    def model_offset(self, off_in):
        off = np.atleast_1d(off_in)
        for ii in xrange(0, self.n_sys['offset']):
            self.model[self.mask['offset'][:, ii]] += off[ii]

    def model_linear(self, m_in):
        m = np.atleast_1d(m_in)
        for ii in xrange(0, self.n_sys['linear']):
            self.model[self.mask['linear'][:, ii]] += m[ii] * self.x0[self.mask['linear'][:, ii]]

    def model_jitter(self, jit_in):
        jit = np.atleast_1d(jit_in)
        for ii in xrange(0, self.n_sys['jitter']):
            self.jitter[self.mask['jitter'][:, ii]] = jit[ii]

    def model_logchi2(self):
        env = 1.0 / (self.e ** 2.0 + self.jitter ** 2.0)
        return -0.5 * (np.sum((self.y - self.model) ** 2 * env - np.log(env)))

    def define_bounds(self, mc):

        for var in self.list_sys:
            if var in self.bounds:
                bounds = self.bounds[var]
            else:
                bounds = self.default_bounds[var]
            for jj in xrange(0, self.n_sys[var]):
                mc.bounds_list.append(bounds)  # bounds for jitter

        mc.variable_list[self.kind] = {}
        mc.variable_list[self.name_ref] = {}

        """adding the systematics variables to the list"""
        for var in self.list_sys:
            mc.variable_list[self.name_ref][var] = np.arange(mc.ndim, mc.ndim + self.n_sys[var], 1)
            mc.ndim += self.n_sys[var]

    def starting_point(self, mc):
        if self.name_ref in self.starts:
            for var in self.starts[self.name_ref]:
                mc.starting_point[mc.variable_list[self.name_ref][var]] = self.starts[self.name_ref][var]

    def print_vars(self, mc, theta):
        for param in self.list_sys:
            id_var = mc.variable_list[self.name_ref][param]
            if np.size(id_var) == 0:
                continue
            if np.size(id_var) == 1:
                mc.pam_names[id_var] = self.name_ref[:-4] + '_' + param
            else:
                for ii in id_var:
                    mc.pam_names[ii] = self.name_ref[:-4] + '_' + param + '_' + repr(ii - id_var[0])

            print self.name_ref, param, ' : ', theta[mc.variable_list[self.name_ref][param]]


class TransitCentralTimes(Dataset):
    def __init__(self, planet_name, input_file):

        self.kind = 'Tcent'
        self.models = ['Tcent']
        # 'RV', 'PHOT', 'ACT'...
        self.name_ref = input_file
        self.planet_name = planet_name

        self.list_sys = ['jitter', 'offset', 'linear']
        self.n_sys = {}
        self.mask = {}
        self.bounds = {}
        self.starts = {}
        self.deltaT = 1.10

        print 'Opening: ', input_file
        self.data = np.atleast_2d(np.loadtxt(input_file))

        self.n_transit = np.asarray(self.data[:, 0], dtype=np.int16)
        self.x = np.asarray(self.data[:, 1], dtype=np.double)
        self.e = np.asarray(self.data[:, 2], dtype=np.double)

        if np.size(self.data[0, :]) > 3:
            print 'TTV jitter found in the dataset'
            self.e = np.sqrt(self.data[:, 2] ** 2 + self.data[:, 3], dtype=np.double)

        self.n = np.size(self.x)
        self.model = np.zeros(self.n, dtype=np.double)

        self.Tref = np.mean(self.x, dtype=np.double)
        self.x0 = self.x - self.Tref

        print 'N = ', self.n
        print

        for var in self.list_sys:
            self.n_sys[var] = 0
            self.mask[var] = np.zeros([self.n, self.n_sys[var]], dtype=bool)

        self.model = np.zeros(self.n, dtype=np.double)
        self.jitter = np.zeros(self.n, dtype=np.double)

        """Default boundaries are defined according to the characteristic of the dataset"""
        self.default_bounds = {'offset': [np.min(self.x), np.max(self.x)],
                               'jitter': [0., 50 * np.max(self.e)],
                               'linear': [-1., 1.]}

    def compute(self, mc, theta):
        # By default, dataset.planet_name == planet_name
        dict_out = mc.pcv.convert(self.planet_name, theta)
        model = np.rint(self.x0 / dict_out['P'] - 1) * dict_out['P'] + \
                kp.kepler_Tcent_T0P(dict_out['P'], dict_out['f'], dict_out['e'], dict_out['o'])
        return model

    def model_logchi2(self):
        # boundaries in Tcent are specific of the dataset and not of a common
        # parameter for different dataset. The check can be internal
        # if np.sum(np.abs(self.x0 - self.model) < self.deltaT) < self.n:
        #    return -np.inf
        env = 1.0 / (self.e ** 2.0)
        return -0.5 * (np.sum((self.x0 - self.model) ** 2 * env - np.log(env)))

    def print_vars(self, mc, theta):
        # period, _, f, e, o = mc.pcv.convert(self.planet_name, theta)
        # model = np.rint(self.x0 / period) * period + kp.kepler_Tcent_T0P(period, f, e, o)

        if self.planet_name in mc.pcv.dynamical:
            dyn_output = mc.pcv.compute_dynamical(mc, theta)
            model = dyn_output[self.name_ref]
        else:
            model = self.compute(mc, theta)

        print 'Tc ', self.planet_name
        for ii in xrange(0, self.n):
            print 'Input Tc: ', self.x0[ii], '  Model Tc: ', model[ii], \
                '  Diff: ', model[ii] - self.x0[ii]


class ModelContainer:
    def __init__(self):
        self.dataset_list = []
        self.n_datasets = 0
        self.scv = SinusoidsCommonVariables()
        self.pcv = PlanetsCommonVariables()
        self.gcv = GaussianProcessCommonVariables()

        # pyde/emcee variables
        self.ngen = 0
        self.nsteps = 0
        self.nburn = 0
        self.npop_multi = 0
        self.nwalkers = 0
        self.thin = 1

        self.model_list = []
        self.bounds_list = []

        self.starting_point = None
        self.starting_point_flag = False
        self.recenter_bounds_flag = True

        self.planet_name = ''

        self.variable_list = {'Common': {}}
        self.bound_list = []
        self.bounds = 0
        self.ndim = 0
        self.pam_names = ''
        self.star_mass_val = 1.0000
        self.star_mass_err = 0.1000

        "Physical parameters are included here, so they are available to all the classes"

        self.G_grav = 6.67398e-11  # Gravitational constant is given in m^3 kg^-1 s^-2
        self.G_ttvfast = 0.000295994511  # Gravitational constant in unknown units from TTVfast
        self.M_sun = 1.98892e30
        self.M_jup = 1.89813e27
        self.M_SJratio = self.M_sun / self.M_jup
        self.M_JEratio = 317.83
        self.M_SEratio = self.M_SJratio * self.M_JEratio
        self.Mu_sun = 132712440018.9
        self.seconds_in_day = 86400
        self.AU_km = 1.4960 * 10 ** 8
        self.AUday2ms = self.AU_km / self.seconds_in_day * 1000.0

    def model_setup(self):
        self.n_datasets = np.size(self.dataset_list)
        for dataset in self.dataset_list:

            if 'sinusoids' in dataset.models:
                self.scv.setup_model_sinusoids(dataset)

            dataset.model_reset()
            for data_model in dataset.models:
                if not (data_model in self.model_list):
                    self.model_list.append(data_model)

    def create_bounds(self):
        # This routine creates the boundary array and at the same time
        # creates a dictionary with the name of the arrays and their
        # positions in bounds/theta array so that they can be accessed
        # without using nested counters

        self.ndim = 0

        for dataset in self.dataset_list:
            dataset.define_bounds(self)

        if 'kepler' in self.model_list:
            self.pcv.define_bounds(self)

        if 'sinusoids' in self.model_list:
            self.scv.define_bounds(self)

        if 'gaussian' in self.model_list:
            self.gcv.define_bounds(self)

        self.bounds = np.asarray(self.bounds_list)

    def create_starting_point(self):

        self.starting_point = np.average(self.bounds, axis=1)

        for dataset in self.dataset_list:
            dataset.starting_point(self)

        if 'kepler' in self.model_list:
            self.pcv.starting_point(self)

        if 'sinusoids' in self.model_list:
            self.scv.starting_point(self)

        if 'gaussian' in self.model_list:
            self.gcv.starting_point(self)

    def check_bounds(self, theta):
        for ii in xrange(0, self.ndim):
            if not (self.bounds[ii, 0] < theta[ii] < self.bounds[ii, 1]):
                return False
        for pl_name in self.pcv.planet_name:
            # if not bool(self.pcv.circular[pl_name]):
            e = self.pcv.variables[pl_name]['e'](theta, self.pcv.fixed, self.pcv.var_list[pl_name]['e'])
            if not self.pcv.bounds[pl_name]['e'][0] <= e < self.pcv.bounds[pl_name]['e'][1]:
                return False

        return True

    def __call__(self, theta):
        if not self.check_bounds(theta):
            return -np.inf
        chi2_out = 0.0

        if 'kepler' in self.model_list:
            for pl_name in self.pcv.planet_name:
                chi2_out += self.pcv.return_priors(pl_name, theta)

        if bool(self.pcv.dynamical):
            """ check if any keyword ahas get the output model from the dynamical tool
            we must do it here because all the planet are involved"""
            dyn_output = self.pcv.compute_dynamical(self, theta)

        if 'sinusoid' in self.model_list:
            chi2_out += giveback_priors(
                self.scv.prior_kind['Prot'], self.scv.prior_pams['Prot'], theta[self.variable_list['Common']['Prot']])

        if 'gaussian' in self.model_list:
            chi2_out += self.gcv.return_priors(theta)

        for dataset in self.dataset_list:
            dataset.model_reset()
            dataset.model_offset(theta[self.variable_list[dataset.name_ref]['offset']])
            dataset.model_jitter(theta[self.variable_list[dataset.name_ref]['jitter']])
            dataset.model_linear(theta[self.variable_list[dataset.name_ref]['linear']])

            if 'kepler' in dataset.models:
                if bool(self.pcv.dynamical):
                    """ we have dynamical computations, so we include them in the model"""
                    dataset.model += dyn_output[dataset.name_ref]
                for pl_name in self.pcv.planet_name:
                    """ we check if there is any planet which model has been obtained by assuming non-intercating
                    keplerians, and then we compute the expected RVs"""
                    pass
                    if pl_name not in self.pcv.dynamical:
                        dataset.model += self.pcv.compute(theta, dataset, pl_name)

            if 'sinusoids' in dataset.models:
                dataset.model += self.scv.compute(self, theta, dataset)

            if 'Tcent' in dataset.models:
                if dataset.planet_name in self.pcv.dynamical:
                    """ we have dynamical computations, so we include them in the model"""
                    dataset.model += dyn_output[dataset.name_ref]
                else:
                    dataset.model += dataset.compute(self, theta)

            # Gaussian Process check MUST be the last one or the program will fail
            if 'gaussian' in dataset.models:
                chi2_out += self.gcv.return_priors(theta, dataset.name_ref)
                chi2_out += self.gcv.lnlk_compute(theta, dataset)
            else:
                chi2_out += dataset.model_logchi2()

        return chi2_out

    def results_resumen(self, theta):
        # Function with two goals:
        # * Unfold and print out the output from theta
        # * give back a parameter name associated to each value in the result array

        self.pam_names = self.ndim * ['']

        for dataset in self.dataset_list:
            dataset.print_vars(self, theta)
            print

        if 'kepler' in self.model_list:
            self.pcv.print_vars(self, theta)
            print

        if 'sinusoids' in self.model_list:
            self.scv.print_vars(self, theta)
            print

        if 'gaussian' in self.model_list:
            self.gcv.print_vars(self, theta)
            print

    def rv_make_model(self, theta, x_range, x_phase):
        # it return the RV model for a single planet, after removing the activity from the RV curve and removing
        # the offsets between the datasets

        model_actv = {}
        model_plan = {}
        model_orbs = {}
        model_dsys = {}

        model_orbs['BJD'] = x_range * 0.0
        model_orbs['pha'] = x_phase * 0.0
        model_plan['BJD'] = {}
        model_plan['pha'] = {}

        if bool(self.pcv.dynamical):
            """Check if the dynamical option has been activated: full RV curve will be computed using
            the dynamical integrator"""
            dyn_output = self.pcv.compute_dynamical(self, theta)
            dyn_output_fullorbit = self.pcv.compute_dynamical(self, theta, full_orbit=(x_range-self.Tref))
            model_orbs['BJD'] += dyn_output_fullorbit['full_orbit']
            print 'WARNING: phase plot generated using non-interacting model!!!'
            print model_orbs['BJD'][:10]
        # computing the orbit for the full dataset
        for pl_name in self.pcv.planet_name:

            dyn_flag = (pl_name in self.pcv.dynamical)
            if dyn_flag:
                dict_pams = self.pcv.kepler_from_dynamical(self, theta, pl_name)
            else:
                dict_pams = self.pcv.convert(pl_name, theta)

            model_plan['BJD'][pl_name] = self.pcv.model_kepler(dict_pams, x_range - self.Tref)
            model_plan['pha'][pl_name] = self.pcv.model_kepler(dict_pams, x_phase * dict_pams['P'])

            model_orbs['pha'] += model_plan['pha'][pl_name]
            if not dyn_flag:
                model_orbs['BJD'] += model_plan['BJD'][pl_name]

        for dataset in self.dataset_list:

            model_actv[dataset.name_ref] = np.zeros(dataset.n)
            model_orbs[dataset.name_ref] = np.zeros(dataset.n)
            model_plan[dataset.name_ref] = {}

            dataset.model_reset()
            dataset.model_offset(theta[self.variable_list[dataset.name_ref]['offset']])
            dataset.model_linear(theta[self.variable_list[dataset.name_ref]['linear']])

            model_dsys[dataset.name_ref] = dataset.model

            if 'kepler' in dataset.models:
                if bool(self.pcv.dynamical):
                    """Dynamical models (with all the interacting planets included in the resulting RVs)"""
                    model_orbs[dataset.name_ref] += dyn_output[dataset.name_ref]

                for pl_name in self.pcv.planet_name:
                    dyn_flag = (pl_name in self.pcv.dynamical)
                    if dyn_flag:
                        dict_pams = self.pcv.kepler_from_dynamical(self, theta, pl_name)
                        model_plan[dataset.name_ref][pl_name] = self.pcv.model_kepler(dict_pams, dataset.x0)
                    else:
                        model_plan[dataset.name_ref][pl_name] = self.pcv.compute(theta, dataset, pl_name)
                        model_orbs[dataset.name_ref] += model_plan[dataset.name_ref][pl_name]

            if 'sinusoids' in dataset.models:
                model_actv[dataset.name_ref] += self.scv.compute(self, theta, dataset)

        return model_dsys, model_plan, model_orbs, model_actv

    # This function recenters the bounds limits for circular variables
    # Also, it extends the range of a variable if the output of PyDE is a fixed number
    def recenter_bounds(self, pop_mean, population):
        ind_list = []

        if 'kepler' in self.model_list:
            for pl_name in self.pcv.planet_name:

                if 'esino' in self.variable_list[pl_name]:
                    esino_list = self.variable_list[pl_name]['esino']
                    ecoso_list = self.variable_list[pl_name]['ecoso']
                    e_pops = population[:, esino_list] ** 2 + population[:, ecoso_list] ** 2
                    o_pops = np.arctan2(population[:, esino_list], population[:, ecoso_list], dtype=np.double)
                    # e_mean = (self.pcv.bounds[pl_name]['e'][0] + self.pcv.bounds[pl_name]['e'][1]) / 2.
                    for ii in xrange(0, self.nwalkers):
                        if not self.pcv.bounds[pl_name]['e'][0] + 0.02 <= e_pops[ii] < \
                                        self.pcv.bounds[pl_name]['e'][1] - 0.02:
                            e_random = np.random.uniform(self.pcv.bounds[pl_name]['e'][0],
                                                         self.pcv.bounds[pl_name]['e'][1])
                            population[ii, esino_list] = np.sqrt(e_random) * np.sin(o_pops[ii])
                            population[ii, ecoso_list] = np.sqrt(e_random) * np.cos(o_pops[ii])

                if 'f' in self.variable_list[pl_name]:
                    ind_list.append(self.variable_list[pl_name]['f'])

                if 'o' in self.variable_list[pl_name]:
                    ind_list.append(self.variable_list[pl_name]['o'])

                if 'lN' in self.variable_list[pl_name]:
                    ind_list.append(self.variable_list[pl_name]['lN'])

        if 'sinusoids' in self.model_list:
            for jj in range(0, self.scv.n_seasons):
                ind_list.extend(self.variable_list[self.scv.season_name[jj] + '_pha'])
            for dataset in self.dataset_list:
                for jj in range(0, self.scv.n_seasons):
                    if dataset.season_flag[jj]:
                        # ind_list.extend(self.variable_list[dataset.planet_name][self.scv.season_name[jj] + '_amp'])
                        if self.scv.use_offset[dataset.kind]:
                            ind_list.append(self.variable_list[dataset.kind][self.scv.season_name[jj] + '_off'])

        if np.size(ind_list) > 0:
            tmp_range = (self.bounds[:, 1] - self.bounds[:, 0]) / 2
            for var_ind in ind_list:
                self.bounds[var_ind, :] = pop_mean[var_ind] + [-tmp_range[var_ind], tmp_range[var_ind]]
                fix_sel = (population[:, var_ind] <= self.bounds[var_ind, 0]) | (
                    population[:, var_ind] >= self.bounds[var_ind, 1])
                population[fix_sel, var_ind] = pop_mean[var_ind]

        for ii in xrange(0, self.ndim):
            if np.amax(population[:, ii]) - np.amin(population[:, ii]) < 10e-7:
                range_restricted = (self.bounds[ii, 1] - self.bounds[ii, 0]) / 100.
                min_bound = np.maximum((pop_mean[ii] - range_restricted / 2.0), self.bounds[ii, 0])
                max_bound = np.minimum((pop_mean[ii] + range_restricted / 2.0), self.bounds[ii, 1])
                population[:, ii] = np.random.uniform(min_bound, max_bound, self.nwalkers)


class ModelContainerPolyChord(ModelContainer):
    def polychord_priors(self, cube):
        theta = (self.bounds[:, 1] - self.bounds[:, 0]) * cube + self.bounds[:, 0]
        print theta.tolist()
        return theta.tolist()

    def polychord_call(self, theta1):
        theta = np.empty(self.ndim)
        for i in xrange(0, self.ndim):
            theta[i] = theta1[i]
        print theta
        phi = [0.0] * self.nDerived
        chi_out = self(theta)
        if chi_out < -0.5e10:
            return -0.5e10, phi
        return chi_out, phi


class ModelContainerMultiNest(ModelContainer):
    def multinest_priors(self, cube, ndim, nparams):
        for i in xrange(0, ndim):
            cube[i] = (self.bounds[i, 1] - self.bounds[i, 0]) * cube[i] + self.bounds[i, 0]

    def multinest_call(self, theta1, ndim, nparams):
        # Workaround for variable selection: if a variable as null index
        # (i.e. it has not been included in the model)
        # the numpy array will give back an empty list, the ctype will give back an error
        theta = np.empty(ndim)
        for i in xrange(0, ndim):
            theta[i] = theta1[i]

        chi_out = self(theta)
        if chi_out < -0.5e10:
            return -0.5e10
        return chi_out


class ModelContainerSPOTPY(ModelContainer):
    def spotpy_priors(self, theta_med):
        pams_priors = []
        for ii in xrange(0, self.ndim):
            pams_priors.append(
                spotpy.parameter.Uniform(
                    'theta' + repr(ii), self.bounds[ii, 0], self.bounds[ii, 1], optguess=theta_med[ii]))

        if 'kepler' in self.model_list:
            for pl_name in self.pcv.planet_name:
                for key in self.pcv.prior_pams[pl_name]:
                    ind = self.variable_list[pl_name][key]
                    pams_priors[ind] = self.assign_priors(self.pcv.prior_kind[pl_name][key], 'theta' + repr(ind),
                                                          self.pcv.prior_pams[pl_name][key], theta_med[ind])

        if 'sinusoids' in self.model_list:
            if 'Prot' in self.scv.prior_pams:
                ind = self.variable_list['Common']['Prot']
                pams_priors[ind] = self.assign_priors(self.scv.prior_kind[key], 'theta' + repr(ind),
                                                      self.scv.prior_pams[key], theta_med[ind])

        if 'gaussian' in self.model_list:
            for key in self.gcv.prior_pams['Common']:
                ind = self.variable_list[key]
                pams_priors[ind] = self.assign_priors(self.gcv.prior_kind[key], 'theta' + repr(ind),
                                                      self.gcv.prior_pams['Common'][key], theta_med[ind])
            for dataset in self.dataset_list:
                for key in self.gcv.prior_pams[dataset.name_ref]:
                    ind = self.variable_list[dataset.name_ref][key]
                    pams_priors[ind] = self.assign_priors(self.pcv.prior_kind[dataset.name_ref][key],
                                                          'theta' + repr(ind),
                                                          self.pcv.prior_pams[dataset.name_ref][key], theta_med[ind])
        return pams_priors

    def assign_priors(self, key_name, prior_name, prior_val, prior_med):
        if key_name == 'Normal' or key_name == 'Gaussian':
            return spotpy.parameter.Normal(prior_name, prior_val[0], prior_val[1], optguess=prior_med)
        if key_name == 'logNormal':
            return spotpy.parameter.logNormal(prior_name, prior_val[0], prior_val[1], optguess=prior_med)
        if key_name == 'Exponential':
            return spotpy.parameter.Exponential(prior_name, prior_val[0], optguess=prior_med)
        return spotpy.parameter.Uniform(prior_name, prior_val[0], prior_val[1], optguess=prior_med)

    def output_concatenated(self):
        output_array = []
        for dataset in self.dataset_list:
            output_array = np.concatenate((output_array, dataset.y), axis=0)
        return output_array

    def __call__(self, theta):
        output_array = []
        for dataset in self.dataset_list:
            dataset.model_reset()
            dataset.model_offset(theta[self.variable_list[dataset.name_ref]['offset']])
            dataset.model_jitter(theta[self.variable_list[dataset.name_ref]['jitter']])
            dataset.model_linear(theta[self.variable_list[dataset.name_ref]['linear']])

            if 'kepler' in dataset.models:
                for pl_name in self.pcv.planet_name:
                    dataset.model += self.pcv.compute(theta, dataset, pl_name)

            if 'sinusoids' in dataset.models:
                dataset.model += self.scv.compute(self, theta, dataset)

            if 'Tcent' in dataset.models:
                dataset.model += dataset.compute(self, theta)

            # Gaussian Process check MUST be the last one or the program will fail
            if 'gaussian' in dataset.models:
                dataset.model += self.gcv.sample_compute(theta, dataset)
            output_array = np.concatenate((output_array, dataset.model), axis=0)

        return output_array


def yaml_parser(file_conf, mc):
    stream = file(file_conf, 'r')
    config_in = yaml.load(stream)

    conf = config_in['Inputs']
    for counter in conf:
        print conf[counter]['Kind'], conf[counter]['File'], conf[counter]['Models']
        mc.dataset_list.append(Dataset(counter, conf[counter]['Kind'], conf[counter]['File'], conf[counter]['Models']))

        if counter == 0:
            mc.Tref = mc.dataset_list[0].Tref
        else:
            mc.dataset_list[counter].common_Tref(mc.Tref)

        if 'Boundaries' in conf[counter]:
            bound_conf = conf[counter]['Boundaries']
            for var in bound_conf:
                mc.dataset_list[counter].bounds[var] = np.asarray(bound_conf[var], dtype=np.double)

        if 'Starts' in conf[counter]:
            mc.starting_point_flag = True
            starts_conf = conf[counter]['Starts']
            for var in bound_conf:
                mc.dataset_list[counter].starts[var] = np.asarray(starts_conf[var], dtype=np.double)

    mc.planet_name = config_in['Output']

    if 'Planets' in config_in:
        conf = config_in['Planets']
        for counter in conf:
            planet_name = 'Planet_' + repr(counter)
            planet_conf = conf[counter]
            mc.pcv.add_planet(planet_name)

            if 'Boundaries' in planet_conf:
                bound_conf = planet_conf['Boundaries']
                for var in bound_conf:
                    mc.pcv.bounds[planet_name][var] = np.asarray(bound_conf[var], dtype=np.double)

            if 'Fixed' in planet_conf:
                fixed_conf = planet_conf['Fixed']
                for var in fixed_conf:
                    mc.pcv.fix_list[planet_name][var] = np.asarray(fixed_conf[var], dtype=np.double)

            if 'Priors' in planet_conf:
                prior_conf = planet_conf['Priors']
                for var in prior_conf:
                    mc.pcv.prior_kind[planet_name][var] = prior_conf[var][0]
                    mc.pcv.prior_pams[planet_name][var] = np.asarray(prior_conf[var][1:], dtype=np.double)
                print mc.pcv.prior_kind
                print mc.pcv.prior_pams

            if 'Starts' in planet_conf:
                mc.starting_point_flag = True
                starts_conf = planet_conf['Starts']
                for var in starts_conf:
                    mc.pcv.starts[planet_name][var] = np.asarray(starts_conf[var], dtype=np.double)

            if 'Orbit' in planet_conf:
                # By default orbits are keplerians
                if planet_conf['Orbit'] == 'circular':
                    mc.pcv.switch_to_circular(planet_name)
                if planet_conf['Orbit'] == 'dynamical':
                    mc.pcv.switch_to_dynamical(planet_name)

            if 'Tcent' in planet_conf:
                mc.dataset_list.append(TransitCentralTimes(planet_name, planet_conf['Tcent']))
                mc.dataset_list[-1].common_Tref(mc.Tref)

            if 'Inclination' in planet_conf:
                mc.pcv.inclination[planet_name] = planet_conf['Inclination']

    if 'Sinusoids' in config_in:
        conf = config_in['Sinusoids']
        mc.scv.Prot_bounds = np.asarray(conf['Prot'], dtype=np.double)
        if 'Priors' in conf:
            mc.scv.prior_kind[var] = conf['Priors'][var][0]
            mc.scv.prior_pams[var] = np.asarray(conf['Priors'][var][1:], dtype=np.double)

        for counter in conf['Seasons']:
            mc.scv.add_season_range(np.asarray(conf['Seasons'][counter][:2], dtype=np.double),
                                    conf['Seasons'][counter][2:])

            # if 'Phase_dataset' in conf:
            #    # Additional activity indicators associated to RVs must the same order and number of sinusoids,
            #    #
            #    mc.scv.phase = np.asarray(conf['Phase_dataset'], dtype=np.int64)

    if 'Gaussian' in config_in:
        conf = config_in['Gaussian']
        for name_ref in conf:
            if name_ref == 'Common':
                dataset_name = 'Common'
            else:
                dataset_name = mc.dataset_list[name_ref].name_ref
            mc.gcv.add_dataset(dataset_name)

            if 'Boundaries' in conf[name_ref]:
                bound_conf = conf[name_ref]['Boundaries']
                for var in bound_conf:
                    mc.gcv.bounds[dataset_name][var] = np.asarray(bound_conf[var], dtype=np.double)

            if 'Fixed' in conf[name_ref]:
                fixed_conf = conf[name_ref]['Fixed']
                for var in fixed_conf:
                    mc.gcv.fix_list[dataset_name][var] = np.asarray(fixed_conf[var], dtype=np.double)

            if 'Priors' in conf[name_ref]:
                prior_conf = conf[name_ref]['Priors']
                for var in prior_conf:
                    mc.gcv.prior_kind[dataset_name][var] = prior_conf[var][0]
                    mc.gcv.prior_pams[dataset_name][var] = np.asarray(prior_conf[var][1:], dtype=np.double)

            if 'Starts' in conf[name_ref]:
                mc.starting_point_flag = True
                starts_conf = conf[name_ref]['Starts']
                for var in starts_conf:
                    mc.gcv.starts[dataset_name][var] = np.asarray(starts_conf[var], dtype=np.double)

    if 'Tref' in config_in:
        mc.Tref = np.asarray(config_in['Tref'])
        for dataset in mc.dataset_list:
            dataset.common_Tref(mc.Tref)

    if 'Ngen' in config_in['emcee']:
        mc.ngen = np.asarray(config_in['emcee']['Ngen'], dtype=np.double)

    if 'Nsteps' in config_in['emcee']:
        mc.nsteps = np.asarray(config_in['emcee']['Nsteps'], dtype=np.int64)

    if 'Nburn' in config_in['emcee']:
        mc.nburn = np.asarray(config_in['emcee']['Nburn'], dtype=np.int64)

    if 'Npop_mult' in config_in['emcee']:
        mc.npop_mult = np.asarray(config_in['emcee']['Npop_mult'], dtype=np.int64)

    if 'Thin' in config_in['emcee']:
        mc.thin = np.asarray(config_in['emcee']['Thin'], dtype=np.int64)

    if 'Recenter_Bounds' in config_in['emcee']:
        # required to avoid a small bug in the code
        # if the dispersion of PyDE walkers around the median value is too broad,
        # then emcee walkers will start outside the bounds, causing an error
        mc.recenter_bounds_flag = config_in['emcee']['Recenter_Bounds']

    if 'Star_Mass' in config_in:
        mc.star_mass_val = np.asarray(config_in['Star_Mass'][0], dtype=np.double)
        mc.star_mass_err = np.asarray(config_in['Star_Mass'][1], dtype=np.double)

    mc.model_setup()
