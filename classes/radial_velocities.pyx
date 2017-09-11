from common import *
from abstract_model import *

"""
New changes: 
    mc.variables  is now called  mc.transformation
    mv.var_list  is now called  mc.variable_index
    
    variable_index is the third argument of transformation
    it identifies which values from theta must be taken to convert the variable_sampler values to the physical parameter
    
    variable_sampler associate the value in theta to their label in the sampler space
"""


class RVkeplerian(AbstractModel):

    model_class = 'rv_keplerian'

    list_pams_common = {
        'P': 'LU',  # Period, log-uniform prior
        'K': 'LU',  # RV semi-amplitude, log-uniform prior
        'f': 'U',  # RV vurve phase, log-uniform prior
        'e': 'U',  # eccentricity, uniform prior - to be fixed
        'o': 'U'}  # argument of pericenter

    list_pams_dataset = {}

    recenter_pams_dataset = {}

    def compute(self, variable_value, dataset):
        #variable_value = self.convert(theta)
        return kp.kepler_RV_T0P(dataset.x0,
                                variable_value['f'],
                                variable_value['P'],
                                variable_value['K'],
                                variable_value['e'],
                                variable_value['o'])


class RVdynamical(AbstractModel):

    model_class = 'rv_dynamical'

    ''' Orbital parameters to be used in the dynamical fit '''
    list_pams_common = {
            'P': 'LU',  # Period in days
            'M': 'LU',  # Mass in Earth masses
            'i': 'U',  # inclination in degrees
            'f': 'U',  # phase - as defined by Malavolta+2016
            'lN': 'U',  # longitude of ascending node
            'e': 'U',  # eccentricity, uniform prior - to be fixed
            'o': 'U'}  # argument of pericenter

    list_pams_dataset = {}

    recenter_pams_dataset = {}


class TransitTimeKeplerian(AbstractModel):

    model_class = 'transit_time_keplerian'

    list_pams_common = {
            'P': 'LU',  # Period, log-uniform prior
            'f': 'U',  # RV vurve phase, log-uniform prior
            'e': 'U',  # eccentricity, uniform prior - to be fixed
            'o': 'U'}  # argument of pericenter

    list_pams_dataset = {}

    recenter_pams_dataset = {}

    def compute(self, variable_value, dataset):
        # By default, dataset.planet_name == planet_name
        #variable_value = self.convert(theta)
        return (np.floor(dataset.x0 / variable_value['P'])) * variable_value['P'] + dataset.Tref + \
                kp.kepler_Tcent_T0P(variable_value['P'], variable_value['f'], variable_value['e'], variable_value['o'])


class TransitTimeDynamical(AbstractModel):

    model_class = 'transit_time_dynamical'

    ''' Orbital parameters to be used in the dynamical fit '''
    list_pams_common = {
            'P': 'LU',  # Period in days
            'M': 'LU',  # Mass in Earth masses
            'i': 'U',  # inclination in degrees
            'f': 'U',  # phase - as defined by Malavolta+2016
            'lN': 'U',  # longitude of ascending node
            'e': 'U',  # eccentricity, uniform prior - to be fixed
            'o': 'U'}  # argument of pericenter

    list_pams_dataset = {}

    recenter_pams_dataset = {}


class DynamicalIntegrator:
    def __init__(self):
        self.model_name = 'dynamical_integrator'
        self.dynamical_integrator = 'TRADES'
        self.dynamical_set = {}

    def prepare(self, mc):
        """
        Prepare the input parameters according to chosen dynamical integrator
        :param mc:
               pcv:
        :return:
        """

        if self.dynamical_integrator == 'TRADES':
            self.prepare_trades(mc)
        if self.dynamical_integrator == 'ttvfast':
            self.prepare_ttvfast(mc)
        return

    def compute(self, mc,  *args, **kwargs):
        """
        Run the appropriate subroutine according to chosen dynamical integrator
        :param mc:
        :return:
        """
        if self.dynamical_integrator == 'TRADES':
            output = self.compute_trades(mc, *args, **kwargs)
        if self.dynamical_integrator == 'ttvfast':
            output = self.compute_ttvfast(mc, *args, **kwargs)
        return output

    def prepare_trades(self, mc):
        """
        :param mc:
        :return:
        """
        dataset_rv = 0
        int_buffer = dict(rv_times=[], t0_times=[], rv_ref=[], t0_ref=[], key_ref={})

        """ Putting all the RV epochs in the same array, flagging in the temporary buffer
            the stored values according to their dataset of origin
        """
        for dataset_name, dataset in mc.dataset_dict.iteritems():
            if dataset.dynamical is False: continue
            if dataset.kind == 'RV':
                int_buffer['rv_times'].extend(dataset.x.tolist())
                int_buffer['rv_ref'].extend(dataset.x * 0.0 + dataset_rv)
                int_buffer['key_ref'][dataset_name] = dataset_rv
                dataset_rv += 1
            elif dataset.kind == 'Tcent':
                int_buffer['t0_times'].extend(dataset.x.tolist())

        """ Creating the flag array after all the RV epochs have been mixed
        """
        self.dynamical_set['data'] = {'selection': {}}
        for dataset_name, dataset in mc.dataset_dict.iteritems():
            if dataset.dynamical is False: continue
            if dataset.kind == 'RV':
                self.dynamical_set['data']['selection'][dataset_name] = \
                    (np.asarray(int_buffer['rv_ref']) == int_buffer['key_ref'][dataset_name])

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
        n_body = 1  # The star is included among the bodies

        period_list = []
        planet_list = []

        for planet_name in mc.dynamical_dict:
            period_list.extend(
                [(mc.planet_commons[planet_name].bounds['P'][1] +
                  mc.planet_commons[planet_name].bounds['P'][0]) / 2.])
            planet_list.extend([planet_name])
        sort_planets = np.argsort(period_list)

        for planet_name in np.asarray(planet_list)[sort_planets]:
            tmp_t0_ntot = [0]
            tmp_t0_flag = [True]
            if planet_name in mc.dynamical_t0_dict:
                tmp_t0_ntot = [mc.dataset_dict[mc.dynamical_t0_dict[planet_name]].n]
                tmp_t0_flag = [True]
            t0_ntot.extend(tmp_t0_ntot)
            t0_flag.extend(tmp_t0_flag)
            self.dynamical_set['data']['plan_ref'][planet_name] = np.asarray(n_body).astype(int)
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

        for planet_name in mc.dynamical_dict:
            plan_n = self.dynamical_set['data']['plan_ref'][planet_name]
            if planet_name in mc.dynamical_t0_dict:
                self.dynamical_set['data']['t0_num'][0:t0_ntot[plan_n], plan_n] = \
                    mc.dataset_dict[mc.dynamical_t0_dict[planet_name]].n_transit[:].astype(int)
                self.dynamical_set['data']['t0_obs'][0:t0_ntot[plan_n], plan_n] = \
                    mc.dataset_dict[mc.dynamical_t0_dict[planet_name]].x
                self.dynamical_set['data']['t0_err'][0:t0_ntot[plan_n], plan_n] = \
                    mc.dataset_dict[mc.dynamical_t0_dict[planet_name]].e

        if self.dynamical_set['fake_t0s']:
            self.dynamical_set['data']['t0_num'][0:3, 1] = np.arange(0, 3, 1, dtype=int)
            self.dynamical_set['data']['t0_obs'][0:3, 1] = np.arange(-1, 2, 1, dtype=np.float64) * 10.0 + mc.Tref
            self.dynamical_set['data']['t0_err'][0:3, 1] = 0.1

        self.dynamical_set['trades']['n_body'] = n_body

        # print   '0  ---> ', self.dynamical_set['trades']['ti_beg']
        # print   '1  ---> ', self.dynamical_set['trades']['ti_ref']
        # print   '2  ---> ', self.dynamical_set['trades']['ti_int']
        # print   '3  ---> ', self.dynamical_set['trades']['n_body']
        # print   '4  ---> ', self.dynamical_set['data']['t0_tot']
        # print   '5  ---> ', self.dynamical_set['data']['t0_num']
        # print   '6  ---> ', self.dynamical_set['data']['t0_obs']
        # print   '7  ---> ', self.dynamical_set['data']['t0_err']

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

    def compute_trades(self, mc, theta, full_orbit=None):
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

        for planet_name in mc.dynamical_dict:
            n_plan = self.dynamical_set['data']['plan_ref'][planet_name]
            dict_pams = mc.common_models[planet_name].convert(theta)
            self.dynamical_set['pams']['M'][n_plan] = dict_pams['M'] / mc.M_SEratio
            self.dynamical_set['pams']['R'][n_plan] = dict_pams['R'] / mc.R_SEratio
            self.dynamical_set['pams']['P'][n_plan] = dict_pams['P']
            self.dynamical_set['pams']['e'][n_plan] = dict_pams['e']
            self.dynamical_set['pams']['o'][n_plan] = dict_pams['o'] * (180. / np.pi)
            self.dynamical_set['pams']['i'][n_plan] = dict_pams['i']
            self.dynamical_set['pams']['lN'][n_plan] = dict_pams['lN'] * (180. / np.pi)
            self.dynamical_set['pams']['mA'][n_plan] = (dict_pams['f'] - dict_pams['o']) * (180. / np.pi)

        # sample_plan[:, convert_out['Tcent']] = mc.Tref + kp.kepler_Tcent_T0P(
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
        #print 'T0_sim: ', t0_sim
        #print 'RV_sim: ', rv_sim
        #t0_sim -= mc.Tref
        for dataset_name, dataset in mc.dataset_dict.iteritems():
            if dataset.dynamical is False: continue
            if dataset.kind == 'RV' and full_orbit is None:
                output[dataset_name] = rv_sim[self.dynamical_set['data']['selection'][dataset_name]]
                # print 'RV out', output[dataset.name_ref]
            elif dataset.kind == 'Tcent' and dataset.planet_name in mc.dynamical_dict:
                n_plan = self.dynamical_set['data']['plan_ref'][dataset.planet_name]
                #print ' T0_sim selected: ', t0_sim[:self.dynamical_set['data']['t0_tot'][n_plan], n_plan]
                output[dataset_name] = t0_sim[:self.dynamical_set['data']['t0_tot'][n_plan], n_plan]
        if full_orbit is not None:
            output['full_orbit'] = rv_sim

        return output

    def prepare_ttvfast(self, mc):

        self.dynamical_set['rv_times'] = []
        self.dynamical_set['data_selection'] = {}

        dataset_rv = 0
        int_buffer = dict(rv_times=[], t0_times=[], rv_ref=[], t0_ref=[], key_ref={})

        for dataset_name, dataset in mc.dataset_dict.iteritems():
            if dataset.dynamical is False: continue
            if dataset.kind == 'RV':
                int_buffer['rv_times'].extend(dataset.x0.tolist())
                int_buffer['rv_ref'].extend(dataset.x0 * 0 + dataset_rv)
                int_buffer['key_ref'][dataset_name] = dataset_rv
                dataset_rv += 1
            elif dataset.kind == 'Tcent':
                int_buffer['t0_times'].extend(dataset.x0.tolist())

        if np.size(int_buffer['t0_times']) == 0:
            int_buffer['t0_times'] = int_buffer['rv_times']

        for dataset_name, dataset in mc.dataset_dict.iteritems():
            if dataset.dynamical is False: continue
            if dataset.kind == 'RV':
                self.dynamical_set['data_selection'][dataset_name] = \
                    (np.asarray(int_buffer['rv_ref']) == int_buffer['key_ref'][dataset_name])

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

        if self.dynamical_set['ttvfast']['t_beg'] > -10.:
            """It means that both rv_minmax[0] and t0_minmax[0] are greater than zero, i.e. that both
            RV and TTV epochs start after Tref """
            self.dynamical_set['ttvfast']['t_beg'] = 0.0000

            ####printself.dynamical_set['ttvfast']['t_beg']

    def compute_ttvfast(self, mc, theta, full_orbit=None):
        """ This function compute the expected TTV and RVs for dynamically interacting planets.
            The user can specify which planets are subject to interactions, e.g. long-period planets can be approximated
            with a Keplerian function"""

        input_flag = 0
        n_plan = 0

        P_min = None

        plan_ref = {}
        params = [mc.G_ttvfast, mc.star_mass[0]]

        for planet_name in mc.dynamical_dict:
            plan_ref[planet_name] = n_plan

            dict_pams = mc.common_models[planet_name].convert(theta)

            mA = (dict_pams['f'] - dict_pams['o']) * (180. / np.pi) + \
                 self.dynamical_set['ttvfast']['t_beg'] / dict_pams['P'] * 360.0000000000

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

        t_step = P_min / 20.
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
        # print 'RV  ------->  ', rv[:10]
        # print self.dynamical_set['len_rv'], rv_meas[:10]
        # print 't_beg', self.dynamical_set['ttvfast']['t_beg']
        # print 't_end', self.dynamical_set['ttvfast']['t_end']
        # print 'Full orbit flag: ', full_orbit
        # print 'Full orbit flag: ', full_orbit
        # print positions[:10,:]

        for dataset_name, dataset in mc.dataset_dict.iteritems():
            if dataset.dynamical is False: continue
            if dataset.kind == 'RV' and full_orbit is None:
                output[dataset_name] = rv_meas[self.dynamical_set['data_selection'][dataset_name]]

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
                # print ref_idf,  np.sum(t0_sel), len(t0_mod), t0_mod

                if np.sum(t0_sel) > np.max(dataset.n_transit) + ref_idf:
                    output[dataset_name] = t0_mod[dataset.n_transit + ref_idf]
                    """With this approach, missing T0s in the input dataset are automatically skipped """
                else:
                    """The output vector contains less T0s than supposed: collision between planets?? """
                    output[dataset_name] = dataset.n_transit * 0.0000

        if full_orbit is not None:
            output['full_orbit'] = rv_meas
        # print rv_meas[:3], self.dynamical_set['rv_times'][:3]
        return output