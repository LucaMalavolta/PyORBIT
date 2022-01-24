from __future__ import print_function
from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.subroutines.results_analysis import get_stellar_parameters

try:
    import ttvfast
except ImportError:
    pass

try:
    from pytrades_lib import pytrades
except ImportError:
    pass

"""
New changes:
    mc.variables  is now called  mc.transformation
    mv.var_list  is now called  mc.variable_index

    variable_index is the third argument of transformation
    it identifies which values from theta must be taken to convert the variable_sampler values to the physical parameter

    variable_sampler associate the value in theta to their label in the sampler spaces
"""


class RVkeplerian(AbstractModel):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'rv_keplerian'

        self.list_pams_common = {
            'P',  # Period
            'K',  # RV semi-amplitude
            'e',  # eccentricity, uniform prior - to be fixed
            'omega'}  # argument of pericenter

        self.use_time_of_transit = False
        self.use_mass_for_planets = False

    def initialize_model(self, mc, **kwargs):

        if mc.common_models[self.planet_ref].use_time_of_transit:
            self.list_pams_common.update(['Tc'])
            self.use_time_of_transit = True
            # Copying the property to the class for faster access
        else:
            self.list_pams_common.update(['mean_long'])

        if mc.common_models[self.planet_ref].use_mass_for_planets:
            self.list_pams_common.update(['M_Me'])
            self.list_pams_common.update(['mass'])
            self.use_mass_for_planets = True
        else:
            self.list_pams_common.update(['K'])

    def compute(self, variable_value, dataset, x0_input=None):

        if self.use_time_of_transit:
            mean_long = kepler_exo.kepler_Tc2phase_Tref(variable_value['P'],
                                                variable_value['Tc'] - dataset.Tref,
                                                variable_value['e'],
                                                variable_value['omega'])
        else:
            mean_long = variable_value['mean_long']

        if self.use_mass_for_planets:

            K = kepler_exo.kepler_K1(variable_value['mass'],
                                     variable_value['M_Me'] / constants.Msear, variable_value['P'], variable_value['i'],
                                     variable_value['e'])
        else:
            K = variable_value['K']

        if x0_input is None:
            return kepler_exo.kepler_RV_T0P(dataset.x0,
                                            mean_long,
                                            variable_value['P'],
                                            K,
                                            variable_value['e'],
                                            variable_value['omega'])
        else:
            return kepler_exo.kepler_RV_T0P(x0_input,
                                            mean_long,
                                            variable_value['P'],
                                            K,
                                            variable_value['e'],
                                            variable_value['omega'])


class RVdynamical(AbstractModel):
    model_class = 'rv_dynamical'

    def __init__(self, *args, **kwargs):
        super(RVdynamical, self).__init__(*args, **kwargs)

        ''' Orbital parameters to be used in the dynamical fit '''
        self.list_pams_common = {
            'P',  # Period in days
            'M_Me',  # Mass in Earth masses
            'Omega',  # longitude of ascending node
            'e',  # eccentricity, uniform prior - to be fixed
            'omega',  # argument of pericenter
            'mass'} #mass of the star (needed for proper dynamical computation and for reversibility)

        self.list_pams_dataset = set()

    def initialize_model(self, mc, **kwargs):

        if mc.common_models[self.planet_ref].use_inclination:
            """ i is the orbital inclination (in degrees) """
            self.list_pams_common.update(['i'])
        else:
            """ b is the impact parameter """
            self.list_pams_common.update(['b'])

            if mc.common_models[self.planet_ref].use_semimajor_axis:
                """ a is the semi-major axis (in units of stellar radii) """
                self.list_pams_common.update(['a_Rs'])
            else:
                """ rho is the density of the star (in solar units) """
                self.list_pams_common.update(['density'])

        if mc.common_models[self.planet_ref].use_time_of_transit:
            self.list_pams_common.update(['Tc'])
        else:
            self.list_pams_common.update(['mean_long'])

    #def compute(self, variable_value, dataset, x0_input=None):
    #    return dataset.external_model

class TransitTimeKeplerian(AbstractModel):
    model_class = 'transit_time_keplerian'

    def __init__(self, *args, **kwargs):
        super(TransitTimeKeplerian, self).__init__(*args, **kwargs)
        self.use_time_of_transit = False

        self.list_pams_common = {'P'}  # Period

        self.list_pams_dataset = set()

    def initialize_model(self, mc, **kwargs):

        if mc.common_models[self.planet_ref].use_time_of_transit:
            self.list_pams_common.update(['Tc'])
            self.use_time_of_transit = True
            # Copying the property to the class for faster access
        else:
            self.list_pams_common.update(['mean_long'])
            self.list_pams_common.update(['e'])
            self.list_pams_common.update(['omega'])
            # mean longitude = argument of pericenter + mean anomaly at Tref

    def compute(self, variable_value, dataset, x0_input=None):

        if self.use_time_of_transit:
            delta_T = variable_value['Tc'] - \
                np.floor((variable_value['Tc'] - dataset.Tref) / variable_value['P']) * variable_value['P']
        else:
            delta_T = dataset.Tref + \
                      kepler_exo.kepler_phase2Tc_Tref(variable_value['P'],
                                                      variable_value['mean_long'],
                                                      variable_value['e'],
                                                      variable_value['omega'])

        if x0_input is None:
            return np.floor(dataset.x0 / variable_value['P']) * variable_value['P'] + delta_T
        else:
            return np.floor(x0_input / variable_value['P']) * variable_value['P'] + delta_T


class TransitTimeDynamical(AbstractModel):
    model_class = 'transit_time_dynamical'

    def __init__(self, *args, **kwargs):
        super(TransitTimeDynamical, self).__init__(*args, **kwargs)

        ''' Orbital parameters to be used in the dynamical fit '''
        self.list_pams_common = {
            'P',    # Period in days
            'M_Me',    # Mass in Earth masses
            'Omega',   # longitude of ascending node
            'e',    # eccentricity, uniform prior - to be fixed
            'R_Rs',    # planet radius (in units of stellar radii)
            'omega',    # argument of pericenter
            'mass'} # mass of the star (needed for proper dynamical computation and for reversibility)

        self.list_pams_dataset = set()

        self.use_semimajor_axis = False
        self.use_inclination = False
        self.use_time_of_transit = False

    def initialize_model(self, mc, **kwargs):

        if mc.common_models[self.planet_ref].use_inclination:
            """ i is the orbital inclination (in degrees) """
            self.list_pams_common.update(['i'])
            self.use_inclination = True
        else:
            """ b is the impact parameter """
            self.list_pams_common.update(['b'])

            if mc.common_models[self.planet_ref].use_semimajor_axis:
                """ a is the semi-major axis (in units of stellar radii) """
                self.list_pams_common.update(['a_Rs'])
                self.use_semimajor_axis = True
            else:
                """ rho is the density of the star (in solar units) """
                self.list_pams_common.update(['density'])

        if mc.common_models[self.planet_ref].use_time_of_transit:
            self.list_pams_common.update(['Tc'])
            self.use_time_of_transit = True
            # Copying the property to the class for faster access
        else:
            self.list_pams_common.update(['mean_long'])
            # mean longitude = argument of pericenter + mean anomaly at Tref


class DynamicalIntegrator:
    def __init__(self):
        self.model_name = 'dynamical_integrator'
        self.dynamical_integrator = 'TRADES'
        self.dynamical_set = {}
        self.n_max_t0 = 0
        self.to_be_initialized = True

    def prepare(self, mc):
        """
        Prepare the input parameters according to chosen dynamical integrator
        :param mc:
               pcv:
        :return:
        """

        if self.dynamical_integrator == 'TRADES':
            try:
                from pytrades_lib import pytrades
            except ImportError:
                print("ERROR: TRADES not installed, this will not work")
                quit()

            self.prepare_trades(mc)

        if self.dynamical_integrator == 'ttvfast':
            try:
                import ttvfast
            except ImportError:
                print("ERROR: ttvfast not installed, this will not work")
                quit()

            self.prepare_ttvfast(mc)
        return

    def compute(self, mc, *args, **kwargs):
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
        for dataset_name, dataset in mc.dataset_dict.items():
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
        for dataset_name, dataset in mc.dataset_dict.items():
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
            period_list.extend([mc.common_models[planet_name].period_average])
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
        self.n_max_t0 = np.max(t0_ntot)

        self.dynamical_set['data']['t0_tot'] = np.asarray(t0_ntot).astype(int)
        self.dynamical_set['data']['t0_flg'] = np.asarray(t0_flag)

        self.dynamical_set['data']['t0_num'] = np.zeros([self.n_max_t0, n_body]).astype(int)
        self.dynamical_set['data']['t0_obs'] = np.zeros([self.n_max_t0, n_body], dtype=np.float64)
        self.dynamical_set['data']['t0_err'] = np.zeros([self.n_max_t0, n_body], dtype=np.float64)

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

        """ When the object is copied, data on RV and other properties are somehow lost
            but the object is still initialized - somehow """

        self.to_be_initialized = False
        print('TRADES parameters', self.dynamical_set['trades'])
        print()

        return

    def compute_trades(self, mc, theta, x_input=None):
        """ This function compute the expected TTV and RVs for dynamically interacting planets.
            The user can specify which planets are subject to interactions, e.g. long-period planets can be approximated
            with a Keplerian function"""

        """ setting up the dictionaries with the orbital parameters required by TRADES,
        """

        if self.to_be_initialized:
            self.prepare_trades(mc)

        self.dynamical_set['pams'] = {
            'M': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64),
            'R': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64),
            'P': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64),
            'e': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64),
            'omega': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64),
            'i': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64),
            'Omega': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64),
            'mA': np.zeros(self.dynamical_set['trades']['n_body'], dtype=np.float64)
        }

        """ Adding star parameters"""
        star_pams = get_stellar_parameters(mc, theta, warnings=False)
        #star_pams = mc.common_models['star_parameters'].convert(theta)

        self.dynamical_set['pams']['mass'][0] = star_pams['mass']
        self.dynamical_set['pams']['radius'][0] = star_pams['radius']

        for planet_name in mc.dynamical_dict:
            n_plan = self.dynamical_set['data']['plan_ref'][planet_name]

            dict_pams = mc.common_models[planet_name].convert(theta)

            if mc.common_models[planet_name].use_inclination:
                self.dynamical_set['pams']['i'][n_plan] = dict_pams['i']
            else:
                if mc.common_models[planet_name].use_semimajor_axis:
                    self.dynamical_set['pams']['i'][n_plan] = \
                        convert_b_to_i(dict_pams['b'],
                                       dict_pams['e'],
                                       dict_pams['omega'],
                                       dict_pams['a_Rs'])
                else:
                    a_temp = convert_rho_to_a(dict_pams['P'], star_pams['density'])
                    self.dynamical_set['pams']['i'][n_plan] = \
                        convert_b_to_i(dict_pams['b'],
                                       dict_pams['e'],
                                       dict_pams['omega'],
                                       a_temp)

            if mc.common_models[planet_name].use_time_of_transit:
                dict_pams['mean_long'] = kepler_exo.kepler_Tc2phase_Tref(dict_pams['P'],
                                                                 dict_pams['Tc'] - mc.Tref,
                                                                 dict_pams['e'],
                                                                 dict_pams['omega'])

            if 'R' in dict_pams:
                """ Converting the radius from Stellar units to Solar units"""
                self.dynamical_set['pams']['R'][n_plan] = dict_pams['R_Rs'] * star_pams['radius']
            else:
                """ Default value: slightly more than 1 Earth radii in Solar units"""
                self.dynamical_set['pams']['R'][n_plan] = 0.02

            self.dynamical_set['pams']['M'][n_plan] = dict_pams['M_Me'] / constants.Msear
            self.dynamical_set['pams']['P'][n_plan] = dict_pams['P']
            self.dynamical_set['pams']['e'][n_plan] = dict_pams['e']
            self.dynamical_set['pams']['omega'][n_plan] = dict_pams['omega']
            self.dynamical_set['pams']['Omega'][n_plan] = dict_pams['Omega']
            self.dynamical_set['pams']['mA'][n_plan] = (dict_pams['mean_long'] - dict_pams['omega'])

        # sample_plan[:, convert_out['Tcent']] = mc.Tref + kepler_exo.kepler_phase2Tc_Tref(
        #    sample_plan[:, convert_out['P']], sample_plan[:, convert_out['mL']],
        #    sample_plan[:, convert_out['e']], sample_plan[:, convert_out['o']])

        """ Extracted from TRADES:
          !!! SUBROUTINE TO RUN TRADES INTEGRATION AND RETURN RV_SIM AND T0_SIM
        !   subroutine kelements_to_data(t_start,t_epoch,step_in,t_int,&
        !     &m_msun,R_rsun,P_day,ecc,argp_deg,mA_deg,inc_deg,lN_deg,&
        !     &t_rv,transit_flag,n_t0,t0_num,& ! input
        !     &rv_sim,t0_sim,& ! output
        !     &n_body,n_rv,n_max_t0) ! dimensions
          subroutine kelements_to_data(t_start,t_epoch,step_in,t_int,&
            &m_msun,R_rsun,P_day,ecc,argp_deg,mA_deg,inc_deg,lN_deg,&
            &t_rv,transit_flag,& ! input
            &rv_sim,t0_sim,& ! output
            &n_body,n_rv,n_max_t0) ! dimensions

            ! INPUT
            ! t_start      == start of the integration
            ! t_epoch      == reference time epoch
            ! step_in      == initial step size of the integration
            ! t_int        == total integration time in days

            ! m_msun       == masses of all the bodies in Msun m_sun(n_body)
            ! R_rsun       == radii of all the bodies in Rsun r_rsun(n_body)
            ! P_day        == periods of all the bodies in days p_day(n_body); p_day(0) = 0
            ! ecc          == eccentricities of all the bodies ecc(n_body); ecc(0) = 0
            ! argp_deg     == argument of pericentre of all the bodies argp_deg(n_body); argp_deg(0) = 0
            ! mA_deg       == mean anomaly of all the bodies mA_deg(n_body); mA_deg(0) = 0
            ! inc_deg      == inclination of all the bodies inc_deg(n_body); inc_deg(0) = 0
            ! lN_deg       == longitude of node of all the bodies lN_deg(n_body); lN_deg(0) = 0

            ! t_rv         == time of the RV datapoints t_rv(n_rv)
            ! transit_flag == logical/boolean vector with which bodies should transit (.true.) or not (.false) transit_flag(n_body); transit_flag(0) = False

            ! OUTPUT
            ! rv_sim       == rv simulated in m/s, same dimension of t_rv
            ! t0_sim       == t0 simulated in days, same dimension of t0_num

            ! DIMENSIONS
            ! n_body       == number of bodies (take into account the star)
            ! n_rv         == number of radial velocities datapoints
            ! n_max_t0     == maxval(n_t0) == maxval of transits available

        """

        if x_input is None:
            rv_sim, t0_sim = pytrades.kelements_to_data(
                self.dynamical_set['trades']['ti_beg'],
                self.dynamical_set['trades']['ti_ref'],
                self.dynamical_set['trades']['i_step'],
                self.dynamical_set['trades']['ti_int'],
                self.dynamical_set['pams']['M'],
                self.dynamical_set['pams']['R'],
                self.dynamical_set['pams']['P'],
                self.dynamical_set['pams']['e'],
                self.dynamical_set['pams']['omega'],
                self.dynamical_set['pams']['mA'],
                self.dynamical_set['pams']['i'],
                self.dynamical_set['pams']['Omega'],
                self.dynamical_set['data']['rv_times'],
                self.dynamical_set['data']['t0_flg'], self.n_max_t0)
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
                self.dynamical_set['pams']['omega'],
                self.dynamical_set['pams']['mA'],
                self.dynamical_set['pams']['i'],
                self.dynamical_set['pams']['Omega'],
                x_input,
                self.dynamical_set['data']['t0_flg'], self.n_max_t0)

        output = {}

        for dataset_name, dataset in mc.dataset_dict.items():
            if dataset.dynamical is False: continue
            if dataset.kind == 'RV':
                if x_input is None:
                    output[dataset_name] = rv_sim[self.dynamical_set['data']['selection'][dataset_name]]
                else:
                    output[dataset_name] = rv_sim
            elif dataset.kind == 'Tcent' and dataset.planet_name in mc.dynamical_dict:
                n_plan = self.dynamical_set['data']['plan_ref'][dataset.planet_name]
                output[dataset_name] = t0_sim[:self.dynamical_set['data']['t0_tot'][n_plan], n_plan]

        return output

    def prepare_ttvfast(self, mc):

        self.dynamical_set['rv_times'] = []
        self.dynamical_set['data_selection'] = {}

        dataset_rv = 0
        int_buffer = dict(rv_times=[], t0_times=[], rv_ref=[], t0_ref=[], key_ref={})

        delta_t = 1000.0

        for dataset_name, dataset in mc.dataset_dict.items():
            if dataset.dynamical is False: continue
            if dataset.kind == 'RV':
                int_buffer['rv_times'].extend(dataset.x0.tolist())
                int_buffer['rv_ref'].extend(dataset.x0 * 0 + dataset_rv)
                int_buffer['key_ref'][dataset_name] = dataset_rv
                dataset_rv += 1
            elif dataset.kind == 'Tcent':
                int_buffer['t0_times'].extend(dataset.x0.tolist())

            delta_t = min(delta_t, np.amin(np.abs(dataset.x0[1:]-dataset.x0[:-1])))

        if np.size(int_buffer['t0_times']) == 0 and np.size(int_buffer['rv_times']) == 0:
            raise ValueError("Error with TTVFAST input: either RVs or Tc arrays must be non-empty")
        elif np.size(int_buffer['t0_times']) == 0:
            int_buffer['t0_times'] = int_buffer['rv_times'].copy()
        elif np.size(int_buffer['rv_times']) == 0:
            int_buffer['rv_times'] = int_buffer['t0_times'].copy()

        for dataset_name, dataset in mc.dataset_dict.items():
            if dataset.dynamical is False: continue
            if dataset.kind == 'RV':
                self.dynamical_set['data_selection'][dataset_name] = \
                    (np.asarray(int_buffer['rv_ref']) == int_buffer['key_ref'][dataset_name])

        self.dynamical_set['rv_times'] = int_buffer['rv_times']
        self.dynamical_set['len_rv'] = np.size(self.dynamical_set['rv_times'])

        self.dynamical_set['min_deltat'] = delta_t

        rv_minmax = [np.amin(int_buffer['rv_times']), np.amax(int_buffer['rv_times'])]
        t0_minmax = [np.amin(int_buffer['t0_times']), np.amax(int_buffer['t0_times'])]

        self.dynamical_set['ttvfast'] = {
            't_beg': np.min([rv_minmax[0], t0_minmax[0]]) - 10,
            't_end': np.max([rv_minmax[1], t0_minmax[1]]) + 10}

        if self.dynamical_set['ttvfast']['t_beg'] > -10.:
            """It means that both rv_minmax[0] and t0_minmax[0] are greater than zero, i.e. that both
            RV and TTV epochs start after Tref """
            self.dynamical_set['ttvfast']['t_beg'] = 0.0000

        self.to_be_initialized = False

    def compute_ttvfast(self, mc, theta, x_input=None):
        """ This function compute the expected TTV and RVs for dynamically interacting planets.
            The user can specify which planets are subject to interactions, e.g. long-period planets can be approximated
            with a Keplerian function"""

        if self.to_be_initialized:
            self.prepare_ttvfast(mc)

        input_flag = 0
        n_plan = 0

        P_min = None

        plan_ref = {}

        star_pams = get_stellar_parameters(mc, theta, warnings=False)
        #star_pams = mc.common_models['star_parameters'].convert(theta)

        # Gravitational constant in G [AU^3/Msun/d^2], stellar mass in Solar units
        params = [constants.Giau, star_pams['mass']]

        for planet_name in mc.dynamical_dict:

            plan_ref[planet_name] = n_plan

            dict_pams = mc.common_models[planet_name].convert(theta)

            if mc.common_models[planet_name].use_inclination:
                i_temp = dict_pams['i']
            else:
                if mc.common_models[planet_name].use_semimajor_axis:
                    i_temp = \
                        convert_b_to_i(dict_pams['b'],
                                       dict_pams['e'],
                                       dict_pams['omega'],
                                       dict_pams['a_Rs'])
                else:
                    a_temp = convert_rho_to_a(dict_pams['P'], star_pams['density'])
                    i_temp = \
                        convert_b_to_i(dict_pams['b'],
                                       dict_pams['e'],
                                       dict_pams['omega'],
                                       a_temp)

            if mc.common_models[planet_name].use_time_of_transit:
                dict_pams['mean_long'] = kepler_exo.kepler_Tc2phase_Tref(dict_pams['P'],
                                                                 dict_pams['Tc'] - mc.Tref,
                                                                 dict_pams['e'],
                                                                 dict_pams['omega'])

            mA = (dict_pams['mean_long'] - dict_pams['omega']) \
                 + self.dynamical_set['ttvfast']['t_beg'] / dict_pams['P'] * 360.0000000000

            params.extend([
                dict_pams['M_Me'] / constants.Msear,  # mass in Solar unit
                dict_pams['P'],
                dict_pams['e'],
                i_temp,
                dict_pams['Omega'],
                dict_pams['omega'],
                mA])

            n_plan += 1
            if P_min is None:
                P_min = dict_pams['P']
            else:
                P_min = np.min(np.asarray([P_min, dict_pams['P']]))

        if x_input is None:
            t_step = min(P_min / 20., self.dynamical_set['min_deltat'])
            pos, rv = ttvfast._ttvfast._ttvfast(params,
                                                t_step,
                                                self.dynamical_set['ttvfast']['t_beg'],
                                                self.dynamical_set['ttvfast']['t_end'],
                                                n_plan, input_flag,
                                                self.dynamical_set['len_rv'],
                                                self.dynamical_set['rv_times'])
        else:
            t_step = min(P_min / 20., np.min(np.abs(x_input[1:]-x_input[:-1])))
            pos, rv = ttvfast._ttvfast._ttvfast(params,
                                                t_step,
                                                self.dynamical_set['ttvfast']['t_beg'],
                                                self.dynamical_set['ttvfast']['t_end'],
                                                n_plan, input_flag,
                                                len(x_input), (x_input - mc.Tref).tolist())

        positions = np.asarray(pos)
        rv_meas = np.asarray(rv) * constants.AU / constants.d2s
        output = {}

        for dataset_name, dataset in mc.dataset_dict.items():
            if dataset.dynamical is False: continue
            if dataset.kind == 'RV':
                if x_input is None:
                    output[dataset_name] = rv_meas[self.dynamical_set['data_selection'][dataset_name]]
                else:
                    output[dataset_name] = rv_meas


            elif dataset.kind == 'Tcent':
                t0_sel = (positions[0][:] == plan_ref[dataset.planet_name])
                t0_mod = positions[2][t0_sel]  # T0 of the transit
                nn_mod = np.asarray(positions[1][t0_sel], dtype=np.int16)  # Number of transit as stored by TTVfast

                if np.size(t0_mod) == 0:
                    output[dataset_name] = dataset.n_transit * 0.0000
                    continue
                """TTVfast transit numbering precedes the one of our dataset by construction, so
                        we need to add an offset to the transit reference number in order to have the
                        right association. We do that by identifying the first TTVfast T0 corresponding to
                        the 0th transit in our dataset"""

                min_idx = np.argmin(np.abs(t0_mod - dataset.x0[0]))
                ref_idf = nn_mod[min_idx]

                if np.sum(t0_sel) > np.max(dataset.n_transit) + ref_idf:
                    output[dataset_name] = t0_mod[dataset.n_transit + ref_idf] + mc.Tref

                    """With this approach, missing T0s in the input dataset are automatically skipped """
                else:
                    """The output vector contains less T0s than supposed: collision between planets?? """
                    output[dataset_name] = dataset.n_transit * 0.0000

        return output
