from __future__ import print_function
from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel

import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.subroutines.results_analysis import get_stellar_parameters


class RVkeplerian(AbstractModel):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'radial_velocities'

        self.list_pams_common = {
            'P',  # Period
            'K',  # RV semi-amplitude
            'e',  # eccentricity, uniform prior - to be fixed
            'omega'}  # argument of pericenter

        self.use_time_inferior_conjunction = False
        self.use_mass = False

    def initialize_model(self, mc, **kwargs):

        if mc.common_models[self.planet_ref].parametrization[:8] == 'Ford2006' \
            and mc.common_models[self.planet_ref].orbit != 'circular':
            self.list_pams_common.discard('e')
            self.list_pams_common.discard('omega')

            self.list_pams_common.update(['e_coso'])
            self.list_pams_common.update(['e_sino'])

        elif mc.common_models[self.planet_ref].parametrization[:8] != 'Standard' \
            and mc.common_models[self.planet_ref].orbit != 'circular':
                # 'Eastman2013' is the standard choice
            self.list_pams_common.discard('e')
            self.list_pams_common.discard('omega')

            self.list_pams_common.update(['sre_coso'])
            self.list_pams_common.update(['sre_sino'])

        if mc.common_models[self.planet_ref].use_time_inferior_conjunction:
            self.list_pams_common.update(['Tc'])
            self.use_time_inferior_conjunction = True
            # Copying the property to the class for faster access
        else:
            self.list_pams_common.update(['mean_long'])

        if mc.common_models[self.planet_ref].use_mass:
            self.list_pams_common.update(['M_Me'])
            self.list_pams_common.update(['mass'])
            self.use_mass = True
        else:
            self.list_pams_common.update(['K'])

    def compute(self, parameter_values, dataset, x0_input=None):

        #if parameter_values['P']<10:
        #    print(parameter_values['P'], parameter_values['K'], parameter_values['e'], parameter_values['omega'])

        if self.use_time_inferior_conjunction:
            mean_long = kepler_exo.kepler_Tc2phase_Tref(parameter_values['P'],
                                                parameter_values['Tc'] - dataset.Tref,
                                                parameter_values['e'],
                                                parameter_values['omega'])
        else:
            mean_long = parameter_values['mean_long']

        if self.use_mass:

            K = kepler_exo.kepler_K1(parameter_values['mass'],
                                     parameter_values['M_Me'] / constants.Msear, parameter_values['P'], parameter_values['i'],
                                     parameter_values['e'])
        else:
            K = parameter_values['K']

        if x0_input is None:
            return kepler_exo.kepler_RV_T0P(dataset.x0,
                                            mean_long,
                                            parameter_values['P'],
                                            K,
                                            parameter_values['e'],
                                            parameter_values['omega'])
        else:
            return kepler_exo.kepler_RV_T0P(x0_input,
                                            mean_long,
                                            parameter_values['P'],
                                            K,
                                            parameter_values['e'],
                                            parameter_values['omega'])


class RVdynamical(AbstractModel):
    model_class = 'radial_velocities'

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

        if mc.common_models[self.planet_ref].use_time_inferior_conjunction:
            self.list_pams_common.update(['Tc'])
        else:
            self.list_pams_common.update(['mean_long'])


class TransitTimeKeplerian(AbstractModel):
    model_class = 'transit_times'

    def __init__(self, *args, **kwargs):
        super(TransitTimeKeplerian, self).__init__(*args, **kwargs)
        self.use_time_inferior_conjunction = False

        self.list_pams_common = {'P'}  # Period

        self.list_pams_dataset = set()

    def initialize_model(self, mc, **kwargs):

        if mc.common_models[self.planet_ref].use_time_inferior_conjunction:
            self.list_pams_common.update(['Tc'])
            self.use_time_inferior_conjunction = True
            # Copying the property to the class for faster access
        else:
            self.list_pams_common.update(['mean_long'])
            self.list_pams_common.update(['e'])
            self.list_pams_common.update(['omega'])
            # mean longitude = argument of pericenter + mean anomaly at Tref

    def compute(self, parameter_values, dataset, x0_input=None):

        if self.use_time_inferior_conjunction:
            delta_T = parameter_values['Tc'] - \
                np.floor((parameter_values['Tc'] - dataset.Tref) / parameter_values['P']) * parameter_values['P']
        else:
            delta_T = dataset.Tref + \
                      kepler_exo.kepler_phase2Tc_Tref(parameter_values['P'],
                                                      parameter_values['mean_long'],
                                                      parameter_values['e'],
                                                      parameter_values['omega'])

        if x0_input is None:
            return np.floor(dataset.x0 / parameter_values['P']) * parameter_values['P'] + delta_T
        else:
            return np.floor(x0_input / parameter_values['P']) * parameter_values['P'] + delta_T


class TransitTimeDynamical(AbstractModel):
    model_class = 'transit_times'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        ''' Orbital parameters to be used in the dynamical fit '''
        self.list_pams_common = {
            'P',     # Period in days
            'M_Me',  # Mass in Earth masses
            'Omega', # longitude of ascending node
            'e',     # eccentricity, uniform prior - to be fixed
            'i',
            'R_Rs',  # planet radius (in units of stellar radii)
            'omega', # argument of pericenter
            'mass'} # mass of the star (needed for proper dynamical computation and for reversibility)

        self.list_pams_dataset = set()

        self.use_semimajor_axis = False
        self.use_inclination = False
        self.use_time_inferior_conjunction = False


    def initialize_model(self, mc, **kwargs):


        if mc.common_models[self.planet_ref].use_inclination:
            """ i is the orbital inclination (in degrees) """
            self.list_pams_common.update(['i'])
        else:
            """ b is the impact parameter """
            self.list_pams_common.update(['b'])

        if mc.common_models[self.planet_ref].use_time_inferior_conjunction:
            self.list_pams_common.update(['Tc'])
        else:
            self.list_pams_common.update(['mean_long'])


class DynamicalIntegrator:
    def __init__(self):
        self.model_name = 'dynamical_integrator'
        self.dynamical_integrator = 'TRADES'

        self.dynamical_set = {}
        self.rv_dataset_idbool = {}
        self.t0_planet_idflag = {}
        self.planet_idflag = {}

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
                from pytrades import pytrades
            except (ModuleNotFoundError,ImportError):
                print("ERROR: TRADES not installed, this will not work")
                quit()

            self.prepare_trades(mc)

        if self.dynamical_integrator == 'ttvfast':
            try:
                import ttvfast
            except (ModuleNotFoundError,ImportError):
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
        int_buffer = dict(rv_time=[], rv_value=[], rv_error=[], rv_ref=[],
                          t0_time=[], t0_error=[], t0_ref=[], key_ref={})


        """ Putting all the RV epochs in the same array, flagging in the temporary buffer
            the stored values according to their dataset of origin
        """
        for dataset_name, dataset in mc.dataset_dict.items():
            if dataset.dynamical is False: continue
            if dataset.kind == 'RV':
                int_buffer['rv_time'].extend(dataset.x.tolist())
                int_buffer['rv_value'].extend(dataset.y.tolist())
                int_buffer['rv_error'].extend(dataset.e.tolist())
                int_buffer['rv_ref'].extend(dataset.x * 0.0 + dataset_rv)
                int_buffer['key_ref'][dataset_name] = dataset_rv
                dataset_rv += 1
            elif dataset.kind == 'transit_time':
                int_buffer['t0_time'].extend(dataset.x.tolist())

        """ Creating the flag array after all the RV epochs have been mixed
        """
        for dataset_name, dataset in mc.dataset_dict.items():
            if dataset.dynamical is False: continue
            if dataset.kind == 'RV':
                self.rv_dataset_idbool[dataset_name] = \
                    (np.asarray(int_buffer['rv_ref']) == int_buffer['key_ref'][dataset_name])

        rv_time = np.float64(int_buffer['rv_time'])
        rv_value = np.float64(int_buffer['rv_value'])
        rv_error = np.float64(int_buffer['rv_error'])

        """ Substituting empty arrays with None to allow minimum determination
            We give for granted that at least one of the two vector is not null, otherwise
            it doesn't make any sense to run the program at all
        """
        if np.size(int_buffer['rv_times']) == 0:
            int_buffer['rv_time'] =  int_buffer['t0_time']

        if np.size(int_buffer['t0_times']) == 0:
            int_buffer['t0_time'] = int_buffer['rv_time']

        rv_minmax = [np.amin(int_buffer['rv_times']), np.amax(int_buffer['rv_times'])]
        t0_minmax = [np.amin(int_buffer['t0_times']), np.amax(int_buffer['t0_times'])]

        self.ti_beg = np.min([rv_minmax[0], t0_minmax[0]]) - 10
        self.ti_end = np.max([rv_minmax[1], t0_minmax[1]]) + 10
        self.ti_int = self.ti_end - self.ti_beg
        self.ti_ref = np.float64(mc.Tref)
        self.i_step = np.float64(1.e-3)

        """First iteration to identify the number of transits
        stored for each planet, including the planets in the dynamical
        simulation but without observed transit
        """
        self.n_body = 1  # The star is included among the bodies
        self.duration_check = 1

        for planet_name in mc.dynamical_dict:
            self.n_body += 1
            self.planet_idflag[planet_name] = self.n_body * 1

        self.t0_flag = np.zeros(self.n_body, dtype=bool)


        # TRADES initialization
        pytrades.args_init(
            self.n_body, 
            self.duration_check,
            t_epoch=self.ti_ref,
            t_start=self.ti_beg,
            t_int=self.ti_int,
        )


        for planet_name, dataset_name in mc.dynamical_t0_dict.items():
            n_plan = self.planet_idflag[planet_name]
            self.t0_flag[n_plan-1] = True

            pytrades.set_t0_dataset(
                n_plan,
                mc.dataset_dict[dataset_name].n_transit[:].astype(int),
                mc.dataset_dict[dataset_name].x,
                mc.dataset_dict[dataset_name].e)

        if len(rv_time)>0:
            pytrades.set_rv_dataset(rv_time, rv_value, rv_error)#, rv_setid=rv_setid, n_rvset=n_rvset)

        self.to_be_initialized = False
        print()

        return

    def compute_trades(self, mc, theta, x_input=None):
        """ This function compute the expected TTV and RVs for dynamically interacting planets.
            The user can specify which planets are subject to interactions, e.g. long-period planets can be approximated
            with a Keplerian function"""

        if self.to_be_initialized:
            self.prepare_trades(mc)

        self.dynamical_pams = {
            'M': np.zeros(self.n_body, dtype=np.float64),
            'R': np.zeros(self.n_body, dtype=np.float64),
            'P': np.zeros(self.n_body, dtype=np.float64),
            'e': np.zeros(self.n_body, dtype=np.float64),
            'omega': np.zeros(self.n_body, dtype=np.float64),
            'i': np.zeros(self.n_body, dtype=np.float64),
            'Omega': np.zeros(self.n_body, dtype=np.float64),
            'mA': np.zeros(self.n_body, dtype=np.float64)
        }

        """ Adding star parameters"""
        star_pams = get_stellar_parameters(mc, theta, warnings=False)
        #star_pams = mc.common_models['star_parameters'].convert(theta)

        self.dynamical_pams['M'][0] = star_pams['mass']
        self.dynamical_pams['R'][0] = star_pams['radius']

        for planet_name in mc.dynamical_dict:
            n_plan = self.planet_idflag[planet_name]

            parameter_values = mc.common_models[planet_name].convert(theta)
            mc.common_models[planet_name].update_parameter_values(parameter_values, self.ti_ref )

            print(parameter_values)

            if 'R_Rs' in parameter_values:
                """ Converting the radius from Stellar units to Solar units"""
                self.dynamical_pams['R'][n_plan] = parameter_values['R_Rs'] * star_pams['radius']
            else:
                """ Default value: slightly more than 1 Earth radii in Solar units"""
                self.dynamical_pams['R'][n_plan] = 0.02

            self.dynamical_pams['M'][n_plan] = parameter_values['M_Me'] / constants.Msear
            self.dynamical_pams['P'][n_plan] = parameter_values['P']
            self.dynamical_pams['e'][n_plan] = parameter_values['e']
            self.dynamical_pams['omega'][n_plan] = parameter_values['omega']
            self.dynamical_pams['Omega'][n_plan] = parameter_values['Omega']
            self.dynamical_pams['mA'][n_plan] = (parameter_values['mean_long'] - parameter_values['omega'])


        if x_input is None:

            rv_sim, body_id_sim, epoch_sim, t0_sim, durations_sim, kep_elem_sim = pytrades.kelements_to_rv_and_t0s(
                self.ti_beg,
                self.ti_ref,
                self.ti_int,
                self.dynamical_pams['M'],
                self.dynamical_pams['R'],
                self.dynamical_pams['P'],
                self.dynamical_pams['e'],
                self.dynamical_pams['omega'],
                self.dynamical_pams['mA'],
                self.dynamical_pams['i'],
                self.dynamical_pams['Omega'],
                self.t0_flag,
            )

        else:

            pytrades.set_rv_dataset(x_input, np.zeros_like(x_input), np.zeros_like(x_input))#, rv_setid=rv_setid, n_rvset=n_rvset)
            self.to_be_initialized = True

            rv_sim, body_id_sim, epoch_sim, t0_sim, durations_sim, kep_elem_sim = pytrades.kelements_to_rv_and_t0s(
                self.ti_beg,
                self.ti_ref,
                self.ti_int,
                self.dynamical_pams['M'],
                self.dynamical_pams['R'],
                self.dynamical_pams['P'],
                self.dynamical_pams['e'],
                self.dynamical_pams['omega'],
                self.dynamical_pams['mA'],
                self.dynamical_pams['i'],
                self.dynamical_pams['Omega'],
                self.t0_flag,
            )
        output = {}

        for dataset_name, dataset in mc.dataset_dict.items():
            if dataset.dynamical is False: continue
            if dataset.kind == 'RV':
                if x_input is None:
                    output[dataset_name] = rv_sim[self.rv_dataset_idbool[dataset_name]]
                else:
                    output[dataset_name] = rv_sim
            elif dataset.kind == 'transit_time' and dataset.planet_name in mc.dynamical_dict:
                n_plan = self.planet_idflag[planet_name]
                data_sel = (body_id_sim==n_plan)
                output[dataset_name] = t0_sim[data_sel]
            elif dataset.kind == 'transit_duration' and dataset.planet_name in mc.dynamical_dict:
                n_plan = self.planet_idflag[planet_name]
                data_sel = (body_id_sim==n_plan)
                output[dataset_name] = durations_sim[data_sel]

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
            elif dataset.kind == 'transit_time':
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
                    a_temp = convert_rho_to_ars(dict_pams['P'], star_pams['density'])
                    i_temp = \
                        convert_b_to_i(dict_pams['b'],
                                       dict_pams['e'],
                                       dict_pams['omega'],
                                       a_temp)

            if mc.common_models[planet_name].use_time_inferior_conjunction:
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


            elif dataset.kind == 'transit_time':
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
