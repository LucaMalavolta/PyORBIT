
from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel
import pyorbit.subroutines.kepler_exo as kepler_exo


try:
    from pytrades import pytrades
except (ModuleNotFoundError,ImportError):
    pass


class AbstractDynamical(object):

    def __init__(self, *args, **kwargs):

        ''' Orbital parameters to be used in the dynamical fit '''
        self.list_pams_common = OrderedSet([
            'P',     # Period in days
            'M_Me',  # Mass in Earth masses
            'Omega', # longitude of ascending node
            'e',     # eccentricity, uniform prior - to be fixed
            'R_Rs',  # planet radius (in units of stellar radii)
            'omega' # argument of pericenter
        ])

        self.list_pams_dataset = OrderedSet()
        self.warning_given = False

    # brainless workaround
    def _initialize_model(self, mc, **kwargs):

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

        if not mc.common_models[self.planet_ref].use_mass:
            print('*** Dynamical modelling requires the mass of the planet as free parameters ***')
            print('*** for efficient exploration of parameter space')
            print('I QUIT')
            quit()

        if mc.common_models[self.planet_ref].use_inclination:
            """ i is the orbital inclination (in degrees) """
            self.list_pams_common.update(['i'])
        else:
            """ b is the impact parameter """
            self.list_pams_common.update(['b'])
            self.list_pams_common.update(['density'])

        if mc.common_models[self.planet_ref].use_semimajor_axis:
            self.list_pams_common.update(['a_Rs'])

        try:
            multivariate_pams = mc.common_models[self.stellar_ref].multivariate_pams
        except AttributeError:
            multivariate_pams = []

        if 'mass' in multivariate_pams and 'radius' in multivariate_pams:
            self.list_pams_common.update(['mass'])
            self.list_pams_common.update(['radius'])

        if mc.common_models[self.stellar_ref].compute_density:
            self.list_pams_common.update(['mass'])
            self.list_pams_common.update(['radius'])
        elif mc.common_models[self.stellar_ref].compute_mass:
            self.list_pams_common.update(['density'])
            self.list_pams_common.update(['radius'])
        elif mc.common_models[self.stellar_ref].compute_radius:
            self.list_pams_common.update(['density'])
            self.list_pams_common.update(['mass'])

        if mc.common_models[self.planet_ref].use_time_inferior_conjunction:
            self.list_pams_common.update(['Tc'])
        else:
            self.list_pams_common.update(['mean_long'])


class RVdynamical(AbstractModel, AbstractDynamical):
    model_class = 'radial_velocities'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

    # brainless workaround
    def initialize_model(self, mc, **kwargs):
        self._initialize_model(mc, **kwargs)

class TransitTimeDynamical(AbstractModel, AbstractDynamical):
    model_class = 'transit_times'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

    # brainless workaround
    def initialize_model(self, mc, **kwargs):
        self._initialize_model(mc, **kwargs)


class PhotoDynamical(AbstractModel, AbstractDynamical):
    model_class = 'transit_times'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

    # brainless workaround
    def initialize_model(self, mc, **kwargs):
        self._initialize_model(mc, **kwargs)

        """ Setting up the limb darkening calculation"""

        self.limb_darkening_model = kwargs['limb_darkening_model']
        self.ld_vars = [0.00] * kwargs['limb_darkening_ncoeff']

        for common_model in self.common_ref:
            if mc.common_models[common_model].model_class == 'limb_darkening':
                ld_parametrization = getattr(mc.common_models[common_model], 'parametrization', 'Standard')

        if ld_parametrization=='Kipping':
            self.ldvars['ld_c1'] = 0
            self.ldvars['ld_c2'] = 1
            self.list_pams_common.update(['ld_q1'])
            self.list_pams_common.update(['ld_q2'])
        else:
            for i_coeff in range(1, kwargs['limb_darkening_ncoeff'] + 1):
                par = 'ld_c' + repr(i_coeff)
                self.ldvars[par] = i_coeff - 1
                self.list_pams_common.update([par])

        if self.limb_darkening_model == 'quadratic':
            self.compute_limb_darkening = self._limb_darkening_quadratic
        else:
            print('ERROR: Selected limb darkening law not implemented')
            quit()

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        supersample_names = ['supersample_factor',
                                'supersample',
                                'supersampling',
                                'oversample_factor',
                                'oversample',
                                'oversampling',
                                'sample_factor',
                                'sample',
                                'sampling'
                                'nsample_factor',
                                'nsample',
                                'nsampling'
                            ]

        sample_factor = 1
        exposure_time = 0.01

        for dict_name in supersample_names:
            if kwargs[dataset.name_ref].get(dict_name, False):
                sample_factor = kwargs[dataset.name_ref][dict_name]
            elif kwargs.get(dict_name, False):
                sample_factor = kwargs[dict_name]

        exptime_names = ['exposure_time',
                            'exposure',
                            'exp_time',
                            'exptime',
                            'obs_duration',
                            'integration',
                        ]

        for dict_name in exptime_names:
            if kwargs[dataset.name_ref].get(dict_name, False):
                exposure_time = kwargs[dataset.name_ref][dict_name]
            elif kwargs.get(dict_name, False):
                exposure_time = kwargs[dict_name]

        self.code_options[dataset.name_ref] = {
            'sample_factor': sample_factor,
            'exp_time': exposure_time / constants.d2s,
        }

    # ????????
    def compute(self, parameter_values, dataset, x0_input=None):
        # ??????
        for par, i_par in self.ldvars.items():
            self.ld_vars[i_par] = parameter_values[par]


    # not sure if needed
    def _limb_darkening_coefficients(self, parameter_values):
        ld_par = np.zeros(2)
        for par, i_par in self.ldvars.items():
            ld_par[i_par] = parameter_values[par]
        return ld_par

    def _limb_darkening_quadratic(self, ld_par, mu):
        return  1 - ld_par[0]*(1. - mu) - ld_par[1]*(1. - mu)**2

class DynamicalIntegrator:
    def __init__(self):
        self.model_name = 'dynamical_integrator'
        self.dynamical_integrator = 'TRADES'

        self.dynamical_set = {}
        self.rv_dataset_idbool = {}
        self.t0_planet_idflag = {}
        self.planet_idflag = {}

        self.to_be_initialized = True

        print('*** Dynamical modelling reuires the use of the stellar mass ***')
        print('This may cause a clash with models requiring stellar density and radius')
        print('As in the case of the RM effect. The use of a multivariate approach is strongly suggested')
        print('You can control the behaviour of mass/radius/density with the specific keywords')
        print('compute_mass, compute_radius, compute_density')
        print()

    def compute(self, mc, theta, *args, **kwargs):
        """
        Run the appropriate subroutine according to chosen dynamical integrator
        :param mc:
        :return:
        """

        if self.dynamical_integrator == 'TRADES':
            output = self.compute_trades(mc, theta, *args, **kwargs)
        #if self.dynamical_integrator == 'ttvfast':
        #    output = self.compute_ttvfast(mc, *args, **kwargs)
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
        if np.size(int_buffer['rv_time']) == 0:
            int_buffer['rv_time'] =  int_buffer['t0_time']

        if np.size(int_buffer['t0_time']) == 0:
            int_buffer['t0_time'] = int_buffer['rv_time']

        rv_minmax = [np.amin(int_buffer['rv_time']), np.amax(int_buffer['rv_time'])]
        t0_minmax = [np.amin(int_buffer['t0_time']), np.amax(int_buffer['t0_time'])]

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
                mc.dataset_dict[dataset_name].n_transit,
                mc.dataset_dict[dataset_name].x,
                mc.dataset_dict[dataset_name].e)

        if len(rv_time)>0:
            pytrades.set_rv_dataset(rv_time, rv_value, rv_error)#, rv_setid=rv_setid, n_rvset=n_rvset)

        self.to_be_initialized = False

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

        star_counter = 0
        for planet_name in mc.dynamical_dict:
            n_plan = self.planet_idflag[planet_name]

            if star_counter == 0:
                star_model = mc.common_models[planet_name].stellar_ref
                star_parameter = mc.common_models[star_model].convert(theta)
                self.dynamical_pams['M'][0] = star_parameter['mass']
                self.dynamical_pams['R'][0] = star_parameter['radius']
                star_counter += 1

            parameter_values = mc.common_models[planet_name].convert(theta)
            parameter_values.update(star_parameter)
            mc.common_models[planet_name].update_parameter_values_for_dynamical(parameter_values, self.ti_ref)

            if 'R_Rs' in parameter_values:
                """ Converting the radius from Stellar units to Solar units"""
                self.dynamical_pams['R'][n_plan-1] = parameter_values['R_Rs'] * star_parameter['radius']
            else:
                """ Default value: slightly more than 1 Earth radii in Solar units"""
                self.dynamical_pams['R'][n_plan-1] = 0.02

            self.dynamical_pams['M'][n_plan-1] = parameter_values['M_Me'] / constants.Msear
            self.dynamical_pams['P'][n_plan-1] = parameter_values['P']
            self.dynamical_pams['e'][n_plan-1] = parameter_values['e']
            self.dynamical_pams['i'][n_plan-1] = parameter_values['i']
            self.dynamical_pams['omega'][n_plan-1] = parameter_values['omega']
            self.dynamical_pams['Omega'][n_plan-1] = parameter_values['Omega']
            self.dynamical_pams['mA'][n_plan-1] = (parameter_values['mean_long'] - parameter_values['omega'])


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

                n_plan = self.planet_idflag[dataset.planet_name]
                data_sel = (body_id_sim==n_plan)
                output[dataset_name] = t0_sim[data_sel]
            elif dataset.kind == 'transit_duration' and dataset.planet_name in mc.dynamical_dict:
                n_plan = self.planet_idflag[planet_name]
                data_sel = (body_id_sim==n_plan)
                output[dataset_name] = durations_sim[data_sel]

        return output
