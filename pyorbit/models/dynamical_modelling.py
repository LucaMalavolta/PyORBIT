
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
            'omega', # argument of pericenter
            #'i',     # inclination in degrees
        ])

        self.list_pams_dataset = OrderedSet()
        self.warning_given = False
        self.dynamical_model = True

    # brainless workaround
    def _prepare_dynamical_parameters(self, mc, **kwargs):

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

        try:
            multivariate_pams = mc.common_models[self.stellar_ref].multivariate_pams
        except AttributeError:
            multivariate_pams = []

        if not mc.common_models[self.planet_ref].use_mass:
            print('*** Dynamical modelling requires the mass of the planet as free parameters ***')
            print('*** for efficient exploration of parameter space')
            print('I QUIT')
            quit()


        try:
            multivariate_pams = mc.common_models[self.stellar_ref].multivariate_pams
        except AttributeError:
            multivariate_pams = []

        if 'mass' in multivariate_pams and 'radius' in multivariate_pams:
            self.list_pams_common.update(['mass'])
            self.list_pams_common.update(['radius'])
        elif mc.common_models[self.stellar_ref].compute_density:
            self.list_pams_common.update(['mass'])
            self.list_pams_common.update(['radius'])
        elif mc.common_models[self.stellar_ref].compute_mass:
            self.list_pams_common.update(['density'])
            self.list_pams_common.update(['radius'])
        elif mc.common_models[self.stellar_ref].compute_radius:
            self.list_pams_common.update(['density'])
            self.list_pams_common.update(['mass'])


        if mc.common_models[self.planet_ref].use_inclination:
            """ i is the orbital inclination (in degrees) """
            self.list_pams_common.update(['i'])
            self.compute_inclination = False
        else:
            """ b is the impact parameter """
            self.list_pams_common.update(['b'])
            self.compute_inclination = True


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
        self._prepare_dynamical_parameters(mc, **kwargs)

class TransitTimeDynamical(AbstractModel, AbstractDynamical):
    model_class = 'transit_times'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

    # brainless workaround
    def initialize_model(self, mc, **kwargs):
        self._prepare_dynamical_parameters(mc, **kwargs)


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

    def compute(self, mc, theta, x_input=None, *args, **kwargs):
        """
        Run the appropriate subroutine according to chosen dynamical integrator
        :param mc:
        :return:
        """

        if self.dynamical_integrator == 'TRADES':
            output = self.compute_trades(mc, theta, x_input, *args, **kwargs)
        #if self.dynamical_integrator == 'ttvfast':
        #    output = self.compute_ttvfast(mc, *args, **kwargs)
        return output

    def prepare_trades(self, mc):
        """
        :param mc:
        :return:
        """
        dataset_rv = 0
        dataset_lc = 0
        int_buffer = dict(rv_time=[], rv_ref=[],
                          lc_time=[], lc_ref=[],
                          t0_time=[], t0_ref=[], key_ref={})


        """ Putting all the RV epochs in the same array, flagging in the temporary buffer
            the stored values according to their dataset of origin
        """
        for dataset_name, dataset in mc.dataset_dict.items():
            if dataset.dynamical is False: continue
            if dataset.kind == 'radial_velocity':
                int_buffer['rv_time'].extend(dataset.x.tolist())
                int_buffer['rv_ref'].extend(dataset.x * 0.0 + dataset_rv)
                int_buffer['key_ref'][dataset_name] = dataset_rv
                dataset_rv += 1
            elif dataset.kind == 'photometry':
                int_buffer['lc_time'].extend(dataset.x.tolist())
                dataset_lc += 1
            elif dataset.kind == 'transit_time':
                int_buffer['t0_time'].extend(dataset.x.tolist())

        """ Creating the flag array after all the RV epochs have been mixed
        """
        for dataset_name, dataset in mc.dataset_dict.items():
            if dataset.dynamical is False: continue
            if dataset.kind == 'radial_velocity':
                self.rv_dataset_idbool[dataset_name] = \
                    (np.asarray(int_buffer['rv_ref']) == int_buffer['key_ref'][dataset_name])

        try:
            rv_minmax = [np.amin(int_buffer['rv_time']), np.amax(int_buffer['rv_time'])]
        except ValueError:
            rv_minmax = [mc.Tref-2, mc.Tref+2]

        try:
            t0_minmax = [np.amin(int_buffer['t0_time']), np.amax(int_buffer['t0_time'])]
        except ValueError:
            t0_minmax = [mc.Tref-2, mc.Tref+2]

        try:
            lc_minmax = [np.amin(int_buffer['lc_time']), np.amax(int_buffer['lc_time'])]
        except ValueError:
            lc_minmax = [mc.Tref-2, mc.Tref+2]

        self.ti_beg = np.min([rv_minmax[0], t0_minmax[0], lc_minmax[0]]) - 10
        self.ti_end = np.max([rv_minmax[1], t0_minmax[1], lc_minmax[1]]) + 10
        self.ti_int = self.ti_end - self.ti_beg
        self.ti_ref = np.float64(mc.Tref)
        self.i_step = np.float64(1.e-3)


        """First iteration to identify the number of transits
        stored for each planet, including the planets in the dynamical
        simulation but without observed transit
        """
        self.n_body = 1  # The star is included among the bodies
        self.duration_check = 1
        self.encounter_check=True
        self.do_hill_check=False
        self.amd_hill_check=False
        self.rv_res_gls=False

        for planet_name in mc.dynamical_dict:
            self.n_body += 1
            self.planet_idflag[planet_name] = self.n_body * 1

        self.rv_epochs_argsort = np.argsort(int_buffer['rv_time'])
        self.rv_epochs = np.asarray(int_buffer['rv_time'], dtype=np.float64)[self.rv_epochs_argsort]

        pytrades.args_init(
            self.n_body, # mandatory
            self.duration_check, # duration_check # mandatory
            t_epoch=self.ti_ref, # not needed here
            t_start=self.ti_beg, # not needed here
            t_int=self.ti_int, # not needed here
            encounter_check=self.encounter_check, # better alway True, we do not want close encounters!
            do_hill_check=self.do_hill_check, # at will, as input option
            amd_hill_check=self.amd_hill_check, # at will, as input option
            rv_res_gls=self.rv_res_gls, # at will, as input option
        )

        self.to_be_initialized = False

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

        return


    def compute_trades(self, mc, theta, x_input):
        """ This function compute the expected TTV and RVs for dynamically interacting planets.
            The user can specify which planets are subject to interactions, e.g. long-period planets can be approximated
            with a Keplerian function"""

        if self.to_be_initialized:
            self.prepare_trades(mc)


        """ Adding star parameters"""

        star_counter = 0
        for planet_name in mc.dynamical_dict:
            n_plan = self.planet_idflag[planet_name]

            if star_counter == 0:
                star_model = mc.common_models[planet_name].stellar_ref
                star_parameters = mc.common_models[star_model].convert(theta)
                self.dynamical_pams['M'][0] = star_parameters['mass']
                self.dynamical_pams['R'][0] = star_parameters['radius']
                star_parameters['density'] = star_parameters['mass']/star_parameters['radius']**3
                star_counter += 1

            parameter_values = mc.common_models[planet_name].convert(theta)
            parameter_values.update(star_parameters)
            mc.common_models[planet_name].update_parameter_values_for_dynamical(parameter_values, self.ti_ref)

            if 'R_Rs' in parameter_values:
                """ Converting the radius from Stellar units to Solar units"""
                self.dynamical_pams['R'][n_plan-1] = parameter_values['R_Rs'] * star_parameters['radius']
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

            (
                time_steps,
                orbits,
                transits,
                durations,
                lambda_rm,
                kep_elem,
                body_flag,
                rv_sim, # {"time": array, "rv": array} the index will match the input t_rv_obs array
                stable,
            ) = pytrades.orbital_parameters_to_transits(
                self.ti_ref,
                self.ti_beg,
                self.ti_int,
                self.dynamical_pams['M'],
                self.dynamical_pams['R'],
                self.dynamical_pams['P'],
                self.dynamical_pams['e'],
                self.dynamical_pams['omega'],
                self.dynamical_pams['mA'],
                self.dynamical_pams['i'],
                self.dynamical_pams['Omega'],
                self.rv_epochs # this can be an empty list [] and it will ignore it, otherwise provide a list time at which compute RV
            )

            rv_sorted = np.zeros_like(rv_sim['rv'], dtype=bool)
            rv_sorted[self.rv_epochs_argsort] = rv_sim['rv']

        else:

            (
                time_steps,
                orbits,
                transits,
                durations,
                lambda_rm,
                kep_elem,
                body_flag,
                rv_sim, # {"time": array, "rv": array} the index will match the input t_rv_obs array
                stable,
            ) = pytrades.orbital_parameters_to_transits(
                self.ti_ref,
                self.ti_beg,
                self.ti_int,
                self.dynamical_pams['M'],
                self.dynamical_pams['R'],
                self.dynamical_pams['P'],
                self.dynamical_pams['e'],
                self.dynamical_pams['omega'],
                self.dynamical_pams['mA'],
                self.dynamical_pams['i'],
                self.dynamical_pams['Omega'],
                x_input # this can be an empty list [] and it will ignore it, otherwise provide a list time at which compute RV
            )
            rv_sorted = rv_sim['rv']

        output = {'stable': stable, 'pass': True}

        #print('PASSED', self.dynamical_pams)


        for dataset_name, dataset in mc.dataset_dict.items():

            if dataset.dynamical is False: continue

            if dataset.kind == 'radial_velocity':
                if x_input is None:
                    output[dataset_name] = rv_sorted[self.rv_dataset_idbool[dataset_name]]
                else:
                    output[dataset_name] = rv_sorted

            elif dataset.kind == 'photometry':

                rp_rs, per, aRs, inc, ecc, w = pytrades.set_transit_parameters(
                    self.dynamical_pams['R'], transits, body_flag, kep_elem)

                output[dataset_name]  = {}

                for planet_name in mc.dynamical_dict:
                    n_plan = self.planet_idflag[planet_name]

                    planet_selection = (body_flag==n_plan)
                    output[dataset_name][planet_name] = {
                        'Rp_Rs': rp_rs[planet_selection],
                        'P': per[planet_selection],
                        'a_Rs': aRs[planet_selection],
                        'i': inc[planet_selection],
                        'e': ecc[planet_selection],
                        'omega': w[planet_selection],
                        'transits': transits[planet_selection],
                        'durations': durations[planet_selection],
                    }

            elif dataset.kind == 'transit_time' and dataset.planet_name in mc.dynamical_dict:

                n_plan = self.planet_idflag[dataset.planet_name]

                transits_planets = transits[body_flag==n_plan]
                epoch_synthetic =  np.rint( ( transits_planets - dataset.x[0] ) / self.dynamical_pams['P'][n_plan-1] ) + dataset.n_transit[0]

                find_epoch = np.isin(epoch_synthetic, dataset.n_transit)
                output[dataset_name] = transits_planets[find_epoch]

                if np.sum(find_epoch) != dataset.n:
                    output[dataset_name] = np.zeros(dataset.n)
                    output['pass'] = False

            elif dataset.kind == 'transit_duration' and dataset.planet_name in mc.dynamical_dict:
                n_plan = self.planet_idflag[dataset.planet_name]

                transits_planets = transits[body_flag==n_plan]
                durations_planets = durations[body_flag==n_plan]
                epoch_synthetic =  np.rint( ( transits_planets - dataset.x[0] ) / self.dynamical_pams['P'][n_plan-1] ) + dataset.n_transit[0]

                find_epoch = np.isin(epoch_synthetic, dataset.n_transit)
                output[dataset_name] = durations_planets[find_epoch]

                if np.sum(find_epoch) != dataset.n:
                    output[dataset_name] = np.zeros(dataset.n)
                    output['pass'] = False


        return output
