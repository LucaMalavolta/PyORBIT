from pyorbit.subroutines.common import np
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit

try:
    from PyAstronomy import modelSuite as PyAstroModelSuite
except ImportError:
    pass

class RossiterMcLaughling_Reloaded(AbstractModel, AbstractTransit):
    model_class = 'rossiter_mclaughlin'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # this calls all constructors up to AbstractModel
        super(AbstractModel, self).__init__(*args, **kwargs)

        self.unitary_model = False

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = {
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'lambda', # Sky-projected angle between stellar rotation axis and normal of orbit plane [deg]
            'R_Rs',  # planet radius (in units of stellar radii)
            'v_sini' # projected rotational velocity of the star
        }

        self.use_stellar_radius = True
        self.use_stellar_period = True
        self.use_stellar_inclination = False

        self.star_grid = {}   # write an empty dictionary


    def initialize_model(self, mc, **kwargs):

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_star_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        #start filling the dictionary with relevant parameters
        self.star_grid['n_grid'] = kwargs.get('star_ngrid', 51 )
        self.star_grid['half_grid'] = int((self.star_grid['n_grid'] - 1) / 2)
        self.star_grid['time_step'] = kwargs.get('time_step', 100 ) # in seconds

        """ Coordinates of the centers of each grid cell (add offset) """
        self.star_grid['xx'] = np.linspace(-1.000000, 1.000000, self.star_grid['n_grid'], dtype=np.double)
        self.star_grid['xc'], self.star_grid['yc'] = np.meshgrid(self.star_grid['xx'], self.star_grid['xx'], indexing='xy')
        # check the Note section of the wiki page of meshgrid
        # https://docs.scipy.org/doc/numpy/reference/generated/numpy.meshgrid.html

        """ Distance of each grid cell from the center of the stellar disk """
        self.star_grid['rc'] = np.sqrt(self.star_grid['xc'] ** 2 + self.star_grid['yc'] ** 2)
        # Must avoid negative numbers inside the square root
        self.star_grid['inside'] = self.star_grid['rc'] < 1.000000
        # Must avoid negative numbers inside the square root
        self.star_grid['outside'] = self.star_grid['rc'] >= 1.000000


        """ Determine the mu angle for each grid cell, as a function of radius. """
        self.star_grid['mu'] = np.zeros([self.star_grid['n_grid'], self.star_grid['n_grid']],dtype=np.double)  # initialization of the matrix with the mu values
        self.star_grid['mu'][self.star_grid['inside']] = np.sqrt(1. - self.star_grid['rc'][self.star_grid['inside']] ** 2)


    def compute(self, parameter_values, dataset, x0_input=None):
        """
        :param parameter_values:
        :param dataset:
        :param x0_input:
        :return:
        """
        #t1_start = process_time()

        

        par_a, par_i = self.retrieve_ai(parameter_values)
        par_tc = self.retrieve_t0(parameter_values, dataset.Tref)
        par_Omega, par_Is = self.retrieve_Omega_Istar(parameter_values)
        par_tc = self.retrieve_t0(parameter_values, dataset.Tref)
        par_lamba = parameter_values['lambda'] * constants.deg2rad

        inclination_rad = par_i * constants.deg2rad
        omega_rad = parameter_values['omega'] * constants.deg2rad


        """ Limb darkening law and coefficients """
        if self.limb_darkening_model == 'quadratic':
            ld_par = np.zeros(2)
            for par, i_par in self.ldvars.items():
                ld_par[i_par] = parameter_values[par]

            self.star_grid['I'] = 1 - ld_par[0]*(1. - self.star_grid['mu']) \
                - ld_par[1]*(1. - self.star_grid['mu'])**2
        else:
            print('ERROR: Selected limb darkening law not implemented')
            quit()

        self.star_grid['I'][self.star_grid['outside']] = 0.000

        """ Intensity normalization"""
        self.star_grid['I0'] = np.sum(self.star_grid['I'])
        self.star_grid['I'] /= self.star_grid['I0']

        self.star_grid['x_ortho'] = self.star_grid['xc'] * np.cos(par_lamba) \
            - self.star_grid['yc'] * np.sin(par_lamba)  # orthogonal distances from the spin-axis
        self.star_grid['y_ortho'] = self.star_grid['xc'] * np.sin(par_lamba) \
            + self.star_grid['yc'] * np.cos(par_lamba)


        star_grid['r_ortho'] = np.sqrt(star_grid['x_ortho'] ** 2 + star_grid['y_ortho'] ** 2)
        star_grid['z_ortho'] = np.zeros([star_grid['n_grid'], star_grid['n_grid']],
                                        dtype=np.double)  # initialization of the matrix
        star_grid['z_ortho'][star_grid['inside']] = np.sqrt(
            1. -star_grid['r_ortho'][star_grid['inside']] ** 2)

        """ rotate the coordinate system around the x_ortho axis by an agle: """
        star_grid['beta'] = (np.pi / 2.) - par_Is * constants.deg2rad



        """ orthogonal distance from the stellar equator """
        ### Equation 7 in Cegla+2016
        star_grid['yp_ortho'] = star_grid['z_ortho'] * np.sin(star_grid['beta']) \
            + star_grid['y_ortho'] * np.cos(star_grid['beta'])

        ### Equation 6 in Cegla+2016
        star_grid['zp_ortho'] = star_grid['z_ortho'] * np.cos(star_grid['beta']) \
            + star_grid['y_ortho'] * np.sin(star_grid['beta'])


        """ stellar rotational velocity for a given position """
        # differential rotation is included considering a sun-like law
        star_grid['v_star'] = star_grid['x_ortho'] * star_dict['vsini'] * (
            1. -star_dict['alpha'] * star_grid['yp_ortho'] ** 2)
        # Null velocity for points outside the stellar surface
        star_grid['v_star'][star_grid['outside']] = 0.0

        """ working arrays for Eq. 1 and 9 of Cegla+2016"""
        star_grid['muI'] = star_grid['I'] * star_grid['mu']
        star_grid['v_starI'] = star_grid['I'] * star_grid['v_star']



        eclipsed_flux = np.zeros_like(bjd)
        mean_mu = np.zeros_like(bjd)
        mean_vstar =  np.zeros_like(bjd)

        for i_obs, bjd_value in enumerate(bjd):
            n_oversampling = int(exptime[i_obs] / star_grid['time_step'])

            """recomputing the oversampling steps to homogeneously cover the
            full integration time """
            if n_oversampling % 2 == 0:
                n_oversampling += 1
                delta_step = exptime[i_obs] / n_oversampling / 86400.

            half_time = exptime[i_obs] / 2 / 86400.

            bjd_oversampling = np.linspace(bjd_value - half_time, bjd_value + half_time, n_oversampling, dtype=np.double)

            true_anomaly, orbital_distance_ratio = kepler_exo.kepler_true_anomaly_orbital_distance(
                bjd_oversampling - dataset.Tref,
                par_tc,
                parameter_values['P'],
                parameter_values['e'],
                omega_rad,
                parameter_values['a_Rs'])

            """ planet position during its orbital motion, in unit of stellar radius"""
            # Following Murray+Correia 2011 , with the argument of the ascending node set to zero.
            # 1) the ascending node coincide with the X axis
            # 2) the reference plance coincide with the plane of the sky
            planet_position = {
                'xp': -orbital_distance_ratio * (np.cos(omega_rad + true_anomaly)),
                'yp': orbital_distance_ratio * (np.sin(omega_rad + true_anomaly) * np.cos(inclination_rad)),
                'zp': orbital_distance_ratio * (np.sin(inclination_rad) * np.sin(omega_rad + true_anomaly))
            }
            # projected distance of the planet's center to the stellar center
            planet_position['rp'] = np.sqrt(planet_position['xp']**2  + planet_position['yp']**2)

            # iterating on the sub-exposures
            I_sum = 0.00
            muI_sum = 0.00
            vstarI_sum = 0.00

            for j, zeta in enumerate(planet_position['zp']):

                if zeta > 0 and planet_position['rp'][j] < 1. + planet_dict['Rp_Rs']:
                    # the planet is in the foreground or inside the stellar disk, continue
                    # adjustment: computation is performed even if only part of the planet is shadowing the star

                    rd = np.sqrt((planet_position['xp'][j] - star_grid['xc']) ** 2 +
                                    (planet_position['yp'][j] - star_grid['yc']) ** 2)

                    """ Seelction of the portion of stars covered by the planet"""
                    sel_eclipsed = (rd <= planet_dict['Rp_Rs']) & star_grid['inside']

                    I_sum +=  np.sum(star_grid['I'][sel_eclipsed])
                    muI_sum +=  np.sum(star_grid['muI'][sel_eclipsed])
                    vstarI_sum += np.sum(star_grid['v_starI'][sel_eclipsed])

            if muI_sum > 0:
                eclipsed_flux[i_obs] = I_sum/n_oversampling
                mean_mu[i_obs] = muI_sum/I_sum
                mean_vstar[i_obs] = vstarI_sum/I_sum














        if self.orbit == 'circular':
            self.rm_ohta.assignValue({"a": par_a,
                            "lambda": parameter_values['lambda']/180.*np.pi,
                            "epsilon": parameter_values['ld_c1'],
                            "P": parameter_values['P'],
                            "T0": par_tc,
                            "i": par_i/180.*np.pi,
                            "Is": par_Is/180.*np.pi,
                            "Omega": par_Omega,
                            "gamma": parameter_values['R_Rs']})
        else:

            if self.use_time_of_transit:
                Tperi  = kepler_exo.kepler_Tc2Tperi_Tref(parameter_values['P'],
                                                         par_tc,
                                                         parameter_values['e'],
                                                         parameter_values['omega'])
            else:
                Tperi  = kepler_exo.kepler_phase2Tperi_Tref(parameter_values['P'],
                                                         parameter_values['mean_long'],
                                                         parameter_values['e'],
                                                         parameter_values['omega'])

            self.rm_ohta.assignValue({"a": par_a,
                "lambda": parameter_values['lambda']/180.*np.pi,
                "epsilon": parameter_values['ld_c1'],
                "P": parameter_values['P'],
                "tau": Tperi,
                "i": par_i/180.*np.pi,
                "w": parameter_values['omega']/180.*np.pi-np.pi,
                "e":parameter_values['e'],
                "Is": par_Is/180.*np.pi,
                "Omega": par_Omega,
                "gamma": parameter_values['R_Rs']})

        if x0_input is None:
            return self.rm_ohta.evaluate(dataset.x0) * parameter_values['radius'] * constants.Rsun * 1000.
        else:
            return self.rm_ohta.evaluate(x0_input) * parameter_values['radius'] * constants.Rsun * 1000.
