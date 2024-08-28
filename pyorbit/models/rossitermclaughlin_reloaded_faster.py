from pyorbit.subroutines.common import np, OrderedSet
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import *

class RossiterMcLaughlin_Reloaded_Faster(AbstractModel, AbstractTransit):
    model_class = 'rossiter_mclaughlin'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # this calls all constructors up to AbstractModel
        super(AbstractModel, self).__init__(*args, **kwargs)

        self.unitary_model = False

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = OrderedSet([
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'lambda', # Sky-projected angle between stellar rotation axis and normal of orbit plane [deg]
            'R_Rs',  # planet radius (in units of stellar radii)
            #'v_sini' # projected rotational velocity of the star
        ])

        self.star_grid = {}   # write an empty dictionary


    def initialize_model(self, mc, **kwargs):

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_star_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        #start filling the dictionary with relevant parameters
        self.star_grid['n_grid'] = kwargs.get('star_ngrid', 51 )
        self.star_grid['half_grid'] = int((self.star_grid['n_grid'] - 1) / 2)
        self.star_grid['time_step'] = kwargs.get('time_step', 145 ) # in seconds

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
        self.update_parameter_values(parameter_values, dataset.Tref)

        for key, key_val in parameter_values.items():
            if np.isnan(key_val):
                return 0.

        ld_par = self._limb_darkening_coefficients(parameter_values)

        istar_rad = parameter_values['i_star'] * constants.deg2rad
        lambda_rad = parameter_values['lambda'] * constants.deg2rad
        inclination_rad = parameter_values['i'] * constants.deg2rad
        omega_rad = parameter_values['omega'] * constants.deg2rad

        #print(parameter_values['v_sini'], parameter_values['veq_star']*np.sin(istar_rad))

        """ Limb darkening law and coefficients """
        star_grid_I = self.compute_limb_darkening(ld_par, self.star_grid['mu'])
        star_grid_I[self.star_grid['outside']] = 0.000

        """ Intensity normalization"""
        star_grid_I /= np.sum(star_grid_I)

        star_grid_x_ortho = self.star_grid['xc'] * np.cos(lambda_rad) \
            - self.star_grid['yc'] * np.sin(lambda_rad)  # orthogonal distances from the spin-axis
        star_grid_y_ortho = self.star_grid['xc'] * np.sin(lambda_rad) \
            + self.star_grid['yc'] * np.cos(lambda_rad)


        star_grid_r_ortho = np.sqrt(star_grid_x_ortho ** 2 + star_grid_y_ortho** 2)
        star_grid_z_ortho = star_grid_r_ortho * 0.  # initialization of the matrix
        star_grid_z_ortho[self.star_grid['inside']] = np.sqrt(
            1. -star_grid_r_ortho[self.star_grid['inside']] ** 2)

        if self.use_differential_rotation:

            """ rotate the coordinate system around the x_ortho axis by an agle: """
            star_grid_beta = (np.pi / 2.) - istar_rad

            """ orthogonal distance from the stellar equator """
            ### Equation 7 in Cegla+2016
            star_grid_yp_ortho = star_grid_z_ortho * np.sin(star_grid_beta) \
                + star_grid_y_ortho * np.cos(star_grid_beta)

            ### Equation 6 in Cegla+2016
            #star_grid_zp_ortho = star_grid_z_ortho * np.cos(star_grid_beta) \
            #    + star_grid_y_ortho * np.sin(star_grid_beta)

            """ stellar rotational velocity for a given position """
            # differential rotation is included considering a sun-like law
            star_grid_v_star = star_grid_x_ortho * parameter_values['v_sini'] * (
                1. -parameter_values['alpha_rotation'] * star_grid_yp_ortho ** 2)
            # Null velocity for points outside the stellar surface
        else:
            star_grid_v_star = star_grid_x_ortho * parameter_values['v_sini']


        convective_c0 = self.retrieve_convective_c0(ld_par, parameter_values)
        star_grid_v_star += convective_c0 + self.retrieve_convective_rv(self.star_grid['mu'], parameter_values)

        star_grid_v_star[self.star_grid['outside']] = 0.0


        """ working arrays for Eq. 1 and 9 of Cegla+2016"""
        star_grid_muI = star_grid_I * self.star_grid['mu']
        star_grid_v_starI = star_grid_I * star_grid_v_star

        if x0_input is None:
            bjd = dataset.x
            exptime = dataset.ancillary['exptime']

        else:
            bjd = x0_input + dataset.Tref
            exptime = np.ones_like(bjd) * np.mean(dataset.ancillary['exptime'])

        #eclipsed_flux = np.zeros_like(bjd)
        mean_mu = np.zeros_like(bjd)
        mean_vstar =  np.zeros_like(bjd)


        for i_obs, bjd_value in enumerate(bjd):
            n_oversampling = int(exptime[i_obs] / self.star_grid['time_step'])

            """recomputing the oversampling steps to homogeneously cover the
            full integration time """
            if n_oversampling % 2 == 0:
                n_oversampling += 1

            half_time = exptime[i_obs] / 2 / 86400.

            bjd_oversampling = np.linspace(bjd_value - half_time, bjd_value + half_time, n_oversampling, dtype=np.double)

            true_anomaly, orbital_distance_ratio = kepler_exo.kepler_true_anomaly_orbital_distance(
                bjd_oversampling - dataset.Tref,
                parameter_values['Tc']-dataset.Tref,
                parameter_values['P'],
                parameter_values['e'],
                parameter_values['omega'],
                parameter_values['a_Rs'])

            """ planet position during its orbital motion, in unit of stellar radius"""
            # Following Murray+Correia 2011 , with the argument of the ascending node set to zero.
            # 1) the ascending node coincide with the X axis
            # 2) the reference plane coincide with the plane of the sky

            planet_position_xp = - orbital_distance_ratio * (np.cos(omega_rad + true_anomaly))
            planet_position_yp = - orbital_distance_ratio * (np.sin(omega_rad + true_anomaly) * np.cos(inclination_rad))
            planet_position_zp = orbital_distance_ratio * (np.sin(inclination_rad) * np.sin(omega_rad + true_anomaly))

            # projected distance of the planet's center to the stellar center
            planet_position_rp = np.sqrt(planet_position_xp**2  + planet_position_yp**2)

            #print('xp', planet_position_xp)
            #print('yp', planet_position_yp)
            #print('rp', planet_position_rp)

            # iterating on the sub-exposures
            I_sum = 0.00
            muI_sum = 0.00
            vstarI_sum = 0.00

            for j, zeta in enumerate(planet_position_zp):

                if zeta > 0 and planet_position_rp[j] < 1. + parameter_values['R_Rs']:
                    # the planet is in the foreground or inside the stellar disk, continue
                    # adjustment: computation is performed even if only part of the planet is shadowing the star

                    rd = np.sqrt((planet_position_xp[j] - self.star_grid['xc']) ** 2 +
                                    (planet_position_yp[j] - self.star_grid['yc']) ** 2)

                    """ Seelction of the portion of stars covered by the planet"""
                    sel_eclipsed = (rd <= parameter_values['R_Rs']) & self.star_grid['inside']

                    I_sum +=  np.sum(star_grid_I[sel_eclipsed])
                    muI_sum +=  np.sum(star_grid_muI[sel_eclipsed])
                    vstarI_sum += np.sum(star_grid_v_starI[sel_eclipsed])

            if muI_sum > 0:
                #eclipsed_flux[i_obs] = I_sum/n_oversampling
                mean_mu[i_obs] = muI_sum/I_sum
                mean_vstar[i_obs] = vstarI_sum/I_sum

        return mean_vstar

