from pyorbit.subroutines.common import np
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import *
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d

try:
    from numba import jit
    @jit(nopython=True)
    def CCF_gauss(x, A, x0, fwhm):
        sigma = fwhm/2.35482004503
        return 1. - A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

    @jit(nopython=True)
    def iter3_CCF_gauss(x, A, fwhm, x0, istar, area, out):
        sigma = fwhm/2.35482004503
        const = (2 * sigma ** 2)

        for i in range(x0.shape[0]):
            for j in range(x0.shape[1]):
                out[i,j, :] = istar[i, j] * (1. - A * np.exp(-(x - x0[i, j]) ** 2 /const)) * area
        return out
except:

    def CCF_gauss(x, A, x0, fwhm):
        sigma = fwhm/2.35482004503
        return 1. - A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

    def iter3_CCF_gauss(x, A, fwhm, x0, istar, area, out):
        sigma = fwhm/2.35482004503
        const = (2 * sigma ** 2)

        for i in range(x0.shape[0]):
            for j in range(x0.shape[1]):
                out[i,j, :] = istar[i, j] * (1. - A * np.exp(-(x - x0[i, j]) ** 2 /const)) * area
        return out


class RossiterMcLaughling_Precise(AbstractModel, AbstractTransit):
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

        self.star_grid = {}   # write an empty dictionary


    def initialize_model(self, mc, **kwargs):

        self.use_differential_rotation = kwargs.get('use_differential_rotation', mc.common_models[self.stellar_ref].use_differential_rotation)

        """ check if the differential rotation should be included in the model"""
        if self.use_differential_rotation:
            self.list_pams_common.discard('v_sini')
            self.list_pams_common.update(['veq_star', 'i_star', 'alpha_rotation'])
            mc.common_models[self.stellar_ref].use_equatorial_velocity =  True
            mc.common_models[self.stellar_ref].use_stellar_inclination =  True

            """ If stellar rotation is provided as a prior or as the outcome of the fit,
                the code will check if the derived stellar radius is consistent with the prior
            """

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_star_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        #start filling the dictionary with relevant parameters
        self.star_grid['n_grid'] = kwargs.get('star_ngrid', 101 )
        self.star_grid['half_grid'] = int((self.star_grid['n_grid'] - 1) / 2)

        """ Coordinates of the centers of each grid cell (add offset) """
        self.star_grid['xx'] = np.linspace(-1.000000, 1.000000, self.star_grid['n_grid'], dtype=np.double)
        self.star_grid['xc'], self.star_grid['yc'] = np.meshgrid(self.star_grid['xx'], self.star_grid['xx'], indexing='xy')
        self.star_grid['area'] = (2. / self.star_grid['n_grid'])**2
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


        #start filling the dictionary with relevant parameters
        self.planet_grid['n_grid'] = kwargs.get('planet_ngrid', 21 )
        self.planet_grid['half_grid'] = int((self.planet_grid['n_grid'] - 1) / 2)
        self.planet_grid['time_step'] = kwargs.get('time_step', 149 ) # in seconds

        """ Coordinates of the centers of each grid cell (add offset) """
        self.planet_grid['xx'] = np.linspace(-1.000000, 1.000000, self.planet_grid['n_grid'], dtype=np.double)
        self.planet_grid['xc'], self.planet_grid['yc'] = np.meshgrid(self.planet_grid['xx'], self.planet_grid['xx'], indexing='xy')
        self.planet_grid['area'] = (2. / self.planet_grid['n_grid'])**2

        # check the Note section of the wiki page of meshgrid
        # https://docs.scipy.org/doc/numpy/reference/generated/numpy.meshgrid.html

        """ Distance of each grid cell from the center of the stellar disk """
        self.planet_grid['rc'] = np.sqrt(self.planet_grid['xc'] ** 2 + self.planet_grid['yc'] ** 2)
        # Must avoid negative numbers inside the square root
        self.planet_grid['inside'] = self.planet_grid['rc'] < 1.000000
        # Must avoid negative numbers inside the square root
        self.planet_grid['outside'] = self.planet_grid['rc'] >= 1.000000

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        self._prepare_dataset_options(mc, dataset, **kwargs)
        #self.batman_models[dataset.name_ref] = \
        #    batman.TransitModel(self.batman_params,
        #                        dataset.x0,
        #                        supersample_factor=self.code_options[dataset.name_ref]['sample_factor'],
        #                        exp_time=self.code_options[dataset.name_ref]['exp_time'],
        #                        nthreads=self.code_options['nthreads'])

        """ Values valid for HARPS-N data"""
        self.ccf_variables = {
            'natural_broadening': 1.5,  # in km/s
            'natural_contrast': 0.5, # relative depth o
            'instrumental_broadening': 1.005, # in km/s
            'rv_min': -20.0,
            'rv_max': 20.0,
            'rv_step': 0.25,
        }

        for dict_name in self.ccf_variables:
            if kwargs[dataset.name_ref].get(dict_name, False):
                self.ccf_variables[dict_name] = kwargs[dataset.name_ref][dict_name]
            elif kwargs.get(dict_name, False):
                self.ccf_variables[dict_name] = kwargs[dict_name]


        self.star_grid['zz'] = np.arange(self.ccf_variables['rv_min'],
                                        self.ccf_variables['rv_max']+self.ccf_variables['rv_step'],
                                        self.ccf_variables['rv_step'],
                                        dtype=np.double)
        self.star_grid['len_zz'] = len(self.star_grid['zz'])
        self.star_grid['rv_step'] = self.ccf_variables['rv_step']

    def compute(self, parameter_values, dataset, x0_input=None):

        """
        :param parameter_values:
        :param dataset:
        :param x0_input:
        :return:
        """
        self.update_parameter_values(parameter_values, dataset.Tref)

        ld_par = self._limb_darkening_coefficients(parameter_values)

        lambda_rad = parameter_values['lambda'] * constants.deg2rad
        inclination_rad = parameter_values['i'] * constants.deg2rad
        omega_rad = parameter_values['omega'] * constants.deg2rad

        sin_lambda = np.sin(lambda_rad)
        cos_lambda = np.cos(lambda_rad)

        if self.use_differential_rotation:
            beta = (np.pi / 2.) - parameter_values['i_star'] * constants.deg2rad
            sin_beta = np.sin(beta)
            cos_beta = np.cos(beta)



        """ Limb darkening law and coefficients """
        star_grid_I = self.compute_limb_darkening(ld_par, self.star_grid['mu'])
        star_grid_I[self.star_grid['outside']] = 0.000

        """ Intensity normalization"""
        star_grid_I /= np.sum(star_grid_I)

        star_grid_x_ortho = self.star_grid['xc'] * cos_lambda - self.star_grid['yc'] * sin_lambda  # orthogonal distances from the spin-axis
        star_grid_y_ortho = self.star_grid['xc'] * sin_lambda + self.star_grid['yc'] * cos_lambda

        star_grid_r2_ortho = star_grid_x_ortho ** 2 + star_grid_y_ortho**2
        star_grid_r2_ortho[self.star_grid['inside']] = 1.
        star_grid_z_ortho = np.sqrt(1. -star_grid_r2_ortho)

        if self.use_differential_rotation:
            istar_rad = parameter_values['i_star'] * constants.deg2rad

            """ rotate the coordinate system around the x_ortho axis by an angle: """
            star_grid_beta = (np.pi / 2.) - istar_rad

            """ orthogonal distance from the stellar equator """
            ### Equation 7 in Cegla+2016
            star_grid_yp_ortho = star_grid_z_ortho * sin_beta + star_grid_y_ortho * cos_beta
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

        """ Addition of the convective contribute"""
        star_grid_v_star += self.retrieve_convective_rv(ld_par, self.star_grid['mu'], parameter_values)

        star_grid_v_star[self.star_grid['outside']] = 0.0

        out_temp = np.empty([self.star_grid['n_grid'], self.star_grid['n_grid'], self.star_grid['len_zz']])
        star_grid_ccf = iter3_CCF_gauss(self.star_grid['zz'],
                            self.ccf_variables['natural_contrast'],
                            self.ccf_variables['natural_broadening'],
                            star_grid_v_star, star_grid_I, self.star_grid['area'], out_temp)

        ccf_total = np.sum(star_grid_ccf[sel_inside,:], axis=0)

        if x0_input is None:
            bjd = dataset.x - dataset.Tref
            exptime = dataset.ancillary['exptime']
            n_vals = dataset.n

        else:
            bjd = x0_input
            exptime = np.ones_like(bjd) * np.mean(dataset.ancillary['exptime'])
            n_vals = len(bjd)

        rv_rml = np.zeros(n_vals)

        p0 = (self.ccf_variables['natural_contrast'], 0.00, self.ccf_variables['instrumental_broadening']/self.star_grid['rv_step'])

        planet_area = self.planet_grid['area'] * parameter_values['R_Rs']**2

        for i_obs, bjd_value in enumerate(bjd):

            n_oversampling = int(exptime[i_obs] / self.planet_grid['time_step'])

            """recomputing the oversampling steps to homogeneously cover the
            full integration time """
            if n_oversampling % 2 == 0:
                n_oversampling += 1

            half_time = exptime[i_obs] / 2 / 86400.

            bjd_oversampling = np.linspace(bjd_value - half_time, bjd_value + half_time, n_oversampling, dtype=np.double)

            true_anomaly, orbital_distance_ratio = kepler_exo.kepler_true_anomaly_orbital_distance(
                bjd_oversampling,
                parameter_values['Tc']-dataset.Tref,
                parameter_values['P'],
                parameter_values['e'],
                parameter_values['omega'],
                parameter_values['a_Rs'])

            """ planet position during its orbital motion, in unit of stellar radius"""
            # Following Murray+Correia 2011 , with the argument of the ascending node set to zero.
            # 1) the ascending node coincide with the X axis
            # 2) the reference plance coincide with the plane of the sky

            planet_position_xp = -orbital_distance_ratio * (np.cos(omega_rad + true_anomaly))
            planet_position_yp = -orbital_distance_ratio * (np.sin(omega_rad + true_anomaly) * np.cos(inclination_rad))
            planet_position_zp = orbital_distance_ratio * (np.sin(inclination_rad) * np.sin(omega_rad + true_anomaly))
            planet_position_rp = np.sqrt(planet_position_xp**2  + planet_position_yp**2)

            # projected distance of the planet's center to the stellar center
            planet_position_rp = np.sqrt(planet_position_xp**2  + planet_position_yp**2)

            if np.amax(planet_position_zp) < 0.:
                continue

            if np.amin(planet_position_rp) > 1.:
                continue


            ccf_out = np.zeros(self.star_grid['len_zz'])



            for s_obs in range(0, n_oversampling):

                xp = self.planet_grid['xc']*parameter_values['R_Rs'] + planet_position_xp[s_obs]
                yp = self.planet_grid['yc']*parameter_values['R_Rs'] + planet_position_yp[s_obs]
                r2 = xp**2 + yp**2
                sel_inside = (r2 < 1.) & (self.planet_grid['inside'])
                r2[~sel_inside] = 1.

                mu_grid =  np.sqrt(1. - r2)
                I_grid = self.compute_limb_darkening(ld_par, mu_grid)

                x_ortho = xp * cos_lambda - yp * sin_lambda  # orthogonal distances from the spin-axis
                y_ortho = xp * sin_lambda + yp * cos_lambda

                r2_ortho = x_ortho ** 2 + y_ortho**2
                r2_ortho[~sel_inside] = 1.
                z_ortho = np.sqrt(1.-r2_ortho)


                if self.use_differential_rotation:

                    """ orthogonal distance from the stellar equator """
                    ### Equation 7 in Cegla+2016
                    yp_ortho = z_ortho * sin_beta + y_ortho * cos_beta

                    """ stellar rotational velocity for a given position """
                    # differential rotation is included considering a sun-like law
                    rv = x_ortho * parameter_values['v_sini'] * (1. -parameter_values['alpha_rotation'] * yp_ortho ** 2)
                    # Null velocity for points outside the stellar surface
                else:
                    rv = x_ortho * parameter_values['v_sini']

                rv += self.retrieve_convective_rv(ld_par, mu_grid, parameter_values)

                out_temp = np.empty([self.planet_grid['n_grid'], self.planet_grid['n_grid'], self.star_grid['len_zz']])
                planet_grid_ccf = iter3_CCF_gauss(self.star_grid['zz'],
                                    self.ccf_variables['natural_contrast'],
                                    self.ccf_variables['natural_broadening'],
                                    rv, I_grid, planet_area, out_temp)

                ccf_out += ccf_total - np.sum(planet_grid_ccf[sel_inside,:], axis=0)

            cont_val = np.amax(ccf_out)

            if cont_val > 0:
                ccf_out /= cont_val

                ccf_broad = gaussian_filter1d(ccf_out, self.ccf_variables['instrumental_broadening']/self.star_grid['rv_step'])
                try:
                    parameters, _ = curve_fit(CCF_gauss, self.star_grid['zz'], ccf_broad, p0=p0, check_finite =False)
                    rv_rml[i_obs] = parameters[1] * 1000.
                except:
                    rv_rml[i_obs] = 0.00

        return rv_rml

