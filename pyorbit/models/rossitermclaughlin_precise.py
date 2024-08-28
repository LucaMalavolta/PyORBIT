from pyorbit.subroutines.common import np, OrderedSet
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
    def iter2_CCF_gauss(x, A, fwhm, x0, istar, out):
        sigma = fwhm/2.35482004503
        const = (2 * sigma ** 2)

        for i in range(x0.shape[0]):
            for j in range(x0.shape[1]):
                out[i,j, :] = istar[i, j] * (1. - A * np.exp(-(x - x0[i, j]) ** 2 /const))
        return out
except:

    def CCF_gauss(x, A, x0, fwhm):
        sigma = fwhm/2.35482004503
        return 1. - A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

    def iter2_CCF_gauss(x, A, fwhm, x0, istar, out):
        sigma = fwhm/2.35482004503
        const = (2 * sigma ** 2)

        for i in range(x0.shape[0]):
            for j in range(x0.shape[1]):
                out[i,j, :] = istar[i, j] * (1. - A * np.exp(-(x - x0[i, j]) ** 2 /const))
        return out

try:
    from lmfit.models import GaussianModel
except (ModuleNotFoundError, ImportError):
    pass

class RossiterMcLaughlin_Precise(AbstractModel, AbstractTransit):
    model_class = 'rossiter_mclaughlin'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # this calls all constructors up to AbstractModel
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            from lmfit.models import GaussianModel
        except (ModuleNotFoundError,ImportError):
            print("ERROR: lmfit not installed, this will not work")
            quit()

        import warnings
        from scipy.optimize import OptimizeWarning
        warnings.filterwarnings("ignore", category=OptimizeWarning)

        self.unitary_model = False

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = OrderedSet([
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'lambda', # Sky-projected angle between stellar rotation axis and normal of orbit plane [deg]
            'R_Rs',  # planet radius (in units of stellar radii)
        ])

        self.star_grid = {}   # write an empty dictionary
        self.model_class = 'rossiter_mclaughlin'


    def initialize_model(self, mc, **kwargs):


        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_star_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        #start filling the dictionary with relevant parameters
        self.star_grid['n_grid'] = kwargs.get('star_ngrid', 101 )
        self.star_grid['half_grid'] = int((self.star_grid['n_grid'] - 1) / 2)
        self.star_grid['time_step'] = kwargs.get('time_step', 149 ) # in seconds

        """ Coordinates of the centers of each grid cell (add offset) """
        self.star_grid['xx'] = np.linspace(-1.000000, 1.000000, self.star_grid['n_grid'], dtype=np.single)
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
        self.star_grid['mu'] = np.zeros([self.star_grid['n_grid'], self.star_grid['n_grid']], dtype=np.single)  # initialization of the matrix with the mu values
        self.star_grid['mu'][self.star_grid['inside']] = np.sqrt(1. - self.star_grid['rc'][self.star_grid['inside']] ** 2)

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
            'instrumental_broadening': 1.108, # in km/s
            'rv_min': -20.0,
            'rv_max': 20.0,
            'rv_step': 0.50,
        }

        # HARPS-N instrumental braodening: FWHM of 3.1 pixels
        # At 530 nm for example, in the center of the CCD,
        # the scale is 0.001415 nm/pixel, and the spectral
        # resolution taking into account the measured spotsize
        # is computed to about R = 124’000.
        # fwhm = 2.355 sigma
        # fwhm = 0.01415 * 3.1 * 299792.458 / 5300 = 2.48 km/s
        # sigma =

        for dict_name in self.ccf_variables:
            if kwargs[dataset.name_ref].get(dict_name, False):
                self.ccf_variables[dict_name] = kwargs[dataset.name_ref][dict_name]
            elif kwargs.get(dict_name, False):
                self.ccf_variables[dict_name] = kwargs[dict_name]

        self.star_grid[dataset.name_ref] = {}
        self.star_grid['zz'] = np.arange(self.ccf_variables['rv_min'],
                                        self.ccf_variables['rv_max']+self.ccf_variables['rv_step'],
                                        self.ccf_variables['rv_step'],
                                        dtype=np.single)
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

        for key, key_val in parameter_values.items():
            if np.isnan(key_val):
                return 0.


        ld_par = self._limb_darkening_coefficients(parameter_values)

        lambda_rad = parameter_values['lambda'] * constants.deg2rad
        inclination_rad = parameter_values['i'] * constants.deg2rad
        omega_rad = parameter_values['omega'] * constants.deg2rad

        """ Limb darkening law and coefficients """
        star_grid_I = self.compute_limb_darkening(ld_par, self.star_grid['mu'])
        star_grid_I[self.star_grid['outside']] = 0.000

        """ Intensity normalization"""
        star_grid_I /= np.sum(star_grid_I, dtype=np.single)

        star_grid_x_ortho = self.star_grid['xc'] * np.cos(lambda_rad, dtype=np.single) \
            - self.star_grid['yc'] * np.sin(lambda_rad, dtype=np.single)  # orthogonal distances from the spin-axis
        star_grid_y_ortho = self.star_grid['xc'] * np.sin(lambda_rad, dtype=np.single) \
            + self.star_grid['yc'] * np.cos(lambda_rad, dtype=np.single)


        star_grid_r_ortho = np.sqrt(star_grid_x_ortho ** 2 + star_grid_y_ortho** 2)
        star_grid_z_ortho = star_grid_r_ortho * 0.  # initialization of the matrix
        star_grid_z_ortho[self.star_grid['inside']] = np.sqrt(
            1. -star_grid_r_ortho[self.star_grid['inside']] ** 2)

        if self.use_differential_rotation:
            istar_rad = parameter_values['i_star'] * constants.deg2rad

            """ rotate the coordinate system around the x_ortho axis by an angle: """
            star_grid_beta = (np.pi / 2.) - istar_rad

            """ orthogonal distance from the stellar equator """
            ### Equation 7 in Cegla+2016
            star_grid_yp_ortho = star_grid_z_ortho * np.sin(star_grid_beta, dtype=np.single) \
                + star_grid_y_ortho * np.cos(star_grid_beta, dtype=np.single)

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
        convective_c0 = self.retrieve_convective_c0(ld_par, parameter_values)
        star_grid_v_star += convective_c0 + self.retrieve_convective_rv(self.star_grid['mu'], parameter_values)

        star_grid_v_star[self.star_grid['outside']] = 0.0

        out_temp = np.empty([self.star_grid['n_grid'], self.star_grid['n_grid'], self.star_grid['len_zz']], dtype=np.single)
        star_grid_ccf = iter2_CCF_gauss(self.star_grid['zz'],
                            self.ccf_variables['natural_contrast'],
                            self.ccf_variables['natural_broadening'],
                            star_grid_v_star, star_grid_I, out_temp)

        if x0_input is None:
            bjd = dataset.x - dataset.Tref
            exptime = dataset.ancillary['exptime']

        else:
            bjd = x0_input
            exptime = np.ones_like(bjd) * np.mean(dataset.ancillary['exptime'])

        rv_rml = np.zeros_like(bjd, dtype=np.single)
        ccf_total = np.sum(star_grid_ccf, axis=(0,1), dtype=np.single)

        p0 = (self.ccf_variables['natural_contrast'], 0.00, self.ccf_variables['instrumental_broadening']/self.star_grid['rv_step'])

        # RV unperturbed CCF
        ccf_broad = gaussian_filter1d(ccf_total/np.amax(ccf_total), self.ccf_variables['instrumental_broadening']/self.star_grid['rv_step'])
        
        gaussian = GaussianModel()
        
        params = gaussian.make_params()
        guess = gaussian.guess(1.-ccf_broad, x=self.star_grid['zz'])
        results = gaussian.fit(1.-ccf_broad, x=self.star_grid['zz'],  params=guess)
        rv_unperturbed = results.params['center'].value * 1000.

        #parameters, _ = curve_fit(CCF_gauss, self.star_grid['zz'], ccf_broad, p0=p0, check_finite =False)
        #rv_unperturbed = parameters[1] * 1000.
        #p0 = (parameters[0], 0.00, parameters[2])

        for i_obs, bjd_value in enumerate(bjd):

            n_oversampling = int(exptime[i_obs] / self.star_grid['time_step'])

            """recomputing the oversampling steps to homogeneously cover the
            full integration time """
            if n_oversampling % 2 == 0:
                n_oversampling += 1

            half_time = exptime[i_obs] / 2 / 86400.

            bjd_oversampling = np.linspace(bjd_value - half_time, bjd_value + half_time, n_oversampling, dtype=np.single)

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

            planet_position_xp = -orbital_distance_ratio * (np.cos(omega_rad + true_anomaly, dtype=np.single))
            planet_position_yp = -orbital_distance_ratio * (np.sin(omega_rad + true_anomaly, dtype=np.single) * np.cos(inclination_rad, dtype=np.single))
            planet_position_zp = orbital_distance_ratio * (np.sin(inclination_rad, dtype=np.single) * np.sin(omega_rad + true_anomaly, dtype=np.single))

            # projected distance of the planet's center to the stellar center
            planet_position_rp = np.sqrt(planet_position_xp**2  + planet_position_yp**2, dtype=np.single)
            ccf_out = np.zeros_like(ccf_total, dtype=np.single)

            for j, zeta in enumerate(planet_position_zp):

                if zeta > 0 and planet_position_rp[j] < 1. + parameter_values['R_Rs']:
                    # the planet is in the foreground or inside the stellar disk, continue
                    # adjustment: computation is performed even if only part of the planet is shadowing the star

                    rd = np.sqrt((planet_position_xp[j] - self.star_grid['xc']) ** 2 +
                                    (planet_position_yp[j] - self.star_grid['yc']) ** 2, dtype=np.single)

                    """ Selection of the portion of stars covered by the planet"""
                    sel_eclipsed = (rd <= parameter_values['R_Rs']) & self.star_grid['inside']

                    ccf_out += ccf_total - np.sum(star_grid_ccf[sel_eclipsed,:], axis=0, dtype=np.single)

            cont_val = np.amax(ccf_out)

            if cont_val > 0:
                ccf_out /= cont_val

                ccf_broad = gaussian_filter1d(ccf_out, self.ccf_variables['instrumental_broadening']/self.star_grid['rv_step'])
                try:
                    #parameters, _ = curve_fit(CCF_gauss, self.star_grid['zz'], ccf_broad, p0=p0, check_finite =False)
                    #rv_rml[i_obs] = parameters[1] * 1000. - rv_unperturbed
                    gaussian = GaussianModel()
        
                    params = gaussian.make_params()
                    guess = gaussian.guess(1.-ccf_broad, x=self.star_grid['zz'])
                    results = gaussian.fit(1.-ccf_broad, x=self.star_grid['zz'],  params=guess)
                    rv_rml[i_obs] = results.params['center'].value * 1000. - rv_unperturbed
                
                except:
                    rv_rml[i_obs] = 0.00
        return rv_rml

