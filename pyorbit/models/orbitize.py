
from os import system

from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel
import pyorbit.subroutines.kepler_exo as kepler_exo
import pyorbit.subroutines.transformations as transformations
import pyorbit.subroutines.constants as constants
import os   

try:
        from orbitize import DATADIR, hipparcos, gaia
        import orbitize.kepler
        import orbitize.lnlike
        from orbitize import read_input, system, priors
except ImportError:
        pass


try:
    from astropy.time import Time
except ImportError:
    pass

class AbstractAstrometry(object):

    def __init__(self, *args, **kwargs):

        ''' Orbital parameters to be used in the astrometric fit '''
        self.list_pams_common = OrderedSet([
            'P',     # Period in days
            'M_Ms',
            'Omega', # longitude of ascending node
            'e',     # eccentricity, uniform prior - to be fixed
            'i',
            'R_Rs',  # planet radius (in units of stellar radii)
            'omega', # argument of pericenter
            #'i',     # inclination in degrees
            'mass', #stellar mass
            'parallax', #stellar parallax
        ])

        self.list_pams_dataset = OrderedSet()
        self.warning_given = False
        self.compute_semimajor_axis = True
        self.compute_scaled_semimajor_axis = True

    def _prepare_astrometry_parameters(self, mc, **kwargs):


        print("*** {0:s} global parameters:".format(self.model_name))

        transformations._prepare_planet_parametrization(self, mc, **kwargs)
        transformations._prepare_planet_scaled_semimajor_axis(self, mc, **kwargs)
        transformations._prepare_planet_semimajor_axis(self, mc, **kwargs)
        transformations._prepare_planet_mass(self, mc, **kwargs)
        transformations._prepare_stellar_mass(self, mc, **kwargs)
        transformations._prepare_planet_inclination(self, mc, **kwargs)
        transformations._prepare_planet_time_inferior_conjunction(self, mc, **kwargs)

    def update_parameter_values(self,
                                parameter_values,
                                input_prepend=''):

        if input_prepend == '':
            prepend = ''
        else:
            prepend = input_prepend + '__'

        print(parameter_values)
        if self.compute_inclination:
            if self.compute_scaled_semimajor_axis:
                parameter_values[prepend+'a_Rs'] = convert_rho_to_ars(parameter_values[prepend+'P'], parameter_values['density'])
            parameter_values[prepend+'i'] = convert_b_to_i(
                parameter_values[prepend+'b'], parameter_values[prepend+'e'], parameter_values[prepend+'omega'], parameter_values[prepend+'a_Rs'])

        if self.compute_semimajor_axis:
            if self.compute_scaled_semimajor_axis:
                parameter_values[prepend+'a_Rs'] = convert_rho_to_ars(parameter_values[prepend+'P'], parameter_values['density'])
            parameter_values[prepend+'a_AU'] = convert_ars_to_a(parameter_values[prepend+'a_Rs'], parameter_values['radius'])

        if self.compute_time_inferior_conjunction:
            parameter_values[prepend+'Tc']= kepler_exo.kepler_compute_deltaTc_from_meanlong(
                parameter_values[prepend+'P'],
                parameter_values[prepend+'mean_long'],
                parameter_values[prepend+'e'],
                parameter_values[prepend+'omega'],
                parameter_values[prepend+'Omega']) + self.Tref

        if self.compute_mean_longitude:
            parameter_values[prepend+'mean_long'] = kepler_exo.kepler_compute_meanlong_from_deltaTc(
                parameter_values[prepend+'P'],
                parameter_values[prepend+'Tc'] - self.Tref,
                parameter_values[prepend+'e'],
                parameter_values[prepend+'omega'],
                parameter_values[prepend+'Omega'])


class Orbitize(AbstractModel, AbstractAstrometry):
    model_class = 'orbitize'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
                from orbitize import DATADIR, hipparcos, gaia
        except ImportError:
                print("    {0:s} WARNING:".format(self.model_name))
                print('        orbitize is not installed. Please install it to use this model.')
                print('        You can install it with pip install orbitize')
                print()


        try:
            from astropy.time import Time
        except ImportError:
            print("    {0:s} WARNING:".format(self.model_name))
            print('        astropy is not installed. Please install it to use this model.')
            print('        You can install it with pip install astropy')
            print()


        print("    {0:s} WARNING:".format(self.model_name))
        print('        Astrometry modelling requires the use of the stellar mass')
        print('        This may cause a clash with models requiring stellar density and radius, e.g., RM modelling')
        print('        The use of a multivariate approach is strongly suggested')
        print('        You can control the behaviour of mass/radius/density with the specific keywords')
        print('        compute_mass, compute_radius, compute_density')
        print()

        self.external_dataset = True
        self.accept_multiple_planets = True

    # brainless workaround
    def initialize_model(self, mc, **kwargs):

        self.iad_filepath = os.path.join('./orbitize_data', "HIP2_HIP004311.d")
        self.gost_filepath = os.path.join('./orbitize_data', "GOST_HIP004311.csv")

        self.hipparcos_lnprob = hipparcos.HipparcosLogProb(self.iad_filepath, 4311, 1)
        self.hgca_lnprob = gaia.HGCALogProb(4311, self.hipparcos_lnprob, self.gost_filepath)

        self.Tref = mc.Tref

    def initialize_model_parameters(self, mc, **kwargs):

        self._prepare_astrometry_parameters(mc, **kwargs)

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        astrometry_filepath = os.path.join('./orbitize_data', "no_data.csv")
        self.data_table = read_input.read_file(astrometry_filepath)

    def compute(self, parameter_values, dataset, x0_input=None):
        print("    {0:s} WARNING:".format(self.model_name), self.parameter_values)
        print('        Astrometry modelling requires the use of the stellar mass')

    def compute_loglikelihood(self, planet_list, dataset):
        print("    {0:s} WARNING:".format(self.model_name))
        print('        loglikelihood not implemented')

        print(self.parameter_values)

        this_system = system.System(
            len(planet_list),
            self.data_table,
            self.parameter_values['mass'],
            self.parameter_values['parallax'],
            tau_ref_epoch=self.Tref - 2400000.5,
            fit_secondary_mass=True,
            gaia=self.hgca_lnprob,
        )

        n_param = len(this_system.labels)
        param_model = np.zeros(n_param) 

        #self.update_parameter_values_for_astrometry(self.parameter_values, self.Tref)

        for i_planet, planet_name in enumerate(planet_list):

            prepend = planet_name + '__'
            
            sau = 1.757
            esino = -0.2299
            ecoso = 0.3395
            inc = 38.7 / 180. * np.pi
            asc = 266.6 / 180. * np.pi
            omega = omega_deg / 180. * np.pi
            msec = 20.0 * constants.Mjups 


            delta_Tc =  kp.kepler_compute_deltaTc_from_deltaTperi(period, delta_Tperi, ecc, omega_deg)

            lam = kp.kepler_compute_meanlong_from_deltaTc(period, delta_Tc, ecc, omega_deg, Omega_deg=180.0) / 180. * np.pi 




            param_model[this_system.param_idx['sma'+repr(i_planet)]] = self.parameter_values[prepend+'a_AU']
            param_model[this_system.param_idx['ecc'+repr(i_planet)]] = self.parameter_values[prepend+'e']
            param_model[this_system.param_idx['inc'+repr(i_planet)]] = self.parameter_values[prepend+'i'] * constants.deg2rad
            param_model[this_system.param_idx['aop'+repr(i_planet)]] = self.parameter_values[prepend+'omega'] * constants.deg2rad
            param_model[this_system.param_idx['pan'+repr(i_planet)]] = self.parameter_values[prepend+'Omega'] * constants.deg2rad
            param_model[this_system.param_idx['tau'+repr(i_planet)]] = self.parameter_values[prepend+'tau']
            param_model[this_system.param_idx['ecc'+repr(i_planet)]] = ecc
            param_model[this_system.param_idx['inc1']] = inc
            param_model[this_system.param_idx['aop1']] = omega
            param_model[this_system.param_idx['pan1']] = asc
            param_model[this_system.param_idx['tau1']] = phase_Tperi
            param_model[this_system.param_idx['plx']] = plx
            param_model[this_system.param_idx['m1']] = msec
        param_model[this_system.param_idx['m0']] = mpri



        chi2_type = 'standard'
        custom_lnlike = None
        like="chi2_lnlike"
        # check if `like` is a string or a function
        if callable(like):
            lnlike = like
        else:
            lnlike = getattr(orbitize.lnlike, like)

        has_corr = np.any(~np.isnan(this_system.data_table["quant12_corr"]))

        # compute the model based on system params
        model, jitter = this_system.compute_model(param_model)

        # fold data/errors to match model output shape. In particualr, quant1/quant2 are interleaved
        data = np.array(
            [this_system.data_table["quant1"], this_system.data_table["quant2"]]
        ).T

        # errors below required for lnlike function below
        errs = np.array(
            [this_system.data_table["quant1_err"], this_system.data_table["quant2_err"]]
        ).T
        # covariances/correlations, if applicable
        # we're doing this check now because the likelihood computation is much faster if we can skip it.
        if has_corr:
            corrs = this_system.data_table["quant12_corr"]
        else:
            corrs = None

        # grab all seppa indices
        seppa_indices = this_system.all_seppa

        # compute lnlike
        lnlikes = lnlike(
            data, errs, corrs, model, jitter, seppa_indices, chi2_type=chi2_type
        )

        # return sum of lnlikes (aka product of likeliehoods)
        lnlikes_sum = np.nansum(lnlikes, axis=(0, 1))

        if custom_lnlike is not None:
            lnlikes_sum += custom_lnlike(param_model)

        if this_system.hipparcos_IAD is not None:
            # compute Ra/Dec predictions at the Hipparcos IAD epochs
            raoff_model, deoff_model, _ = this_system.compute_all_orbits(
                param_model, epochs=this_system.hipparcos_IAD.epochs_mjd
            )

            (
                raoff_model_hip_epoch,
                deoff_model_hip_epoch,
                _,
            ) = this_system.compute_all_orbits(
                param_model, epochs=Time([1991.25], format="decimalyear").mjd
            )

            # subtract off position of star at reference Hipparcos epoch
            raoff_model[:, 0, :] -= raoff_model_hip_epoch[:, 0, :]
            deoff_model[:, 0, :] -= deoff_model_hip_epoch[:, 0, :]

            # select body 0 raoff/deoff predictions & feed into Hip IAD lnlike fn
            lnlikes_sum += this_system.hipparcos_IAD.compute_lnlike(
                raoff_model[:, 0, :],
                deoff_model[:, 0, :],
                param_model,
                this_system.param_idx,
            )

        if this_system.gaia is not None:
            gaiahip_epochs = Time(
                np.append(
                    this_system.gaia.hipparcos_epoch, this_system.gaia.gaia_epoch
                ),
                format="decimalyear",
            ).mjd

            # compute Ra/Dec predictions at the Gaia epoch
            raoff_model, deoff_model, _ = this_system.compute_all_orbits(
                param_model, epochs=gaiahip_epochs
            )

            # select body 0 raoff/deoff predictions & feed into Gaia module lnlike fn
            lnlikes_sum += this_system.gaia.compute_lnlike(
                raoff_model[:, 0, :],
                deoff_model[:, 0, :],
                param_model,
                this_system.param_idx,
            )
        return lnlikes_sum
