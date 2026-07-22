
from os import system

from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_astrometry import AbstractAstrometry

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


class Orbitize(AbstractModel, AbstractAstrometry):
    model_class = 'orbitize'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        ''' Orbital parameters to be used in the astrometric fit '''
        self.list_pams_common = OrderedSet([
            'P',     # Period in days
            #'M_Me', # planet mass in Earth masses
            'Omega', # longitude of ascending node
            'e',     # eccentricity, uniform prior - to be fixed
            'i',
            'omega', # argument of pericenter
            #'i',     # inclination in degrees
            'mass', #stellar mass
            'parallax', #stellar parallax
        ])

        self.list_pams_dataset = OrderedSet()


        try:
                from orbitize import DATADIR, hipparcos, gaia
        except ImportError:
                print("    {0:s} WARNING:".format(self.model_name))
                print('        orbitize is not installed. Please install it to use this model.')
                print('        You can install it with pip install orbitize')
                print()
                quit()

        try:
            from astropy.time import Time
        except ImportError:
            print("    {0:s} WARNING:".format(self.model_name))
            print('        astropy is not installed. Please install it to use this model.')
            print('        You can install it with pip install astropy')
            print()
            quit()

        self.external_dataset = True
        self.accept_multiple_planets = True
        self.force_model_dataset_initialization = True

    def print_warning(self):
        print("{0:s} WARNING:".format(self.model_name))
        print('    Astrometric orbit modelling requires the use of the stellar mass')
        print('    This may cause a clash with models requiring stellar density and radius, e.g., RM modelling')
        print('    The use of a multivariate approach is strongly suggested')
        print('    You can control the behaviour of mass/radius/density with the specific keywords')
        print('    compute_mass, compute_radius, compute_density')
        print()

    def initialize_model(self, mc, **kwargs):


        if kwargs.get('iad_filepath', False):
            self.iad_filepath = kwargs.get('iad_filepath')
            self.hipparcos_lnprob = hipparcos.HipparcosLogProb(self.iad_filepath, 4311, 1)
        else:
            self.hipparcos_lnprob = None

        if kwargs.get('gost_filepath', False):
            if not kwargs.get('iad_filepath', False):
                print("UNRECOVERABLE ERROR model {0:s} :".format(self.model_name))
                print('To use HGCA fit, Hipparcos IAD data is required. Please provide the iad_filepath keyword')
                raise ValueError()

            hipparcos_ID = kwargs.get('hipparcos_ID', False)
            if not hipparcos_ID:
                print("UNRECOVERABLE ERROR model {0:s} :".format(self.model_name))
                print('To use HGCA fit, Hipparcos ID is required. Please provide the hipparcos_ID keyword')
                raise ValueError()

            self.gost_filepath = kwargs.get('gost_filepath')
            self.hgca_lnprob = gaia.HGCALogProb(4311, self.hipparcos_lnprob, self.gost_filepath)
            self.hipparcos_lnprob = None
        else:
            self.hgca_lnprob = None

        self.Tref = mc.Tref

    def initialize_model_parameters(self, mc, **kwargs):

        self._prepare_astrometry_parameters(mc, **kwargs)

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        self.data_table = read_input.read_file(dataset.input_file)

    def compute(self, parameter_values, dataset, x0_input=None):
        pass

    def compute_loglikelihood(self, planet_list, dataset):

        orbitize_tref = 58849.0 # MJD

        this_system = system.System(
            len(planet_list),
            self.data_table,
            self.parameter_values['mass'],
            self.parameter_values['parallax'],
            tau_ref_epoch=orbitize_tref,
            fit_secondary_mass=True,
            hipparcos_IAD = self.hipparcos_lnprob,
            gaia=self.hgca_lnprob,
        )

        n_param = len(this_system.labels)
        param_model = np.zeros(n_param) 

        #self.update_parameter_values_for_astrometry(self.parameter_values, self.Tref)

        for i0_planet, planet_name in enumerate(planet_list):
            i_planet = i0_planet + 1 
            prepend = planet_name + '__'

            Tperi_MJD =  self.parameter_values[prepend+'Tperi'] + self.Tref - 2400000.5
            tau = orbitize.basis.tp_to_tau(Tperi_MJD, orbitize_tref, self.parameter_values[prepend+'P'])

            param_model[this_system.param_idx['sma'+repr(i_planet)]] = self.parameter_values[prepend+'a_AU']
            param_model[this_system.param_idx['ecc'+repr(i_planet)]] = self.parameter_values[prepend+'e']
            param_model[this_system.param_idx['inc'+repr(i_planet)]] = self.parameter_values[prepend+'i'] * constants.deg2rad
            param_model[this_system.param_idx['aop'+repr(i_planet)]] = self.parameter_values[prepend+'omega'] * constants.deg2rad
            param_model[this_system.param_idx['pan'+repr(i_planet)]] = self.parameter_values[prepend+'Omega'] * constants.deg2rad
            param_model[this_system.param_idx['tau'+repr(i_planet)]] = tau
            param_model[this_system.param_idx['m'+repr(i_planet)]] = self.parameter_values[prepend+'M_Me'] * constants.Mears

        param_model[this_system.param_idx['plx']] = self.parameter_values['parallax']
        param_model[this_system.param_idx['m0']] =  self.parameter_values['mass']

        #print("param_model = ", param_model)
        #print("this_system.param_idx = ", this_system.param_idx)

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
