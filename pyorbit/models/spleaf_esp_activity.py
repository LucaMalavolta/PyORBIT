from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
from pyorbit.keywords_definitions import *

from scipy.linalg import cho_factor, cho_solve, lapack, LinAlgError
from scipy import matrix, spatial
import sys

__all__ = ['SPLEAF_ESP']


try:
    from spleaf import cov as spleaf_cov
    from spleaf import term as spleaf_term
except (ModuleNotFoundError,ImportError):
    pass


class SPLEAF_ESP(AbstractModel):
    ''' Three parameters out of four are the same for all the datasets, since they are related to
    the properties of the physical process rather than the observed effects on a dataset
     From Grunblatt+2015, Affer+2016
     - theta: is usually related to the rotation period of the star( or one of its harmonics);
     - lambda: is the correlation decay timescale, and it can be related to the lifetime of the active regions.
     - omega: is the length scale of the periodic component, and can be linked to the size evolution of the active regions;
     - h: represents the amplitude of the correlations '''

    default_common = 'activity'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'spleaf_esp'

        self.internal_likelihood = True

        self.list_pams_common = OrderedSet([
            'Prot',  # Rotational period of the star
            'Pdec',  # Decay timescale of activity
            'Oamp',  # Granulation of activity
        ])
        self.list_pams_dataset = OrderedSet([
            'Hamp'  # Amplitude of the signal in the covariance matrix
        ])

        try:
            from spleaf import cov as spleaf_cov
            from spleaf import term as spleaf_term
        except (ModuleNotFoundError,ImportError):
            print("ERROR: S+LEAF package not installed, this will not work")
            quit()

        self.n_harmonics = 4


    def initialize_model(self, mc,  **kwargs):

        self.n_harmonics = kwargs.get('n_harmonics', self.n_harmonics)
        print(' S+LEAF model, number of harmonics:', self.n_harmonics)
        print()

        if kwargs.get('hyperparameters_condition', False):
            self.hyper_condition = self._hypercond_01
        else:
            self.hyper_condition = self._hypercond_00

        if kwargs.get('rotation_decay_condition', False):
            self.rotdec_condition = self._hypercond_02
        else:
            self.rotdec_condition = self._hypercond_00

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'activity':
                self.use_stellar_rotation_period = getattr(mc.common_models[common_ref], 'use_stellar_rotation_period', False)
                break

        for keyword in keywords_stellar_rotation:
            self.use_stellar_rotation_period = kwargs.get(keyword, self.use_stellar_rotation_period)

        if self.use_stellar_rotation_period:
            self.list_pams_common.update(['rotation_period'])
            self.list_pams_common.discard('Prot')


    def lnlk_compute(self, parameter_values, dataset):

        if self.use_stellar_rotation_period:
            parameter_values['Prot'] = parameter_values['rotation_period']

        if not self.hyper_condition(parameter_values):
            return -np.inf
        if not self.rotdec_condition(parameter_values):
            return -np.inf

        """ I'm creating the kernel here has """
        D = spleaf_cov.Cov(dataset.x0,
            err=spleaf_term.Error(dataset.e ** 2.0 + dataset.jitter ** 2.0),
            GP=spleaf_term.ESPKernel(parameter_values['Hamp'],
                                    parameter_values['Prot'],
                                    parameter_values['Pdec'],
                                    parameter_values['Oamp'],
                                    nharm=self.n_harmonics))



        return D.loglike(dataset.residuals)


    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        if self.use_stellar_rotation_period:
            parameter_values['Prot'] = parameter_values['rotation_period']

        """ I'm creating the kernel here has """
        D = spleaf_cov.Cov(dataset.x0,
            err=spleaf_term.Error(dataset.e ** 2.0 + dataset.jitter ** 2.0),
            GP=spleaf_term.ESPKernel(parameter_values['Hamp'],
                                    parameter_values['Prot'],
                                    parameter_values['Pdec'],
                                    parameter_values['Oamp'],
                                    nharm=self.n_harmonics))

        if x0_input is None:
            t_predict = dataset.x0
        else:
            t_predict = x0_input

        mu, var = D.conditional(dataset.residuals, t_predict, calc_cov='diag')

        if return_variance:
            return mu, np.sqrt(var)
        else:
            return mu

    @staticmethod
    def _hypercond_00(parameter_values):
        #Condition from Rajpaul 2017, Rajpaul+2021
        return True

    @staticmethod
    def _hypercond_01(parameter_values):
        # Condition from Rajpaul 2017, Rajpaul+2021
        # Taking into account that Pdec^2 = 2*lambda_2^2
        return parameter_values['Pdec']**2 > (3. / 2. / np.pi) * parameter_values['Oamp']**2 * parameter_values['Prot']**2

    @staticmethod
    def _hypercond_02(parameter_values):
        #Condition on Rotation period and decay timescale
        return parameter_values['Pdec'] > 2. * parameter_values['Prot']
