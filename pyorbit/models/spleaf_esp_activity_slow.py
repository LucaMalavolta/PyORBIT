from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_gaussian_processes import AbstractGaussianProcesses
from pyorbit.keywords_definitions import *

from scipy.linalg import cho_factor, cho_solve, lapack, LinAlgError
from scipy import matrix, spatial
import sys

__all__ = ['SPLEAF_ESP_slow']


try:
    from spleaf import cov as spleaf_cov
    from spleaf import term as spleaf_term
except (ModuleNotFoundError,ImportError):
    pass


class SPLEAF_ESP_slow(AbstractModel, AbstractGaussianProcesses):
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
        super(AbstractModel, self).__init__(*args, **kwargs)

        self.model_class = 'gaussian_process'

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
        print(self.model_name,  ' S+LEAF model, number of harmonics:', self.n_harmonics)
        print()

        self._prepare_hyperparameter_conditions(mc, **kwargs)
        self._prepare_rotation_replacement(mc, **kwargs)
        self._prepare_decay_replacement(mc, **kwargs)

    def lnlk_compute(self, parameter_values, dataset):

        self.update_parameter_values(parameter_values)
        
        pass_conditions = self.check_hyperparameter_values(parameter_values)
        if not pass_conditions:
            return pass_conditions

        """ I'm creating the kernel here has """
        D = spleaf_cov.Cov(dataset.x0,
            err=spleaf_term.Error(np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)),
            GP=spleaf_term.ESPKernel(parameter_values['Hamp'],
                                    parameter_values['Prot'],
                                    parameter_values['Pdec'],
                                    parameter_values['Oamp'],
                                    nharm=self.n_harmonics))



        return D.loglike(dataset.residuals)


    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        self.update_parameter_values(parameter_values)

        """ I'm creating the kernel here has """
        D = spleaf_cov.Cov(dataset.x0,
            err=spleaf_term.Error(np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)),
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
