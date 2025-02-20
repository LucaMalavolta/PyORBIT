from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_gaussian_processes import AbstractGaussianProcesses
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


class SPLEAF_ESP(AbstractModel, AbstractGaussianProcesses):
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

        self._jitter_mask = {}
        self._sorting_mask = {}
        self._n_jitter = {}
        self.D_spleaf = {}

    def initialize_model(self, mc,  **kwargs):

        self.n_harmonics = kwargs.get('n_harmonics', self.n_harmonics)
        print(self.model_name,  ' S+LEAF model, number of harmonics:', self.n_harmonics)
        print()

        self._prepare_hyperparameter_conditions(mc, **kwargs)
        self._prepare_rotation_replacement(mc, **kwargs)
        self._prepare_decay_replacement(mc, **kwargs)

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        """ when reloading the .p files, the object is not reinitialized, so we have to skip the
        incremental addition of datasets if they are already present  """
        if dataset.name_ref in self._n_jitter:
            return

        self._sorting_mask[dataset.name_ref] = np.argsort(dataset.x0)
        temp_sorting = self._sorting_mask[dataset.name_ref]

        self._jitter_mask[dataset.name_ref] = []
        self._n_jitter[dataset.name_ref] = 0

        temp_mask = np.arange(0, dataset.n, 1, dtype=int)
        for var in dataset.list_pams:
            if dataset.variable_expanded[var] != 'jitter':
                continue
            temp_jitmask = dataset.mask[var][temp_sorting]

            self._jitter_mask[dataset.name_ref].append(temp_mask[temp_jitmask])
            self._n_jitter[dataset.name_ref] += 1

        parameter_values = {
            'Hamp': 10.0,
            'Prot': 30.0,
            'Pdec': 100.0,
            'Oamp': 0.35
        }

        self._reset_kernel(parameter_values, dataset, temp_sorting)


    def lnlk_compute(self, parameter_values, dataset):

        self.update_parameter_values(parameter_values)
        
        pass_conditions = self.check_hyperparameter_values(parameter_values)
        if not pass_conditions:
            return pass_conditions

        try:
            temp_sorting = self._sorting_mask[dataset.name_ref]
        except:
            temp_sorting = np.argsort(dataset.x0)

        """
        Randomly reset the kernel with a probability of 0.1%
        To prevent memory allocations issues I suspect are happening
        """
        random_selector = np.random.randint(1000)
        if random_selector == 50:
            try:
                self._reset_kernel(parameter_values, dataset, temp_sorting)
            except TypeError:
                self._reset_kernel(parameter_values, dataset)


        jitter_values = np.zeros(self._n_jitter[dataset.name_ref])
        for n_jit in range(0, self._n_jitter[dataset.name_ref]):
            jitter_values[n_jit] = dataset.jitter[self._jitter_mask[dataset.name_ref][n_jit]][0]

        input_param = np.concatenate(([parameter_values['Hamp'],
                                    parameter_values['Prot'],
                                    parameter_values['Pdec'],
                                    parameter_values['Oamp']],
                                    jitter_values))


        self.D_spleaf[dataset.name_ref].set_param(input_param, self.D_spleaf[dataset.name_ref].param)
        return  self.D_spleaf[dataset.name_ref].loglike(dataset.residuals[temp_sorting])


    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):
        
        self.update_parameter_values(parameter_values)
        
        if x0_input is None:
            t_predict = dataset.x0
            sorting_predict = np.argsort(dataset.x0)
        else:
            t_predict = x0_input
            sorting_predict = np.argsort(x0_input)

        temp_sorting = np.argsort(dataset.x0)

        """ I'm creating the kernel here """
        D = spleaf_cov.Cov(dataset.x0[temp_sorting],
            err=spleaf_term.Error(np.sqrt(dataset.e[temp_sorting] ** 2.0 + dataset.jitter[temp_sorting] ** 2.0)),
            GP=spleaf_term.ESPKernel(parameter_values['Hamp'],
                                    parameter_values['Prot'],
                                    parameter_values['Pdec'],
                                    parameter_values['Oamp'],
                                    nharm=self.n_harmonics))


        mu = np.empty_like(t_predict)
        var = np.empty_like(t_predict)
        mu[sorting_predict], var[sorting_predict] = D.conditional(dataset.residuals[temp_sorting], t_predict[sorting_predict], calc_cov='diag')

        if return_variance:
            return mu, np.sqrt(var)
        else:
            return mu

    def _reset_kernel(self, parameter_values, dataset, argsorting=None):

        if argsorting is None:
            argsorting = np.argsort(dataset.x0)

        kwargs = {
            'err': spleaf_term.Error(dataset.e[argsorting]),
            'GP': spleaf_term.ESPKernel(parameter_values['Hamp'],
                                    parameter_values['Prot'],
                                    parameter_values['Pdec'],
                                    parameter_values['Oamp'],
                                    nharm=self.n_harmonics)
        }
        for n_jit in range(0, self._n_jitter[dataset.name_ref]):
            kwargs['jitter_'+repr(n_jit)] = spleaf_term.InstrumentJitter(self._jitter_mask[dataset.name_ref][n_jit],
                                                                            dataset.jitter[self._jitter_mask[dataset.name_ref][n_jit]][0])

        self.D_spleaf[dataset.name_ref] = spleaf_cov.Cov(dataset.x0[argsorting], **kwargs)
