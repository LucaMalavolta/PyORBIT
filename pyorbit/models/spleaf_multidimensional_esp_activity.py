from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_gaussian_processes import AbstractGaussianProcesses
from pyorbit.keywords_definitions import *

from scipy.linalg import cho_factor, cho_solve, lapack, LinAlgError
from scipy import matrix, spatial
import sys

__all__ = ['SPLEAF_Multidimensional_ESP']


try:
    from spleaf import cov as spleaf_cov
    from spleaf import term as spleaf_term
except (ModuleNotFoundError,ImportError):
    pass


class SPLEAF_Multidimensional_ESP(AbstractModel, AbstractGaussianProcesses):
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

        self.model_class = 'multidimensional_gaussian_process'

        self.internal_likelihood = True
        self.delayed_lnlk_computation = True

        self.list_pams_common = OrderedSet([
            'Prot',  # Rotational period of the star
            'Pdec',  # Decay timescale of activity
            'Oamp',  # Granulation of activity
        ])
        self.list_pams_dataset = OrderedSet([
            'rot_amp', # Amplitude of the covariance matrix
            'con_amp' # Amplitude of the first derivative of the covariance matrix
        ])

        try:
            from spleaf import cov as spleaf_cov
            from spleaf import term as spleaf_term
        except (ModuleNotFoundError,ImportError):
            print("ERROR: S+LEAF package not installed, this will not work")
            quit()


        self.internal_parameter_values = None

        self._dataset_x0 = []
        self._dataset_label = []
        self._dataset_err = []
        self._dataset_names = {}

        self._jitter_mask = []
        self._dataset_jitmask = []

        self._dataset_nindex = {}
        self._dataset_njitter = {}

        self._dataset_temporary_jitmask = {}

        self.internal_coeff_prime = None
        self.internal_coeff_deriv = None

        self._dataset_ej2 = None
        self._dataset_res = None

        self._added_datasets = 0
        self._added_jitters = 0
        self.n_harmonics = 4

        self.D_spleaf = None
        self.D_param = None
        self.input_param = None

        ### NEW Addded in PyORBIT v10.10
        self._dataset_rvflag_dict = {}

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
        if dataset.name_ref in self._dataset_nindex:
            return

        ## NEW Addded in PyORBIT v10.10
        if dataset.kind == 'RV':
            self._dataset_rvflag_dict[dataset.name_ref] = True
        else:
            self._dataset_rvflag_dict[dataset.name_ref] = False

        self._dataset_x0.append(dataset.x0)
        self._dataset_err.append(dataset.e)
        self._dataset_nindex[dataset.name_ref] = self._added_datasets

        self._dataset_njitter[dataset.name_ref] = []


        self._dataset_temporary_jitmask[dataset.name_ref] = []

        for var in dataset.list_pams:
            if dataset.variable_expanded[var] != 'jitter':
                continue
            self._dataset_temporary_jitmask[dataset.name_ref].append(dataset.mask[var])
            self._dataset_njitter[dataset.name_ref].append(self._added_jitters)
            self._added_jitters += 1

        self.spleaf_time, self.spleaf_res, self.spleaf_err, self.spleaf_series_index = \
            spleaf_cov.merge_series(self._dataset_x0, self._dataset_err, self._dataset_err)


        self._jitter_mask = []
        self._dataset_jitmask = []
        for dataset_name, temporary_jitmask in self._dataset_temporary_jitmask.items():
            d_ind = self._dataset_nindex[dataset_name]

            temporary_mask = self.spleaf_series_index[d_ind]
            for jit_list in temporary_jitmask:
                self._dataset_jitmask.append(jit_list)
                self._jitter_mask.append(temporary_mask[jit_list])

        self._added_datasets += 1

        d_ind = self._dataset_nindex[dataset.name_ref]

        self.spleaf_res[self.spleaf_series_index[d_ind]] = dataset.residuals

        self.internal_parameter_values = {
            'Prot': 30.0,
            'Pdec': 100.0,
            'Oamp': 0.35
        }

        self.internal_jitter = np.ones(self._added_jitters)
        self.internal_coeff_prime = np.ones(self._added_datasets)
        self.internal_coeff_deriv = np.ones(self._added_datasets)

        self._reset_kernel()

        self._set_derivative_option(mc, dataset, **kwargs) 

        return

    def add_internal_dataset(self, parameter_values, dataset):

        self.update_parameter_values(parameter_values)
        self.internal_parameter_values = parameter_values

        d_ind = self._dataset_nindex[dataset.name_ref]
        j_ind = self._dataset_njitter[dataset.name_ref]

        self.spleaf_res[self.spleaf_series_index[d_ind]] = dataset.residuals

        for j_jit in j_ind:
            self.internal_jitter[j_jit] =  dataset.jitter[self._dataset_jitmask[j_jit]][0]

        self.internal_coeff_prime[d_ind] = parameter_values['con_amp']
        self.internal_coeff_deriv[d_ind] = parameter_values['rot_amp']


    def lnlk_compute(self):

        pass_conditions = self.check_hyperparameter_values(self.internal_parameter_values)
        if not np.isfinite(pass_conditions):
            return pass_conditions

        """
        Randomly reset the kernel with a probability of 0.1%
        To prevent memory allocations issues I suspect are happening
        """
        random_selector = np.random.randint(1000)
        if random_selector == 50:
            self._reset_kernel()


        input_param = np.concatenate(([self.internal_parameter_values['Prot'],
                        self.internal_parameter_values['Pdec'],
                        self.internal_parameter_values['Oamp']],
                        self.internal_coeff_prime,
                        self.internal_coeff_deriv,
                        self.internal_jitter))
        self.D_spleaf.set_param(input_param, self.D_param)

        return self.D_spleaf.loglike(self.spleaf_res)


    ## NEW Addded in PyORBIT v10.10
    def lnlk_rvonly_compute(self):
        """
        Randomly reset the kernel with a probability of 0.1%
        To prevent memory allocations issues I suspect are happening
        """
        random_selector = np.random.randint(1000)
        if random_selector == 50:
            self._reset_kernel()


        input_param = np.concatenate(([self.internal_parameter_values['Prot'],
                        self.internal_parameter_values['Pdec'],
                        self.internal_parameter_values['Oamp']],
                        self.internal_coeff_prime,
                        self.internal_coeff_deriv,
                        self.internal_jitter))

        self.D_spleaf.set_param(input_param, self.D_param)
        #self.D_spleaf.kernel['GP'].set_conditional_coef(series_id=None)
        #gp_model  = self.D_spleaf.conditional(self.spleaf_res, self.spleaf_time)
        #residuals = (self.spleaf_res - gp_model)
        env = np.zeros_like(self.spleaf_err)

        rv_loglike = 0.0

        for dataset_name, d_ind in self._dataset_nindex.items():
            if not self._dataset_rvflag_dict[dataset_name]:
                continue

            j_ind = self._dataset_njitter[dataset_name]
            for j_jit in j_ind:
                env[self._jitter_mask[j_jit]] = 1. / (self.spleaf_err[self._jitter_mask[j_jit]]**2 +  self.internal_jitter[j_jit] **2)

            dataset_selection = self.spleaf_series_index[d_ind]
            n = len(dataset_selection)

            self.D_spleaf.kernel['GP'].set_conditional_coef(series_id=d_ind)
            gp_model = self.D_spleaf.conditional(self.spleaf_res, self.spleaf_time[dataset_selection])

            residuals = (self.spleaf_res[dataset_selection] - gp_model)
            rv_loglike += -0.5 * (n * np.log(2 * np.pi) + np.sum(residuals** 2 * env[dataset_selection] - np.log(env[dataset_selection])))

            #rv_loglike += -0.5 * (n * np.log(2 * np.pi) + np.sum(residuals[dataset_selection] ** 2 * env[dataset_selection] - np.log(env[dataset_selection])))

        return rv_loglike



    def sample_predict(self, dataset, x0_input=None, return_covariance=False, return_variance=False):

        input_param = np.concatenate(([self.internal_parameter_values['Prot'],
                        self.internal_parameter_values['Pdec'],
                        self.internal_parameter_values['Oamp']],
                        self.internal_coeff_prime,
                        self.internal_coeff_deriv,
                        self.internal_jitter))

        self.D_spleaf.set_param(input_param, self.D_param)
        """ I'm creating the kernel here has """

        d_ind = self._dataset_nindex[dataset.name_ref]

        self.D_spleaf.kernel['GP'].set_conditional_coef(series_id=d_ind)


        if x0_input is None:
            t_predict = dataset.x0
            sorting_predict = np.argsort(dataset.x0)
        else:
            t_predict = x0_input
            sorting_predict = np.argsort(x0_input)

        mu = np.empty_like(t_predict)
        var = np.empty_like(t_predict)
        mu[sorting_predict], var[sorting_predict]  = self.D_spleaf.conditional(self.spleaf_res, t_predict[sorting_predict], calc_cov='diag')

        if return_variance:
            return mu, np.sqrt(var)
        else:
            return mu

    def _reset_kernel(self):

        kwargs = {
            'err': spleaf_term.Error(self.spleaf_err),
            'GP': spleaf_term.MultiSeriesKernel(spleaf_term.ESPKernel(1.,
                                                                    self.internal_parameter_values['Prot'],
                                                                    self.internal_parameter_values['Pdec'],
                                                                    self.internal_parameter_values['Oamp'],
                                                                    nharm=self.n_harmonics),
                                            self.spleaf_series_index,
                                            self.internal_coeff_prime,
                                            self.internal_coeff_deriv)
        }
        for n_jit in range(0, self._added_jitters):
            kwargs['jitter_'+repr(n_jit)] = spleaf_term.InstrumentJitter(self._jitter_mask[n_jit], self.internal_jitter[n_jit])

        self.D_spleaf = spleaf_cov.Cov(self.spleaf_time, **kwargs)
        self.D_param = self.D_spleaf.param[1:]
