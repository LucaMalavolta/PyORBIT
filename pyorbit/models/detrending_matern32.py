from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
import pyorbit.subroutines.constants as constants
from numpy.polynomial import polynomial
from numpy.lib.recfunctions import append_fields
from scipy.interpolate import interp1d

class Detrending_Matern32(AbstractModel):

    default_common = 'detrending'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        try:
            from george import kernels, GP
            from george.metrics import Metric
        except:
            print("ERROR: celerite2 not installed, this will not work")
            quit()

        import warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning)


        self.model_class = 'detrending_celerite2_matern32'
        self.internal_likelihood = True

        self.local_model = True

        self.list_pams_common = set()
        self.list_pams_dataset = {
            'det_m32_sigma',  # sigma
        }
        self.interpolated = {}
        self.gp_ndim = None
        self.gp_rvector = {}
        self.gp_metric_index = {}

        self.gp_kernel = {}
        self.gp = {}

    def initialize_model(self, mc, **kwargs):

        self.local_model = kwargs.get('local_model', self.local_model)

        self.common_model = ~self.local_model
        for keyword in ['use_common_parameters', 'common_parameters', 'use_common_model', 'common_model']:
            self.common_model = kwargs.get(keyword, self.common_model)
        self.local_model = ~self.common_model

        if self.common_model:
            pam = 'det_m32_sigma'
            self.list_pams_common.update([pam])
            self.list_pams_dataset.discard(pam)

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == self.default_common:
                self.common_detrending = common_ref
                break

        self.gp_ndim= len(kwargs['detrending_variables'])


    def initialize_model_dataset(self, mc, dataset, **kwargs):

        self.gp_rvector[dataset.name_ref] = np.empty([dataset.n, self.gp_ndim])
        self.gp_metric_index[dataset.name_ref] = {}

        """ For plotting purposes only """
        self.interpolated[dataset.name_ref] = {}

        """ starting from one as the first variable is the amplitude of the GP"""
        index = 1

        for data_name in kwargs['detrending_variables']:

            if data_name not in dataset.ancillary.dtype.names:
                dataset.ancillary = append_fields(dataset.ancillary, data_name, np.zeros(dataset.n))

            par_original = 'det_m32_rho'
            par_addition = 'det_' + data_name + '_m32_rho'

            if self.local_model:
                self.transfer_parameter_properties(mc, dataset, par_original, par_addition, keywords=kwargs, dataset_pam=True)
            else:
                self.transfer_parameter_properties(mc, dataset, par_original, par_addition, keywords=kwargs, common_pam=True)

            self.gp_metric_index[dataset.name_ref][par_addition] = index
            self.gp_rvector[dataset.name_ref][:,index] = dataset.ancillary['data_name']

            index += 1

            """ For plotting purposes only """
            self.interpolated[dataset.name_ref][data_name]=interp1d(
            dataset.x0,
            dataset.ancillary[data_name],
            bounds_error=False,
            fill_value=(np.amin(dataset.ancillary[data_name]), np.amax(dataset.ancillary[data_name])))

        self.define_kernel(dataset)

    def define_kernel(self, dataset):
        # random initialization

        metric = [0.]*self.gp_ndim
        #self.gp_kernel[dataset.name_ref] = 1.0 * kernels.Matern32Kernel(metric=metric, ndim=self.gp_ndim)
        kernel = 1.0 * kernels.Matern32Kernel(metric=metric, ndim=self.gp_ndim)
        self.gp[dataset.name_ref] = GP(kernel)

        """ I've decided to add the jitter in quadrature instead of using a constant kernel to allow the use of
        different / selective jitter within the dataset
        """
        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)

        self.gp[dataset.name_ref].compute(dataset.x0, env)

        return


    def lnlk_compute(self, parameter_values, dataset):

        gp_pams = [0.]*self.gp_ndim

        gp_pams[0] = np.log(parameter_values['det_m32_sigma']/self.gp_ndim)
        for data_name, pam_index in self.gp_metric_index[dataset.name_ref].items():
            par_name = 'det_' + data_name + '_m32_rho'
            gp_pams[pam_index] = np.log(parameter_values[par_name]**2)

        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].compute(dataset.x0, env)
        return self.gp[dataset.name_ref].log_likelihood(dataset.residuals, quiet=True)

    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        gp_pams = [0.]*self.gp_ndim

        gp_pams[0] = np.log(parameter_values['det_m32_sigma']/self.gp_ndim)
        for data_name, pam_index in self.gp_metric_index[dataset.name_ref].items():
            par_name = 'det_' + data_name + '_m32_rho'
            gp_pams[pam_index] = np.log(gp_pams[par_name]**2)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)
        if x0_input is None:
            return self.gp[dataset.name_ref].predict(dataset.residuals, dataset.x0, return_cov=return_covariance, return_var=return_variance)
        else:
            return self.gp[dataset.name_ref].predict(dataset.residuals, x0_input, return_cov=return_covariance, return_var=return_variance)

    def sample_conditional(self, parameter_values, dataset, x0_input=None):

        gp_pams = [0.]*self.gp_ndim

        gp_pams[0] = np.log(parameter_values['det_m32_sigma']/self.gp_ndim)
        for data_name, pam_index in self.gp_metric_index[dataset.name_ref].items():
            par_name = 'det_' + data_name + '_m32_rho'
            gp_pams[pam_index] = np.log(gp_pams[par_name]**2)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(dataset.x0, env)
        if x0_input is None:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, dataset.x0)
        else:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, x0_input)
