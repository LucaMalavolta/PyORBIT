from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
import pyorbit.subroutines.constants as constants
from numpy.polynomial import polynomial
from numpy.lib.recfunctions import append_fields
from scipy.interpolate import interp1d
from pyorbit.keywords_definitions import *

try:
    import george
except:
    pass

class Detrending_Matern32(AbstractModel):

    default_common = 'detrending'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        try:
            import george
        except (ModuleNotFoundError,ImportError):
            print("ERROR: george not installed, this will not work")
            quit()

        import warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning)


        self.model_class = 'detrending_celerite2_matern32'
        self.internal_likelihood = True

        self.local_model = True

        self.list_pams_common = OrderedSet()
        self.list_pams_dataset = OrderedSet([
            'det_m32_sigma',  # sigma
        ])
        self.interpolated = {}
        self.gp_ndim = None
        self.gp_rvector = {}
        self.gp_metric_index = {}

        self.standardize = False

        self.gp_kernel = {}
        self.gp = {}

    def initialize_model(self, mc, **kwargs):

        self.local_model = kwargs.get('local_model', self.local_model)
        self.standardize = kwargs.get('standardize', self.standardize)

        self.common_model = ~self.local_model
        for keyword in keywords_detrending_common_parameters:
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

        index = 0

        for data_name in kwargs['detrending_variables']:

            if data_name not in dataset.ancillary.dtype.names:
                self.gp_rvector[dataset.name_ref][:,index] = np.zeros(dataset.n)
            else:

                if self.standardize:
                    mean = np.mean(dataset.ancillary[data_name])
                    std = np.std(dataset.ancillary[data_name])
                    self.gp_rvector[dataset.name_ref][:,index] = (dataset.ancillary[data_name]-mean)/std
                    print('  {0:s} {1:s} mean:{2:f} std:{3:f}'.format(dataset.name_ref, data_name, mean, std))
                else:
                    self.gp_rvector[dataset.name_ref][:,index] = dataset.ancillary[data_name]

            par_original = 'det_m32_rho'
            par_addition = 'det_' + data_name + '_m32_rho'

            if self.local_model:
                self.transfer_parameter_properties(mc, dataset, par_original, par_addition, keywords=kwargs, dataset_pam=True)
            else:
                self.transfer_parameter_properties(mc, dataset, par_original, par_addition, keywords=kwargs, common_pam=True)

            self.gp_metric_index[dataset.name_ref][data_name] = index

            self.interpolated[dataset.name_ref][data_name]=interp1d(
            dataset.x0,
            self.gp_rvector[dataset.name_ref][:,index],
            bounds_error=False,
            fill_value=(np.amin(self.gp_rvector[dataset.name_ref][:,index]), np.amax(self.gp_rvector[dataset.name_ref][:,index])))

            index += 1

        self.define_kernel(dataset)

    def define_kernel(self, dataset):
        # random initialization

        metric = [1.]*self.gp_ndim
        #self.gp_kernel[dataset.name_ref] = 1.0 * kernels.Matern32Kernel(metric=metric, ndim=self.gp_ndim)
        kernel = 0.1 * george.kernels.Matern32Kernel(metric=metric, ndim=self.gp_ndim)
        self.gp[dataset.name_ref] = george.GP(kernel)

        """ I've decided to add the jitter in quadrature instead of using a constant kernel to allow the use of
        different / selective jitter within the dataset
        """
        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)

        self.gp[dataset.name_ref].compute(self.gp_rvector[dataset.name_ref], env)

        return


    def lnlk_compute(self, parameter_values, dataset):

        gp_pams = [0.]*(self.gp_ndim+1)

        gp_pams[0] = np.log(parameter_values['det_m32_sigma']**2/self.gp_ndim)
        for data_name, pam_index in self.gp_metric_index[dataset.name_ref].items():
            par_name = 'det_' + data_name + '_m32_rho'
            gp_pams[pam_index+1] = np.log(parameter_values[par_name]**2)

        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].compute(self.gp_rvector[dataset.name_ref], env)
        return self.gp[dataset.name_ref].log_likelihood(dataset.residuals, quiet=True)

    def sample_predict(self, parameter_values, dataset, x0_input=None, return_covariance=False, return_variance=False):

        gp_pams = [0.]*(self.gp_ndim+1)

        if x0_input is None:
            gp_rvector = np.empty([dataset.n, self.gp_ndim])
        else:
            gp_rvector = np.empty([len(x0_input), self.gp_ndim])

        gp_pams[0] = np.log(parameter_values['det_m32_sigma']**2/self.gp_ndim)
        for data_name, pam_index in self.gp_metric_index[dataset.name_ref].items():
            par_name = 'det_' + data_name + '_m32_rho'
            gp_pams[pam_index+1] = np.log(parameter_values[par_name]**2)

            if x0_input is None:
                gp_rvector[:,pam_index] = self.gp_rvector[dataset.name_ref][:,pam_index]*1.
            else:
                gp_rvector[:,pam_index] = self.interpolated[dataset.name_ref][data_name](x0_input)


        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(self.gp_rvector[dataset.name_ref], env)
        return self.gp[dataset.name_ref].predict(dataset.residuals, gp_rvector, return_cov=return_covariance, return_var=return_variance)

    def sample_conditional(self, parameter_values, dataset, x0_input=None):

        gp_pams = [0.]*(self.gp_ndim+1)

        if x0_input is None:
            gp_rvector = np.empty([dataset.n, self.gp_ndim])
        else:
            gp_rvector = np.empty([len(x0_input), self.gp_ndim])

        gp_pams[0] = np.log(parameter_values['det_m32_sigma']**2/self.gp_ndim)
        for data_name, pam_index in self.gp_metric_index[dataset.name_ref].items():
            par_name = 'det_' + data_name + '_m32_rho'
            gp_pams[pam_index+1] = np.log(parameter_values[par_name]**2)

            if x0_input is None:
                gp_rvector[:,pam_index] = self.gp_rvector[dataset.name_ref][:,pam_index]*1.
            else:
                gp_rvector[:,pam_index] = self.interpolated[dataset.name_ref][data_name](x0_input)

        env = np.sqrt(dataset.e ** 2.0 + dataset.jitter ** 2.0)
        self.gp[dataset.name_ref].set_parameter_vector(gp_pams)
        self.gp[dataset.name_ref].compute(self.gp_rvector[dataset.name_ref], env)
        return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, gp_rvector)

