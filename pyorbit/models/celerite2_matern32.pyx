from pyorbit.classes.common import np, dummy_import_4args
from pyorbit.models.abstract_model import AbstractModel

try:
    import celerite2
    #import autograd.numpy as np
    from celerite2 import terms
except ImportError:
    import numpy as np
    terms = dummy_import_4args


class Celerite2_Matern32(AbstractModel):

    r"""A term that approximates a Matern-3/2 function

    Args:
       log_sigma (float) – The log of the parameter σ.
       log_rho (float) – The log of the parameter ρ.
       eps (Optional[float]) – The value of the parameter ϵ. (default: 0.01)
    """

    internal_likelihood = True

    model_class = 'celerite2_matern32'

    list_pams_common = {}

    list_pams_dataset = {
        'matern32_log10_sigma',  # sigma
        'matern32_log10_rho', # rho
        }

    recenter_pams_dataset = {}

    n_pams = 2

    def __init__(self, *args, **kwargs):
        super(Celerite2_Matern32, self).__init__(*args, **kwargs)
        self.gp = {}

        try:
            import celerite
        except:
            print("ERROR: celerite not installed, this will not work")
            quit()


    def setup_dataset(self, mc, dataset, **kwargs):
        self.define_kernel(dataset)
        return

    def define_kernel(self, dataset):
        # random initialization

        kernel = terms.Matern32Term(0.0, 0.0)
        self.gp[dataset.name_ref] = celerite2.GaussianProcess(kernel, mean=0.0)

        """ I've decided to add the jitter in quadrature instead of using a constant kernel to allow the use of
        different / selective jitter within the dataset
        """
        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag)
        return

    def lnlk_compute(self, variable_value, dataset):
        """ 2 steps:
        In celerite2 the old function "set_parameter_vector" has been removed
        and the kernel has to be defined every time
        """

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = terms.SHOTerm(
            10**variable_value['matern32_log10_sigma'],
            10**variable_value['matern32_log10_rho'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        return self.gp[dataset.name_ref].log_likelihood(dataset.residuals)

    def sample_predict(self, variable_value, dataset, x0_input=None, return_covariance=False, return_variance=False):

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = terms.SHOTerm(
            10**variable_value['matern32_log10_sigma'],
            10**variable_value['matern32_log10_rho'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        if x0_input is None:
            return self.gp[dataset.name_ref].predict(dataset.residuals, dataset.x0, return_cov=return_covariance, return_var=return_variance)
        else:
            return self.gp[dataset.name_ref].predict(dataset.residuals, x0_input, return_cov=return_covariance, return_var=return_variance)

    def sample_conditional(self, variable_value, dataset,  x0_input=None):

        self.gp[dataset.name_ref].mean = 0.
        self.gp[dataset.name_ref].kernel = terms.SHOTerm(
            10**variable_value['matern32_log10_sigma'],
            10**variable_value['matern32_log10_rho'])

        diag = dataset.e ** 2.0 + dataset.jitter ** 2.0
        self.gp[dataset.name_ref].compute(dataset.x0, diag=diag, quiet=True)

        if x0_input is None:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, dataset.x0)
        else:
            return self.gp[dataset.name_ref].sample_conditional(dataset.residuals, x0_input)
