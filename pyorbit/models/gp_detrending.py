from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
import pyorbit.subroutines.constants as constants
from numpy.polynomial import polynomial
from numpy.lib.recfunctions import append_fields
from scipy.interpolate import interp1d

class GPDetrending(AbstractModel):

    default_common = 'detrending'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        try:
            import celerite2
        except:
            print("ERROR: celerite2 not installed, this will not work")
            quit()

        import warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning)


        self.model_class = 'detrending'
        self.internal_likelihood = True

        self.local_model = True

        self.list_pams_common = set()
        self.list_pams_dataset = {
            'det_m32_sigma',  # sigma
        }
        self.interpolated = {}

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

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        """ For plotting purposes only """
        self.interpolated[dataset.name_ref] = {}

        for data_name in kwargs['detrending_variables']:

            if data_name not in dataset.ancillary.dtype.names:
                dataset.ancillary = append_fields(dataset.ancillary, data_name, np.zeros(dataset.n))
                self.pams_order[data_name] = 0
                continue

            for i_order in range(1, self.pams_order[data_name]+1):

                par_original = 'det_m32_rho'+repr(i_order)
                par_addition = 'det_' + data_name + '_c'+repr(i_order)

                if self.local_model:
                    self.transfer_parameter_properties(mc, dataset, par_original, par_addition, keywords=kwargs, dataset_pam=True)
                else:
                    self.transfer_parameter_properties(mc, dataset, par_original, par_addition, keywords=kwargs, common_pam=True)

            self.interpolated[dataset.name_ref][data_name]=interp1d(
            dataset.x0,
            dataset.ancillary[data_name],
            bounds_error=False,
            fill_value=(np.amin(dataset.ancillary[data_name]), np.amax(dataset.ancillary[data_name])))


    def compute(self, parameter_values, dataset, x0_input=None):

        if x0_input is None:

            trend = np.ones(dataset.n) * parameter_values.get('det_c0', self.baseline_value)

            for data_name, data_order in self.pams_order.items():

                coeff = np.zeros(data_order+1)
                for i_order in range(1, data_order+1):
                    par_addition = 'det_' + data_name + '_c'+repr(i_order)
                    coeff[i_order] = parameter_values[par_addition]

                trend += polynomial.polyval(dataset.ancillary[data_name]-parameter_values['x_zero_'+data_name], coeff)
                #if dataset.name_ref == 'GROND_r':
                #    print(self.model_name, data_name, trend[20:23], coeff)
                #    print()
            if self.exponential_detrending:
                if self.natural_base:
                    return np.exp(trend)
                else:
                    return 10**(trend)
            else:
                return trend

        else:

            trend = np.ones_like(x0_input) * parameter_values.get('det_c0', self.baseline_value)

            for data_name, data_order in self.pams_order.items():

                coeff = np.zeros(data_order+1)
                for i_order in range(1, data_order+1):
                    par_addition = 'det_' + data_name + '_c'+repr(i_order)
                    coeff[i_order] = parameter_values[par_addition]

                trend += polynomial.polyval(self.interpolated[dataset.name_ref][data_name](x0_input)-parameter_values['x_zero_'+data_name], coeff)

            if self.exponential_detrending:
                if self.natural_base:
                    return np.exp(trend)
                else:
                    return 10**(trend)
            else:
                return trend


class PolynomialDetrending(Detrending):

    def __init__(self, *args, **kwargs):
        AbstractModel.__init__(self, *args, **kwargs)
        Detrending.__init__(self, *args, **kwargs)
        self.exponential_detrending = False

class ExponentialDetrending(Detrending):

    def __init__(self, *args, **kwargs):
        AbstractModel.__init__(self, *args, **kwargs)
        Detrending.__init__(self, *args, **kwargs)
        self.exponential_detrending = True