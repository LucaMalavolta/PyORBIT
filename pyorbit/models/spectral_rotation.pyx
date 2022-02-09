from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
from numpy.polynomial import polynomial

try:
    from PyAstronomy.pyasl import fastRotBroad as PyAstroFastRotBroad
except:
    pass


class SpectralRotation(AbstractModel):

    default_common = 'star_parameters'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        try:
            from PyAstronomy.pyasl import fastRotBroad as PyAstroFastRotBroad
        except ImportError:
            print("ERROR: PyAstronomy not installed, this will not work")
            quit()

        self.model_class = 'spectral_rotation'
        self.unitary_model = True
        self.normalization_model = False

        self.ldvars = {}
        self.ld_ncoeff = 2
        self.parametrization = 'Standard'

        self.list_pams_common = {
            'line_contrast',
            'line_fwhm',
            'v_sini',
        }

        self.list_pams_dataset = {
            'rv_center',
        }

    def initialize_model(self, mc, **kwargs):

        self.reference_wavelength = kwargs.get('reference_wavelength', 5500.)
        self.baseline_RV = kwargs.get('baseline_RV', True)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

    def compute(self, variable_value, dataset, x0_input=None):
        """
        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        if x0_input is None:

            sigma = variable_value['line_fwhm']/constants.sigma2FWHM
            solar_flux = - (variable_value['line_contrast']/100.) \
                * np.exp(-(dataset.x - variable_value['rv_center'])**2
                        / (2 * sigma**2))

            if self.baseline_RV:
                wave_array = self.reference_wavelength * \
                    (1. + 1000. * dataset.x / constants.speed)
            else:
                wave_array = dataset.x

            return PyAstroFastRotBroad(wave_array,
                                   solar_flux,
                                   variable_value['ld_c1'],
                                   variable_value['v_sini'],
                                   effWvl =self.reference_wavelength)

        else:
            return x0_input*0.

    def _prepare_limb_darkening_coefficients(self, mc, **kwargs):
        """ Setting up the limb darkening calculation"""

        self.limb_darkening_model = kwargs['limb_darkening_model']
        self.ld_vars = [0.00] * kwargs['limb_darkening_ncoeff']
        for i_coeff in range(1, kwargs['limb_darkening_ncoeff'] + 1):
            var = 'ld_c' + repr(i_coeff)
            self.ldvars[var] = i_coeff - 1
            self.list_pams_common.update([var])


class SubsetSpectralRotation(AbstractModel):

    default_common = 'star_parameters'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        try:
            from PyAstronomy.pyasl import rotBroad as PyAstroRotBroad
        except ImportError:
            print("ERROR: PyAstronomy not installed, this will not work")
            quit()

        self.model_class = 'subset_spectral_rotation'
        self.unitary_model = True
        self.normalization_model = False

        self.ldvars = {}
        self.ld_ncoeff = 2
        self.parametrization = 'Standard'

        self.list_pams_common = {
            'line_contrast',
            'line_fwhm',
            'v_sini',
        }

        self.list_pams_dataset = set()

    def initialize_model(self, mc, **kwargs):

        self.reference_wavelength = kwargs.get('reference_wavelength', 5500.)
        self.baseline_RV = kwargs.get('baseline_RV', True)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        if not dataset.submodel_flag:
            return

        for i_sub in range(0, dataset.submodel_flag):

            var_original = 'rv_center'
            var_subset = 'rv_center_sub'+repr(i_sub)
            self._subset_transfer_priors(mc, dataset, var_original, var_subset)

    def compute(self, variable_value, dataset, x0_input=None):
        """
        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """

        if x0_input is None:

            sigma = variable_value['line_fwhm']/constants.sigma2FWHM

            if self.baseline_RV:
                wave_array = self.reference_wavelength * \
                    (1. + 1000. * dataset.x / constants.speed)
            else:
                wave_array = dataset.x

            y_output = np.zeros(dataset.n)

            for i_sub in range(0, dataset.submodel_flag):
                var = 'rv_center_sub'+repr(i_sub)
                sel_data = (dataset.submodel_id == i_sub)

                solar_flux = -(variable_value['line_contrast']/100.) \
                    * np.exp(-(dataset.x[sel_data] - variable_value[var])**2 / (2 * sigma**2))

                y_output[sel_data] = PyAstroFastRotBroad(wave_array[sel_data],
                                                     solar_flux,
                                                     variable_value['ld_c1'],
                                                     variable_value['v_sini'],
                                                     effWvl=self.reference_wavelength)
            return y_output
        else:
            return x0_input*0.

    def _prepare_limb_darkening_coefficients(self, mc, **kwargs):
        """ Setting up the limb darkening calculation"""

        self.limb_darkening_model = kwargs['limb_darkening_model']
        self.ld_vars = [0.00] * kwargs['limb_darkening_ncoeff']
        for i_coeff in range(1, kwargs['limb_darkening_ncoeff'] + 1):
            var = 'ld_c' + repr(i_coeff)
            self.ldvars[var] = i_coeff - 1
            self.list_pams_common.update([var])


class SubsetSpectralRotationPolynomial(AbstractModel):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        try:
            from PyAstronomy.pyasl import rotBroad as PyAstroRotBroad
        except ImportError:
            print("ERROR: PyAstronomy not installed, this will not work")
            quit()

        self.model_class = 'subset_spectral_rotation_polynomial'
        self.unitary_model = True
        self.normalization_model = False

        self.ldvars = {}
        self.ld_ncoeff = 2
        self.parametrization = 'Standard'

        self.list_pams_common = {
            'line_contrast',
            'line_fwhm',
            'v_sini',
        }

        self.list_pams_dataset = set()

        self.order = 1
        self.starting_order = 1

    def initialize_model(self, mc, **kwargs):

        self.order = kwargs.get('polynomial_order', 6)
        for i_order in range(self.starting_order, self.order+1):
            var = 'poly_c'+repr(i_order)
            self.list_pams_common.update([var])

        self.reference_wavelength = kwargs.get('reference_wavelength', 5500.)
        self.baseline_RV = kwargs.get('baseline_RV', True)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        if not dataset.submodel_flag:
            return

        for i_sub in range(0, dataset.submodel_flag):

            var_original = 'rv_center'
            var_subset = 'rv_center_sub'+repr(i_sub)
            self._subset_transfer_priors(mc, dataset, var_original, var_subset)

    def compute(self, variable_value, dataset, x0_input=None):
        """
        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """

        coeff = np.zeros(self.order+1)
        for i_order in range(self.starting_order, self.order+1):
            var = 'poly_c'+repr(i_order)
            coeff[i_order] = variable_value[var]

        if x0_input is None:

            sigma = variable_value['line_fwhm']/constants.sigma2FWHM

            if self.baseline_RV:
                wave_array = self.reference_wavelength * \
                    (1. + 1000. * dataset.x / constants.speed)
            else:
                wave_array = dataset.x

            y_output = np.zeros(dataset.n)

            for i_sub in range(0, dataset.submodel_flag):
                var = 'rv_center_sub'+repr(i_sub)
                sel_data = (dataset.submodel_id == i_sub)

                solar_flux = -(variable_value['line_contrast']/100.) \
                    * np.exp(-(dataset.x[sel_data] - variable_value[var])**2 / (2 * sigma**2))

                y_output[sel_data] = PyAstroFastRotBroad(wave_array[sel_data],
                                                     solar_flux,
                                                     variable_value['ld_c1'],
                                                     variable_value['v_sini'],
                                                     effWvl=self.reference_wavelength)
                y_output[sel_data] += polynomial.polyval(dataset.x[sel_data]-variable_value[var], coeff)

            return y_output

        else:
            return x0_input*0.


    def _prepare_limb_darkening_coefficients(self, mc, **kwargs):
        """ Setting up the limb darkening calculation"""

        self.limb_darkening_model = kwargs['limb_darkening_model']
        self.ld_vars = [0.00] * kwargs['limb_darkening_ncoeff']
        for i_coeff in range(1, kwargs['limb_darkening_ncoeff'] + 1):
            var = 'ld_c' + repr(i_coeff)
            self.ldvars[var] = i_coeff - 1
            self.list_pams_common.update([var])
