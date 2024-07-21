from pyorbit.subroutines.common import *
from pyorbit.subroutines.constants import *
from pyorbit.models.abstract_model import *
from pyorbit.models.abstract_model import *
from numpy.polynomial import polynomial

class Sinusoid(AbstractModel):

    default_common = 'sinusoid'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'sinusoid'

        self.list_pams_common = OrderedSet([
            'sine_period'  # Period, log-uniform prior
            'sine_amp',
            'sine_phase'
        ])

    def compute(self, parameter_values, dataset, x0_input=None):
        if x0_input is None:
            return parameter_values['sine_amp'] * np.sin(dataset.x0/parameter_values['sine_period'] - parameter_values['sine_phase']*deg2rad )
        else:
            return parameter_values['sine_amp'] * np.sin(x0_input/parameter_values['sine_period'] - parameter_values['sine_phase']*deg2rad )



class LocalSinusoid(AbstractModel):

    default_common = 'sinusoid'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'sinusoid'

        self.list_pams_dataset = OrderedSet([
            'sine_period'  # Period, log-uniform prior
            'sine_amp',
            'sine_phase'
        ])

    def compute(self, parameter_values, dataset, x0_input=None):
        if x0_input is None:
            return parameter_values['sine_amp'] * np.sin(dataset.x0/parameter_values['sine_period'] - parameter_values['sine_phase']*deg2rad )
        else:
            return parameter_values['sine_amp'] * np.sin(x0_input/parameter_values['sine_period'] - parameter_values['sine_phase']*deg2rad )


class SinusoidCommonPeriod(AbstractModel):

    default_common = 'sinusoid'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'sinusoid_common_period'

        self.list_pams_common = OrderedSet([
            'sine_period'  # Period, log-uniform prior
        ])

        self.list_pams_dataset = OrderedSet([
            'sine_amp',
            'sine_phase'
        ])

    def compute(self, parameter_values, dataset, x0_input=None):
        if x0_input is None:
            return parameter_values['sine_amp'] * np.sin(dataset.x0/parameter_values['sine_period'] - parameter_values['sine_phase']*deg2rad )

            #return kepler_exo.kepler_RV_T0P(dataset.x0,
            #                                parameter_values['f'],
            #                                parameter_values['P'],
            #                                parameter_values['K'],
            #                                0.00,
            #                                np.pi/2.)
        else:
            return parameter_values['sine_amp'] * np.sin(x0_input/parameter_values['sine_period'] - parameter_values['sine_phase']*deg2rad )

            #return kepler_exo.kepler_RV_T0P(x0_input,
            #                                parameter_values['f'],
            #                                parameter_values['P'],
            #                                parameter_values['K'],
            #                                0.00,
            #                                np.pi / 2.)



class SinusoidPolynomialModulation(AbstractModel):

    default_common = 'sinusoid'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'sinusoid_polynomial_modulation'

        self.list_pams_common = OrderedSet([
            'sine_period',  # Period, log-uniform prior
            'sine_phase',
            'x_zero'
        ])

        self.list_pams_dataset = OrderedSet([
            'sine_amp',
        ])

        self.order = 1
        self.starting_order = 0

        """
        The x-intercept must be defined within the interval of at least one dataset,
        otherwise there will be a degeneracy between the offset parameter and the coefficients
        of the polynomial
        """
        self.x_zero = None
        self.common_poly_ref = None

        self.time_interval = 1.0000000
        self.time_offset = False
        self.count_dataset = 0
        self.reference_dataset = None
        self.reference_kind = None

    def initialize_model(self, mc, **kwargs):

        self.order = kwargs.get('order', 1)

        """ The user may decide to include the 0th order anyway - be aware of correlations with dataset offset!"""
        if kwargs.get('include_zero_point', False):
            self.starting_order = 0

        """ The user may decide to compute the polynomial parameters over a different time interval
            useful for long-term with very slow variations over a single day
        """
        self.time_interval = kwargs.get('time_interval', 1.000000000)

        self.time_offset = kwargs.get('time_offset', False)
        if self.time_offset:
            self.list_pams_dataset.update(['x_offset'])

            try:
                self.reference_dataset = kwargs.get('reference_dataset')
            except:
                self.reference_kind = kwargs.get('reference_kind', 'RV')


        for i_order in range(self.starting_order, self.order+1):
            par = 'poly_c'+repr(i_order)
            self.list_pams_common.update([par])

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'sinusoid':
                self.common_poly_ref = common_ref
                break

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        try:
            mc.common_models[self.common_poly_ref].fix_list['x_zero'] = np.asarray([kwargs['x_zero'], 0.0000], dtype=np.double)
        except (KeyError, ValueError):
            if np.amin(dataset.x) < mc.Tref < np.amax(dataset.x):
                self.x_zero = mc.Tref
            elif not self.x_zero:
                self.x_zero = np.average(dataset.x)
            mc.common_models[self.common_poly_ref].fix_list['x_zero'] = np.asarray([self.x_zero, 0.0000])

        if dataset.name_ref == self.reference_dataset or dataset.kind == self.reference_kind:
            self.fix_list[dataset.name_ref]['x_offset'] = np.asarray([0.0000, 0.0000])

        mc.common_models[self.common_poly_ref].fix_list['poly_c0'] = np.asarray([1.000000, 0.0000])
        #mc.common_models[self.common_poly_ref].fix_list['poly_c1'] = np.asarray([1.000000, 0.0000])


    def compute(self, parameter_values, dataset, x0_input=None):

        coeff = np.zeros(self.order+1)

        if 'x_offset' in parameter_values:
            x_offset = parameter_values['x_offset']
        else:
            x_offset = 0

        for i_order in range(self.starting_order, self.order+1):
            par = 'poly_c'+repr(i_order)
            coeff[i_order] = parameter_values[par]

        if x0_input is None:
            return parameter_values['sine_amp'] \
                * np.sin(dataset.x0/parameter_values['sine_period'] - parameter_values['sine_phase']*deg2rad ) \
                * polynomial.polyval((dataset.x-parameter_values['x_zero']-x_offset)/self.time_interval, coeff)

        else:
            return parameter_values['sine_amp'] \
                * np.sin(x0_input/parameter_values['sine_period'] - parameter_values['sine_phase']*deg2rad ) \
                * polynomial.polyval((x0_input+dataset.Tref-parameter_values['x_zero']-x_offset)/self.time_interval, coeff)
