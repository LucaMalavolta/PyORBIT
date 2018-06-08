from ..classes.common import *
from abstract_common import *
from abstract_model import *

"""
New changes:
    mc.variables  is now called  mc.transformation
    mv.var_list  is now called  mc.variable_index

    variable_index is the third argument of transformation
    it identifies which values from theta must be taken to convert the variable_sampler values to the physical parameter

    variable_sampler associate the value in theta to their label in the sampler space
"""


class SinusoidCommonPeriod(AbstractModel):

    model_class = 'sinusoid_common_period'

    list_pams_common = {
        'P': 'LU'  # Period, log-uniform prior
    }

    list_pams_dataset = {
        'K': 'LU',  # RV semi-amplitude, log-uniform prior
        'f': 'U'  # RV vurve phase, log-uniform prior
    }

    recenter_pams_dataset = {'f'}

    def compute(self, variable_value, dataset, x0_input=None):
        if x0_input is None:
            return kepler_exo.kepler_RV_T0P(dataset.x0,
                                            variable_value['f'],
                                            variable_value['P'],
                                            variable_value['K'],
                                            0.00,
                                            np.pi/2.)
        else:
            return kepler_exo.kepler_RV_T0P(x0_input,
                                            variable_value['f'],
                                            variable_value['P'],
                                            variable_value['K'],
                                            0.00,
                                            np.pi / 2.)

