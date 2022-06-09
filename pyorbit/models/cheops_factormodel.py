from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
import pyorbit.subroutines.constants as constants

from scipy.interpolate import interp1d
from numpy.lib.recfunctions import append_fields

class CheopsFactorModel(AbstractModel):

    default_common = 'cheops_modelling'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'cheops_factormodel'
        self.unitary_model = False
        self.normalization_model = True

        self.list_pams_dataset = {
            "scale_factor",
            "dfdbg",
            "dfdcontam",
            "dfdsmear",
            "dfdx",
            "dfdy",
            "d2fdx2",
            "d2fdxdy",
            "d2fdy2",
        }

        self.roll_angle_parameters = {
            "dfdsinphi",
            "dfdcosphi",
            "dfdcos2phi",
            "dfdsin2phi",
            "dfdcos3phi",
            "dfdsin3phi",
        }

        self.cheops_diagnostics = {
            'roll_angle': 'None',
            'ramp': 'max', # ???
            'smear': 'max',
            'deltaT': 'None',
            'xoff': 'range',
            'yoff': 'range',
            'bg': 'max',
            'contam': 'max',
            'sinphi': 'None',
            'cosphi': 'None'
        }

        self.cheops_instrumental = {}
        self.cheops_interpolated = {}

        self.retrieve_roll_angle = None

    def initialize_model(self, mc, **kwargs):

        if kwargs.get('fit_roll_angle', True):

            for var in self.roll_angle_parameters:
                self.list_pams_dataset.update([var])

            self.retrieve_roll_angle = self._retrieve_roll_angle_mod01
        else:
            self.retrieve_roll_angle = self._retrieve_roll_angle_mod00

        if kwargs.get('fit_constant_trend', False):
            self.retrieve_trend = self._retrieve_trend_mod00
        elif kwargs.get('fit_linear_trend', False):
            self.list_pams_dataset.update(['dfdt'])
            self.retrieve_trend = self._retrieve_trend_mod01
        else:
            self.list_pams_dataset.update(['dfdt'])
            self.list_pams_dataset.update(['d2fdt2'])
            self.retrieve_trend = self._retrieve_trend_mod02

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        self.cheops_instrumental[dataset.name_ref] = {}
        self.cheops_interpolated[dataset.name_ref] = {}

        if 'sinphi' not in dataset.ancillary.dtype.names:
            dataset.ancillary = append_fields(dataset.ancillary, 'sinphi', np.sin(dataset.ancillary["roll_angle"]/180.*np.pi))
        if 'cosphi' not in dataset.ancillary.dtype.names:
            dataset.ancillary = append_fields(dataset.ancillary, 'cosphi', np.cos(dataset.ancillary["roll_angle"]/180.*np.pi))

        for diag_name, diag_scale in self.cheops_diagnostics.items():
            if diag_name in dataset.ancillary.dtype.names:
                if diag_scale == 'max':
                    self.cheops_instrumental[dataset.name_ref][diag_name] = \
                        (dataset.ancillary[diag_name] - np.amin(dataset.ancillary[diag_name])) \
                            / np.ptp(dataset.ancillary[diag_name])
                elif diag_scale == 'range':
                    self.cheops_instrumental[dataset.name_ref][diag_name] = \
                        ( 2 * dataset.ancillary[diag_name]
                        - (np.amin(dataset.ancillary[diag_name]) + np.amax(dataset.ancillary[diag_name]))) \
                            / np.ptp(dataset.ancillary[diag_name])
                else:
                    self.cheops_instrumental[dataset.name_ref][diag_name] = dataset.ancillary[diag_name]
                print()
            else:
                self.cheops_instrumental[dataset.name_ref][diag_name] = np.zeros_like(dataset.ancillary['time'])

            self.cheops_interpolated[dataset.name_ref][diag_name] =interp1d(
                dataset.ancillary['time'],
                self.cheops_instrumental[dataset.name_ref][diag_name],
                  bounds_error=False,
                  fill_value=(self.cheops_instrumental[dataset.name_ref][diag_name][0],self.cheops_instrumental[dataset.name_ref][diag_name][-1]))


    def compute(self, variable_value, dataset, x0_input=None):

        if x0_input is None:

            trend = variable_value['dfdbg']* self.cheops_instrumental[dataset.name_ref]['bg']

            trend += variable_value['dfdcontam']* self.cheops_instrumental[dataset.name_ref]['contam']

            trend += variable_value['dfdsmear']* self.cheops_instrumental[dataset.name_ref]['smear']

            trend += self.cheops_instrumental[dataset.name_ref]['ramp'] * self.cheops_instrumental[dataset.name_ref]['deltaT'] /1e6

            trend += variable_value['dfdx']* self.cheops_instrumental[dataset.name_ref]['xoff'] \
                + variable_value['d2fdx2']*self.cheops_instrumental[dataset.name_ref]['xoff']**2

            trend += variable_value['dfdy']* self.cheops_instrumental[dataset.name_ref]['yoff'] \
                + variable_value['d2fdy2'] * self.cheops_instrumental[dataset.name_ref]['yoff']**2
            trend += variable_value['d2fdxdy'] \
                * self.cheops_instrumental[dataset.name_ref]['xoff'] \
                * self.cheops_instrumental[dataset.name_ref]['yoff']

            trend += self.retrieve_trend(variable_value, dataset.x)
            trend += self.retrieve_roll_angle(variable_value, self.cheops_instrumental[dataset.name_ref]['sinphi'], self.cheops_instrumental[dataset.name_ref]['cosphi'])
            #if self.fit_roll_angle:
            #    trend += variable_value['dfdsinphi']*self.cheops_instrumental[dataset.name_ref]['sinphi'] \
            #        + variable_value['dfdcosphi']*self.cheops_instrumental[dataset.name_ref]['cosphi']
            #    trend += variable_value['dfdsin2phi']*(2*self.cheops_instrumental[dataset.name_ref]['sinphi']*self.cheops_instrumental[dataset.name_ref]['cosphi'])
            #    trend += variable_value['dfdcos2phi']*(2*self.cheops_instrumental[dataset.name_ref]['cosphi']**2 - 1)
            #    trend += variable_value['dfdsin3phi']*(3*self.cheops_instrumental[dataset.name_ref]['sinphi'] - 4* self.cheops_instrumental[dataset.name_ref]['sinphi']**3)
            #    trend += variable_value['dfdcos3phi']*(4*self.cheops_instrumental[dataset.name_ref]['cosphi']**3 - 3*self.cheops_instrumental[dataset.name_ref]['cosphi'])

        else:
            t = x0_input + dataset.Tref
            trend = variable_value['dfdbg']* self.cheops_interpolated[dataset.name_ref]['bg'](t)

            trend += variable_value['dfdcontam']* self.cheops_interpolated[dataset.name_ref]['contam'](t)

            trend += variable_value['dfdsmear']* self.cheops_interpolated[dataset.name_ref]['smear'](t)

            trend += self.cheops_interpolated[dataset.name_ref]['ramp'](t) * self.cheops_interpolated[dataset.name_ref]['deltaT'](t) /1e6

            trend += variable_value['dfdx']* self.cheops_interpolated[dataset.name_ref]['xoff'](t) + variable_value['d2fdx2']*self.cheops_interpolated[dataset.name_ref]['xoff'](t)**2

            trend += variable_value['dfdy']* self.cheops_interpolated[dataset.name_ref]['yoff'](t) + variable_value['d2fdy2'] * self.cheops_interpolated[dataset.name_ref]['yoff'](t)**2
            trend += variable_value['d2fdxdy'] * self.cheops_interpolated[dataset.name_ref]['xoff'](t) * self.cheops_interpolated[dataset.name_ref]['yoff'](t)

            trend += self.retrieve_trend(variable_value, t)
            trend += self.retrieve_roll_angle(variable_value, self.cheops_interpolated[dataset.name_ref]['sinphi'](t), self.cheops_interpolated[dataset.name_ref]['cosphi'](t))

        return trend * variable_value['dfdbg']

    @staticmethod
    def _retrieve_roll_angle_mod00(variable_value, sinphi, cosphi):
        # internal_dictionary = self.cheops_instrumental[dataset.name_ref]
        return 0.000

    @staticmethod
    def _retrieve_roll_angle_mod01(variable_value, sinphi, cosphi):
        # internal_dictionary = self.cheops_instrumental[dataset.name_ref]
        return variable_value['dfdsinphi']*sinphi \
                    + variable_value['dfdcosphi']*cosphi \
                    + variable_value['dfdsin2phi']*(2*sinphi*cosphi) \
                    + variable_value['dfdcos2phi']*(2*cosphi**2 - 1) \
                    + variable_value['dfdsin3phi']*(3*sinphi - 4* sinphi**3) \
                    + variable_value['dfdcos3phi']*(4*cosphi**3 - 3*cosphi)

    @staticmethod
    def _retrieve_trend_mod00(variable_value, t):
        return 0.000

    @staticmethod
    def _retrieve_trend_mod01(variable_value, t):
        dt = t - np.median(t)
        return variable_value['dfdt'] *dt

    @staticmethod
    def _retrieve_trend_mod02(variable_value, t):
        dt = t - np.median(t)
        return variable_value['dfdt'] *dt +  variable_value['d2fdt2'] *dt**2
