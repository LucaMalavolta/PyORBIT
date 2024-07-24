from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
import pyorbit.subroutines.constants as constants

from scipy.interpolate import interp1d
from numpy.lib.recfunctions import append_fields

class CheopsDetrending(AbstractModel):

    default_common = 'cheops_modelling'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'cheops_detrending'
        self.unitary_model = True
        self.normalization_model = False

        self.list_pams_dataset = OrderedSet([
            "dfdbg",
            "dfdcontam",
            "dfdsmear",
            "dfdx",
            "dfdy",
            "d2fdx2",
            "d2fdxdy",
            "d2fdy2",
        ])

        self.fit_roll_angle = True

        self.roll_angle_parameters = OrderedSet([
            "dfdsinphi",
            "dfdcosphi",
            "dfdcos2phi",
            "dfdsin2phi",
            "dfdcos3phi",
            "dfdsin3phi",
        ])

        self.cheops_diagnostics = {
            'roll_angle': {
                'scale': 'None',
                'pams': []
                },
            'ramp': {
                'scale': 'max',
                'pams': []
                },  # ???
            'smear': {
                'scale': 'max',
                'pams': ['dfdsmear']
                },
            'deltaT': {
                'scale': 'None',
                'pams': []
                },
            'xoff': {
                'scale': 'range',
                'pams': ['dfdx', 'd2fdx2','d2fdxdy']
                },
            'yoff': {
                'scale': 'range',
                'pams': ['dfdy', 'd2fdy2','d2fdxdy']
                },
            'bg': {
                'scale': 'max',
                'pams': ['dfdbg']
                },
            'contam': {
                'scale': 'max',
                'pams': ['dfdcontam']
                },
            'sinphi': {
                'scale': 'None',
                'pams': ['dfdsinphi','dfdsin2phi','dfdsin3phi']
                },
            'cosphi': {
                'scale': 'None',
                'pams': ['dfdcosphi','dfdsin2phi','dfdcos2phi','dfdcos3phi' ]
                }
        }


        self.cheops_instrumental = {}
        self.cheops_interpolated = {}

        self.retrieve_roll_angle = None

    def initialize_model(self, mc, **kwargs):

        if kwargs.get('fit_roll_angle', True):

            for par in self.roll_angle_parameters:
                self.list_pams_dataset.update([par])

            self.retrieve_roll_angle = self._retrieve_roll_angle_mod01
        else:
            self.retrieve_roll_angle = self._retrieve_roll_angle_mod00

        if kwargs.get('fit_quadratic_trend', False):
            self.list_pams_dataset.update(['dfdt'])
            self.list_pams_dataset.update(['d2fdt2'])

            self.retrieve_trend = self._retrieve_trend_mod02

        elif kwargs.get('fit_linear_trend', False):
            self.list_pams_dataset.update(['dfdt'])
            self.retrieve_trend = self._retrieve_trend_mod01
        else:
            self.retrieve_trend = self._retrieve_trend_mod00

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        self.cheops_instrumental[dataset.name_ref] = {}
        self.cheops_interpolated[dataset.name_ref] = {}

        if 'sinphi' not in dataset.ancillary.dtype.names:
            dataset.ancillary = append_fields(dataset.ancillary, 'sinphi', np.sin(dataset.ancillary["roll_angle"]/180.*np.pi))
        if 'cosphi' not in dataset.ancillary.dtype.names:
            dataset.ancillary = append_fields(dataset.ancillary, 'cosphi', np.cos(dataset.ancillary["roll_angle"]/180.*np.pi))

        for diag_name, diag_dict in self.cheops_diagnostics.items():
            if diag_name in dataset.ancillary.dtype.names:
                if diag_dict['scale'] == 'max':
                    self.cheops_instrumental[dataset.name_ref][diag_name] = \
                        (dataset.ancillary[diag_name] - np.amin(dataset.ancillary[diag_name])) \
                            / np.ptp(dataset.ancillary[diag_name])
                elif diag_dict['scale'] == 'range':
                    self.cheops_instrumental[dataset.name_ref][diag_name] = \
                        ( 2 * dataset.ancillary[diag_name]
                        - (np.amin(dataset.ancillary[diag_name]) + np.amax(dataset.ancillary[diag_name]))) \
                            / np.ptp(dataset.ancillary[diag_name])
                else:
                    self.cheops_instrumental[dataset.name_ref][diag_name] = dataset.ancillary[diag_name]
                print()
            else:
                for pams_fixed  in diag_dict['pams']:
                    self.fix_list[dataset.name_ref][pams_fixed] = np.asarray([0.000000, 0.0000])
                self.cheops_instrumental[dataset.name_ref][diag_name] = np.zeros_like(dataset.ancillary['time'])

            self.cheops_interpolated[dataset.name_ref][diag_name] =interp1d(
                dataset.ancillary['time'],
                self.cheops_instrumental[dataset.name_ref][diag_name],
                  bounds_error=False,
                  fill_value=(self.cheops_instrumental[dataset.name_ref][diag_name][0],self.cheops_instrumental[dataset.name_ref][diag_name][-1]))


    def compute(self, parameter_values, dataset, x0_input=None):

        if x0_input is None:

            trend = parameter_values['dfdbg']* self.cheops_instrumental[dataset.name_ref]['bg']

            trend += parameter_values['dfdcontam']* self.cheops_instrumental[dataset.name_ref]['contam']

            trend += parameter_values['dfdsmear']* self.cheops_instrumental[dataset.name_ref]['smear']

            trend += self.cheops_instrumental[dataset.name_ref]['ramp'] * self.cheops_instrumental[dataset.name_ref]['deltaT'] /1e6

            trend += parameter_values['dfdx']* self.cheops_instrumental[dataset.name_ref]['xoff'] \
                + parameter_values['d2fdx2']*self.cheops_instrumental[dataset.name_ref]['xoff']**2

            trend += parameter_values['dfdy']* self.cheops_instrumental[dataset.name_ref]['yoff'] \
                + parameter_values['d2fdy2'] * self.cheops_instrumental[dataset.name_ref]['yoff']**2

            trend += parameter_values['d2fdxdy'] \
                * self.cheops_instrumental[dataset.name_ref]['xoff'] \
                * self.cheops_instrumental[dataset.name_ref]['yoff']

            trend += self.retrieve_trend(parameter_values, dataset.x)
            trend += self.retrieve_roll_angle(parameter_values, self.cheops_instrumental[dataset.name_ref]['sinphi'], self.cheops_instrumental[dataset.name_ref]['cosphi'])
            #if self.fit_roll_angle:
            #    trend += parameter_values['dfdsinphi']*self.cheops_instrumental[dataset.name_ref]['sinphi'] \
            #        + parameter_values['dfdcosphi']*self.cheops_instrumental[dataset.name_ref]['cosphi']
            #    trend += parameter_values['dfdsin2phi']*(2*self.cheops_instrumental[dataset.name_ref]['sinphi']*self.cheops_instrumental[dataset.name_ref]['cosphi'])
            #    trend += parameter_values['dfdcos2phi']*(2*self.cheops_instrumental[dataset.name_ref]['cosphi']**2 - 1)
            #    trend += parameter_values['dfdsin3phi']*(3*self.cheops_instrumental[dataset.name_ref]['sinphi'] - 4* self.cheops_instrumental[dataset.name_ref]['sinphi']**3)
            #    trend += parameter_values['dfdcos3phi']*(4*self.cheops_instrumental[dataset.name_ref]['cosphi']**3 - 3*self.cheops_instrumental[dataset.name_ref]['cosphi'])

        else:
            t = x0_input + dataset.Tref
            trend = parameter_values['dfdbg']* self.cheops_interpolated[dataset.name_ref]['bg'](t)

            trend += parameter_values['dfdcontam']* self.cheops_interpolated[dataset.name_ref]['contam'](t)

            trend += parameter_values['dfdsmear']* self.cheops_interpolated[dataset.name_ref]['smear'](t)

            trend += self.cheops_interpolated[dataset.name_ref]['ramp'](t) * self.cheops_interpolated[dataset.name_ref]['deltaT'](t) /1e6

            trend += parameter_values['dfdx']* self.cheops_interpolated[dataset.name_ref]['xoff'](t) + parameter_values['d2fdx2']*self.cheops_interpolated[dataset.name_ref]['xoff'](t)**2

            trend += parameter_values['dfdy']* self.cheops_interpolated[dataset.name_ref]['yoff'](t) + parameter_values['d2fdy2'] * self.cheops_interpolated[dataset.name_ref]['yoff'](t)**2
            trend += parameter_values['d2fdxdy'] * self.cheops_interpolated[dataset.name_ref]['xoff'](t) * self.cheops_interpolated[dataset.name_ref]['yoff'](t)

            trend += self.retrieve_trend(parameter_values, t)
            trend += self.retrieve_roll_angle(parameter_values, self.cheops_interpolated[dataset.name_ref]['sinphi'](t), self.cheops_interpolated[dataset.name_ref]['cosphi'](t))

        return trend

    @staticmethod
    def _retrieve_roll_angle_mod00(parameter_values, sinphi, cosphi):
        # internal_dictionary = self.cheops_instrumental[dataset.name_ref]
        return 0.000

    @staticmethod
    def _retrieve_roll_angle_mod01(parameter_values, sinphi, cosphi):
        # internal_dictionary = self.cheops_instrumental[dataset.name_ref]
        return parameter_values['dfdsinphi']*sinphi \
                    + parameter_values['dfdcosphi']*cosphi \
                    + parameter_values['dfdsin2phi']*(2*sinphi*cosphi) \
                    + parameter_values['dfdcos2phi']*(2*cosphi**2 - 1) \
                    + parameter_values['dfdsin3phi']*(3*sinphi - 4* sinphi**3) \
                    + parameter_values['dfdcos3phi']*(4*cosphi**3 - 3*cosphi)

    @staticmethod
    def _retrieve_trend_mod00(parameter_values, t):
        return 0.000

    @staticmethod
    def _retrieve_trend_mod01(parameter_values, t):
        dt = t - np.median(t)
        return parameter_values['dfdt'] *dt

    @staticmethod
    def _retrieve_trend_mod02(parameter_values, t):
        dt = t - np.median(t)
        return parameter_values['dfdt'] *dt +  parameter_values['d2fdt2'] *dt**2
