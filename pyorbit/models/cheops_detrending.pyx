from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
import pyorbit.subroutines.constants as constants

from scipy.interpolate import interp1d
from numpy.lib.recfunctions import append_fields

class CheopsDetrending(AbstractModel):

    default_common = 'cheops_detrending'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'cheops_detrending'
        self.unitary_model = True
        self.normalization_model = False

        self.list_pams_dataset = {
            "dfdt",
            "d2fdt2",
            "dfdbg",
            "dfdcontam",
            "dfdsmear",
            "dfdx",
            "dfdy",
            "d2fdx2",
            "d2fdxdy",
            "d2fdy2",
        }

        self.fit_roll_angle = True

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

    def initialize_model(self, mc, **kwargs):

        self.fit_roll_angle = kwargs.get('fit_roll_angle', True)
        if self.fit_roll_angle:

            for var in self.roll_angle_parameters:
                self.list_pams_dataset.update([var])

        print(self.list_pams_dataset)
    def initialize_model_dataset(self, mc, dataset, **kwargs):

        self.cheops_instrumental[dataset.name_ref] = {}
        self.cheops_interpolated[dataset.name_ref] = {}

        if 'sinphi' not in dataset.ancillary.dtype.names:
            dataset.ancillary = append_fields(dataset.ancillary, 'sinphi', np.sin(dataset.ancillary["roll_angle"]/180.*np.pi))
        if 'cosphi' not in dataset.ancillary.dtype.names:
            dataset.ancillary = append_fields(dataset.ancillary, 'cosphi', np.cos(dataset.ancillary["roll_angle"]/180.*np.pi))

        for diag_name, diag_scale in self.cheops_diagnostics.items():
            if diag_name in dataset.ancillary.dtype.names:
                if diag_scale is 'max':
                    self.cheops_instrumental[dataset.name_ref][diag_name] = \
                        (dataset.ancillary[diag_name] - np.amin(dataset.ancillary[diag_name])) \
                            / np.ptp(dataset.ancillary[diag_name])
                elif diag_scale is 'range':
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

            dt = dataset.x - np.median(dataset.x)
            # average model is zero - it is added to the light curve 
            trend = variable_value['dfdt'] *dt +  variable_value['d2fdt2'] *dt**2 

            trend += variable_value['dfdbg']* self.cheops_instrumental[dataset.name_ref]['bg']

            trend += variable_value['dfdcontam']* self.cheops_instrumental[dataset.name_ref]['contam']

            trend += variable_value['dfdsmear']* self.cheops_instrumental[dataset.name_ref]['smear']

            trend += self.cheops_instrumental[dataset.name_ref]['smear'] * self.cheops_instrumental[dataset.name_ref]['deltaT'] /1e6

            trend += variable_value['dfdx']* self.cheops_instrumental[dataset.name_ref]['xoff'] + variable_value['d2fdx2']*self.cheops_instrumental[dataset.name_ref]['xoff']**2

            trend += variable_value['dfdy']* self.cheops_instrumental[dataset.name_ref]['yoff'] + variable_value['d2fdy2'] * self.cheops_instrumental[dataset.name_ref]['yoff']**2
            trend += variable_value['d2fdxdy'] * self.cheops_instrumental[dataset.name_ref]['xoff'] * self.cheops_instrumental[dataset.name_ref]['yoff']

            if self.fit_roll_angle:
                trend += variable_value['dfdsinphi']*self.cheops_instrumental[dataset.name_ref]['sinphi'] + variable_value['dfdcosphi']*self.cheops_instrumental[dataset.name_ref]['cosphi']
                trend += variable_value['dfdsin2phi']*(2*self.cheops_instrumental[dataset.name_ref]['sinphi']*self.cheops_instrumental[dataset.name_ref]['cosphi'])
                trend += variable_value['dfdcos2phi']*(2*self.cheops_instrumental[dataset.name_ref]['cosphi']**2 - 1)
                trend += variable_value['dfdsin3phi']*(3*self.cheops_instrumental[dataset.name_ref]['sinphi'] - 4* self.cheops_instrumental[dataset.name_ref]['sinphi']**3)
                trend += variable_value['dfdcos3phi']*(4*self.cheops_instrumental[dataset.name_ref]['cosphi']**3 - 3*self.cheops_instrumental[dataset.name_ref]['cosphi'])
            return trend

        else:

            t = x0_input + dataset.Tref
            dt = t - np.median(dataset.x)
            # average model is zero - it is added to the light curve
            trend = variable_value['dfdt'] *dt +  variable_value['d2fdt2'] *dt**2 

            trend += variable_value['dfdbg']* self.cheops_interpolated[dataset.name_ref]['bg'](t)

            trend += variable_value['dfdcontam']* self.cheops_interpolated[dataset.name_ref]['contam'](t)

            trend += variable_value['dfdsmear']* self.cheops_interpolated[dataset.name_ref]['smear'](t)

            trend += self.cheops_interpolated[dataset.name_ref]['smear'](t) * self.cheops_interpolated[dataset.name_ref]['deltaT'](t) /1e6

            trend += variable_value['dfdx']* self.cheops_interpolated[dataset.name_ref]['xoff'](t) + variable_value['d2fdx2']*self.cheops_interpolated[dataset.name_ref]['xoff'](t)**2

            trend += variable_value['dfdy']* self.cheops_interpolated[dataset.name_ref]['yoff'](t) + variable_value['d2fdy2'] * self.cheops_interpolated[dataset.name_ref]['yoff'](t)**2
            trend += variable_value['d2fdxdy'] * self.cheops_interpolated[dataset.name_ref]['xoff'](t) * self.cheops_interpolated[dataset.name_ref]['yoff'](t)

            if self.fit_roll_angle:
                trend += variable_value['dfdsinphi']*self.cheops_interpolated[dataset.name_ref]['sinphi'](t) + variable_value['dfdcosphi']*self.cheops_interpolated[dataset.name_ref]['cosphi'](t)
                trend += variable_value['dfdsin2phi']*(2*self.cheops_interpolated[dataset.name_ref]['sinphi'](t)*self.cheops_interpolated[dataset.name_ref]['cosphi'](t))
                trend += variable_value['dfdcos2phi']*(2*self.cheops_interpolated[dataset.name_ref]['cosphi'](t)**2 - 1)
                trend += variable_value['dfdsin3phi']*(3*self.cheops_interpolated[dataset.name_ref]['sinphi'](t) - 4* self.cheops_interpolated[dataset.name_ref]['sinphi'](t)**3)
            return trend

