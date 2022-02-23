from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
import pyorbit.subroutines.constants as constants

try:
    import pycheops
    from lmfit import Parameter, Parameters
except ImportError:
    pass


class CheopsDetrending(AbstractModel):

    default_common = 'cheops_detrending'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'cheops_detrending'
        self.unitary_model = False
        self.normalization_model = True

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
            "dfdsinphi",
            "dfdcosphi",
            "dfdcos2phi",
            "dfdsin2phi",
            "dfdcos3phi",
            "dfdsin3phi",
        }

        self.fit_roll_angle = True

        self.parameters_copy = self.list_pams_dataset.copy()

        ## serve??????
        self.common_cheops_detrending = None
        self.pycheops_params ={}
        self.detrend_model = {}


        self.cheops_diagnostics = ['ramp', etc etc]

        self.cheops_instrumental = {}
        self.cheops_interpolated = {}

        #self.FM = pycheops.models.FactorModel
        #self.make_interp = pycheops.dataset._make_interp

    def initialize_model(self, mc, **kwargs):

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'cheops_detrending':
                self.common_cheops_detrending = common_ref
                break


        self.fit_roll_angle = kwargs.get('fit_roll_angle', True)
        if not self.fit_roll_angle:
            for p in self.parameters_copy:
                if 'phi' in p:
                    self.list_pams_dataset.remove(p)







    def initialize_model_dataset(self, mc, dataset, **kwargs):


        self.cheops_instrumental[dataset.name_ref] = {}
        self.cheops_interpolated[dataset.name_ref] = {}

        for cheops_diagnosi

        if 'ramp' in dataset.ancillary:
            self.cheops_instrumental[dataset.name_ref]['ramp'] = dataset.ancillary['ramp']
        else:
            self.cheops_instrumental[dataset.name_ref]['ramp'] = np.zeros_like(dataset.ancillary['time'])


        #self.cheops_interpolate

        #        ramp = np.zeros_like(lc["time"])
        #        
        #        phi = np.array(lc["roll_angle"])/180.*np.pi
        #        dx     = make_interp(lc["time"], lc["xoff"], scale="range")
        #        dy     = make_interp(lc["time"], lc["yoff"], scale="range")
        #        sinphi = make_interp(lc["time"], np.sin(phi))
        #        cosphi = make_interp(lc["time"], np.cos(phi))
        #        bg     = make_interp(lc["time"], lc["bg"], scale="max")
        #        contam = make_interp(lc["time"], lc["contam"], scale="max")
        #        smear  = make_interp(lc["time"], lc["smear"], scale="max")
        #        deltaT = make_interp(lc["time"], lc["deltaT"])








        for var in self.parameters_copy:
            if var in self.list_pams_dataset: continue

            self.fix_list[dataset.name_ref][var] = np.asarray([
                mc.common_models[self.common_cheops_detrending].default_fixed[var],
                0.000
            ])





        self.pycheops_params[dataset.name_ref] = Parameters()

        # DEFAULT VALUES: use them all
        default_value, default_vary, default_min, default_max = 0.2, True, 0.0, 1.0
        for p in self.parameters_copy:
            if p not in self.list_pams_dataset:
                self.pycheops_params[dataset.name_ref][p] = \
                    Parameter(p,
                        value= 0.0,
                        vary = False,
                        min = -1.,
                        max =  1,
                    )
            else:
                self.pycheops_params[dataset.name_ref][p] = \
                    Parameter(p,
                        value=mc.common_models[self.common_cheops_detrending].default_fixed[p],
                        vary = True,
                        min  = mc.common_models[self.common_cheops_detrending].default_bounds[p][0],
                        max  = mc.common_models[self.common_cheops_detrending].default_bounds[p][1],
                    )

        # tooks functions from pycheops
        #FM = pycheops.models.FactorModel
        #make_interp = pycheops.dataset._make_interp

        """ Function that set the detrend model of cheops data

        :dataset.ancillary: lc data with time, flux, flux_err, roll_angle, xoff, yoff, bg, contam, smear (optional), deltaT (optional)
        :type lc: Dict or numpy.struct_array

        : detrend_model: the detrending model to use
        :rtype detrend_model: FactorModel function
        """
        print(kwargs)
        print(self.fix_list)
        print(self.fixed)
        #self.detrend_model[dataset.name_ref]=  kwargs[dataset.name_ref]['pycheops_detrend_model']

    def compute(self, variable_value, dataset, x0_input=None):

        for p in self.list_pams_dataset:
            self.pycheops_params[dataset.name_ref][p].value = variable_value[p]

        print(variable_value)
        quit()
        return 5.
        #if x0_input is None:
        #    return self.detrend_model[dataset.name_ref].eval(self.pycheops_params[dataset.name_ref], t=dataset.x) # full detrend
        #else:
        #    return self.detrend_model[dataset.name_ref].eval(self.pycheops_params[dataset.name_ref], t=x0_input + dataset.Tref) # full detrend



    ### from pycheops
    def _make_interp(t,x,scale=None):
        if scale is None:
            z = x
        elif np.ptp(x) == 0:
            z = np.zeros_like(x)
        elif scale == 'max':
            z = (x-min(x))/np.ptp(x) 
        elif scale == 'range':
            z = (2*x-(x.min()+x.max()))/np.ptp(x)
        else:
            raise ValueError('scale must be None, max or range')
        return interp1d(t,z,bounds_error=False, fill_value=(z[0],z[-1]))




def broken():

    time = np.array(lc["time"])
    phi = np.array(lc["roll_angle"])*180./np.pi
    try:
        smear = lc["smear"]
    except KeyError:
        smear = np.zeros_like(time)
    try:
        deltaT = lc["deltaT"]
    except KeyError:
        deltaT = np.zeros_like(time)
        
    from scipy.interpolate import interp1d

    def compute_model(variable_values, t, lc):

                ramp = np.zeros_like(lc["time"])
                
                phi = np.array(lc["roll_angle"])/180.*np.pi
                dx     = make_interp(lc["time"], lc["xoff"], scale="range")
                dy     = make_interp(lc["time"], lc["yoff"], scale="range")
                sinphi = make_interp(lc["time"], np.sin(phi))
                cosphi = make_interp(lc["time"], np.cos(phi))
                bg     = make_interp(lc["time"], lc["bg"], scale="max")
                contam = make_interp(lc["time"], lc["contam"], scale="max")
                smear  = make_interp(lc["time"], lc["smear"], scale="max")
                deltaT = make_interp(lc["time"], lc["deltaT"])
                
                
                dt = t - np.median(t)
                trend = 1. +  variable_values['dfdt'] *dt +  variable_values['d2fdt2'] *dt**2 

                trend += variable_values['dfdbg']* bg(t)

                trend += variable_values['dfdcontam']* contam(t)

                trend += variable_values['dfdsmear']* smear(t)

                #trend += ramp*deltaT(t)/1e6

                trend += variable_values['dfdx']* dx(t) + variable_values['d2fdx2']* dx(t)**2

                trend += variable_values['dfdy']* dy(t) + variable_values['d2fdy2'] * dy(t)**2
                trend += variable_values['d2fdxdy'] * dx(t) * dy(t)
                sinphit = sinphi(t)
                cosphit = cosphi(t)
                trend += variable_values['dfdsinphi']*sinphit + variable_values['dfdcosphi']*cosphit
                trend += variable_values['dfdsin2phi']*(2*sinphit*cosphit)
                trend += variable_values['dfdcos2phi']*(2*cosphit**2 - 1)
                trend += variable_values['dfdsin3phi']*(3*sinphit - 4* sinphit**3)
                trend += variable_values['dfdcos3phi']*(4*cosphit**3 - 3*cosphit)

                return trend

    detrend_flux = compute_model(detrending_params, lc['time'], lc)