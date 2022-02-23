from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
import pyorbit.subroutines.constants as constants

try:
    import pycheops
    from lmfit import Parameter, Parameters
except ImportError:
    pass


FM = pycheops.models.FactorModel
make_interp = pycheops.dataset._make_interp


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
        self.common_cheops_detrending = None
        self.pycheops_params ={}
        self.detrend_model = {}

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

        time = np.array(dataset.ancillary["time"])
        phi = np.array(dataset.ancillary["roll_angle"])*constants.deg2rad
        try:
            smear = dataset.ancillary["smear"]
        except KeyError:
            smear = np.zeros_like(time)
        try:
            deltaT = dataset.ancillary["deltaT"]
        except KeyError:
            deltaT = np.zeros_like(time)

        self.detrend_model[dataset.name_ref]= FM(
            dx     = make_interp(time, dataset.ancillary["xoff"], scale="range"),
            dy     = make_interp(time, dataset.ancillary["yoff"], scale="range"),
            sinphi = make_interp(time, np.sin(phi)),
            cosphi = make_interp(time, np.cos(phi)),
            bg     = make_interp(time, dataset.ancillary["bg"], scale="max"),
            contam = make_interp(time, dataset.ancillary["contam"], scale="max"),
            smear  = make_interp(time, smear, scale="max"),
            deltaT = make_interp(time, deltaT)
        )


    def compute(self, variable_value, dataset, x0_input=None):

        for p in self.list_pams_dataset:
            self.pycheops_params[dataset.name_ref][p].value = variable_value[p]

        if x0_input is None:
            return self.detrend_model[dataset.name_ref].eval(self.pycheops_params[dataset.name_ref], t=dataset.x) # full detrend
        else:
            return self.detrend_model[dataset.name_ref].eval(self.pycheops_params[dataset.name_ref], t=x0_input + dataset.Tref) # full detrend

