from pyorbit.subroutines.common import np, OrderedSet, sys, os
import pyorbit.subroutines.constants as constants
import pyorbit.subroutines.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit

try:
    import pyarome
except (ModuleNotFoundError,ImportError):
    pass

class RossiterMcLaughlin_Pyarome(AbstractModel, AbstractTransit):
    model_class = 'rossiter_mclaughlin'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)  # this calls all constructors up to AbstractModel
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            import pyarome
            #https://github.com/andres-jordan/PyARoME
        except (ModuleNotFoundError,ImportError):
            print("ERROR: arome/pyarome not installed, this will not work")
            quit()

        self.unitary_model = False

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = OrderedSet([
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'lambda', # Sky-projected angle between stellar rotation axis and normal of orbit plane [deg]
            'R_Rs',  # planet radius (in units of stellar radii)
        ])

        self.arome_parameters = {}
        self.model_class = 'rossiter_mclaughlin'

    def initialize_model(self, mc, **kwargs):

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_star_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        print('    Note: assumption of quadratic limb darkening for RML computation')

        self.arome_parameters['macroturbulence'] =  kwargs.get('macroturbulence', 1.0)
        self.arome_parameters['instrumental_broadening'] =  kwargs.get('instrumental_broadening', 1.0)
        self.arome_parameters['measurement_technique'] =  kwargs.get('measurement_technique', 'ccf')
        if self.arome_parameters['measurement_technique'] not in ['ccf', 'iodine']:
            raise ValueError('{0:s} error: measurement_technique must be either "ccf" or "iodine"'.format(self.model_name))

        try:
            self.arome_parameters['measured_ccf_width'] =  kwargs.get('measured_ccf_width') / constants.sigma2FWHM
        except (KeyError, TypeError):
            raise ValueError('{0:s} error: must provide measured_ccf_width (FWHM) in km/s'.format(self.model_name))

        print("    {0:s} global parameters:".format(self.model_name))
        for key, value in self.arome_parameters.items():
            print(f"        {key}: {value}")


    def compute(self, parameter_values, dataset, x0_input=None):
        """
        :param parameter_values:
        :param dataset:
        :param x0_input:
        :return:
        """
        #t1_start = process_time()

        self.update_parameter_values(parameter_values, dataset.Tref)

        for key, key_val in parameter_values.items():
            if np.isnan(key_val):
                return 0.

        ld_par = self._limb_darkening_coefficients(parameter_values)


        if x0_input is not None:
            x0 = x0_input
        else:
            x0 = dataset.x0

        try:
            rv_model = pyarome.RM(parameter_values['lambda'],
                            x0,
                            parameter_values['P'],
                            parameter_values['i'],
                            parameter_values['e'],
                            parameter_values['omega'],
                            parameter_values['Tc']-dataset.Tref,
                            parameter_values['a_Rs'],
                            ld_par[0],
                            ld_par[1],
                            self.arome_parameters['instrumental_broadening'],
                            parameter_values['v_sini'],
                            self.arome_parameters['measured_ccf_width'],
                            self.arome_parameters['macroturbulence'],
                            6,
                            parameter_values['R_Rs'])
        except RuntimeError:
            return -np.inf

        #print('RV model', rv_model[2], rv_model[0]*1000.)

        if rv_model[2]!=0:
            return -np.inf

        if self.arome_parameters['measurement_technique'] == 'ccf':
            return rv_model[0]*1000.
        elif self.arome_parameters['measurement_technique'] == 'iodine':
            return rv_model[1]*1000.
