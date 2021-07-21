
from pyorbit.classes.common import np, convert_rho_to_a, convert_b_to_i
import pyorbit.classes.constants as constants
import pyorbit.classes.kepler_exo as kepler_exo
from pyorbit.models.abstract_model import AbstractModel

#from time import process_time

try:
    from pytransit import QuadraticModel
except ImportError:
    pass


class PyTransit_Transit(AbstractModel):
    model_class = 'transit'
    unitary_model = True

    default_bounds = {}
    default_spaces = {}
    default_priors = {}

    recenter_pams_dataset = {}

    def __init__(self, *args, **kwargs):

        super(PyTransit_Transit, self).__init__(*args, **kwargs)

        try:
            from pytransit import QuadraticModel
            from pytransit import RoadRunnerModel
            from pytransit import QPower2Model
        except ImportError:
            print("ERROR: PyTransit not installed, this will not work")
            quit()

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = {
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'o',  # argument of pericenter (in radians)
            'R',  # planet radius (in units of stellar radii)
        }
        self.list_pams_dataset = {}

        self.pytransit_ldvars = {}
        self.ld_ncoeff = 2
        self.parametrization = 'Standard'

        self.use_semimajor_axis = False
        self.use_inclination = False
        self.use_time_of_transit = False
        #self.nthreads = 1

        #self.transittype = 'primary'
        #self.batman_params = None
        self.pytransit_models = {}
        self.pytransit_plot = {}
        self.pytransit_options = {}

        self.limb_darkening_model = None

        self.retrieve_ai = None
        self.retrieve_t0 = None

    def initialize_model(self, mc, **kwargs):

        if mc.common_models[self.planet_ref].use_semimajor_axis:
            """ a is the semi-major axis (in units of stellar radii) """
            self.list_pams_common.update({'a': None})
            self.use_semimajor_axis = True
        else:
            """ rho is the density of the star (in solar units) """
            self.list_pams_common.update({'rho': None})

        if mc.common_models[self.planet_ref].use_inclination:
            """ i is the orbital inclination (in degrees) """
            self.list_pams_common.update({'i': None})
            self.use_inclination = True
        else:
            """ b is the impact parameter """
            self.list_pams_common.update({'b': None})

        if mc.common_models[self.planet_ref].use_time_of_transit:
            self.list_pams_common.update({'Tc': None})
            self.use_time_of_transit = True
            # Copying the property to the class for faster access
        else:
            self.list_pams_common.update({'f': None})
            # mean longitude = argument of pericenter + mean anomaly at Tref

        """ The appropriate function for variable conversion is stored internally
        """
        if self.use_semimajor_axis and self.use_inclination:
            self.retrieve_ai = self._internal_transformation_mod03
        elif self.use_semimajor_axis:
            self.retrieve_ai = self._internal_transformation_mod02
        elif self.use_inclination:
            self.retrieve_ai = self._internal_transformation_mod01
        else:
            self.retrieve_ai = self._internal_transformation_mod00

        if self.use_time_of_transit:
            self.retrieve_t0 = self._internal_transformation_mod04
        else:
            self.retrieve_t0 = self._internal_transformation_mod05

        """ Setting up the limb darkening calculation"""

        self.limb_darkening_model = kwargs['limb_darkening_model']
        self.ld_vars = [0.00] * kwargs['limb_darkening_ncoeff']
        for i_coeff in range(1, kwargs['limb_darkening_ncoeff'] + 1):
            var = 'ld_c' + repr(i_coeff)
            self.pytransit_ldvars[var] = i_coeff - 1
            self.list_pams_common.update({var: None})

    def setup_dataset(self, mc, dataset, **kwargs):

        supersample_names = ['supersample_factor',
                             'supersample',
                             'supersampling',
                             'oversample_factor',
                             'oversample',
                             'oversampling',
                             'sample_factor',
                             'sample',
                             'sampling'
                             'nsample_factor',
                             'nsample',
                             'nsampling'
                             ]

        sample_factor = 1
        exposure_time = 0.01

        for dict_name in supersample_names:
            if kwargs[dataset.name_ref].get(dict_name, False):
                sample_factor = kwargs[dataset.name_ref][dict_name]
            elif kwargs.get(dict_name, False):
                sample_factor = kwargs[dict_name]

        exptime_names = ['exposure_time',
                         'exposure',
                         'exp_time',
                         'exptime',
                         'obs_duration',
                         'integration',
                         ]

        for dict_name in exptime_names:
            if kwargs[dataset.name_ref].get(dict_name, False):
                exposure_time = kwargs[dataset.name_ref][dict_name]
            elif kwargs.get(dict_name, False):
                exposure_time = kwargs[dict_name]

        self.pytransit_options[dataset.name_ref] = {
            'sample_factor': sample_factor,
            'exp_time': exposure_time,
        }

        if self.limb_darkening_model == 'quadratic':
            self.pytransit_models[dataset.name_ref] = QuadraticModel()
            self.pytransit_plot[dataset.name_ref] = QuadraticModel()

        self.pytransit_models[dataset.name_ref].set_data(
            dataset.x0, exptimes=exposure_time, nsamples=sample_factor)

    def _internal_transformation_mod00(self, variable_value):
        """ this function transforms b and rho to i and a  """
        a = convert_rho_to_a(variable_value['P'], variable_value['rho'])
        i = convert_b_to_i(
            variable_value['b'], variable_value['e'], variable_value['o'], a)
        return a, i

    def compute(self, variable_value, dataset, x0_input=None):
        """
        :param variable_value:
        :param dataset:
        :param x0_input:
        :return:
        """
        #t1_start = process_time()

        pams_a, pams_i = self.retrieve_ai(variable_value)
        pams_t0 = self.retrieve_t0(variable_value, dataset.Tref)

        for var, i_var in self.pytransit_ldvars.items():
            self.ld_vars[i_var] = variable_value[var]

        if x0_input is None:
            ##model = self.batman_models[dataset.name_ref].light_curve(self.batman_params) - 1.
            ##t1_stop = process_time()
            ##
            ##print("Elapsed time:", t1_stop-t1_start)
            # return model
            return self.pytransit_models[dataset.name_ref].evaluate_ps(
                variable_value['R'],
                self.ld_vars,
                pams_t0, variable_value['P'], pams_a, pams_i, variable_value['e'], variable_value['o']) - 1.

        else:
            self.pytransit_plot[dataset.name_ref].set_data(x0_input,
                                                           exptimes=self.pytransit_options[dataset.name_ref]['exp_time'],
                                                           nsamples=self.pytransit_options[dataset.name_ref]['sample_factor'])

            return self.pytransit_plot[dataset.name_ref].evaluate_ps(
                variable_value['R'],
                self.ld_vars,
                pams_t0, variable_value['P'], pams_a, pams_i, variable_value['e'], variable_value['o']) - 1.


    """ function for internal transformation of variables, to avoid if calls"""
    def _internal_transformation_mod01(self, variable_value):
        """ this function transforms b to i"""
        i = convert_b_to_i(
            variable_value['b'], variable_value['e'], variable_value['o'], variable_value['a'])
        return variable_value['a'], i

    def _internal_transformation_mod02(self, variable_value):
        """ this function transforms rho to a  """
        a = convert_rho_to_a(variable_value['P'], variable_value['rho'])
        return a, variable_value['i']

    def _internal_transformation_mod03(self, variable_value):
        """ no transformation needed  """
        return variable_value['a'], variable_value['i']

    def _internal_transformation_mod04(self, variable_value, Tref):
        """ this function transforms Tc into Tc- Tref t"""
        return variable_value['Tc'] - Tref

    def _internal_transformation_mod05(self, variable_value, Tref):
        """ this function transforms Tc into Tc- Tref t"""
        return kepler_exo.kepler_phase2Tc_Tref(variable_value['P'],
                                               variable_value['f'],
                                               variable_value['e'],
                                               variable_value['o'])
