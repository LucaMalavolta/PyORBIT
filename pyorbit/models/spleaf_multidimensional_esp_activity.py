from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
from pyorbit.keywords_definitions import *

from scipy.linalg import cho_factor, cho_solve, lapack, LinAlgError
from scipy import matrix, spatial
import sys

__all__ = ['SPLEAF_Multidimensional_ESP']


try:
    from spleaf import cov as spleaf_cov
    from spleaf import term as spleaf_term
except (ModuleNotFoundError,ImportError):
    pass


class SPLEAF_Multidimensional_ESP(AbstractModel):
    ''' Three parameters out of four are the same for all the datasets, since they are related to
    the properties of the physical process rather than the observed effects on a dataset
     From Grunblatt+2015, Affer+2016
     - theta: is usually related to the rotation period of the star( or one of its harmonics);
     - lambda: is the correlation decay timescale, and it can be related to the lifetime of the active regions.
     - omega: is the length scale of the periodic component, and can be linked to the size evolution of the active regions;
     - h: represents the amplitude of the correlations '''

    default_common = 'activity'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model_class = 'spleaf_multidimensional_esp'

        self.internal_likelihood = True
        self.delayed_lnlk_computation = True

        self.list_pams_common = OrderedSet([
            'Prot',  # Rotational period of the star
            'Pdec',  # Decay timescale of activity
            'Oamp',  # Granulation of activity
        ])
        self.list_pams_dataset = OrderedSet([
            'rot_amp', # Amplitude of the covariance matrix
            'con_amp' # Amplitude of the first derivative of the covariance matrix
        ])

        try:
            from spleaf import cov as spleaf_cov
            from spleaf import term as spleaf_term
        except (ModuleNotFoundError,ImportError):
            print("ERROR: S+LEAF package not installed, this will not work")
            quit()


        self.internal_parameter_values = None
        #self._dist_t1 = None
        #self._dist_t2 = None
        #self._added_datasets = 0
        #self.dataset_ordering = {}
        #self.inds_cache = None

        self._dataset_x0 = []
        self._dataset_label = []
        self._dataset_err = []
        self._dataset_names = {}

        self._jitter_mask = []
        self._dataset_jitmask = []

        self._dataset_nindex = {}
        self._dataset_njitter = {}

        self._dataset_temporary_jitmask = {}


        #self.use_derivative_dict = {}

        self.internal_coeff_prime = None
        self.internal_coeff_deriv = None

        self._dataset_ej2 = None
        self._dataset_res = None

        self._added_datasets = 0
        self._added_jitters = 0
        self.n_harmonics = 4

        self.D_spleaf = None
        self.D_param = None
        self.input_param = None


        #self.pi2 = np.pi * np.pi


    def initialize_model(self, mc,  **kwargs):

        self.n_harmonics = kwargs.get('n_harmonics', self.n_harmonics)
        print(' S+LEAF model, number of harmonics:', self.n_harmonics)
        print()

        if kwargs.get('hyperparameters_condition', False):
            self.hyper_condition = self._hypercond_01
        else:
            self.hyper_condition = self._hypercond_00

        if kwargs.get('rotation_decay_condition', False):
            self.rotdec_condition = self._hypercond_02
        else:
            self.rotdec_condition = self._hypercond_00

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'activity':
                self.use_stellar_rotation_period = getattr(mc.common_models[common_ref], 'use_stellar_rotation_period', False)
                break

        for keyword in keywords_stellar_rotation:
            self.use_stellar_rotation_period = kwargs.get(keyword, self.use_stellar_rotation_period)

        if self.use_stellar_rotation_period:
            self.list_pams_common.update(['rotation_period'])
            self.list_pams_common.discard('Prot')


        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'activity':
                self.use_stellar_activity_decay = getattr(mc.common_models[common_ref], 'use_stellar_activity_decay', False)
                break

        for keyword in keywords_stellar_activity_decay:
            self.use_stellar_activity_decay = kwargs.get(keyword, self.use_stellar_activity_decay)

        if self.use_stellar_activity_decay:
            self.list_pams_common.update(['activity_decay'])
            self.list_pams_common.discard('Pdec')

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        """ when reloading the .p files, the object is not reinitialized, so we have to skip the
        incremental addition of datasets if they are already present  """
        if dataset.name_ref in self._dataset_names:
            return

        self._dataset_x0.append(dataset.x0)
        self._dataset_err.append(dataset.e)
        self._dataset_nindex[dataset.name_ref] = self._added_datasets

        self._dataset_njitter[dataset.name_ref] = []


        self._dataset_temporary_jitmask[dataset.name_ref] = []

        for var in dataset.list_pams:
            if dataset.variable_expanded[var] != 'jitter':
                continue
            self._dataset_temporary_jitmask[dataset.name_ref].append(dataset.mask[var])
            self._dataset_njitter[dataset.name_ref].append(self._added_jitters)
            self._added_jitters += 1

        self.spleaf_time, self.spleaf_res, self.spleaf_err, self.spleaf_series_index = \
            spleaf_cov.merge_series(self._dataset_x0, self._dataset_err, self._dataset_err)


        self._jitter_mask = []
        self._dataset_jitmask = []
        for dataset_name, temporary_jitmask in self._dataset_temporary_jitmask.items():
            d_ind = self._dataset_nindex[dataset_name]
            
            temporary_mask = self.spleaf_series_index[d_ind]
            for jit_list in temporary_jitmask:
                self._dataset_jitmask.append(jit_list)
                self._jitter_mask.append(temporary_mask[jit_list])


        self._added_datasets += 1


        d_ind = self._dataset_nindex[dataset.name_ref]
        j_ind = self._dataset_njitter[dataset.name_ref]

        self.spleaf_res[self.spleaf_series_index[d_ind]] = dataset.residuals

        self.internal_parameter_values = {
            'Prot': 30.0,
            'Pdec': 100.0,
            'Oamp': 0.35
        }

        self.internal_jitter = np.ones(self._added_jitters)
        self.internal_coeff_prime = np.ones(self._added_datasets)
        self.internal_coeff_deriv = np.ones(self._added_datasets)

        self._reset_kernel()

        if 'derivative'in kwargs:
            use_derivative = kwargs['derivative'].get(dataset.name_ref, False)
        elif dataset.name_ref in kwargs:
            use_derivative = kwargs[dataset.name_ref].get('derivative', False)
        else:
            if dataset.kind == 'H-alpha' or \
                dataset.kind == 'S_index' or \
                dataset.kind == 'Ca_HK' or \
                dataset.kind == 'FWHM':
                    use_derivative = False
            else:
                use_derivative = True

        if not use_derivative:
            self.fix_list[dataset.name_ref] = {'rot_amp': [0., 0.]}

        return

    def add_internal_dataset(self, parameter_values, dataset):

        if self.use_stellar_rotation_period:
            parameter_values['Prot'] = parameter_values['rotation_period']

        if self.use_stellar_activity_decay:
            parameter_values['Pdec'] = parameter_values['activity_decay']

        self.internal_parameter_values = parameter_values

        d_ind = self._dataset_nindex[dataset.name_ref]
        j_ind = self._dataset_njitter[dataset.name_ref]

        self.spleaf_res[self.spleaf_series_index[d_ind]] = dataset.residuals

        for j_jit in j_ind:
            self.internal_jitter[j_jit] =  dataset.jitter[self._dataset_jitmask[j_jit]][0]

        self.internal_coeff_prime[d_ind] = parameter_values['con_amp']
        self.internal_coeff_deriv[d_ind] = parameter_values['rot_amp']


    def lnlk_compute(self):


        if not self.hyper_condition(self.internal_parameter_values):
            return -np.inf
        if not self.rotdec_condition(self.internal_parameter_values):
            return -np.inf
    
        input_param = np.concatenate(([self.internal_parameter_values['Prot'],
                        self.internal_parameter_values['Pdec'],
                        self.internal_parameter_values['Oamp']],
                        self.internal_coeff_prime,
                        self.internal_coeff_deriv,
                        self.internal_jitter))
        self.D_spleaf.set_param(input_param, self.D_param)

        return self.D_spleaf.loglike(self.spleaf_res)


    def sample_predict(self, dataset, x0_input=None, return_covariance=False, return_variance=False):


        input_param = np.concatenate(([self.internal_parameter_values['Prot'],
                        self.internal_parameter_values['Pdec'],
                        self.internal_parameter_values['Oamp']],
                        self.internal_coeff_prime,
                        self.internal_coeff_deriv,
                        self.internal_jitter))

        self.D_spleaf.set_param(input_param, self.D_param)
        """ I'm creating the kernel here has """

        d_ind = self._dataset_nindex[dataset.name_ref]

        self.D_spleaf.kernel['GP'].set_conditional_coef(series_id=d_ind)


        if x0_input is None:
            t_predict = dataset.x0
        else:
            t_predict = x0_input

        mu, var = self.D_spleaf.conditional(self.spleaf_res, t_predict, calc_cov='diag')

        if return_variance:
            return mu, np.sqrt(var)
        else:
            return mu

    def _reset_kernel(self):

        kwargs = {
            'err': spleaf_term.Error(self.spleaf_err),
            'GP': spleaf_term.MultiSeriesKernel(spleaf_term.ESPKernel(1.,
                                                                    self.internal_parameter_values['Prot'],
                                                                    self.internal_parameter_values['Pdec'],
                                                                    self.internal_parameter_values['Oamp'],
                                                                    nharm=self.n_harmonics),
                                            self.spleaf_series_index,
                                            self.internal_coeff_prime,
                                            self.internal_coeff_deriv)
        }
        for n_jit in range(0, self._added_jitters):
            kwargs['jitter_'+repr(n_jit)] = spleaf_term.InstrumentJitter(self._jitter_mask[n_jit], self.internal_jitter[n_jit])

        self.D_spleaf = spleaf_cov.Cov(self.spleaf_time, **kwargs)
        self.D_param = self.D_spleaf.param[1:]
        #print(self.D_param)

    @staticmethod
    def _hypercond_00(parameter_values):
        #Condition from Rajpaul 2017, Rajpaul+2021
        return True

    @staticmethod
    def _hypercond_01(parameter_values):
        # Condition from Rajpaul 2017, Rajpaul+2021
        # Taking into account that Pdec^2 = 2*lambda_2^2
        return parameter_values['Pdec']**2 > (3. / 2. / np.pi) * parameter_values['Oamp']**2 * parameter_values['Prot']**2

    @staticmethod
    def _hypercond_02(parameter_values):
        #Condition on Rotation period and decay timescale
        return parameter_values['Pdec'] > 2. * parameter_values['Prot']

