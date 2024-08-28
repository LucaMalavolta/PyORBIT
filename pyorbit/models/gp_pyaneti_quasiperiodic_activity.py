from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import *
from pyorbit.keywords_definitions import *

from scipy.linalg import cho_factor, cho_solve, lapack, LinAlgError
from scipy import matrix, spatial

try:
    import pyaneti
except:
    pass

class GP_Pyaneti_QuasiPeriodicActivity(AbstractModel):
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

        try:
            import pyaneti
        except (ModuleNotFoundError,ImportError):
            print("ERROR: pyaneti not installed, this will not work")
            quit()

        self.model_class = 'gp_pyaneti_quasiperiodic_activity'

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


        self.internal_parameter_values = None
        self._dist_t1 = None
        self._dist_t2 = None
        self._added_datasets = 0
        self.dataset_ordering = {}
        self.inds_cache = None

        self._original_x0 = []
        self._dataset_x0 = []
        self._dataset_e2 = []
        self._dataset_names = {}

        self._dataset_nindex = []

        self.use_derivative_dict = {}

        self.internal_coefficients = []

        self._dataset_ej2 = []
        self._dataset_res = []

        self._added_datasets = 0
        self._n_cov_matrix = 0

        self.pi2 = np.pi * np.pi


    def initialize_model(self, mc,  **kwargs):
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

        self._dataset_nindex.append([self._n_cov_matrix,
                                     self._n_cov_matrix+dataset.n])
        self.internal_coefficients.append([0.00, 0.00])

        self._original_x0.append(dataset.x0)

        self._dataset_x0 = np.append(self._dataset_x0, dataset.x0)
        self._dataset_e2 = np.append(self._dataset_e2, dataset.e**2)

        self._dataset_names[dataset.name_ref] = self._added_datasets
        self._n_cov_matrix += dataset.n
        self._added_datasets += 1

        self._dataset_ej2 = self._dataset_e2 * 1.
        self._dataset_res = self._dataset_e2 * 0.
        self._fake_jitter = np.zeros(self._added_datasets)


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

    ## WHICH ONE SHOULD I KEEP???
    #def common_initialization_with_dataset(self, dataset):
    #    self.define_kernel(dataset)
    #    return

    #def define_kernel(self):

    #    # Prot, Pdec, Oamp
    #    return

    def add_internal_dataset(self, parameter_values, dataset):

        if self.use_stellar_rotation_period:
            parameter_values['Prot'] = parameter_values['rotation_period']

        if self.use_stellar_activity_decay:
            parameter_values['Pdec'] = parameter_values['activity_decay']

        self.internal_parameter_values = parameter_values

        d_ind = self._dataset_names[dataset.name_ref]
        d_nstart, d_nend = self._dataset_nindex[d_ind]

        self._dataset_ej2[d_nstart:d_nend] = self._dataset_e2[d_nstart:d_nend] + dataset.jitter**2.0
        self._dataset_res[d_nstart:d_nend] = dataset.residuals

        self.internal_coefficients[d_ind] = [parameter_values['con_amp'], parameter_values['rot_amp']]

    def lnlk_compute(self):
        if not self.hyper_condition(self.internal_parameter_values):
            return -np.inf
        if not self.rotdec_condition(self.internal_parameter_values):
            return -np.inf

        self._gp_pams = np.empty(2*self._added_datasets + 3)

        for l_dataset in range(0, self._added_datasets):
            self._gp_pams[2*l_dataset], self._gp_pams[2*l_dataset+1], = self.internal_coefficients[l_dataset]

        self._gp_pams[-3] = self.internal_parameter_values['Pdec']
        self._gp_pams[-2] = self.internal_parameter_values['Oamp']
        self._gp_pams[-1] = self.internal_parameter_values['Prot']

        kernel_name = 'MQ' + repr(self._added_datasets)
        fake_ljitter = np.zeros(self._n_cov_matrix)
        fake_jitter = np.zeros(1)


        output = pyaneti.nll_gp(self._gp_pams,
                                kernel_name,
                                self._dataset_x0,
                                self._dataset_res,
                                np.sqrt(self._dataset_ej2),
                                fake_jitter, fake_ljitter)
        #print(output)
        #cov_matrix = self._compute_cov_matrix()

        return output[0]

        #cov_matrix = self._compute_cov_matrix(add_diagonal_errors=True)
        #chi2 = np.dot(_3res,np.matmul(inv_M,_3res))
        #
        # try:
        #    alpha = cho_solve(cho_factor(cov_matrix), self._3res)
        #    (s, d) = np.linalg.slogdet(cov_matrix)
        #
        #    return -0.5 * (self.n * np.log(2 * np.pi) +
        #                   np.dot(self._3res, alpha) + d)
        # except:
        #    return -np.inf

    def sample_predict(self, dataset, x0_input=None, return_covariance=False, return_variance=False):

        dataset_index = self._dataset_names[dataset.name_ref]

        if x0_input is None:
            faster_computation = False
            t_predict = dataset.x0
            l_nstart, l_nend = self._dataset_nindex[dataset_index]
        else:
            faster_computation = True
            t_predict = x0_input
            l_nstart, l_nend = len(x0_input)*dataset_index, len(x0_input)*(dataset_index+1)

        kernel_name = 'MQ' + repr(self._added_datasets)
        cov_matrix = pyaneti.covfunc(kernel_name,self._gp_pams,self._dataset_x0,self._dataset_x0)

        if faster_computation:
            Ks = self._compute_cov_Ks(t_predict)
        else:
            Ks = cov_matrix - np.diag(self._dataset_ej2)

        alpha = cho_solve(cho_factor(cov_matrix), self._dataset_res)
        mu = np.dot(Ks, alpha).flatten()
        (s, d) = np.linalg.slogdet(cov_matrix)

        B = cho_solve(cho_factor(cov_matrix), Ks.T)
        KsB_dot_diag = np.diag(np.dot(Ks, B))

        B = None
        Ks = None
        if faster_computation:
            Kss = self._compute_cov_diag(t_predict)
        else:
            Kss = np.diag(cov_matrix)
        cov_matrix = None

        std = np.sqrt(np.array(Kss - KsB_dot_diag).flatten())

        Kss = None

        if return_covariance:
            print('Covariance matrix output not implemented - ERROR')
            quit()

        if return_variance:
            return mu[l_nstart:l_nend], std[l_nstart:l_nend]
        else:
            return mu[l_nstart:l_nend]

    def sample_conditional(self, dataset, x0_input=None):
        val, std = self.sample_predict(dataset, x0_input)
        return val

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
