from pyorbit.subroutines.common import *
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_gaussian_processes import AbstractGaussianProcesses
from pyorbit.keywords_definitions import *

from scipy.linalg import cho_factor, cho_solve, lapack, LinAlgError
from scipy import matrix, spatial

try:
    import pyaneti
except:
    pass

class GP_Pyaneti_QuasiPeriodicActivity(AbstractModel, AbstractGaussianProcesses):
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
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            import pyaneti
        except (ModuleNotFoundError,ImportError):
            print("ERROR: pyaneti not installed, this will not work")
            quit()

        self.model_class = 'multidimensional_gaussian_process'

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

        ### NEW Addded in PyORBIT v10.10
        self._dataset_rvflag_dict = {}

    def initialize_model(self, mc,  **kwargs):

        self._prepare_hyperparameter_conditions(mc, **kwargs)
        self._prepare_rotation_replacement(mc, **kwargs)
        self._prepare_decay_replacement(mc, **kwargs)

    def initialize_model_dataset(self, mc, dataset, **kwargs):

        """ when reloading the .p files, the object is not reinitialized, so we have to skip the
        incremental addition of datasets if they are already present  """
        if dataset.name_ref in self._dataset_names:
            return

        ## NEW Addded in PyORBIT v10.10
        if dataset.kind == 'RV':
            self._dataset_rvflag_dict[dataset.name_ref] = True
        else:
            self._dataset_rvflag_dict[dataset.name_ref] = False

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

        self._set_derivative_option(mc, dataset, **kwargs) 

        return


    def add_internal_dataset(self, parameter_values, dataset):

        self.update_parameter_values(parameter_values)

        self.internal_parameter_values = parameter_values

        d_ind = self._dataset_names[dataset.name_ref]
        d_nstart, d_nend = self._dataset_nindex[d_ind]

        self._dataset_ej2[d_nstart:d_nend] = self._dataset_e2[d_nstart:d_nend] + dataset.jitter**2.0
        self._dataset_res[d_nstart:d_nend] = dataset.residuals

        self.internal_coefficients[d_ind] = [parameter_values['con_amp'], parameter_values['rot_amp']]

    def lnlk_compute(self):

        pass_conditions = self.check_hyperparameter_values(self.internal_parameter_values)
        if not pass_conditions:
            return pass_conditions
        
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

        return output[0]


    ## NEW Addded in PyORBIT v10.10
    def lnlk_rvonly_compute(self):

        rv_loglike = 0.0

        kernel_name = 'MQ' + repr(self._added_datasets)
        cov_matrix = pyaneti.covfunc(kernel_name,self._gp_pams,self._dataset_x0,self._dataset_x0)

        Ks = cov_matrix - np.diag(self._dataset_ej2)

        alpha = cho_solve(cho_factor(cov_matrix), self._dataset_res)
        mu = np.dot(Ks, alpha).flatten()

        for dataset_name, d_ind in self._dataset_names.items():

            l_nstart, l_nend = self._dataset_nindex[d_ind]

            env = 1. / self._dataset_ej2[l_nstart:l_nend]
            residuals = (self._dataset_res[l_nstart:l_nend] - mu[l_nstart:l_nend])
            n = len(residuals)

            rv_loglike += -0.5 * (n * np.log(2 * np.pi) + np.sum(residuals ** 2 * env - np.log(env)))

        return rv_loglike



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
