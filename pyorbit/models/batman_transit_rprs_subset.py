from pyorbit.subroutines.common import np, OrderedSet
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit

try:
    import batman
except (ModuleNotFoundError,ImportError):
    pass


class Batman_Transit_RpRs_Subset(AbstractModel, AbstractTransit):

    def __init__(self, *args, **kwargs):
        # this calls all constructors up to AbstractModel
        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            import batman
        except (ModuleNotFoundError,ImportError):
            print("ERROR: batman not installed, this will not work")
            quit()

        self.list_pams_common = OrderedSet([
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'R_Rs',  # planet radius (in units of stellar radii)
        ])

        """ Model-specifc parameters, not declared in the abstract class """
        self.batman_params = None
        self.batman_models = {}
        self.code_options = {}

        self.dataset_x0 ={}
        self.dataset_Tref = {}

    def initialize_model(self, mc, **kwargs):

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        self.code_options['nthreads'] = kwargs.get('nthreads', 1)
        try:
            import multiprocessing
            if self.code_options['nthreads'] > multiprocessing.cpu_count():
                print('Batman nthreads automatically lowered to the maximum CPU count')
                self.code_options['nthreads'] = multiprocessing.cpu_count()
        except:
            self.code_options['nthreads'] = 1

        #if not batman.openmp.detect():
        #    print('OpenMP not supported, batman nthreads automatically lowered to 1')
        #    self.code_options['nthreads'] = 1

        self.batman_params = batman.TransitParams()

        """ Initialization with random transit parameters"""
        self.batman_params.t0 = 0.  # time of inferior conjunction
        self.batman_params.per = 1.  # orbital period
        # planet radius (in units of stellar radii)
        self.batman_params.rp = 0.1
        # semi-major axis (in units of stellar radii)
        self.batman_params.a = 15.
        self.batman_params.inc = 87.  # orbital inclination (in degrees)
        self.batman_params.ecc = 0.  # eccentricity
        self.batman_params.w = 90.  # longitude of periastron (in degrees)

        """ Setting up the limb darkening calculation"""

        self.batman_params.limb_dark = kwargs['limb_darkening_model']
        self.batman_params.u = np.ones(kwargs['limb_darkening_ncoeff'],
                                       dtype=np.double) * 0.1  # limb darkening coefficients

        """ And now we remove the planetary radius from the common parameters, and add it back as a dataset-specific parameter """

        self.list_pams_common.discard('R_Rs')
        # self.list_pams_dataset.update(['Tc'])

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        """ Reading some code-specific keywords from the configuration file"""
        self._prepare_dataset_options(mc, dataset, **kwargs)

        for i_sub in range(0, dataset.submodel_flag):

            par_original = 'R_Rs'
            par_subset = 'R_Rs_'+repr(i_sub)

            self.transfer_parameter_properties(mc, dataset, par_original, par_subset, dataset_pam=True)

            sub_dataset = dataset.x[(dataset.submodel_id == i_sub)]

            if kwargs[dataset.name_ref].get('boundaries', False):
                par_update = kwargs[dataset.name_ref]['boundaries'].get(par_subset, [0.00001, 0.5])
            elif kwargs.get('boundaries', False):
                par_update = kwargs['boundaries'].get(par_subset, [0.00001, 0.5])
            else:
                par_update = [0.00001, 0.5]

            self.bounds[dataset.name_ref].update({par_subset: par_update})

            self.batman_models[dataset.name_ref + '_'+repr(i_sub)] = \
                batman.TransitModel(self.batman_params,
                                    sub_dataset,
                                    supersample_factor=self.code_options[dataset.name_ref]['sample_factor'],
                                    exp_time=self.code_options[dataset.name_ref]['exp_time'],
                                    nthreads=self.code_options['nthreads'])

    def compute(self, parameter_values, dataset, x0_input=None):
        """
        :param parameter_values:
        :param dataset:
        :param x0_input:
        :return:
        """

        self.update_parameter_values(parameter_values, dataset.Tref)

        for key, key_val in parameter_values.items():
            if np.isnan(key_val):
                return 0.

        self.batman_params.a = parameter_values['a_Rs']
        self.batman_params.inc = parameter_values['i']
        self.batman_params.t0 = parameter_values['Tc'] - dataset.Tref

        self.batman_params.per = parameter_values['P']  # orbital period
        # planet radius (in units of stellar radii)
        self.batman_params.ecc = parameter_values['e']  # eccentricity
        # longitude of periastron (in degrees)
        self.batman_params.w = parameter_values['omega']

        """
        print 'a    ', self.batman_params.a
        print 'inc  ', self.batman_params.inc
        print 't0   ', self.batman_params.t0
        print 'per  ', self.batman_params.per
        print 'rp   ', self.batman_params.rp
        print 'ecc  ', self.batman_params.ecc
        print 'w    ', self.batman_params.w
        print 'u    ', self.batman_params.u
        """
        for par, i_par in self.ldvars.items():
            self.batman_params.u[i_par] = parameter_values[par]

        """
        From the batman manual:
        Reinitializing the model is by far the slowest component of batman,
        because it calculates the optimal step size
        for the integration starting from a very small value.
        -> However, we estimated the optimal step size from random parameters,
        so at some point we'll need to
        reinitialize the model so that the correct step size is computed.
        """

        random_selector = np.random.randint(1000)

        if x0_input is None:
            y_output = np.zeros(dataset.n)
        else:
            y_output = x0_input * 0.

        for i_sub in range(0,dataset.submodel_flag):

            par_subset = 'R_Rs_'+repr(i_sub)
            self.batman_params.rp = parameter_values[par_subset]


            #if not self.use_inclination:
            #    if parameter_values['b'] > 1. + par_subset :
            #        return y_output


            if x0_input is None:
                sel_data = (dataset.submodel_id==i_sub)

                if random_selector == 50:
                    self.batman_models[dataset.name_ref + '_'+repr(i_sub)] = \
                        batman.TransitModel(self.batman_params,
                                            dataset.x0[sel_data],
                                            supersample_factor=self.code_options[
                                                dataset.name_ref]['sample_factor'],
                                            exp_time=self.code_options[dataset.name_ref]['exp_time'],
                                            nthreads=self.code_options['nthreads'])

                y_output[sel_data] = self.batman_models[dataset.name_ref+ '_'+repr(i_sub)].light_curve(self.batman_params) - 1.


            else:
                original_dataset = dataset.x0[(dataset.submodel_id==i_sub)]
                sel_data = (x0_input >= np.amin(original_dataset)) &  (x0_input <= np.amax(original_dataset))

                temporary_model = batman.TransitModel(self.batman_params,
                                            x0_input[sel_data],
                                            supersample_factor=self.code_options[
                                                dataset.name_ref]['sample_factor'],
                                            exp_time=self.code_options[dataset.name_ref]['exp_time'],
                                            nthreads=self.code_options['nthreads'])

                y_output[sel_data] = temporary_model.light_curve(self.batman_params) - 1.

        return y_output
