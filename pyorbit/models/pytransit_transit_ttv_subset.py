
from pyorbit.subroutines.common import np, constants, OrderedSet
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit

try:
    from pytransit import QuadraticModel
    from pytransit import RoadRunnerModel
    from pytransit import QPower2Model
except ImportError:
    pass


class PyTransit_Transit_TTV_Subset(AbstractModel, AbstractTransit):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            from pytransit import QuadraticModel
            from pytransit import RoadRunnerModel
            from pytransit import QPower2Model
        except ImportError:
            print("ERROR: PyTransit not installed, this will not work")
            quit()

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = OrderedSet([
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'R_Rs',  # planet radius (in units of stellar radii)
        ])

        self.pytransit_models = {}
        self.pytransit_plot = {}


        self.Tc_number = {}
        self.Tc_names = {}
        self.Tc_arrays = {}


    def initialize_model(self, mc, **kwargs):
        """ Force the use of the time of inferior conjunction"""
        mc.common_models[self.planet_ref].use_time_inferior_conjunction = True

        self.use_roadrunner = kwargs.get('use_roadrunner', True)
        if self.use_roadrunner:
            print('Using RoadRunner Model from PyTransit')

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        self.list_pams_common.discard('Tc')

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        """ Reading some code-specific keywords from the configuration file"""
        self._prepare_dataset_options(mc, dataset, **kwargs)

        #TODO remove in version 11
        try:
            self.start_flag = dataset.submodel_minflag
            self.end_flag = dataset.submodel_maxflag
        except AttributeError:
            self.start_flag = 0
            self.end_flag = dataset.submodel_flag

        subset_flag = np.zeros_like(dataset.x, dtype=int) + self.end_flag  + 5

        self.Tc_number[dataset.name_ref] = []
        self.Tc_names[dataset.name_ref] = []

        transit_index = 0
        for i_sub in range(self.start_flag, self.end_flag):

            par_original = 'Tc'
            par_subset = 'Tc_'+repr(i_sub)

            if np.amin(np.abs(dataset.submodel_id-i_sub)) > 0.5: continue

            self.Tc_names[dataset.name_ref].append(par_subset)
            self.Tc_number[dataset.name_ref].append(i_sub)

            subset_flag[(dataset.submodel_id == i_sub)] = transit_index

            sub_dataset = dataset.x[(dataset.submodel_id == i_sub)]

            if kwargs[dataset.name_ref].get('boundaries', False):
                par_update = kwargs[dataset.name_ref]['boundaries'].get(
                    par_subset, [min(sub_dataset), max(sub_dataset)])
            elif kwargs.get('boundaries', False):
                par_update = kwargs['boundaries'].get(par_subset, [min(sub_dataset), max(sub_dataset)])
            else:
                par_update = [min(sub_dataset), max(sub_dataset)]

            if self.use_shared_ttvs:
                self.transfer_parameter_properties(mc, dataset, par_original, par_subset, keywords=kwargs, common_pam=True)
                mc.common_models[self.planet_ref].bounds.update({par_subset: par_update})

            else:
                self.transfer_parameter_properties(mc, dataset, par_original, par_subset, keywords=kwargs, dataset_pam=True)
                self.bounds[dataset.name_ref].update({par_subset: par_update})

            if self.use_roadrunner:
                self.pytransit_models[dataset.name_ref+ '_'+repr(i_sub)] = RoadRunnerModel(self.limb_darkening_model)
                self.pytransit_plot[dataset.name_ref+ '_'+repr(i_sub)] = RoadRunnerModel(self.limb_darkening_model)
            elif self.limb_darkening_model == 'quadratic':
                self.pytransit_models[dataset.name_ref+ '_'+repr(i_sub)] = QuadraticModel()
                self.pytransit_plot[dataset.name_ref+ '_'+repr(i_sub)] = QuadraticModel()


            self.pytransit_models[dataset.name_ref+ '_'+repr(i_sub)].set_data(sub_dataset-dataset.Tref,
                                                            exptimes=self.code_options[dataset.name_ref]['exp_time'],
                                                            nsamples=self.code_options[dataset.name_ref]['sample_factor'])


            transit_index += 1

        transit_id = np.arange(0, transit_index, dtype=int)
        #print(subset_flag)
        #print(transit_id)

        #if self.use_roadrunner:
        #    self.pytransit_models[dataset.name_ref] = RoadRunnerModel(self.limb_darkening_model)
        #    self.pytransit_plot[dataset.name_ref] = RoadRunnerModel(self.limb_darkening_model)
        #elif self.limb_darkening_model == 'quadratic':
        #    self.pytransit_models[dataset.name_ref] = QuadraticModel()
        #    self.pytransit_plot[dataset.name_ref] = QuadraticModel()

        #self.pytransit_models[dataset.name_ref].set_data(dataset.x0,
        #                                                    lcids=subset_flag, epids=transit_id,
        #                                                    exptimes=self.code_options[dataset.name_ref]['exp_time'],
        #                                                    nsamples=self.code_options[dataset.name_ref]['sample_factor'])


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

        for par, i_par in self.ldvars.items():
            self.ld_vars[i_par] = parameter_values[par]

        if x0_input is None:
            y_output = np.zeros(dataset.n)
        else:
            y_output = x0_input * 0.

        for i_tc, n_tc in enumerate(self.Tc_names[dataset.name_ref]):

            i_sub = self.Tc_number[dataset.name_ref][i_tc]
            par_subset = n_tc
            Tc = parameter_values[par_subset] - dataset.Tref

            if x0_input is None:
                sel_data = (dataset.submodel_id==i_sub)
                y_output[sel_data] = self.pytransit_models[dataset.name_ref+ '_'+repr(i_sub)].evaluate(
                    parameter_values['R_Rs'],
                    self.ld_vars,
                    Tc,
                    parameter_values['P'],
                    parameter_values['a_Rs'],
                    parameter_values['i'],
                    parameter_values['e'],
                    parameter_values['omega'] * constants.deg2rad) - 1.

            else:
                original_dataset = dataset.x0[(dataset.submodel_id==i_sub)]
                sel_data = (x0_input >= np.amin(original_dataset)) &  (x0_input <= np.amax(original_dataset))

                self.pytransit_plot[dataset.name_ref+ '_'+repr(i_sub)].set_data(x0_input[sel_data],
                                                                exptimes=self.code_options[dataset.name_ref]['exp_time'],
                                                                nsamples=self.code_options[dataset.name_ref]['sample_factor'])

                y_output[sel_data] = self.pytransit_plot[dataset.name_ref+ '_'+repr(i_sub)].evaluate(
                    parameter_values['R_Rs'],
                    self.ld_vars,
                    Tc,
                    parameter_values['P'],
                    parameter_values['a_Rs'],
                    parameter_values['i'],
                    parameter_values['e'],
                    parameter_values['omega'] * constants.deg2rad) - 1.


        return y_output


        """ 
        Tc_array = []
        for n_tc in self.Tc_names[dataset.name_ref]:
            Tc_array.append(parameter_values[n_tc] - dataset.Tref)

        print('aaaaaaaa: ', Tc_array)

        if x0_input is None:
            return self.pytransit_models[dataset.name_ref].evaluate(
                parameter_values['R_Rs'],
                self.ld_vars,
                Tc_array,
                parameter_values['P'],
                parameter_values['a_Rs'],
                parameter_values['i'],
                parameter_values['e'],
                parameter_values['omega'] * constants.deg2rad) - 1.

        else:
            self.pytransit_plot[dataset.name_ref].set_data(x0_input,
                                                            exptimes=self.code_options[dataset.name_ref]['exp_time'],
                                                            nsamples=self.code_options[dataset.name_ref]['sample_factor'])

            return self.pytransit_plot[dataset.name_ref].evaluate(
                parameter_values['R_Rs'],
                self.ld_vars,
                Tc_array,
                parameter_values['P'],
                parameter_values['a_Rs'],
                parameter_values['i'],
                parameter_values['e'],
                parameter_values['omega'] * constants.deg2rad) - 1.
        """
