
from pyorbit.subroutines.common import np, constants, OrderedSet
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit

try:
    from pytransit import QuadraticModel
    from pytransit import RoadRunnerModel
    from pytransit import QPower2Model
except (ModuleNotFoundError,ImportError):
    pass


class PyTransit_Transit_TTV_Ancillary(AbstractModel, AbstractTransit):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            from pytransit import QuadraticModel
            from pytransit import RoadRunnerModel
            from pytransit import QPower2Model
        except (ModuleNotFoundError,ImportError):
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
        self.subset_selection = {}


    def initialize_model(self, mc, **kwargs):
        """ Force the use of the time of inferior conjunction"""
        mc.common_models[self.planet_ref].use_time_inferior_conjunction = True

        self.use_roadrunner = kwargs.get('use_roadrunner', True)
        if self.use_roadrunner:
            print('Using RoadRunner Model from PyTransit')

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        self.tc_flag_name = mc.common_models[self.planet_ref].tc_flag

        self.list_pams_common.discard('Tc')

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        """ Reading some code-specific keywords from the configuration file"""
        self._prepare_dataset_options(mc, dataset, **kwargs)

        submodel_id = dataset.ancillary[self.tc_flag_name]
        flag_sel = (submodel_id >= -0.5)
        self.start_flag = np.int64(np.amin(submodel_id[flag_sel]))
        self.end_flag = np.int64(np.amax(submodel_id[flag_sel])) + 1

        subset_flag = np.zeros_like(dataset.x, dtype=int) - 1

        self.Tc_number[dataset.name_ref] = []
        self.Tc_names[dataset.name_ref] = []

        transit_index = 0
        for i_sub in range(self.start_flag, self.end_flag):

            par_original = 'Tc'
            par_subset = 'Tc_'+repr(i_sub)

            if np.amin(np.abs(submodel_id-i_sub)) > 0.5: continue

            self.Tc_names[dataset.name_ref].append(par_subset)
            self.Tc_number[dataset.name_ref].append(i_sub)

            subset_flag[(submodel_id == i_sub)] = transit_index

            sub_dataset = dataset.x[(submodel_id == i_sub)]

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

            transit_index += 1

        transit_id = np.arange(0, transit_index, dtype=int)

        if self.use_roadrunner:
            self.pytransit_models[dataset.name_ref] = RoadRunnerModel(self.limb_darkening_model)
            self.pytransit_plot[dataset.name_ref] = RoadRunnerModel(self.limb_darkening_model)
        elif self.limb_darkening_model == 'quadratic':
            self.pytransit_models[dataset.name_ref] = QuadraticModel()
            self.pytransit_plot[dataset.name_ref] = QuadraticModel()

        if self.code_options[dataset.name_ref]['sample_factor'] == 1:
            self.code_options[dataset.name_ref]['exp_time'] = 0.

        exptimes= np.ones(transit_index) * self.code_options[dataset.name_ref]['exp_time']
        nsamples= np.ones(transit_index) * self.code_options[dataset.name_ref]['sample_factor']
        
        self.subset_selection[dataset.name_ref] = (subset_flag >= 0)

        self.pytransit_models[dataset.name_ref].set_data(dataset.x0[self.subset_selection[dataset.name_ref]],
                                                            lcids=subset_flag[self.subset_selection[dataset.name_ref]],
                                                            epids=transit_id,
                                                            exptimes=exptimes,
                                                            nsamples=nsamples)


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
        
        Tc_array = []
        for n_tc in self.Tc_names[dataset.name_ref]:
            Tc_array.append(parameter_values[n_tc] - dataset.Tref)

        if x0_input is None:
            y_output[self.subset_selection[dataset.name_ref]] = self.pytransit_models[dataset.name_ref].evaluate(
                parameter_values['R_Rs'],
                self.ld_vars,
                Tc_array,
                parameter_values['P'],
                parameter_values['a_Rs'],
                parameter_values['i'] * constants.deg2rad,
                parameter_values['e'],
                parameter_values['omega'] * constants.deg2rad) - 1.
            
        else:
            subset_flag = np.zeros_like(x0_input, dtype=int) - 1
            transit_id = np.arange(0, len(Tc_array), dtype=int)

            submodel_id = dataset.ancillary[self.tc_flag_name]

            for i_tc, n_tc in enumerate(self.Tc_number[dataset.name_ref]):
                original_dataset = dataset.x0[(submodel_id==n_tc)]
                sel_data = (x0_input >= np.amin(original_dataset)) &  (x0_input <= np.amax(original_dataset))
                subset_flag[sel_data] = i_tc
            
            subset_selection = (subset_flag >= 0)

            self.pytransit_plot[dataset.name_ref].set_data(x0_input[subset_selection],
                                                            lcids=subset_flag[subset_selection],
                                                            epids=transit_id,
                                                            exptimes=[self.code_options[dataset.name_ref]['exp_time']]*len(Tc_array),
                                                            nsamples=[self.code_options[dataset.name_ref]['sample_factor']]*len(Tc_array))

            y_output[subset_selection] = self.pytransit_plot[dataset.name_ref].evaluate(
                parameter_values['R_Rs'],
                self.ld_vars,
                Tc_array,
                parameter_values['P'],
                parameter_values['a_Rs'],
                parameter_values['i'] * constants.deg2rad,
                parameter_values['e'],
                parameter_values['omega'] * constants.deg2rad) - 1.
        
        return y_output
