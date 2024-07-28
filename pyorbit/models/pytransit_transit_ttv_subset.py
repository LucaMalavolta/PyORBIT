
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
        self.Tc_names = {}
        self.Tc_arrays = {}


    def initialize_model(self, mc, **kwargs):
        """ Force the use of the time of inferior conjunction"""
        mc.common_models[self.planet_ref].use_time_inferior_conjunction = True

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

        subset_flag = np.zeros_like(dataset.x, dtype=int) - 1

        self.Tc_number = int(self.end_flag)
        self.Tc_names[dataset.name_ref] = []

        transit_index = 0
        for i_sub in range(self.start_flag, self.end_flag):
            par_original = 'Tc'
            par_subset = 'Tc_'+repr(i_sub)

            if np.amin(np.abs(dataset.submodel_id-i_sub)) > 0.5: continue

            self.Tc_names[dataset.name_ref].append(par_subset)
            subset_flag[(dataset.submodel_id == i_sub)] = transit_index

            self.transfer_parameter_properties(mc, dataset, par_original, par_subset, dataset_pam=True)

            sub_dataset = dataset.x[(dataset.submodel_id == i_sub)]

            if kwargs[dataset.name_ref].get('boundaries', False):
                par_update = kwargs[dataset.name_ref]['boundaries'].get(
                    par_subset, [min(sub_dataset), max(sub_dataset)])
            elif kwargs.get('boundaries', False):
                par_update = kwargs['boundaries'].get(par_subset, [min(sub_dataset), max(sub_dataset)])
            else:
                par_update = [min(sub_dataset), max(sub_dataset)]

            self.bounds[dataset.name_ref].update({par_subset: par_update})

            transit_index += 1

        transit_id = np.arange(0, transit_index, dtype=int)

        if self.limb_darkening_model == 'quadratic':
            self.pytransit_models[dataset.name_ref] = QuadraticModel()
            self.pytransit_plot[dataset.name_ref] = QuadraticModel()

        self.pytransit_models[dataset.name_ref].set_data(dataset.x0,
                                                            lcids=subset_flag, epids=transit_id,
                                                            exptimes=self.code_options[dataset.name_ref]['exp_time'],
                                                            nsamples=self.code_options[dataset.name_ref]['sample_factor'])


    def compute(self, parameter_values, dataset, x0_input=None):
        """
        :param parameter_values:
        :param dataset:
        :param x0_input:
        :return:
        """

        self.update_parameter_values(parameter_values, dataset.Tref)

        for par, i_par in self.ldvars.items():
            self.ld_vars[i_par] = parameter_values[par]

        Tc_array = []
        for n_tc in self.Tc_names[dataset.name_ref]:
            Tc_array.append(parameter_values[n_tc] - dataset.Tref)

        if x0_input is None:
            return self.pytransit_models[dataset.name_ref].evaluate_ps(
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
