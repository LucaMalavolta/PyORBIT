
from pyorbit.subroutines.common import np, constants
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit

try:
    from pytransit import QuadraticModel
    from pytransit import RoadRunnerModel
    from pytransit import QPower2Model
except ImportError:
    pass


class PyTransit_Transit_TTV_Ancillary(AbstractModel, AbstractTransit):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            from pytransit import QuadraticModel
        except ImportError:
            print("ERROR: PyTransit not installed, this will not work")
            quit()

        # Must be moved here because it will updated depending on the selected limb darkening
        self.list_pams_common = {
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'R_Rs',  # planet radius (in units of stellar radii)
        }

        self.pytransit_models = {}
        self.pytransit_plot = {}
        self.Tc_names = {}
        self.Tc_array = {}

    def initialize_model(self, mc, **kwargs):
        """ Force the use of the central time of transit"""
        self.use_time_of_transit = True

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        self.list_pams_common.discard('Tc')

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        """ Reading some code-specific keywords from the configuration file"""
        self._prepare_dataset_options(mc, dataset, **kwargs)

        planet_selection = (dataset.ancillary_str['planet'] == self.planet_ref )
        self.Tc_number = int(np.sum(planet_selection))

        """ extracting the central time of transit and the transit duration
        for the planet under analysis"""
        Tc_selected = dataset.ancillary['transit_time'][planet_selection]
        Td_selected = dataset.ancillary['duration'][planet_selection]
        try:
            transit_id = dataset.ancillary['transit_id'][planet_selection]
        except:
            if self.use_shared_ttvs:
                print('ERROR: transit_id should be provided in the ancillary file')
                print('     : when flag:use_shared_ttvs is activated ')
                quit()
            transit_id = np.arange(0, self.Tc_number, dtype=np.int64)

        self.Tc_names[dataset.name_ref] = []
        self.Tc_array[dataset.name_ref] = np.zeros(self.Tc_number, dtype=np.double)

        for i_sub in range(0, self.Tc_number):
            par_original = 'Tc'

            id_sub = transit_id[i_sub]
            par_subset = 'Tc_'+repr(id_sub)
            self.Tc_names[dataset.name_ref].append(par_subset)

            par_update = [Tc_selected[i_sub]-Td_selected[i_sub]/2., Tc_selected[i_sub]+Td_selected[i_sub]/2.]

            if self.use_shared_ttvs:
                self.transfer_parameter_properties(mc, dataset, par_original, par_subset, common_pam=True)
                mc.common_models[self.planet_ref].bounds.update({par_subset: par_update})

            else:
                self.transfer_parameter_properties(mc, dataset, par_original, par_subset, dataset_pam=True)
                self.bounds[dataset.name_ref].update({par_subset: par_update})

        if self.limb_darkening_model == 'quadratic':
            self.pytransit_models[dataset.name_ref] = QuadraticModel()
            self.pytransit_plot[dataset.name_ref] = QuadraticModel()

        self.pytransit_models[dataset.name_ref].set_data(dataset.x0, lcids=dataset.submodel_id.astype(int), epids=transit_id.astype(int),
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

        for i_tc, n_tc in enumerate(self.Tc_names[dataset.name_ref]):
            self.Tc_array[dataset.name_ref][i_tc] = parameter_values[n_tc] - dataset.Tref

        if x0_input is None:
            return self.pytransit_models[dataset.name_ref].evaluate_ps(
                parameter_values['R_Rs'],
                self.ld_vars,
                self.Tc_array[dataset.name_ref],
                parameter_values['P'],
                parameter_values['a_Rs'],
                parameter_values['i'],
                parameter_values['e'],
                parameter_values['omega'] * constants.deg2rad) - 1.

        else:
            self.pytransit_plot[dataset.name_ref].set_data(x0_input,
                                                            exptimes=self.code_options[dataset.name_ref]['exp_time'],
                                                            nsamples=self.code_options[dataset.name_ref]['sample_factor'])

            return self.pytransit_plot[dataset.name_ref].evaluate_ps(
                parameter_values['R_Rs'],
                self.ld_vars,
                self.Tc_array[dataset.name_ref],
                parameter_values['P'],
                parameter_values['a_Rs'],
                parameter_values['i'],
                parameter_values['e'],
                parameter_values['omega'] * constants.deg2rad) - 1.