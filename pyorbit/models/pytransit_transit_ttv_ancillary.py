
from pyorbit.subroutines.common import np, constants, OrderedSet
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
        self.list_pams_common = OrderedSet([
            'P',  # Period, log-uniform prior
            'e',  # eccentricity, uniform prior
            'omega',  # argument of pericenter (in radians)
            'R_Rs',  # planet radius (in units of stellar radii)
        ])

        self.pytransit_models = {}
        self.pytransit_plot = {}
        self.Tc_names = {}
        self.Tc_array = {}

        self.subset_flag = {}

    def initialize_model(self, mc, **kwargs):
        """ Force the use of the time of inferior conjunction"""
        mc.common_models[self.planet_ref].use_time_inferior_conjunction = True

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        self.list_pams_common.discard('Tc')

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        """ Reading some code-specific keywords from the configuration file"""
        self._prepare_dataset_options(mc, dataset, **kwargs)

        self.subset_flag[dataset.name_ref] = np.zeros_like(dataset.x) - 1.

        planet_selection = (dataset.ancillary_str['planet'] == self.planet_ref )
        number_of_transit_times = int(np.sum(planet_selection))

        self.Tc_names[dataset.name_ref] = []
        self.Tc_names[dataset.name_ref] = []


        """ extracting the central time of transit and the transit duration
        for the planet under analysis"""

        transit_time = dataset.ancillary['transit_time'][planet_selection]
        transit_id = int(dataset.ancillary['transit_id'][planet_selection])
        transit_window = dataset.ancillary['transit_window'][planet_selection]

        for i_sub in range(0, number_of_transit_times):

            min_distance = np.abs(dataset.x - transit_time[i_sub]).min()
            if min_distance > transit_window[i_sub]: continue

            par_original = 'Tc'

            id_sub = transit_id[i_sub]
            par_subset = 'Tc_'+repr(id_sub)
            self.Tc_names[dataset.name_ref].append(par_subset)

            par_update = [transit_time[i_sub]-transit_window[i_sub]/2., transit_time[i_sub]+transit_window[i_sub]/2.]

            if self.use_shared_ttvs:
                self.transfer_parameter_properties(mc, dataset, par_original, par_subset, common_pam=True)
                mc.common_models[self.planet_ref].bounds.update({par_subset: par_update})

            else:
                self.transfer_parameter_properties(mc, dataset, par_original, par_subset, dataset_pam=True)
                self.bounds[dataset.name_ref].update({par_subset: par_update})

            times_sel = (np.abs(dataset.x - transit_time[i_sub]) < transit_window[i_sub]/2.)
            self.subset_flag[dataset.name_ref][time_sel] = i_sub


        print()

        quit()


        try:
            Tc_selected = dataset.ancillary['transit_time'][planet_selection]
            try:
                Td_selected = dataset.ancillary['transit_window'][planet_selection]
            except:
                Td_selected = dataset.ancillary['duration'][planet_selection]

            try:
                transit_id = dataset.ancillary['transit_id'][planet_selection]
            except:
                if self.use_shared_ttvs:
                    print('ERROR: transit_id should be provided in the ancillary file')
                    print('     : when flag:use_shared_ttvs is activated ')
                    quit()
                transit_id = np.arange(0, self.Tc_number, dtype=np.int64)
        except IndexError:
            print('AAAAAAAAAAAAAAAAAAAAAAAAAAA')
            quit()

        self.Tc_names[dataset.name_ref] = []
        self.Tc_array[dataset.name_ref] = np.zeros(self.Tc_number, dtype=np.int64)

        for i_sub in range(0, self.Tc_number):
            par_original = 'Tc'

            id_sub = transit_id[i_sub]
            par_subset = 'Tc_'+repr(id_sub)
            self.Tc_names[dataset.name_ref].append(par_subset)

            print(par_subset)
            print(self.Tc_number)
            print('BBBBBBBBBBBBBBBBBB')
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

        quit()

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
