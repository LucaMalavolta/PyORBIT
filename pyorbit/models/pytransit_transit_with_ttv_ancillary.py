
from pyorbit.subroutines.common import np, constants
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit

try:
    from pytransit import QuadraticModel
    from pytransit import RoadRunnerModel
    from pytransit import QPower2Model
except ImportError:
    pass


class PyTransit_Transit_With_TTV_Ancillary(AbstractModel, AbstractTransit):

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

        """ Dataset-specific time of transit boundaries are stored here"""
        self.transit_time_boundaries = {}


    def initialize_model(self, mc, **kwargs):
        """ Force the use of the central time of transit"""
        self.use_time_of_transit = True

        self._prepare_planetary_parameters(mc, **kwargs)
        self._prepare_limb_darkening_coefficients(mc, **kwargs)

        self.list_pams_common.discard('Tc')

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        """ Reading some code-specific keywords from the configuration file"""
        self._prepare_dataset_options(mc, dataset, **kwargs)

        column_selection = self.ancillary_str_index['planet']
        planet_selection = (dataset.ancillary_str[:,column_selection] == self.planet_ref )
        self.Tc_number = int(np.sum(planet_selection))

        """ extracting the central time of transit and the transit duration 
        for the planet under analysis"""
        Tc_selected = dataset.ancillary['transit_time'][planet_selection]
        Td_selected = dataset.ancillary['duration'][planet_selection]
        try:
            transit_id = dataset.ancillary['transit_id'][planet_selection]
        except:
            if self.use_shared_ttvs: 
                print('transit_id should be provided in the ancillary file')
            transit_id = np.arange(0, self.Tc_number, dtype=np.int64)

        for i_sub in range(0, self.Tc_number):
            par_original = 'Tc'

            id_sub = transit_id[i_sub]
            par_subset = 'Tc_'+repr(id_sub)

            self._subset_transfer_priors(mc, dataset, par_original, par_subset)
            par_update = [Tc_selected[i_sub]-Tc_selected[i_sub]/2., 
            Tc_selected[i_sub]+Tc_selected[i_sub]/2.]

            self.bounds[dataset.name_ref].update({par_subset: par_update})

        #!NEW 
        #!NEW 
        #!NEW 
        #!NEW 
        #!NEW 
        #!NEW 
        #!NEW 
        #!NEW 
        #!NEW 
        #!NEW 


        if self.limb_darkening_model == 'quadratic':
            self.pytransit_models[dataset.name_ref] = QuadraticModel()
            self.pytransit_plot[dataset.name_ref] = QuadraticModel()

        self.pytransit_models[dataset.name_ref].set_data(dataset.x0,
                                                            exptimes=self.code_options[dataset.name_ref]['exp_time'],
                                                            nsamples=self.code_options[dataset.name_ref]['sample_factor'])

    def define_special_parameter_properties(self,
                                            ndim,
                                            output_lists,
                                            dataset_name,
                                            par):

        if par == 'Tc' and (par not in self.bounds[dataset_name]):
            self.bounds[dataset_name][par] = self.code_options[dataset_name]['Tc_boundaries']
        return ndim, output_lists, False

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

        if x0_input is None:
            return self.pytransit_models[dataset.name_ref].evaluate_ps(
                parameter_values['R_Rs'],
                self.ld_vars,
                parameter_values['Tc'] - dataset.Tref,
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
                parameter_values['Tc'] - dataset.Tref,
                parameter_values['P'],
                parameter_values['a_Rs'],
                parameter_values['i'],
                parameter_values['e'],
                parameter_values['omega'] * constants.deg2rad) - 1.