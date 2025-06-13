
from pyorbit.subroutines.common import np, constants, OrderedSet
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.abstract_transit import AbstractTransit

try:
    from pytransit import QuadraticModel
    from pytransit import RoadRunnerModel
    from pytransit import QPower2Model
except (ModuleNotFoundError,ImportError):
    pass


class PyTransit_Dynamical(AbstractModel, AbstractTransit):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)
        super(AbstractModel, self).__init__(*args, **kwargs)

        try:
            from pytransit import QuadraticModel
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

    def initialize_model(self, mc, **kwargs):


        self.use_roadrunner = kwargs.get('use_roadrunner', True)
        if self.use_roadrunner:
            print('Using RoadRunner Model from PyTransit')

        """ Planetary parameters initialization is taken care of by the Dynamical integration model"""
        #self._prepare_planetary_parameters(mc, **kwargs)


        self._prepare_limb_darkening_coefficients(mc, **kwargs)

    def initialize_model_dataset(self, mc, dataset, **kwargs):
        self._prepare_dataset_options(mc, dataset, **kwargs)

        if self.use_roadrunner:
            self.pytransit_models[dataset.name_ref] = RoadRunnerModel(self.limb_darkening_model)
            self.pytransit_plot[dataset.name_ref] = RoadRunnerModel(self.limb_darkening_model)
        elif self.limb_darkening_model == 'quadratic':
            self.pytransit_models[dataset.name_ref] = QuadraticModel()
            self.pytransit_plot[dataset.name_ref] = QuadraticModel()

        self.pytransit_models[dataset.name_ref].set_data(dataset.x0,
                                                            exptimes=self.code_options[dataset.name_ref]['exp_time'],
                                                            nsamples=self.code_options[dataset.name_ref]['sample_factor'])


    def compute_dynamical(self, dataset,  parameter_values, transits, durations,
                          rp_rs, per, aRs, inc, ecc, w,
                          x0_input=None):

        ### TODO this one will run inside the dynamical model
        if x0_input is None:

            t = dataset.x
            t = np.asarray(dataset.x) # for convenience
        else:
            t = np.asarray(x0_input + dataset.Tref) # add the reference time to the input time

        f0 = np.zeros_like(t) # create a model at 0.0 for all time points

        """ compute half duration in days for all the transits of all the planets"""
        half_dur_d = durations * constants.min2day

        """ select transit times of all planets in the time range """
        tra_in_t = np.logical_and(transits >= t.min(), transits <= t.max())
        n_tra = np.sum(tra_in_t)


        self.update_parameter_values(parameter_values, dataset.Tref)

        if parameter_values['i'] == 0.0:
            return 0.

        for par, i_par in self.ldvars.items():
            self.ld_vars[i_par] = parameter_values[par]



        # select partial transits of all planets in the time range
        tra_dur_in_t = np.logical_and(
            transits - half_dur_d >= t.min(),
            transits + half_dur_d <= t.max(),
        )
        n_dur = np.sum(tra_dur_in_t)
        # number of events based on the max between n_tra and n_dur
        if n_tra >= n_dur:
            n = n_tra
            sel_tra = tra_in_t
        else:
            n = n_dur
            sel_tra = tra_dur_in_t

        # compute transit model only if full or partial transits have been found
        if n > 0:
            # select transits and parameters
            tra_sel = np.atleast_1d(transits[sel_tra])
            rp_rs_sel = np.atleast_1d(rp_rs[sel_tra])
            per_sel = np.atleast_1d(per[sel_tra])
            aRs_sel = np.atleast_1d(aRs[sel_tra])
            inc_sel = np.atleast_1d(inc[sel_tra])
            ecc_sel = np.atleast_1d(ecc[sel_tra])
            w_sel = np.atleast_1d(w[sel_tra])

            flux_ = []
            for itra, tra in enumerate(tra_sel): # loop in the selected transits (independent of the body)
                ff = f0.copy()
                sel_t = np.logical_and(
                    t >= tra - 0.5 * per[itra],
                    t <= tra + 0.5 * per[itra],
                ) # select portion of the light curve centered on the transit time that cover a full period
                self.pytransit_models[dataset.name_ref].set_data(t[sel_t],
                                                            exptimes=self.code_options[dataset.name_ref]['exp_time'],
                                                            nsamples=self.code_options[dataset.name_ref]['sample_factor']) # set the pytransit time data with oversampling if needed
                ff[sel_t] = tm.evaluate(
                    k=rp_rs_sel[itra],
                    ldc=self.ld_vars,
                    t0=tra,
                    p=per_sel[itra],
                    a=aRs_sel[itra],
                    i=inc_sel[itra],
                    e=ecc_sel[itra],
                    w=w_sel[itra],
                ) # compute the model and associate it only for the selected portion close to the transit
                flux_.append(ff) # append it
            f2d = np.atleast_2d(flux_)
            flux = np.sum(f2d - 1.0, axis=0) # in one step it removes 1, sum flux for each time point
        else: # set to 0.0 the model flux if there are not transits (full or partials) in this photometry
            flux = f0

        return flux



    def compute(self, parameter_values, dataset, x0_input=None):
        """
        :param parameter_values:
        :param dataset:
        :param x0_input:
        :return:
        """

        return dataset.external_model

        self.update_parameter_values(parameter_values, dataset.Tref)

        if parameter_values['i'] == 0.0:
            return 0.

        for par, i_par in self.ldvars.items():
            self.ld_vars[i_par] = parameter_values[par]


        if x0_input is None:
            return self.pytransit_models[dataset.name_ref].evaluate(
                parameter_values['R_Rs'],
                self.ld_vars,
                parameter_values['Tc'] - dataset.Tref,
                parameter_values['P'],
                parameter_values['a_Rs'],
                parameter_values['i']* constants.deg2rad,
                parameter_values['e'],
                parameter_values['omega'] * constants.deg2rad) - 1.

        else:
            self.pytransit_plot[dataset.name_ref].set_data(x0_input,
                                                            exptimes=self.code_options[dataset.name_ref]['exp_time'],
                                                            nsamples=self.code_options[dataset.name_ref]['sample_factor'])

            return self.pytransit_plot[dataset.name_ref].evaluate(
                parameter_values['R_Rs'],
                self.ld_vars,
                parameter_values['Tc'] - dataset.Tref,
                parameter_values['P'],
                parameter_values['a_Rs'],
                parameter_values['i']* constants.deg2rad,
                parameter_values['e'],
                parameter_values['omega'] * constants.deg2rad) - 1.
