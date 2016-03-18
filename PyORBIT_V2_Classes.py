import numpy as np
import kepler_exo as kp

# A class for common variables (number of planets, number of sinusoids...)
# Info to unpack the variables inside emcee must be included here
# Physical effects must be included here

class PlanetsCommonVariables:
    def __init__(self):
        self.n_planets = 0
        self.name_ref = []

        self.n_orbpams = {}
        self.logP_bounds = {}
        self.logK_bounds = {}
        self.e_bounds = {}

    def add_planet(self, name_ref, n_orbs):
        self.n_planets += 1
        self.name_ref.append(name_ref)
        self.n_orbpams[name_ref] = n_orbs

    # def convert_params(self,in_pams):
    #    out_pams = np.zeros(5,dtype=np.double)
    #    out_pams[0] = 10.0 ** (in_pams[0]) #P
    #    out_pams[1] = 10.0 ** (in_pams[1]) #K
    #    out_pams[2] = in_pams[2]
    #    if np.size(in_pams)==5:
    #        out_pams[3] = in_pams[3] ** 2 + in_pams[4] ** 2
    #        out_pams[4] = np.arctan2(in_pams[3], in_pams[4])
    #    return out_pams
    def compute(self, mc, theta, dataset, planet_name):
        in_pams = theta[mc.variable_list[planet_name]['kepler_pams']]
        if self.n_orbpams[planet_name] == 5:
            conv_pams = self.convert_params(in_pams[0], in_pams[1], in_pams[2], in_pams[3], in_pams[4])
        else:
            conv_pams = self.convert_params(in_pams[0], in_pams[1], in_pams[2], 0.00, 0.00)
        return self.model_kepler(conv_pams, dataset.x0)

    def print_vars(self, mc, theta):
        for planet_name in self.name_ref:
            in_pams = theta[mc.variable_list[planet_name]['kepler_pams']]
            if self.n_orbpams[planet_name] == 5:
                for param in ['logP', 'logK', 'phase', 'esino', 'ecoso']:
                    mc.pam_names[mc.variable_list[planet_name][param]] = planet_name + '_' + param
                conv_pams = self.convert_params(in_pams[0], in_pams[1], in_pams[2], in_pams[3], in_pams[4])
            else:
                for param in ['logP', 'logK', 'phase']:
                    mc.pam_names[mc.variable_list[planet_name][param]] = planet_name + '_' + param
                conv_pams = self.convert_params(in_pams[0], in_pams[1], in_pams[2], 0.00, 0.00)

            print planet_name, ' vars: ', np.asarray(theta[mc.variable_list[planet_name]['kepler_pams']])
            print planet_name, ' pams: ', conv_pams[:]

    def define_bounds(self, mc):
        for planet_name in self.name_ref:
            mc.bounds_list.append(self.logP_bounds[planet_name][:])
            mc.bounds_list.append(self.logK_bounds[planet_name][:])
            mc.bounds_list.append([0, 2 * np.pi])
            if self.n_orbpams[planet_name] == 5:
                mc.bounds_list.append([-1.0, 1.0])
                mc.bounds_list.append([-1.0, 1.0])

            mc.variable_list[planet_name] = {}
            if self.n_orbpams[planet_name] == 5:
                mc.variable_list[planet_name]['logP'] = mc.ndim
                mc.variable_list[planet_name]['logK'] = mc.ndim + 1
                mc.variable_list[planet_name]['phase'] = mc.ndim + 2
                mc.variable_list[planet_name]['esino'] = mc.ndim + 3
                mc.variable_list[planet_name]['ecoso'] = mc.ndim + 4
                mc.variable_list[planet_name]['kepler_pams'] = np.arange(mc.ndim, mc.ndim + 5, 1)
                mc.ndim += 5
            else:
                mc.variable_list[planet_name]['logP'] = mc.ndim
                mc.variable_list[planet_name]['logK'] = mc.ndim + 1
                mc.variable_list[planet_name]['phase'] = mc.ndim + 2
                mc.variable_list[planet_name]['kepler_pams'] = np.arange(mc.ndim, mc.ndim + 3, 1)
                mc.ndim += 3


    @staticmethod
    def model_kepler(orbit_pams, x0):
        P, K, phase, e, omega = orbit_pams
        rv_out = kp.kepler_RV_T0P(x0, phase, P, K, e, omega)
        return rv_out

    @staticmethod
    def convert_params(logP, logK, phase, esino, ecoso):
        P = np.exp2(logP)
        K = np.exp2(logK)
        e = np.square(ecoso) + np.square(esino)
        o = np.arctan2(esino, ecoso)
        return P, K, phase, e, o


class SinusoidsCommonVariables:
    def __init__(self):

        self.n_datasets = 0

        self.offset_coherence = True  # Same phase offset across seasons
        self.offset_common_id = -1
        self.offset_reference_name = ''
        self.use_offset = {}

        self.phase_coherence = False
        self.phase_sincro = False

        self.n_pha = 0

        self.period_sel = False
        self.n_periods = 1
        self.period_name = ['Period_0']
        self.period_list = [0., 5000000.0]
        self.period_range = np.asarray(self.period_list, dtype=np.double)

        # self.Prot_bounded = False
        # self.pha_bounded = False

        self.Prot_bounds = [0., 30.0]
        self.pof_bounds = [0., 1.0]
        self.pha_bounds = [0., 1.0]

        self.phase_list = []

    def add_period_range(self, range_input, phase_input):
        # Reset default values at first call
        if not self.period_sel:
            self.n_periods = 0
            self.period_list = []
            self.period_name = []
            self.period_sel = True

        self.period_name.append('Period_' + `self.n_periods`)
        self.period_list.append(range_input)
        self.period_range = np.asarray(self.period_list, dtype=np.double)

        self.phase_list.append(phase_input)
        self.phase = np.asarray(self.phase_list, dtype=np.int64)
        self.n_pha = np.amax(self.phase, axis=1)  # maximum value for each period
        self.n_pha_max = np.amax(self.n_pha)  # maximum value for each period

        self.n_periods += 1
        return

    # def add_phase_offset(self, dataset, period_name):
    #    # Check if the dataset needs an offset with respect to the
    #    # reference dataset (simply the first one)
    #    add_bound = False
    #    try:
    #        value = self.use_offset[dataset.kind][period_name + '_off']
    #    except KeyError:
    #        # Key is not present
    #        if self.offset_skip_first:
    #            self.use_offset[dataset.kind][period_name] = False
    #            self.offset_skip_first = False
    #        else:
    #            self.use_offset[dataset.kind][period_name] = True
    #            add_bound = True
    #    return add_bound

    def setup_model_sinusoids(self, dataset):
        dataset.n_amp = self.phase[:, dataset.ind]
        dataset.p_mask = np.zeros([self.n_periods, dataset.n], dtype=bool)
        dataset.n_periods = 0
        dataset.periods_flag = np.zeros(self.n_periods, dtype=bool)

        # If this is the first dataset, it is assumed as the reference one for the offset
        # otherwise, the offset is applied
        try:
            value = self.use_offset[dataset.kind]
        except KeyError:
            if self.offset_reference_name == '':
                self.offset_reference_name = dataset.kind
                self.use_offset[dataset.kind] = False
            else:
                self.use_offset[dataset.kind] = not self.phase_sincro #True

        for ii in xrange(0, self.n_periods):
            p_sel = (self.period_range[ii, 0] < dataset.x) & (dataset.x < self.period_range[ii, 1])
            if np.sum(p_sel) == 0:
                'No data points within the specified range'
            else:
                dataset.periods_flag[ii] = True
                dataset.p_mask[ii, :] = p_sel[:]
                dataset.n_periods += 1

        if dataset.kind == 'RV': dataset.sinamp_bounds = np.asarray([0., 100.], dtype=np.double)
        if dataset.kind == 'Phot': dataset.sinamp_bounds = np.asarray([0., 0.5], dtype=np.double)
        if dataset.kind == 'FWHM': dataset.sinamp_bounds = np.asarray([0., 2000.], dtype=np.double)
        if dataset.kind == 'BIS': dataset.sinamp_bounds = np.asarray([0., 60.], dtype=np.double)
        if dataset.kind == 'Act': dataset.sinamp_bounds = np.asarray([0., 10.], dtype=np.double)
        return

    def compute(self, mc, theta, dataset):
        # MC = Model_Container object
        # Prot and pha_in could be brought out from the llop, but I would not work
        # for plaet-only analysis
        Prot = theta[mc.variable_list['Prot']]
        pha_in = np.zeros([self.n_periods, self.n_pha_max], dtype=np.double)
        off_in = np.zeros([self.n_periods], dtype=np.double)
        amp_in = np.zeros([self.n_periods, self.n_pha_max], dtype=np.double)
        for jj in range(0, self.n_periods):
            pha_in[jj, :self.n_pha[jj]] = theta[mc.variable_list[self.period_name[jj] + '_pha']]

            if dataset.periods_flag[jj]:
                if self.use_offset[dataset.kind]:
                    off_in[:] = theta[mc.variable_list[dataset.kind][self.period_name[jj] + '_off']]
                amp_in[jj, :dataset.n_amp[jj]] = \
                        theta[mc.variable_list[dataset.name_ref][self.period_name[jj] + '_amp']]

        return self.model_sinusoids(dataset, Prot, amp_in, pha_in, off_in)

    def print_vars(self, mc, theta):
        # Prot and pha_in could be brought out from the llop, but I would not work
        # for plaet-only analysis
        mc.pam_names[mc.variable_list['Prot']] = 'Prot'
        print 'Prot ', theta[mc.variable_list['Prot']]

        for jj in range(0, self.n_periods):
            id_var = mc.variable_list[self.period_name[jj] + '_pha']
            if np.size(id_var) == 0:
                continue
            if np.size(id_var) == 1:
                mc.pam_names[id_var] = self.period_name[jj] + '_pha'
            else:
                for ii in id_var:
                    mc.pam_names[ii] = self.period_name[jj] + '_' + repr(ii - id_var[0]) + '_pha'

            print self.period_name[jj], '_pha', theta[id_var]

        for dataset in mc.dataset_list:
            for jj in range(0, self.n_periods):
                if dataset.periods_flag[jj]:
                    if self.use_offset[dataset.kind]:
                        id_var = mc.variable_list[dataset.kind][self.period_name[jj] + '_off']
                        if np.size(id_var) == 0:
                            continue
                        if np.size(id_var) == 1:
                            mc.pam_names[id_var] = dataset.kind + '_' + self.period_name[jj] + '_off'
                        else:
                            for ii in id_var:
                                mc.pam_names[ii] = dataset.kind + '_' + self.period_name[jj] + \
                                                     '_' + repr(ii - id_var[0]) + '_off'

                        print dataset.name_ref, dataset.kind, self.period_name[jj], '_off', \
                            theta[mc.variable_list[dataset.kind][self.period_name[jj] + '_off']]

                    id_var = mc.variable_list[dataset.name_ref][self.period_name[jj] + '_amp']
                    if np.size(id_var) == 0:
                        continue
                    if np.size(id_var) == 1:
                        mc.pam_names[id_var] = dataset.name_ref[:-4] + '_' + self.period_name[jj] + '_amp'
                    else:
                        for ii in id_var:
                            mc.pam_names[ii] = dataset.name_ref[:-4] + '_' + self.period_name[jj] + \
                                                 '_' + repr(ii - id_var[0]) + '_amp'

                    print dataset.name_ref, self.period_name[jj], '_amp',\
                        theta[mc.variable_list[dataset.name_ref][self.period_name[jj] + '_amp']]

    def define_bounds(self, mc):
        mc.bounds_list.append(self.Prot_bounds[:])
        mc.variable_list['Prot'] = mc.ndim
        mc.ndim += 1

        for jj in range(0, self.n_periods):
            for kk in range(0, self.n_pha[jj]):
                mc.bounds_list.append(self.pha_bounds[:])
            mc.variable_list[self.period_name[jj] + '_pha'] = np.arange(mc.ndim,
                mc.ndim + self.n_pha[jj], 1)
            mc.ndim += self.n_pha[jj]

        for dataset in mc.dataset_list:

            # two nested case:
            # 1) dataset.kind has an offset or not (if it is the reference offset)
            # 2) the offset is the same for every season or not
            # WARNING: the offset is defined for each DATASET.KIND and not for each DATASET.NAME_REF
            # since the offset is a phyisical effect (not an instrumental one)

            for jj in range(0, self.n_periods):

                if dataset.periods_flag[jj]:

                    if self.use_offset[dataset.kind]:
                        if (not self.offset_coherence) or \
                                (self.offset_coherence and self.offset_common_id < 0):
                            mc.bounds_list.append(self.pof_bounds[:])
                            mc.variable_list[dataset.kind][self.period_name[jj] + '_off'] = mc.ndim
                            self.offset_common_id = mc.ndim
                            mc.ndim += 1
                        else:
                            mc.variable_list[dataset.kind][self.period_name[jj] + '_off'] = \
                                self.offset_common_id

                    for kk in xrange(0, dataset.n_amp[jj]):
                        mc.bounds_list.append(dataset.sinamp_bounds)

                    mc.variable_list[dataset.name_ref][self.period_name[jj] + '_amp'] = \
                        np.arange(mc.ndim, mc.ndim + dataset.n_amp[jj], 1)
                    mc.ndim += dataset.n_amp[jj]

    @staticmethod
    def model_sinusoids(dataset, p_rot, amp, pha, off):
        # np.size(SinAmp)==np.sum(n_amp) ???? What if an activity dataset is missing?
        # cycle for np.sum(rva_mask[:, jj]) <= 0 ??

        # Prot= Rotational periodo of the star
        # n_amp = number of sinusoids in the fit
        # amp = amplitudes
        # sin = phase of each sinusoid
        # off = overall offset of the sinusoids
        # sel = if the model has to be restricted to a specific temporal range
        model = np.zeros(dataset.n, dtype=np.double)
        xph = (dataset.x0 / p_rot) % 1
        har = np.arange(1, np.size(amp, axis=1) + 1, 1., dtype=np.double)

        for ii in xrange(0, dataset.n_periods):
            for jj in xrange(0, np.size(pha[ii, :])):
                model += dataset.p_mask[ii, :] * amp[ii, jj] * np.sin(
                        (har[jj] * xph + pha[ii, jj] + off[ii]) * 2. * np.pi)
        return model

class Dataset:
    def __init__(self, ind, kind, input_file, models):
        self.ind = ind
        self.kind = kind
        # 'RV', 'PHOT', 'ACT'...
        self.models = models
        self.name_ref = input_file
        print 'Opening: ', input_file
        self.data = np.loadtxt(input_file)
        self.x = np.asarray(self.data[:, 0], dtype=np.double)
        self.y = np.asarray(self.data[:, 1], dtype=np.double)
        self.e = np.asarray(self.data[:, 2], dtype=np.double)
        self.j = np.asarray(self.data[:, 3], dtype=np.double)
        # fit for different RV jitters
        self.o = np.asarray(self.data[:, 4], dtype=np.double)
        self.l = np.asarray(self.data[:, 5], dtype=np.double)
        # use different offsets for the data
        # off must start from zero
        # -1 values for self.j and self.l mean that these

        # self.a = np.asarray(self.data[:, 5], dtype=np.double)
        # Flag for activity fit: we can choose to fit the same RV amplitude
        # of the activity signal for all the datasets, use dfferent values
        # for each dataset or exclude some of the datasets

        # Model for RV systematics
        self.n_o = np.max(self.o.astype(np.int64)) + 1
        self.n_j = np.max(self.j.astype(np.int64)) + 1
        self.n_l = np.max(self.l.astype(np.int64)) + 1
        # self.n_a = np.max(self.a.astype(np.int64)) + 1

        self.n = np.size(self.x)
        print 'N = ', self.n
        print 'N jitter = ', self.n_j
        print 'N offset = ', self.n_o
        print 'N linear = ', self.n_l
        # print 'N activ. = ', self.n_a
        print

        self.Tref = np.mean(self.x, dtype=np.double)
        self.x0 = self.x - self.Tref

        self.o_mask = np.zeros([self.n, self.n_o], dtype=bool)
        for ii in xrange(0, self.n_o):
            self.o_mask[(abs(self.o - ii) < 0.1), ii] = True

        self.j_mask = np.zeros([self.n, self.n_j], dtype=bool)
        for ii in xrange(0, self.n_j):
            self.j_mask[(abs(self.j - ii) < 0.1), ii] = True

        self.l_mask = np.zeros([self.n, self.n_l], dtype=bool)
        for ii in xrange(0, self.n_l):
            self.l_mask[(abs(self.l - ii) < 0.1), ii] = True

        self.model = np.zeros(self.n, dtype=np.double)
        self.jitter = np.zeros(self.n, dtype=np.double)

        self.y_model = self.y * 0.00

    def common_Tref(self, Tref_in):
        self.Tref = Tref_in
        self.x0 = self.x - self.Tref
        return

    def model_reset(self):
        self.model[:] = 0.0
        self.jitter[:] = 0.0
        return

    def model_offset(self, off_in):
        off = np.atleast_1d(off_in)
        for ii in xrange(0, self.n_o):
            self.model[self.o_mask[:, ii]] += off[ii]

    def model_linear(self, m_in):
        m = np.atleast_1d(m_in)
        for ii in xrange(0, self.n_l):
            self.model[self.l_mask[:, ii]] += m[ii] * self.x0[self.l_mask[:, ii]]

    def model_jitter(self, jit_in):
        jit = np.atleast_1d(jit_in)
        for ii in xrange(0, self.n_j):
            self.jitter[self.j_mask[:, ii]] = jit[ii]

    def model_logchi2(self):
        env = 1.0 / (self.e ** 2.0 + self.jitter ** 2.0)
        return -0.5 * (np.sum((self.y - self.model) ** 2 * env - np.log(env)))

    def define_bounds(self, mc):
        for jj in xrange(0, self.n_j):
            mc.bounds_list.append([0., 50 * np.max(self.e)])  # bounds for jitter
        for jj in xrange(0, self.n_o):
            mc.bounds_list.append([np.min(self.y), np.max(self.y)])
        for jj in xrange(0, self.n_l):
            mc.bounds_list.append([-1., 1.])

        mc.variable_list[self.kind] = {}
        mc.variable_list[self.name_ref] = {}

        mc.variable_list[self.name_ref]['jitter'] = np.arange(mc.ndim, mc.ndim + self.n_j, 1)
        mc.ndim += self.n_j
        mc.variable_list[self.name_ref]['offset'] = np.arange(mc.ndim, mc.ndim + self.n_o, 1)
        mc.ndim += self.n_o
        mc.variable_list[self.name_ref]['linear'] = np.arange(mc.ndim, mc.ndim + self.n_l, 1)
        mc.ndim += self.n_l

class ModelContainer:
    def __init__(self):
        self.dataset_list = []
        self.n_datasets = 0
        self.scv = SinusoidsCommonVariables()
        self.pcv = PlanetsCommonVariables()

        # pyde/emcee variables
        self.ngen = 0
        self.nsteps = 0
        self.nburn = 0
        self.npop_multi = 0
        self.nwalkers = 0
        self.thin = 0

        self.model_list = []

        self.recenter_bounds_flag = True

        self.planet_name = ''

        self.variable_list = {}
        self.bound_list = []
        self.bounds = 0
        self.ndim = 0
        self.pam_names = ''
        self.star_mass_val = 1.0000
        self.star_mass_err = 0.1000

    def model_setup(self):
        self.n_datasets = np.size(self.dataset_list)
        for dataset in self.dataset_list:
            # if 'kepler' in dataset.models:

            if 'sinusoids' in dataset.models:
                self.scv.setup_model_sinusoids(dataset)

            dataset.model_reset()
            for data_model in dataset.models:
                if not (data_model in self.model_list):
                    self.model_list.append(data_model)

    def create_bounds(self):
        # This routine creates the boundary array and at the same time
        # creates a dictionary with the name of the arrays and their
        # positions in bounds/theta array so that they can be accessed
        # without using nested counters

        self.bounds_list = []
        self.ndim = 0

        for dataset in self.dataset_list:
            dataset.define_bounds(self)

        if 'kepler' in self.model_list:
            self.pcv.define_bounds(self)

        if 'sinusoids' in self.model_list:
            self.scv.define_bounds(self)

        self.bounds = np.asarray(self.bounds_list)

    def check_bounds(self, theta):
        # former lnprior(theta)
        for ii in xrange(0, self.ndim):
            if not (self.bounds[ii, 0] < theta[ii] < self.bounds[ii, 1]):
                return False
        for planet_name in self.pcv.name_ref:
            if self.pcv.n_orbpams[planet_name] == 5:
                e = theta[self.variable_list[planet_name]['esino']] ** 2 + \
                    theta[self.variable_list[planet_name]['ecoso']] ** 2
                if not self.pcv.e_bounds[planet_name][0] <= e < self.pcv.e_bounds[planet_name][1]:
                    return False

        return True

    def __call__(self, theta):
        if not self.check_bounds(theta):
            return -np.inf
        chi2_out = 0.0

        for dataset in self.dataset_list:
            dataset.model_reset()
            dataset.model_offset(theta[self.variable_list[dataset.name_ref]['offset']])
            dataset.model_jitter(theta[self.variable_list[dataset.name_ref]['jitter']])
            dataset.model_linear(theta[self.variable_list[dataset.name_ref]['linear']])

            if 'kepler' in dataset.models:
                for planet_name in self.pcv.name_ref:
                    dataset.model += self.pcv.compute(self, theta, dataset, planet_name)
                    #in_pams = theta[self.variable_list[planet_name]['kepler_pams']]
                    #if self.pcv.n_orbpams[planet_name] == 5:
                    #    conv_pams = self.pcv.convert_params(in_pams[0], in_pams[1], in_pams[2], in_pams[3], in_pams[4])
                    #else:
                    #    conv_pams = self.pcv.convert_params(in_pams[0], in_pams[1], in_pams[2], 0.00, 0.00)
                    #dataset.model += self.pcv.model_kepler(conv_pams, dataset.x0)

            if 'sinusoids' in dataset.models:
                #dataset.model += self.sinusoids_model(theta,dataset)
                dataset.model += self.scv.compute(self, theta, dataset)

            chi2_out += dataset.model_logchi2()
        return chi2_out

    # def sinusoids_model(self,theta,dataset):
    #     # Prot and pha_in could be brought out from the llop, but I would not work
    #     # for plaet-only analysis
    #     Prot = theta[self.variable_list['Prot']]
    #     pha_in = np.zeros([self.scv.n_periods, self.scv.n_pha_max], dtype=np.double)
    #     off_in = np.zeros([self.scv.n_periods], dtype=np.double)
    #     amp_in = np.zeros([self.scv.n_periods, self.scv.n_pha_max], dtype=np.double)
    #     for jj in range(0, self.scv.n_periods):
    #         pha_in[jj, :self.scv.n_pha[jj]] = theta[self.variable_list[self.scv.period_name[jj] + '_pha']]
    #
    #         if dataset.periods_flag[jj]:
    #             if self.scv.use_offset[dataset.kind]:
    #                 off_in[:] = theta[self.variable_list[dataset.kind][self.scv.period_name[jj] + '_off']]
    #             amp_in[jj, :dataset.n_amp[jj]] = \
    #                     theta[self.variable_list[dataset.name_ref][self.scv.period_name[jj] + '_amp']]
    #
    #     return self.scv.model_sinusoids(dataset, Prot, amp_in, pha_in, off_in)
    #
    # def kepler_model(self, theta, dataset, planet_name):
    #     in_pams = theta[self.variable_list[planet_name]['kepler_pams']]
    #     if self.pcv.n_orbpams[planet_name] == 5:
    #         conv_pams = self.pcv.convert_params(in_pams[0], in_pams[1], in_pams[2], in_pams[3], in_pams[4])
    #     else:
    #         conv_pams = self.pcv.convert_params(in_pams[0], in_pams[1], in_pams[2], 0.00, 0.00)
    #     return self.pcv.model_kepler(conv_pams, dataset.x0)

    def results_resumen(self, theta):
        # Function with two goals:
        # * Unfold eand print out the output from theta
        # * give back a parameter name associated to each value in the result array

        self.pam_names= self.ndim * ['']

        for dataset in self.dataset_list:
            for param in ['offset', 'jitter', 'linear']:
                id_var = self.variable_list[dataset.name_ref][param]
                if np.size(id_var) == 0:
                    continue
                if np.size(id_var) == 1:
                    self.pam_names[id_var] = dataset.name_ref[:-4] + '_' + param
                else:
                    for ii in id_var:
                        self.pam_names[ii] = dataset.name_ref[:-4] + '_' + param + '_' + repr(ii - id_var[0])

                print dataset.name_ref, param, ' : ', theta[self.variable_list[dataset.name_ref][param]]

        if 'kepler' in self.model_list:
            self.pcv.print_vars(self, theta)

        if 'sinusoids' in self.model_list:
            self.scv.print_vars(self, theta)

    # it return the RV model for a single planet, after removing the activity from the RV curve and removing the offsets
    # between the datasets
    def rv_1planet_model(self, theta, planet_ref, dataset):

        dataset.model_reset()
        dataset.model_offset(theta[self.variable_list[dataset.name_ref]['offset']])
        dataset.model_linear(theta[self.variable_list[dataset.name_ref]['linear']])

        if 'kepler' in dataset.models:
            for planet_name in self.pcv.name_ref:
                dataset.model += self.pcv.model(self, theta, dataset, planet_name)

        if 'sinusoids' in dataset.models:
            dataset.model += self.scv.model(self, theta, dataset)

        return dataset.y - dataset.model, self.pvc.model(self, theta, dataset, planet_ref)

    # This function recenters the bounds limits for circular variables
    # Also, it extends the range of a variable if the output of PyDE is a fixed number
    def recenter_bounds(self, pop_mean, population):
        ind_list = []

        if 'kepler' in self.model_list:
            for planet_name in self.pcv.name_ref:
                ind_list.append(self.variable_list[planet_name]['phase'])

                if self.pcv.n_orbpams[planet_name] == 5:
                    esino_list =  self.variable_list[planet_name]['esino']
                    ecoso_list =  self.variable_list[planet_name]['ecoso']

                    e_pops = population[:, esino_list]**2 + population[:, ecoso_list]**2
                    e_mean = (self.pcv.e_bounds[planet_name][0]+self.pcv.e_bounds[planet_name][1])/2.
                    for ii in xrange(0, self.nwalkers):
                        if not self.pcv.e_bounds[planet_name][0] <= e_pops[ii] < self.pcv.e_bounds[planet_name][1]:
                            population[ii, esino_list] = np.sqrt(e_mean)
                            population[ii, ecoso_list] = np.sqrt(e_mean)

        if 'sinusoids' in self.model_list:
            for jj in range(0, self.scv.n_periods):
                ind_list.extend(self.variable_list[self.scv.period_name[jj] + '_pha'])
            for dataset in self.dataset_list:
                for jj in range(0, self.scv.n_periods):
                    if dataset.periods_flag[jj]:
                        # ind_list.extend(self.variable_list[dataset.name_ref][self.scv.period_name[jj] + '_amp'])
                        if self.scv.use_offset[dataset.kind]:
                            ind_list.append(self.variable_list[dataset.kind][self.scv.period_name[jj] + '_off'])

        if np.size(ind_list) > 0:
            tmp_range = (self.bounds[:, 1]-self.bounds[:, 0])/2
            for var_ind in ind_list:
                self.bounds[var_ind, :] = pop_mean[var_ind] + [-tmp_range[var_ind], tmp_range[var_ind]]
                fix_sel = (population[:, var_ind] <= self.bounds[var_ind, 0]) | (
                    population[:, var_ind] >= self.bounds[var_ind, 1])
                population[fix_sel, var_ind] = pop_mean[var_ind]

        for ii in xrange(0, self.ndim):
            if np.amax(population[:, ii])-np.amin(population[:, ii]) < 10e-7 :
                range_restricted = (self.bounds[ii, 1]-self.bounds[ii, 0])/100.
                min_bound = np.maximum((pop_mean[ii]-range_restricted/2.0), self.bounds[ii,0])
                max_bound = np.minimum((pop_mean[ii]+range_restricted/2.0), self.bounds[ii,1])
                population[:, ii] = np.random.uniform(min_bound, max_bound, self.nwalkers)


def get_pyorbit_input(file_conf, mc):
    file_in = open(file_conf, 'r')
    row = file_in.readlines()

    pl_row = 0
    pr_row = 0
    in_row = 0

    for line in row:
        info = line.split()
        if line.find("Input_" + `in_row` + " ") > -1:
            input_file = info[1]
            input_kind = info[2]
            input_model = info[3:]
            mc.dataset_list.append(Dataset(in_row, input_kind, input_file, input_model))
            if in_row > 0:
                mc.dataset_list[in_row].common_Tref(mc.dataset_list[0].Tref)
            mc.Tref = mc.dataset_list[0].Tref
            in_row += 1

        if line.find("Output") > -1:
            mc.planet_name = info[1]
        # if line.find("Nplanets") > -1:
        #    n_planets = np.asarray(info[1], dtype=np.int64)
        if line.find("Planet_" + `pl_row` + " ") > -1:
            mc.pcv.add_planet(info[0], np.asarray(info[1], dtype=np.int64))

            p_min = np.asarray(info[2], dtype=np.double)
            p_max = np.asarray(info[3], dtype=np.double)
            mc.pcv.logP_bounds[info[0]] = np.asarray([np.log2(p_min), np.log2(p_max)])

            k_min = np.asarray(info[4], dtype=np.double)
            k_max = np.asarray(info[5], dtype=np.double)
            mc.pcv.logK_bounds[info[0]] = np.asarray([np.log2(k_min), np.log2(k_max)])

            if len(info) > 6:
                mc.pcv.e_bounds[info[0]] = np.asarray([info[6], info[7]], dtype=np.double)
            else:
                mc.pcv.e_bounds[info[0]] = np.asarray([0., 1.], dtype=np.double)
            pl_row += 1
            # if line.find("LongTerm") > -1:
            # Two possible drift are possible: one caused by the instrument
            # (and this one is coded in the data through the rv_l flag)
            # and one due to the star, and common to all the datasets
            # Use this option in the second case
            # info= line.split()
            # n_trends = np.asarray(info[2], dtype=np.int64)
        if line.find("Prot") > -1:
            # info= line.split()
            p_min = np.asarray(info[2], dtype=np.double)
            p_max = np.asarray(info[3], dtype=np.double)
            # kind_sin = np.asarray(info[1], dtype=np.int64)

            mc.scv.Prot_bounds = [p_min, p_max]

        # The number of chuncks for the activity fit
        # if line.find("Nperiods") > -1:
        #    scv.n_periods = np.asarray(info[1], dtype=np.int64)

        # Each chunk has a different rotational period, still within the boundaries
        # defined in  Prot, and it is defined as an interval in BJD
        if line.find("Period_" + `pr_row` + " ") > -1:
            P0_tmp = np.asarray(info[1], dtype=np.double)
            P1_tmp = np.asarray(info[2], dtype=np.double)
            mc.scv.add_period_range([P0_tmp, P1_tmp], info[3:])
            # print period_range[pr_row,1]
            # print period_range[pr_row ][0]
            # print period_range[pr_row ][1]

            pr_row += 1

        if line.find("Ngen") > -1:
            mc.ngen = np.asarray(info[1], dtype=np.int64)
        if line.find("Nsteps") > -1:
            mc.nsteps = np.asarray(info[1], dtype=np.int64)
        if line.find("Nburn") > -1:
            mc.nburn = np.asarray(info[1], dtype=np.int64)
        if line.find("Npop_mult") > -1:
            mc.npop_mult = np.asarray(info[1], dtype=np.int64)
        if line.find("Thin") > -1:
            mc.thin = np.asarray(info[1], dtype=np.int64)

        if line.find("Phase_dataset") > -1:
            # Number of sinusoids that must be used for each dataset for the activity fit
            # If Period_n are not specified in the input file
            # So, if we are using two dataset where the first one has greater precision
            # we can specifiy to use 3 sinusoids for the first dataset and only two
            # for the second one
            # Additional activity indicators associated to RVs must the same order and number of sinusoids,
            #
            # scv.n_spec = np.asarray(info[1], dtype=np.int64)
            mc.scv.phase = np.asarray(info[1:], dtype=np.int64)

        if line.find("Recenter_Bounds") > -1:
            # required to avoid a small bug in the code
            # if the dispersion of PyDE walkers around the median value is too broad,
            # then emcee walkers will start outside the bounds, causing an error
            if np.asarray(info[1], dtype=np.int64) < 1:
                mc.recenter_bounds_flag = False
        if line.find("Star_Mass") > -1:
            mc.star_mass_val = np.asarray(info[1], dtype=np.double)
            mc.star_mass_err = np.asarray(info[2], dtype=np.double)
    #if mc.scv.n_periods != pr_row:
    #    print 'PROBLEM'

    mc.model_setup()
