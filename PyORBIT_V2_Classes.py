import numpy as np
import kepler_exo as kp
import yaml
import george

def get_var_log(var, fix, i):
    return np.log2(var[i], dtype=np.double)

def get_var_exp(var, fix, i):
    return np.exp2(var[i], dtype=np.double)


def get_var_val(var, fix, i):
    return var[i]


def get_var_sqr(var, fix, i):
    return var[i]**2


def get_var_inv(var, fix, i):
    return 1./var[i]


def get_var_2inv(var, fix, i):
    return 1./(2*var[i])


def get_fix_log(var, fix, i):
    return np.log2(fix[i], dtype=np.double)


def get_fix_exp(var, fix, i):
    return np.exp2(fix[i], dtype=np.double)


def get_fix_val(var, fix, i):
    return fix[i]


def get_fix_sqr(var, fix, i):
    return fix[i]**2


def get_fix_inv(var, fix, i):
    return 1./fix[i]


def get_fix_2inv(var, fix, i):
    return 1./(2*fix[i])


def get_2var_e(var, fix, i):
    ecoso = var[i[0]]
    esino = var[i[1]]
    return np.square(ecoso, dtype=np.double) + np.square(esino, dtype=np.double)


def get_2var_o(var, fix, i):
    ecoso = var[i[0]]
    esino = var[i[1]]
    return np.arctan2(esino, ecoso, dtype=np.double)


# A class for common variables (number of planets, number of sinusoids...)
# Info to unpack the variables inside emcee must be included here
# Physical effects must be included here
class PlanetsCommonVariables:
    def __init__(self):
        self.n_planets = 0
        self.name_ref = []

        self.bounds = {}

        self.variables = {}
        self.var_list = {}
        self.fix_list = {}

        self.circular = {}

        self.fixed = []
        self.nfix = 0

        self.list_pams = ['P', 'K', 'f', 'e', 'o']

    def add_planet(self, name_ref):
        self.n_planets += 1
        self.name_ref.append(name_ref)

        self.fix_list[name_ref] = {}
        self.var_list[name_ref] = {}
        self.variables[name_ref] = {}

        self.bounds[name_ref] = {}
        self.bounds[name_ref]['logP'] = [0.0, 7.0]
        self.bounds[name_ref]['logK'] = [0.0, 7.0]
        self.bounds[name_ref]['phase'] = [0.0, 2 * np.pi]
        self.bounds[name_ref]['ecoso'] = [0.0, 1.0]
        self.bounds[name_ref]['esino'] = [0.0, 1.0]
        self.bounds[name_ref]['e'] = [0.0, 1.0]
        self.bounds[name_ref]['o'] = [0.0, 2 * np.pi]

    # def convert_params(self,in_pams):
    #    out_pams = np.zeros(5,dtype=np.double)
    #    out_pams[0] = 10.0 ** (in_pams[0]) #P
    #    out_pams[1] = 10.0 ** (in_pams[1]) #K
    #    out_pams[2] = in_pams[2]
    #    if np.size(in_pams)==5:
    #        out_pams[3] = in_pams[3] ** 2 + in_pams[4] ** 2
    #        out_pams[4] = np.arctan2(in_pams[3], in_pams[4])
    #    return out_pams

    def define_bounds(self, mc):

        for planet_name in self.name_ref:
            mc.variable_list[planet_name] = {}

            ndim_buffer = mc.ndim
            if 'P' in self.fix_list[planet_name]:
                self.variables[planet_name]['P'] = get_fix_val
                self.var_list[planet_name]['P'] = self.nfix
                self.fixed.append(self.fix_list[planet_name]['P'])
                self.nfix += 1
            else:
                self.variables[planet_name]['P'] = get_var_exp
                self.var_list[planet_name]['P'] = mc.ndim
                mc.variable_list[planet_name]['logP'] = mc.ndim
                mc.bounds_list.append(self.bounds[planet_name]['logP'])
                mc.ndim += 1

            if 'K' in self.fix_list[planet_name]:
                self.variables[planet_name]['K'] = get_fix_val
                self.var_list[planet_name]['K'] = self.nfix
                self.fixed.append(self.fix_list[planet_name]['K'])
                self.nfix += 1
            else:
                self.variables[planet_name]['K'] = get_var_exp
                self.var_list[planet_name]['K'] = mc.ndim
                mc.variable_list[planet_name]['logK'] = mc.ndim
                mc.bounds_list.append(self.bounds[planet_name]['logK'])
                mc.ndim += 1

            if 'f' in self.fix_list[planet_name]:
                self.variables[planet_name]['f'] = get_fix_val
                self.var_list[planet_name]['f'] = self.nfix
                self.fixed.append(self.fix_list[planet_name]['f'])
                self.nfix += 1
            else:
                self.variables[planet_name]['f'] = get_var_val
                self.var_list[planet_name]['f'] = mc.ndim
                mc.variable_list[planet_name]['phase'] = mc.ndim
                mc.bounds_list.append(self.bounds[planet_name]['phase'])
                mc.ndim += 1

            if 'e' in self.fix_list[planet_name] and 'o' in self.fix_list[planet_name]:
                self.variables[planet_name]['e'] = get_fix_val
                self.var_list[planet_name]['e'] = self.nfix
                self.fixed.append(self.fix_list[planet_name]['e'])
                self.nfix += 1
                self.variables[planet_name]['o'] = get_fix_val
                self.var_list[planet_name]['o'] = self.nfix
                self.fixed.append(self.fix_list[planet_name]['o'])
                self.nfix += 1
            elif 'e' in self.fix_list[planet_name]:
                self.variables[planet_name]['e'] = get_fix_val
                self.var_list[planet_name]['e'] = self.nfix
                self.fixed.append(self.fix_list[planet_name]['e'])
                self.nfix += 1
                self.variables[planet_name]['o'] = get_var_val
                self.var_list[planet_name]['o'] = mc.ndim
                mc.variable_list[planet_name]['o'] = mc.ndim
                mc.bounds_list.append(self.bounds[planet_name]['o'])
                mc.ndim += 1
            elif 'o' in self.fix_list[planet_name]:
                self.variables[planet_name]['o'] = get_fix_val
                self.var_list[planet_name]['o'] = self.nfix
                self.fixed.append(self.fix_list[planet_name]['o'])
                self.nfix += 1
                self.variables[planet_name]['e'] = get_var_val
                self.var_list[planet_name]['e'] = mc.ndim
                mc.variable_list[planet_name]['e'] = mc.ndim
                mc.bounds_list.append(self.bounds[planet_name]['e'])
                mc.ndim += 1
            else:
                self.variables[planet_name]['e'] = get_2var_e
                self.var_list[planet_name]['e'] = [mc.ndim, mc.ndim + 1]
                self.variables[planet_name]['o'] = get_2var_o
                self.var_list[planet_name]['o'] = [mc.ndim, mc.ndim + 1]
                mc.variable_list[planet_name]['ecoso'] = mc.ndim
                mc.variable_list[planet_name]['esino'] = mc.ndim + 1
                mc.bounds_list.append(self.bounds[planet_name]['ecoso'])
                mc.bounds_list.append(self.bounds[planet_name]['esino'])
                mc.ndim += 2
            mc.variable_list[planet_name]['kepler_pams'] = np.arange(ndim_buffer, mc.ndim, 1)

    def convert(self, planet_name, theta):
        list_out = []
        for name in self.list_pams:
            list_out.append(self.variables[planet_name][name](theta, self.fixed, self.var_list[planet_name][name]))
        return list_out

        # return self.variables[planet_name]['P'](theta, fixed, ), self.variables[planet_name]['K'](theta, fixed, 2), \
        #       self.variables[planet_name]['f'](theta, fixed, 2), self.variables[planet_name]['e'](theta, fixed, 2), \
        #       self.variables[planet_name]['o'](theta, fixed, 2)

    def compute(self, theta, dataset, planet_name):
        pams = self.convert(planet_name, theta)
        return self.model_kepler(pams, dataset.x0)

    def print_vars(self, mc, theta):
        for planet_name in self.name_ref:
            out_list = self.convert(planet_name, theta)

            for var_name in mc.variable_list[planet_name]:
                if var_name != 'kepler_pams':
                    mc.pam_names[mc.variable_list[planet_name][var_name]] = var_name
            print planet_name, ' vars: ', np.asarray(theta[mc.variable_list[planet_name]['kepler_pams']])
            print planet_name, ' pams: ', out_list[:]

    @staticmethod
    def model_kepler(orbit_pams, x0):
        P, K, phase, e, omega = orbit_pams
        rv_out = kp.kepler_RV_T0P(x0, phase, P, K, e, omega)
        return rv_out

    #@staticmethod
    #def convert_params(logP, logK, phase, esino, ecoso):
    #    P = np.exp2(logP)
    #    K = np.exp2(logK)
    #    e = np.square(ecoso) + np.square(esino)
    #    o = np.arctan2(esino, ecoso)
    #    return P, K, phase, e, o


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

        self.season_sel = False
        self.n_seasons = 1
        self.season_name = ['Season_0']
        self.season_list = [0., 5000000.0]
        self.season_range = np.asarray(self.season_list, dtype=np.double)

        # self.Prot_bounded = False
        # self.pha_bounded = False

        self.Prot_bounds = [0., 30.0]
        self.pof_bounds = [0., 1.0]
        self.pha_bounds = [0., 1.0]

        self.phase_list = []

    def add_season_range(self, range_input, phase_input):
        # Reset default values at first call
        if not self.season_sel:
            self.n_seasons = 0
            self.season_list = []
            self.season_name = []
            self.season_sel = True

        self.season_name.append('Season_' + repr(self.n_seasons))
        self.season_list.append(range_input)
        self.season_range = np.asarray(self.season_list, dtype=np.double)

        self.phase_list.append(phase_input)
        self.phase = np.asarray(self.phase_list, dtype=np.int64)

        self.n_pha = np.amax(self.phase, axis=1)  # maximum value for each period
        self.n_pha_max = np.amax(self.n_pha)  # maximum value for each period

        self.n_seasons += 1
        return

    # def add_phase_offset(self, dataset, season_name):
    #    # Check if the dataset needs an offset with respect to the
    #    # reference dataset (simply the first one)
    #    add_bound = False
    #    try:
    #        value = self.use_offset[dataset.kind][season_name + '_off']
    #    except KeyError:
    #        # Key is not present
    #        if self.offset_skip_first:
    #            self.use_offset[dataset.kind][season_name] = False
    #            self.offset_skip_first = False
    #        else:
    #            self.use_offset[dataset.kind][season_name] = True
    #            add_bound = True
    #    return add_bound

    def setup_model_sinusoids(self, dataset):
        dataset.n_amp = self.phase[:, dataset.ind]
        dataset.p_mask = np.zeros([self.n_seasons, dataset.n], dtype=bool)
        dataset.n_seasons = 0
        dataset.season_flag = np.zeros(self.n_seasons, dtype=bool)

        # If this is the first dataset, it is assumed as the reference one for the offset
        # otherwise, the offset is applied
        try:
            value = self.use_offset[dataset.kind]
        except KeyError:
            if self.offset_reference_name == '':
                self.offset_reference_name = dataset.kind
                self.use_offset[dataset.kind] = False
            else:
                self.use_offset[dataset.kind] = not self.phase_sincro  # True

        for ii in xrange(0, self.n_seasons):
            p_sel = (self.season_range[ii, 0] < dataset.x) & (dataset.x < self.season_range[ii, 1])
            if np.sum(p_sel) == 0:
                'No data points within the specified range'
            else:
                dataset.season_flag[ii] = True
                dataset.p_mask[ii, :] = p_sel[:]
                dataset.n_seasons += 1

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
        pha_in = np.zeros([self.n_seasons, self.n_pha_max], dtype=np.double)
        off_in = np.zeros([self.n_seasons], dtype=np.double)
        amp_in = np.zeros([self.n_seasons, self.n_pha_max], dtype=np.double)
        for jj in range(0, self.n_seasons):
            pha_in[jj, :self.n_pha[jj]] = theta[mc.variable_list[self.season_name[jj] + '_pha']]

            if dataset.season_flag[jj]:
                if self.use_offset[dataset.kind]:
                    off_in[:] = theta[mc.variable_list[dataset.kind][self.season_name[jj] + '_off']]
                amp_in[jj, :dataset.n_amp[jj]] = \
                    theta[mc.variable_list[dataset.name_ref][self.season_name[jj] + '_amp']]

        return self.model_sinusoids(dataset, Prot, amp_in, pha_in, off_in)

    def print_vars(self, mc, theta):
        # Prot and pha_in could be brought out from the llop, but I would not work
        # for plaet-only analysis
        mc.pam_names[mc.variable_list['Prot']] = 'Prot'
        print 'Prot ', theta[mc.variable_list['Prot']]

        for jj in range(0, self.n_seasons):
            id_var = mc.variable_list[self.season_name[jj] + '_pha']
            if np.size(id_var) == 0:
                continue
            if np.size(id_var) == 1:
                mc.pam_names[id_var] = self.season_name[jj] + '_pha'
            else:
                for ii in id_var:
                    mc.pam_names[ii] = self.season_name[jj] + '_' + repr(ii - id_var[0]) + '_pha'

            print self.season_name[jj], '_pha', theta[id_var]

        for dataset in mc.dataset_list:
            for jj in range(0, self.n_seasons):
                if dataset.season_flag[jj]:
                    if self.use_offset[dataset.kind]:
                        id_var = mc.variable_list[dataset.kind][self.season_name[jj] + '_off']
                        if np.size(id_var) == 0:
                            continue
                        if np.size(id_var) == 1:
                            mc.pam_names[id_var] = dataset.kind + '_' + self.season_name[jj] + '_off'
                        else:
                            for ii in id_var:
                                mc.pam_names[ii] = dataset.kind + '_' + self.season_name[jj] + \
                                                   '_' + repr(ii - id_var[0]) + '_off'

                        print dataset.name_ref, dataset.kind, self.season_name[jj], '_off', \
                            theta[mc.variable_list[dataset.kind][self.season_name[jj] + '_off']]

                    id_var = mc.variable_list[dataset.name_ref][self.season_name[jj] + '_amp']
                    if np.size(id_var) == 0:
                        continue
                    if np.size(id_var) == 1:
                        mc.pam_names[id_var] = dataset.name_ref[:-4] + '_' + self.season_name[jj] + '_amp'
                    else:
                        for ii in id_var:
                            mc.pam_names[ii] = dataset.name_ref[:-4] + '_' + self.season_name[jj] + \
                                               '_' + repr(ii - id_var[0]) + '_amp'

                    print dataset.name_ref, self.season_name[jj], '_amp', \
                        theta[mc.variable_list[dataset.name_ref][self.season_name[jj] + '_amp']]

    def define_bounds(self, mc):
        mc.bounds_list.append(self.Prot_bounds[:])
        mc.variable_list['Prot'] = mc.ndim
        mc.ndim += 1

        for jj in range(0, self.n_seasons):
            for kk in range(0, self.n_pha[jj]):
                mc.bounds_list.append(self.pha_bounds[:])
            mc.variable_list[self.season_name[jj] + '_pha'] = np.arange(mc.ndim,
                                                                        mc.ndim + self.n_pha[jj], 1)
            mc.ndim += self.n_pha[jj]

        for dataset in mc.dataset_list:

            # two nested case:
            # 1) dataset.kind has an offset or not (if it is the reference offset)
            # 2) the offset is the same for every season or not
            # WARNING: the offset is defined for each DATASET.KIND and not for each DATASET.NAME_REF
            # since the offset is a physical effect (not an instrumental one)

            for jj in range(0, self.n_seasons):

                if dataset.season_flag[jj]:

                    if self.use_offset[dataset.kind]:
                        if (not self.offset_coherence) or \
                                (self.offset_coherence and self.offset_common_id < 0):
                            mc.bounds_list.append(self.pof_bounds[:])
                            mc.variable_list[dataset.kind][self.season_name[jj] + '_off'] = mc.ndim
                            self.offset_common_id = mc.ndim
                            mc.ndim += 1
                        else:
                            mc.variable_list[dataset.kind][self.season_name[jj] + '_off'] = \
                                self.offset_common_id

                    for kk in xrange(0, dataset.n_amp[jj]):
                        mc.bounds_list.append(dataset.sinamp_bounds)

                    mc.variable_list[dataset.name_ref][self.season_name[jj] + '_amp'] = \
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

        for ii in xrange(0, dataset.n_seasons):
            for jj in xrange(0, np.size(pha[ii, :])):
                model += dataset.p_mask[ii, :] * amp[ii, jj] * np.sin(
                    (har[jj] * xph + pha[ii, jj] + off[ii]) * 2. * np.pi)
        return model


class GaussianProcessCommonVariables:
    def __init__(self):
        self.fix_list = {}
        self.bounds = {}
        self.var_list = {}
        self.fix_list = {}
        self.variables = {}

        self.fixed = []
        self.nfix = 0

        # Three parameters out of four are the same for all the datasets, since they are
        # related
        # From Affer et al. 2016
        # - theta: is usually related to the rotation period of the star( or one of its harmonics);
        # - lambda: is the correlation decay timescale, and it can be related to the lifetime of the active regions.
        # - omega: is the length scale of the periodic component, and can be linked to the size evolution of the active regions;
        # - h: represents the amplitude of the correlations;

        self.list_pams_human = ['Prot', 'Pdec', 'Oamp', 'Hamp']
        self.list_pams_george = ['lnProt', 'Pds2', 'inv2O', 'H2']

        # There si a 1-1 correspondence ("biunivoca") between variables, so we use this trick
        self.list_pams_corr = {'Prot': 'lnProt', 'Pdec': 'Pds2', 'Oamp': 'inv2O', 'Hamp': 'H2'}

        self.list_pams_common = ['lnProt', 'Pds2', 'inv2O']
        self.list_pams_dataset = ['H2']

        return

    def convert_val(self, name, input_pam):
        # conversion of the GP parameters into physically meaningful ones
        if name == self.list_pams_george[0]:
            return np.exp(input_pam)
        if name == self.list_pams_human[0]:
            return np.log(input_pam)
        if name == self.list_pams_george[1]:
            return np.sqrt(input_pam)
        if name == self.list_pams_human[1]:
            return input_pam**2
        if name == self.list_pams_george[2]:
            return np.sqrt(1./(2.0*input_pam))
        if name == self.list_pams_human[2]:
            return 1./(2*input_pam**2)
        if name == self.list_pams_george[3]:
            return np.sqrt(input_pam)
        if name == self.list_pams_human[3]:
            return input_pam*input_pam

    def add_dataset(self, name_ref='Common'):
        # name_ref : name of the dataset, otherwise 'common'

        if name_ref == 'Common':
            for name in self.list_pams_common:
                self.bounds[name] = [0.1, 20]
        else:
            self.bounds[name_ref] = {}
            self.variables[name_ref] = {}
            self.fix_list[name_ref] = {}
            self.var_list[name_ref] = {}
            for name in self.list_pams_common:
                self.bounds[name_ref][name]= [0.00001, 20]

    def define_bounds(self, mc):

        for var in self.list_pams_common:
            if var in self.fix_list:
                self.variables[var] = get_fix_val
                self.var_list[var] = self.nfix
                self.fixed.append(self.fix_list[var])
                self.nfix += 1
            else:
                self.variables[var] = get_var_val
                self.var_list[var] = mc.ndim
                mc.variable_list[var] = mc.ndim
                mc.bounds_list.append(self.bounds[var])
                mc.ndim += 1

        for dataset in mc.dataset_list:
            if 'gaussian' in dataset.models:
                for var in self.list_pams_dataset:
                    if var in self.fix_list[dataset.name_ref]:
                        self.variables[dataset.name_ref][var] = get_fix_val
                        self.var_list[dataset.name_ref][var] = self.nfix
                        self.fixed.append(self.fix_list[dataset.name_ref][var])
                        self.nfix += 1
                    else:
                        self.variables[dataset.name_ref][var] = get_var_val
                        self.var_list[dataset.name_ref][var] = mc.ndim
                        mc.variable_list[dataset.name_ref][var] = mc.ndim
                        mc.bounds_list.append(self.bounds[dataset.name_ref][var])
                        mc.ndim += 1

    def obtain(self, dataset_name, theta):
        list_out = []
        for name in self.list_pams_common:
            list_out.append(self.variables[name](theta, self.fixed, self.var_list[name]))
        for name in self.list_pams_dataset:
            list_out.append(self.variables[dataset_name][name](theta, self.fixed, self.var_list[dataset_name][name]))
        return list_out

    def lnlk_compute(self, mc, theta, dataset):

        gp_pams = self.obtain(dataset.name_ref, theta)
        # gp_pams[0] = ln_theta = ln_Period -> ExpSine2Kernel(gamma, ln_period)
        # gp_pams[1] = metric = r^2 = lambda**2  -> ExpSquaredKernel(metric=r^2)
        # gp_pams[2] = Gamma =  1/ (2 omega**2) -> ExpSine2Kernel(gamma, ln_period)
        # gp_pams[3] = h^2 -> h^2 * ExpSquaredKernel * ExpSine2Kernel
        kernel = gp_pams[3] * george.kernels.ExpSquaredKernel(metric=gp_pams[1]) * \
            george.kernels.ExpSine2Kernel(gamma=gp_pams[2], ln_period=gp_pams[0])

        gp = george.GP(kernel)
        gp.compute(dataset.x0, dataset.e)

        return gp.lnlikelihood(dataset.y - dataset.model, quiet=True)

    def print_vars(self, mc, theta):

        for ii in xrange(0,np.size(self.list_pams_common)):
            name = self.list_pams_common[ii]
            var  = self.variables[name](theta, self.fixed, self.var_list[name])
            print 'GaussianProcess ', self.list_pams_common[ii], var, '     ',\
                self.list_pams_human[ii], '  ', self.convert_val(name, var)

        for dataset in mc.dataset_list:
            for ii in xrange(0, np.size(self.list_pams_dataset)):
                name = self.list_pams_dataset[ii]
                ii_add = ii + np.size(self.list_pams_common)
                var = self.variables[dataset.name_ref][name](theta, self.fixed, self.var_list[dataset.name_ref][name])
                print 'GaussianProcess ', dataset.name_ref, '   ',\
                    self.list_pams_dataset[ii], var, '     ', \
                    self.list_pams_human[ii+ii_add], '  ', self.convert_val(name, var)


class Dataset:
    def __init__(self, ind, kind, input_file, models):
        self.ind = ind
        self.kind = kind
        # 'RV', 'PHOT', 'ACT'...
        self.models = models
        self.name_ref = input_file

        print 'Opening: ', input_file
        self.data = np.loadtxt(input_file)

        n_cols = np.size(self.data, axis=1)

        self.x = np.asarray(self.data[:, 0], dtype=np.double)
        self.y = np.asarray(self.data[:, 1], dtype=np.double)
        self.e = np.asarray(self.data[:, 2], dtype=np.double)

        self.n = np.size(self.x)

        if n_cols > 3:
            self.j = np.asarray(self.data[:, 3], dtype=np.double)
        else:
            self.j = np.zeros(self.n, dtype=np.double)
        # fit for different dataset jitters

        if n_cols > 4:
            self.o = np.asarray(self.data[:, 4], dtype=np.double)
        else:
            self.o = np.zeros(self.n, dtype=np.double)

        if n_cols > 5:
            self.l = np.asarray(self.data[:, 5], dtype=np.double)
        else:
            self.l = np.zeros(self.n, dtype=np.double) -1



        # use different offsets for the data
        # off must start from zero
        # -1 values for self.j and self.l mean that these will not be used

        # self.a = np.asarray(self.data[:, 5], dtype=np.double)
        # Flag for activity fit: we can choose to fit the same RV amplitude
        # of the activity signal for all the datasets, use dfferent values
        # for each dataset or exclude some of the datasets

        # Model for RV systematics
        self.n_o = np.max(self.o.astype(np.int64)) + 1
        self.n_j = np.max(self.j.astype(np.int64)) + 1
        self.n_l = np.max(self.l.astype(np.int64)) + 1
        # self.n_a = np.max(self.a.astype(np.int64)) + 1

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

    def print_vars(self, mc, theta):
        for param in ['offset', 'jitter', 'linear']:
            id_var = mc.variable_list[self.name_ref][param]
            if np.size(id_var) == 0:
                continue
            if np.size(id_var) == 1:
                mc.pam_names[id_var] = self.name_ref[:-4] + '_' + param
            else:
                for ii in id_var:
                    mc.pam_names[ii] = self.name_ref[:-4] + '_' + param + '_' + repr(ii - id_var[0])

            print self.name_ref, param, ' : ', theta[mc.variable_list[self.name_ref][param]]


# Transit times cannot be trated as regular datasets since they miss the independent value - they are just
#a collection of epochs. So a different Class to treat them is created. Moreover, each dataset is specific
# to a given planet, while the other datasets
# class TransitCentralTimes(Dataset):
#     def __init__(self, planet_name, input_file):
#
#         self.kind = 'Tcent'
#         self.models = ['Tcent']
#         # 'RV', 'PHOT', 'ACT'...
#         self.name_ref = input_file
#         self.planet_name = planet_name
#
#         self.deltaT = 1.10
#
#         print 'Opening: ', input_file
#         self.data = np.atleast_2d(np.loadtxt(input_file))
#
#         self.x = np.asarray(self.data[:, 0], dtype=np.double)
#         self.e = np.asarray(self.data[:, 1], dtype=np.double)
#
#         self.n = np.size(self.x)
#         self.model = np.zeros(self.n, dtype=np.double)
#
#         n_cols = np.size(self.data, axis=1)
#         if n_cols > 2:
#             self.j = np.asarray(self.data[:, 2], dtype=np.double)
#         else:
#             self.j = np.zeros(self.n, dtype=np.double)
#         # fit for different time measurement jitters (different instruments for the same planet)
#
#         self.n_j = np.max(self.j.astype(np.int64)) + 1
#
#         self.Tref = np.mean(self.x, dtype=np.double)
#         self.x0 = self.x - self.Tref
#
#         print 'N = ', self.n
#         print 'N jitter = ', self.n_j
#         print
#
#         self.j_mask = np.zeros([self.n, self.n_j], dtype=bool)
#         for ii in xrange(0, self.n_j):
#             self.j_mask[(abs(self.j - ii) < 0.1), ii] = True
#
#         self.n_o = 0
#         self.o_mask = np.zeros([self.n, self.n_o], dtype=bool)
#
#         self.n_l = 0
#         self.l_mask = np.zeros([self.n, self.n_l], dtype=bool)
#
#         self.model = np.zeros(self.n, dtype=np.double)
#         self.jitter = np.zeros(self.n, dtype=np.double)
#
#     def compute(self, mc, theta):
#         # By default, dataset.name_ref == planet_name
#         period, _, f, e, o = mc.pcv.convert(self.planet_name, theta)
#         model = np.rint(self.x0 / period) * period + kp.kepler_Tcent_T0P(period, f, e, o)
#         #print self.x0, model, self.x0-model, period, f, e, o
#         return model
#
#     def model_logchi2(self):
#         # boundaries in Tcent are specific of the dataset and not of a common
#         # parameter for different dataset. The check can be internal
#         #if np.sum(np.abs(self.x0 - self.model) < self.deltaT) < self.n:
#         #    return -np.inf
#         env = 1.0 / (self.e ** 2.0 + self.jitter ** 2.0)
#         return -0.5 * (np.sum((self.x0 - self.model) ** 2 * env - np.log(env)))

#Without jitter
class TransitCentralTimes(Dataset):
    def __init__(self, planet_name, input_file):

        self.kind = 'Tcent'
        self.models = ['Tcent']
        # 'RV', 'PHOT', 'ACT'...
        self.name_ref = input_file
        self.planet_name = planet_name

        self.deltaT = 1.10

        print 'Opening: ', input_file
        self.data = np.atleast_2d(np.loadtxt(input_file))

        self.x = np.asarray(self.data[:, 0], dtype=np.double)
        self.e = np.asarray(self.data[:, 1], dtype=np.double)

        self.n = np.size(self.x)
        self.model = np.zeros(self.n, dtype=np.double)

        self.Tref = np.mean(self.x, dtype=np.double)
        self.x0 = self.x - self.Tref

        print 'N = ', self.n
        print

        self.n_j = 0
        self.j_mask = np.zeros([self.n, self.n_j], dtype=bool)

        self.n_o = 0
        self.o_mask = np.zeros([self.n, self.n_o], dtype=bool)

        self.n_l = 0
        self.l_mask = np.zeros([self.n, self.n_l], dtype=bool)

        self.model = np.zeros(self.n, dtype=np.double)
        self.jitter = np.zeros(self.n, dtype=np.double)

    def compute(self, mc, theta):
        # By default, dataset.name_ref == planet_name
        period, _, f, e, o = mc.pcv.convert(self.planet_name, theta)
        model = np.rint(self.x0 / period) * period + kp.kepler_Tcent_T0P(period, f, e, o)
        return model

    def model_logchi2(self):
        # boundaries in Tcent are specific of the dataset and not of a common
        # parameter for different dataset. The check can be internal
        # if np.sum(np.abs(self.x0 - self.model) < self.deltaT) < self.n:
        #    return -np.inf
        env = 1.0 / (self.e ** 2.0)
        return -0.5 * (np.sum((self.x0 - self.model) ** 2 * env - np.log(env)))

    def print_vars(self, mc, theta):
        period, _, f, e, o = mc.pcv.convert(self.planet_name, theta)
        model = np.rint(self.x0 / period) * period + kp.kepler_Tcent_T0P(period, f, e, o)

        print 'Tc ', self.planet_name
        for ii in xrange(0,self.n):
                print 'Input Tc: ',  self.x0[ii] , '  Model Tc: ', model[ii],\
                    '  Diff: ', model[ii]-self.x0[ii]


class ModelContainer:
    def __init__(self):
        self.dataset_list = []
        self.n_datasets = 0
        self.scv = SinusoidsCommonVariables()
        self.pcv = PlanetsCommonVariables()
        self.gcv = GaussianProcessCommonVariables()

        # pyde/emcee variables
        self.ngen = 0
        self.nsteps = 0
        self.nburn = 0
        self.npop_multi = 0
        self.nwalkers = 0
        self.thin = 1

        self.model_list = []
        self.bounds_list = []

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

        self.ndim = 0

        for dataset in self.dataset_list:
            dataset.define_bounds(self)

        if 'kepler' in self.model_list:
            self.pcv.define_bounds(self)

        if 'sinusoids' in self.model_list:
            self.scv.define_bounds(self)

        if 'gaussian' in self.model_list:
            self.gcv.define_bounds(self)

        self.bounds = np.asarray(self.bounds_list)

    def check_bounds(self, theta):
        for ii in xrange(0, self.ndim):
            if not (self.bounds[ii, 0] < theta[ii] < self.bounds[ii, 1]):
                return False
        for planet_name in self.pcv.name_ref:
            e = self.pcv.variables[planet_name]['e'](theta, self.pcv.fixed, self.pcv.var_list[planet_name]['e'])
            if not self.pcv.bounds[planet_name]['e'][0] <= e < self.pcv.bounds[planet_name]['e'][1]:
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
                    dataset.model += self.pcv.compute(theta, dataset, planet_name)

            if 'sinusoids' in dataset.models:
                dataset.model += self.scv.compute(self, theta, dataset)

            if 'Tcent' in dataset.models:
                dataset.model += dataset.compute(self, theta)

            # Gaussian Process check MUST be the last one or the program will fail
            if 'gaussian' in dataset.models:
                chi2_out += self.gcv.lnlk_compute(self, theta, dataset)
            else:
                chi2_out += dataset.model_logchi2()

        return chi2_out

    def pymultinest_priors(self, theta, ndim, nparams):
        theta = (self.bounds[:, 1]-self.bounds[:, 0])*theta
        return theta

    def pymultinest_call(self, theta, ndim, nparams, lnew):
        if not self.check_bounds(theta):
            return -0.5e10

        chi2_out = 0.0

        for dataset in self.dataset_list:
            dataset.model_reset()
            dataset.model_offset(theta[self.variable_list[dataset.name_ref]['offset']])
            dataset.model_jitter(theta[self.variable_list[dataset.name_ref]['jitter']])
            dataset.model_linear(theta[self.variable_list[dataset.name_ref]['linear']])

            if 'kepler' in dataset.models:
                for planet_name in self.pcv.name_ref:
                    dataset.model += self.pcv.compute(theta, dataset, planet_name)

            if 'sinusoids' in dataset.models:
                dataset.model += self.scv.compute(self, theta, dataset)

            if 'Tcent' in dataset.models:
                dataset.model += dataset.compute(self, theta)

            # Gaussian Process check MUST be the last one or the program will fail
            if 'gaussian' in dataset.models:
                chi2_out += self.gcv.lnlk_compute(self, theta, dataset)
            else:
                chi2_out += dataset.model_logchi2()

        return chi2_out

    def results_resumen(self, theta):
        # Function with two goals:
        # * Unfold and print out the output from theta
        # * give back a parameter name associated to each value in the result array

        self.pam_names = self.ndim * ['']

        for dataset in self.dataset_list:
            dataset.print_vars(self, theta)
            print

        if 'kepler' in self.model_list:
            self.pcv.print_vars(self, theta)
            print

        if 'sinusoids' in self.model_list:
            self.scv.print_vars(self, theta)
            print

        if 'gaussian' in self.model_list:
            self.gcv.print_vars(self, theta)
            print

    def rv_make_model(self, theta, x_range, x_phase):
        # it return the RV model for a single planet, after removing the activity from the RV curve and removing
        # the offsets between the datasets

        model_actv = {}
        model_plan = {}
        model_orbs = {}
        model_dsys = {}

        model_orbs['BJD'] = x_range*0.0
        model_orbs['pha'] = x_phase*0.0
        model_plan['BJD'] = {}
        model_plan['pha'] = {}

        # computing the orbit for the full dataset
        for planet_name in self.pcv.name_ref:
            pams = self.pcv.convert(planet_name, theta)
            model_plan['BJD'][planet_name] = self.pcv.model_kepler(pams, x_range-self.Tref)
            model_orbs['BJD'] += model_plan['BJD'][planet_name]
            model_plan['pha'][planet_name] = self.pcv.model_kepler(pams, x_phase*pams[0])
            model_orbs['pha'] += model_plan['pha'][planet_name]

        for dataset in self.dataset_list:

            model_actv[dataset.name_ref] = np.zeros(dataset.n)
            model_orbs[dataset.name_ref] = np.zeros(dataset.n)
            model_plan[dataset.name_ref] = {}

            dataset.model_reset()
            dataset.model_offset(theta[self.variable_list[dataset.name_ref]['offset']])
            dataset.model_linear(theta[self.variable_list[dataset.name_ref]['linear']])

            model_dsys[dataset.name_ref] = dataset.model

            if 'kepler' in dataset.models:
                for planet_name in self.pcv.name_ref:
                    model_plan[dataset.name_ref][planet_name] = self.pcv.compute(theta, dataset, planet_name)
                    model_orbs[dataset.name_ref] += model_plan[dataset.name_ref][planet_name]

            if 'sinusoids' in dataset.models:
                model_actv[dataset.name_ref] += self.scv.compute(self, theta, dataset)

        return model_dsys, model_plan, model_orbs, model_actv

    # This function recenters the bounds limits for circular variables
    # Also, it extends the range of a variable if the output of PyDE is a fixed number
    def recenter_bounds(self, pop_mean, population):
        ind_list = []

        if 'kepler' in self.model_list:
            for planet_name in self.pcv.name_ref:

                if 'esino' in self.variable_list[planet_name]:
                    esino_list = self.variable_list[planet_name]['esino']
                    ecoso_list = self.variable_list[planet_name]['ecoso']
                    e_pops = population[:, esino_list] ** 2 + population[:, ecoso_list] ** 2
                    o_pops = np.arctan2(population[:, esino_list], population[:, ecoso_list], dtype=np.double)
                    #e_mean = (self.pcv.bounds[planet_name]['e'][0] + self.pcv.bounds[planet_name]['e'][1]) / 2.
                    for ii in xrange(0, self.nwalkers):
                        if not self.pcv.bounds[planet_name]['e'][0]+0.02 <= e_pops[ii] < \
                                self.pcv.bounds[planet_name]['e'][1] - 0.02:
                            e_random = np.random.uniform(self.pcv.bounds[planet_name]['e'][0],
                                                         self.pcv.bounds[planet_name]['e'][1])
                            population[ii, esino_list] = np.sqrt(e_random) * np.sin(o_pops[ii])
                            population[ii, ecoso_list] = np.sqrt(e_random) * np.cos(o_pops[ii])

                if 'phase' in self.variable_list[planet_name]:
                    ind_list.append(self.variable_list[planet_name]['phase'])

                if 'o' in self.variable_list[planet_name]:
                    ind_list.append(self.variable_list[planet_name]['o'])

        if 'sinusoids' in self.model_list:
            for jj in range(0, self.scv.n_seasons):
                ind_list.extend(self.variable_list[self.scv.season_name[jj] + '_pha'])
            for dataset in self.dataset_list:
                for jj in range(0, self.scv.n_seasons):
                    if dataset.season_flag[jj]:
                        # ind_list.extend(self.variable_list[dataset.name_ref][self.scv.season_name[jj] + '_amp'])
                        if self.scv.use_offset[dataset.kind]:
                            ind_list.append(self.variable_list[dataset.kind][self.scv.season_name[jj] + '_off'])

        if np.size(ind_list) > 0:
            tmp_range = (self.bounds[:, 1] - self.bounds[:, 0]) / 2
            for var_ind in ind_list:
                self.bounds[var_ind, :] = pop_mean[var_ind] + [-tmp_range[var_ind], tmp_range[var_ind]]
                fix_sel = (population[:, var_ind] <= self.bounds[var_ind, 0]) | (
                    population[:, var_ind] >= self.bounds[var_ind, 1])
                population[fix_sel, var_ind] = pop_mean[var_ind]

        for ii in xrange(0, self.ndim):
            if np.amax(population[:, ii]) - np.amin(population[:, ii]) < 10e-7:
                range_restricted = (self.bounds[ii, 1] - self.bounds[ii, 0]) / 100.
                min_bound = np.maximum((pop_mean[ii] - range_restricted / 2.0), self.bounds[ii, 0])
                max_bound = np.minimum((pop_mean[ii] + range_restricted / 2.0), self.bounds[ii, 1])
                population[:, ii] = np.random.uniform(min_bound, max_bound, self.nwalkers)


def yaml_parser(file_conf, mc):
    stream = file(file_conf, 'r')
    config_in = yaml.load(stream)

    conf = config_in['Inputs']
    for counter in conf:
        print conf[counter]['Kind'], conf[counter]['File'], conf[counter]['Models']
        mc.dataset_list.append(Dataset(counter, conf[counter]['Kind'], conf[counter]['File'], conf[counter]['Models']))

        if counter == 0:
            mc.Tref = mc.dataset_list[0].Tref
        else:
            mc.dataset_list[counter].common_Tref(mc.Tref)

    mc.planet_name = config_in['Output']

    if 'Planets' in config_in:
        conf = config_in['Planets']
        for counter in conf:
            planet_name = 'Planet_' + repr(counter)
            planet_conf = conf[counter]
            mc.pcv.add_planet(planet_name)

            if 'Boundaries' in planet_conf:
                bound_conf = planet_conf['Boundaries']
                for var in bound_conf:
                    if var == 'P' or var == 'K':
                        mc.pcv.bounds[planet_name]['log'+var] = np.log2(np.asarray(bound_conf[var], dtype=np.double))
                    else:
                        mc.pcv.bounds[planet_name][var] = np.asarray(bound_conf[var], dtype=np.double)

            if 'Fixed' in planet_conf:
                fixed_conf = planet_conf['Fixed']
                for var in fixed_conf:
                    mc.pcv.fix_list[planet_name][var] = np.asarray(fixed_conf[var], dtype=np.double)

            if 'Orbit' in planet_conf and planet_conf['Orbit'] == 'circular':
                mc.pcv.circular[planet_name] = True
                mc.pcv.fix_list[planet_name]['e'] = 0.00000
                mc.pcv.fix_list[planet_name]['o'] = 0.00000

            if 'Tcent' in planet_conf:
                mc.dataset_list.append(TransitCentralTimes(planet_name, planet_conf['Tcent']))
                mc.dataset_list[-1].common_Tref(mc.Tref)

    if 'Sinusoids' in config_in:
        conf = config_in['Sinusoids']
        mc.scv.Prot_bounds = np.asarray(conf['Prot'], dtype=np.double)
        for counter in conf['Seasons']:
            mc.scv.add_season_range(np.asarray(conf['Seasons'][counter][:2], dtype=np.double),
                                    conf['Seasons'][counter][2:])

            # if 'Phase_dataset' in conf:
            #    # Additional activity indicators associated to RVs must the same order and number of sinusoids,
            #    #
            #    mc.scv.phase = np.asarray(conf['Phase_dataset'], dtype=np.int64)

    if 'Gaussian' in config_in:
        conf = config_in['Gaussian']
        for name_ref in conf:
            if name_ref == 'Common':
                mc.gcv.add_dataset(name_ref)
                for var in conf[name_ref] :
                    if var in mc.gcv.list_pams_george:
                        var_gp = var
                        converted = np.asarray(conf[name_ref][var], dtype=np.double)
                    else:
                        var_gp = mc.gcv.list_pams_corr[var]
                        converted = mc.gcv.convert_val(var, np.asarray(conf[name_ref][var], dtype=np.double))

                    if np.size(conf[name_ref][var]) == 1:
                        mc.gcv.fix_list[var_gp] = converted
                    else:
                        mc.gcv.bounds[var_gp] = np.sort(converted)
            else:
                #num_ref = np.asarray(name_ref, dtype=np.double)
                # dataset_name = mc.dataset_list[num_ref].name_ref
                dataset_name = mc.dataset_list[name_ref].name_ref
                mc.gcv.add_dataset(dataset_name)

                for var in conf[name_ref]:
                    if var in mc.gcv.list_pams_george:
                        var_gp = var
                        converted = np.asarray(conf[name_ref][var], dtype=np.double)
                    else:
                        var_gp = mc.gcv.list_pams_corr[var]
                        converted = mc.gcv.convert_val(var, np.asarray(conf[name_ref][var], dtype=np.double))

                    if np.size(conf[name_ref][var]) == 1:
                        mc.gcv.fix_list[dataset_name][var_gp] = converted
                    else:
                        mc.gcv.bounds[dataset_name][var_gp] = np.sort(converted)

    if 'Tref' in config_in:
        mc.Tref = np.asarray(config_in['Tref'])
        for dataset in mc.dataset_list:
            dataset.common_Tref(mc.Tref)

    if 'Ngen' in config_in['emcee']:
        mc.ngen = np.asarray(config_in['emcee']['Ngen'], dtype=np.double)

    if 'Nsteps' in config_in['emcee']:
        mc.nsteps = np.asarray(config_in['emcee']['Nsteps'], dtype=np.int64)

    if 'Nburn' in config_in['emcee']:
        mc.nburn = np.asarray(config_in['emcee']['Nburn'], dtype=np.int64)

    if 'Npop_mult' in config_in['emcee']:
        mc.npop_mult = np.asarray(config_in['emcee']['Npop_mult'], dtype=np.int64)

    if 'Thin' in config_in['emcee']:
        mc.thin = np.asarray(config_in['emcee']['Thin'], dtype=np.int64)

    if 'Recenter_Bounds' in config_in['emcee']:
        # required to avoid a small bug in the code
        # if the dispersion of PyDE walkers around the median value is too broad,
        # then emcee walkers will start outside the bounds, causing an error
        mc.recenter_bounds_flag = config_in['emcee']['Recenter_Bounds']

    if 'Star_Mass' in config_in:
        mc.star_mass_val = np.asarray(config_in['Star_Mass'][0], dtype=np.double)
        mc.star_mass_err = np.asarray(config_in['Star_Mass'][1], dtype=np.double)

    mc.model_setup()


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
        #    scv.n_seasons = np.asarray(info[1], dtype=np.int64)

        # Each chunk has a different rotational period, still within the boundaries
        # defined in  Prot, and it is defined as an interval in BJD
        if line.find("Period_" + `pr_row` + " ") > -1:
            P0_tmp = np.asarray(info[1], dtype=np.double)
            P1_tmp = np.asarray(info[2], dtype=np.double)
            mc.scv.add_season_range([P0_tmp, P1_tmp], info[3:])
            # print season_range[pr_row,1]
            # print season_range[pr_row ][0]
            # print season_range[pr_row ][1]

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
    # if mc.scv.n_seasons != pr_row:
    #    print 'PROBLEM'

    mc.model_setup()
