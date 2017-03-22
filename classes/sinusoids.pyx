from common import *


class SinusoidsCommonVariables:
    ''' This is among the first classes I've been writing in Python, and differently from other
    classes in this suite it has not updated (only bug-fixed) while my Python skilled were improving
    So you may find many sub-optimal approaches
    '''

    def __init__(self):

        self.n_datasets = 0

        self.offset_coherence = True  # Same phase offset across seasons
        self.offset_common_id = -1
        self.offset_reference_name = ''
        self.use_offset = {}

        self.phase_coherence = False
        self.phase_sincro = False

        self.n_pha = 0

        # No 'Season' selection in the configuration file, so we make believe the code that all the data
        # is contained within one season
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

        self.prior_kind = {}
        self.prior_pams = {}

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

    def model_setup(self, dataset):
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
        Prot = theta[mc.variable_list['Common']['Prot']]
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

    def initialize(self, mc, theta):
        mc.pam_names[mc.variable_list['Prot']] = 'Prot'

        for jj in range(0, self.n_seasons):
            id_var = mc.variable_list[self.season_name[jj] + '_pha']
            if np.size(id_var) == 0:
                continue
            if np.size(id_var) == 1:
                mc.pam_names[id_var] = self.season_name[jj] + '_pha'
            else:
                for ii in id_var:
                    mc.pam_names[ii] = self.season_name[jj] + '_' + repr(ii - id_var[0]) + '_pha'

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

                    id_var = mc.variable_list[dataset.name_ref][self.season_name[jj] + '_amp']
                    if np.size(id_var) == 0:
                        continue
                    if np.size(id_var) == 1:
                        mc.pam_names[id_var] = dataset.name_ref[:-4] + '_' + self.season_name[jj] + '_amp'
                    else:
                        for ii in id_var:
                            mc.pam_names[ii] = dataset.name_ref[:-4] + '_' + self.season_name[jj] + \
                                               '_' + repr(ii - id_var[0]) + '_amp'

    def print_vars(self, mc, theta):
        # Prot and pha_in could be brought out from the llop, but I would not work
        # for planet-only analysis

        print 'Prot ', theta[mc.variable_list['Prot']]

        for jj in range(0, self.n_seasons):
            id_var = mc.variable_list[self.season_name[jj] + '_pha']
            print self.season_name[jj], '_pha', theta[id_var]

        for dataset in mc.dataset_list:
            for jj in range(0, self.n_seasons):
                if dataset.season_flag[jj]:
                    if self.use_offset[dataset.kind]:
                        print dataset.name_ref, dataset.kind, self.season_name[jj], '_off', \
                                theta[mc.variable_list[dataset.kind][self.season_name[jj] + '_off']]

                    print dataset.name_ref, self.season_name[jj], '_amp', \
                            theta[mc.variable_list[dataset.name_ref][self.season_name[jj] + '_amp']]
        print

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

    def starting_point(self, mc):
        pass
        ### MUST BE FIXED ###

    def model_sinusoids(self, dataset, p_rot, amp, pha, off):
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

        for ii in xrange(0, self.n_seasons):
            for jj in xrange(0, np.size(pha[ii, :])):
                model += dataset.p_mask[ii, :] * amp[ii, jj] * np.sin(
                    (har[jj] * xph + pha[ii, jj] + off[ii]) * 2. * np.pi)
        return model
