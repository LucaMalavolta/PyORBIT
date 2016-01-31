import numpy as np
import kepler_exo as kp

# A class for common variables (number of planets, number of sinusoids...)
# Info to unpack the variables inside emcee must be included here
# Physical effects must be included here

class Planets_Common_Variables:
    def __init__(self):
        self.n_planets = 0
        self.name_ref = []

        self.n_orbpams = {}
        self.logP_bounds={}
        self.logK_bounds={}
        self.e_bounds={}

    def add_planet(self,name_ref,n_orbs):
        self.n_planets += 1
        self.name_ref.append(name_ref)
        self.n_orbpams[name_ref]= n_orbs

    #def convert_params(self,in_pams):
    #    out_pams = np.zeros(5,dtype=np.double)
    #    out_pams[0] = 10.0 ** (in_pams[0]) #P
    #    out_pams[1] = 10.0 ** (in_pams[1]) #K
    #    out_pams[2] = in_pams[2]
    #    if np.size(in_pams)==5:
    #        out_pams[3] = in_pams[3] ** 2 + in_pams[4] ** 2
    #        out_pams[4] = np.arctan2(in_pams[3], in_pams[4])
    #    return out_pams

    def convert_params(self,logP,logK,phase,esino,ecoso):
        P = np.exp2(logP)
        K = np.exp2(logK)
        e = np.square(ecoso) + np.square(esino)
        o = np.arctan2(esino, ecoso)
        return P, K, phase, e, o

    def model_kepler(self,orbit_pams,x0):
        P, K, phase, e, omega = orbit_pams
        rv_out = kp.kepler_RV_T0P(x0, phase, P, K, e, omega)
        return rv_out


class Sinusoids_Common_Variables:
    def __init__(self):

        self.n_datasets = 0

        self.offset_coherence = True #Same phase offset across seasons
        self.phase_coherence = False
        self.phase_sincro = False

        self.n_pha     = 0

        self.period_sel = False
        self.n_periods = 1
        self.period_name =['Period_0']
        self.period_list = [0., 5000000.0]
        self.period_range = np.asarray(self.period_list,dtype=np.double)

        #self.Prot_bounded = False
        #self.pha_bounded = False

        self.Prot_bounds = [0., 30.0]
        self.pof_bounds = [0., 1.0]
        self.pha_bounds = [0., 1.0]

        self.phase_list = []

        self.use_offset = {}
        self.offset_skip_first = True

    def add_period_range(self,range_input,phase_input):
        #Reset default values at first call
        if self.period_sel==False:
            self.n_periods = 0
            self.period_list=[]
            self.period_name=[]
            self.period_sel=True

        self.period_name.append('Period_'+`self.n_periods`)
        self.period_list.append(range_input)
        self.period_range=np.asarray(self.period_list,dtype=np.double)

        self.phase_list.append(phase_input)
        self.phase = np.asarray(self.phase_list,dtype=np.int64)
        self.n_pha = np.amax(self.phase,axis=1) #maximum value for each period
        self.n_pha_max = np.amax(self.n_pha) #maximum value for each period
        self.n_periods += 1
        return

    def add_phase_offset(self,dataset):
        add_bound = False
        try:
            value = self.use_offset[dataset.kind]
        except KeyError:
            # Key is not present
            if self.offset_skip_first:
                self.use_offset[dataset.kind]=False
                self.offset_skip_first = False
            else:
                self.use_offset[dataset.kind]=True
                add_bound = True
        return add_bound

    def setup_model_sinusoids(self,dataset):
        dataset.n_amp  = self.phase[:,dataset.ind]
        dataset.p_mask = np.zeros([self.n_periods,dataset.n], dtype=bool)
        dataset.n_periods = 0
        dataset.periods_flag = np.zeros(self.n_periods,dtype=bool)

        for ii in xrange(0, self.n_periods):
            p_sel = (self.period_range[ii,0] < dataset.x) & (dataset.x < self.period_range[ii,1])
            if (np.sum(p_sel)==0) :
                'No data points within the specified range'
                break
            else:
                dataset.periods_flag[ii] = True
                dataset.p_mask[ii,:] = p_sel[:]
                dataset.n_periods += 1

        if (dataset.kind=='RV'  ): dataset.sinamp_bounds = np.asarray([0., 100.], dtype=np.double)
        if (dataset.kind=='Phot'): dataset.sinamp_bounds = np.asarray([0., 0.5], dtype=np.double)
        if (dataset.kind=='FWHM'): dataset.sinamp_bounds = np.asarray([0., 2000.], dtype=np.double)
        if (dataset.kind=='BIS' ): dataset.sinamp_bounds = np.asarray([0., 60.], dtype=np.double)
        return

    def model_sinusoids(self,dataset,Prot,amp,pha,off):
        #np.size(SinAmp)==np.sum(n_amp) ???? What if an activity dataset is missing?
        #cycle for np.sum(rva_mask[:, jj]) <= 0 ??

        #Prot= Rotational periodo of the star
        #n_amp = number of sinusoids in the fit
        #amp = amplitudes
        #sin = phase of each sinusoid
        #off = overall offset of the sinusoids
        #sel = if the model has to be restricted to a specific temporal range
        model = np.zeros(dataset.n,dtype=np.double)
        xph = (dataset.x0 / Prot) % 1
        har = np.arange(1,np.size(amp,axis=1)+1,1.,dtype=np.double)

        for ii in xrange(0,dataset.n_periods):
            for jj in xrange (0,np.size(pha[ii,:])):
                model +=  dataset.p_mask[ii,:] * amp[ii,jj] * np.sin((har[jj] * xph + pha[ii,jj] + off[ii])* 2. * np.pi)
        return model


class Dataset:
    def __init__(self,ind,kind,input_file,models):
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
        ## fit for different RV jitters
        self.o = np.asarray(self.data[:, 4], dtype=np.double)
        self.l = np.asarray(self.data[:, 5], dtype=np.double)
        ## use different offsets for the data
        ## off must start from zero
        ## -1 values for self.j and self.l mean that these

        # self.a = np.asarray(self.data[:, 5], dtype=np.double)
        # Flag for activity fit: we can choose to fit the same RV amplitude
        # of the activity signal for all the datasets, use dfferent values
        #for each dataset or exclude some of the datasets

        ## Model for RV systematics
        self.n_o = np.max(self.o.astype(np.int64)) + 1
        self.n_j = np.max(self.j.astype(np.int64)) + 1
        self.n_l = np.max(self.l.astype(np.int64)) + 1
        #self.n_a = np.max(self.a.astype(np.int64)) + 1

        self.n = np.size(self.x)
        print 'N = ', self.n
        print 'N jitter = ', self.n_j
        print 'N offset = ', self.n_o
        print 'N linear = ', self.n_l
        #print 'N activ. = ', self.n_a
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

    def common_Tref(self,Tref_in):
        self.Tref = Tref_in
        self.x0 = self.x - self.Tref
        return

    def model_reset(self):
        self.model =np.zeros(self.n,dtype=np.double)
        self.jitter=np.zeros(self.n,dtype=np.double)
        return

    def model_offset(self,off_in):
        off=np.atleast_1d(off_in)
        for ii in xrange(0,self.n_o):
            self.model[self.o_mask[:,ii]]+= off[ii]

    def model_linear(self,m_in):
        m=np.atleast_1d(m_in)
        for ii in xrange(0,self.n_l):
            self.model[self.l_mask[:,ii]]+= m[ii]*sel.x0[self.l_mask[:,ii]]

    def model_jitter(self,jit_in):
        jit=np.atleast_1d(jit_in)
        for ii in xrange(0, self.n_j):
            self.jitter[self.j_mask[:, ii]] = jit[ii]

    def model_logchi2(self):
        env = 1.0 / (self.e ** 2.0 + self.jitter ** 2.0)
        return -0.5 * ( np.sum((self.y - self.model) ** 2 * env - np.log(env)))


class Model_Container:
    def __init__(self):
        self.dataset_list = []
        self.scv = Sinusoids_Common_Variables()
        self.pcv = Planets_Common_Variables()
        self.ndim = 0

        #pyde/emcee variables
        self.ngen = 0
        self.nsteps = 0
        self.nburn  = 0
        self.npop_multi = 0
        self.nwalkers = 0
        self.thin = 0

        self.recenter_bounds_flag = True

    def Model_setup(self):
        self.n_datasets = np.size(self.dataset_list)
        self.model_list=[]
        for dataset in self.dataset_list:
            self.scv.setup_model_sinusoids(dataset)
            dataset.model_reset()
            for data_model in dataset.models:
                if not (data_model in self.model_list):
                    self.model_list.append(data_model)

    def create_bounds(self):
        #This routine creates the boundary array and at the same time
        #creates a dictionary with the name of the arrays and their
        #positions in bounds/theta array so that they can be accessed
        #without using nested counters
        bounds_list=[]

        self.variable_list={}
        var_count=0

        for dataset in self.dataset_list:
            for jj in xrange(0,dataset.n_j):
                bounds_list.append([0., 50 * np.max(dataset.e)]) #bounds for jitter
            for jj in xrange(0,dataset.n_o):
                bounds_list.append([np.min(dataset.y), np.max(dataset.y)])
            for jj in xrange(0,dataset.n_l):
                bounds_list.append([-1.,1.])

            self.variable_list[dataset.name_ref]={}
            self.variable_list[dataset.name_ref]['jitter'] = xrange(var_count,var_count+dataset.n_j)
            var_count += dataset.n_j
            self.variable_list[dataset.name_ref]['offset'] = xrange(var_count,var_count+dataset.n_o)
            var_count += dataset.n_o
            self.variable_list[dataset.name_ref]['linear'] = xrange(var_count,var_count+dataset.n_l)
            var_count += dataset.n_l

        if 'kepler' in self.model_list:
            for planet_name in  self.pcv.name_ref:
                bounds_list.append(self.pcv.logP_bounds[planet_name][:])
                bounds_list.append(self.pcv.logK_bounds[planet_name][:])
                bounds_list.append([0,2*np.pi])
                if (self.pcv.n_orbpams[planet_name] == 5):
                    bounds_list.append([-1.0, 1.0])
                    bounds_list.append([-1.0, 1.0])

                self.variable_list[planet_name]={}
                if (self.pcv.n_orbpams[planet_name] == 5):
                    self.variable_list[planet_name]['logP'] =var_count
                    self.variable_list[planet_name]['logK'] =var_count+1
                    self.variable_list[planet_name]['phase']=var_count+2
                    self.variable_list[planet_name]['esino']=var_count+3
                    self.variable_list[planet_name]['ecoso']=var_count+4
                    self.variable_list[planet_name]['kepler_pams'] = np.arange(var_count,var_count+5,1)
                    var_count += 5
                else:
                    self.variable_list[planet_name]['logP'] =var_count
                    self.variable_list[planet_name]['logK'] =var_count+1
                    self.variable_list[planet_name]['phase']=var_count+2
                    self.variable_list[planet_name]['kepler_pams'] = np.arange(var_count,var_count+3,1)
                    var_count += 3

        if 'sinusoids' in self.model_list:
            bounds_list.append(self.scv.Prot_bounds[:])
            self.variable_list['Prot']=var_count
            var_count += 1

            for jj in range(0,self.scv.n_periods):
                for kk in range(0, self.scv.n_pha[jj]):
                    bounds_list.append(self.scv.pha_bounds[:])
                self.variable_list[self.scv.period_name[jj]+'_pha']= xrange(var_count,var_count+self.scv.n_pha[jj])
                var_count += self.scv.n_pha[jj]

            for dataset in self.dataset_list:

                #two nested case:
                # 1) dataset has an offset or not (if it is the reference offset)
                # 2) the offset is the same for every season or not
                add_bound = self.scv.add_phase_offset(dataset)
                if (add_bound and self.scv.offset_coherence):
                    bounds_list.append(self.scv.pof_bounds[:])
                    for jj in range(0,self.scv.n_periods):
                        if (dataset.periods_flag[jj]):
                            self.variable_list[dataset.name_ref][self.scv.period_name[jj]+'_off']= var_count
                    var_count += 1
                if (add_bound and (not self.scv.offset_coherence)):
                    for jj in range(0,self.scv.n_periods):
                        if (dataset.periods_flag[jj]):
                            bounds_list.append(self.scv.pof_bounds[:])
                            self.variable_list[dataset.name_ref][self.scv.period_name[jj]+'_off']= var_count
                            var_count += 1

                for jj in range(0,self.scv.n_periods):
                    if (dataset.periods_flag[jj]):

                        for kk in xrange(0,dataset.n_amp[jj]):
                            bounds_list.append(dataset.sinamp_bounds)

                        self.variable_list[dataset.name_ref][self.scv.period_name[jj]+'_amp']= \
                          xrange(var_count,var_count+dataset.n_amp[jj])
                        var_count += dataset.n_amp[jj]

        self.bounds = np.asarray(bounds_list)
        self.ndim = var_count

    def check_bounds(self,theta):
        #former lnprior(theta)
        for ii in xrange(0,self.ndim):
            if not (self.bounds[ii,0] < theta[ii] < self.bounds[ii,1]):
                return False
        for planet_name in self.pcv.name_ref:
            e = theta[self.variable_list[planet_name]['esino']]**2 + theta[self.variable_list[planet_name]['ecoso']]**2
            #omega = np.arctan2(theta[self.variable_list[self.pcv.name_ref[pp]]['ecoso']], theta[self.variable_list[self.pcv.name_ref[pp]]['esino']])
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
                    in_pams=theta[self.variable_list[planet_name]['kepler_pams']]
                    if (self.pcv.n_orbpams[planet_name]==5):
                        conv_pams = self.pcv.convert_params(in_pams[0],in_pams[1],in_pams[2],in_pams[3],in_pams[4])
                    else:
                        conv_pams = self.pcv.convert_params(in_pams[0],in_pams[1],in_pams[2],0.00,0.00)
                    dataset.model += self.pcv.model_kepler(conv_pams,dataset.x0)

            if 'sinusoids' in dataset.models:
                #Prot and pha_in could be brought out from the llop, but I would not work
                #for plaet-only analysis
                Prot   = theta[self.variable_list['Prot']]
                pha_in = np.zeros([self.scv.n_periods,self.scv.n_pha_max],dtype=np.double)
                off_in = np.zeros([self.scv.n_periods],dtype=np.double)
                amp_in = np.zeros([self.scv.n_periods,self.scv.n_pha_max],dtype=np.double)
                for jj in range(0,self.scv.n_periods):
                    pha_in[jj,:self.scv.n_pha[jj]] = theta[self.variable_list[self.scv.period_name[jj]+'_pha']]


                    if (dataset.periods_flag[jj]):
                        if (self.scv.use_offset[dataset.kind]):
                            off_in[:] = theta[self.variable_list[dataset.name_ref][self.scv.period_name[jj]+'_off']]
                        amp_in[jj,:dataset.n_amp[jj]] = \
                            theta[self.variable_list[dataset.name_ref][self.scv.period_name[jj]+'_amp']]

                dataset.model += self.scv.model_sinusoids(dataset,Prot,amp_in,pha_in,off_in)

            chi2_out += dataset.model_logchi2()

        return chi2_out

    def results_resumen(self, theta):

        for dataset in self.dataset_list:
            print dataset.name_ref, ' offset: ', theta[self.variable_list[dataset.name_ref]['offset']]
            print dataset.name_ref, ' jitter: ', theta[self.variable_list[dataset.name_ref]['jitter']]
            print dataset.name_ref, ' linear: ', theta[self.variable_list[dataset.name_ref]['linear']]

        if 'kepler' in self.model_list:

            for planet_name in self.pcv.name_ref:
                print planet_name, theta[self.variable_list[planet_name]['kepler_pams']]
                in_pams=theta[self.variable_list[planet_name]['kepler_pams']]
                if (self.pcv.n_orbpams[planet_name]==5):
                    conv_pams = self.pcv.convert_params(in_pams[0],in_pams[1],in_pams[2].in_pams[3],in_pams[4])
                else:
                    conv_pams = self.pcv.convert_params(in_pams[0],in_pams[1],in_pams[2],0.00,0.00)
                print planet_name,' pams: ',conv_pams[:]

        if 'sinusoids' in self.model_list:
            #Prot and pha_in could be brought out from the llop, but I would not work
            #for plaet-only analysis
            print 'Prot ', theta[self.variable_list['Prot']]
            for jj in range(0,self.scv.n_periods):
                print self.scv.period_name[jj],'_pha', theta[self.variable_list[self.scv.period_name[jj]+'_pha']]

            for dataset in self.dataset_list:
                for jj in range(0,self.scv.n_periods):
                    if (dataset.periods_flag[jj]):
                        if (self.scv.use_offset[dataset.kind]):
                            print dataset.name_ref,self.scv.period_name[jj],'_off',theta[self.variable_list[dataset.name_ref][self.scv.period_name[jj]+'_off']]
                        print dataset.name_ref,self.scv.period_name[jj],'_amp',theta[self.variable_list[dataset.name_ref][self.scv.period_name[jj]+'_amp']]

    def recenter_bounds(self,pop_mean,population):
        ind_list = []

        if 'kepler' in self.model_list:
            for planet_name in self.pcv.name_ref:
                for name_var in ['esino','ecoso']:
                    ind_list.extend(self.variable_list[planet_name][name_var])

        if 'sinusoids' in self.model_list:
            for jj in range(0,self.scv.n_periods):
                var_ind = self.variable_list[self.scv.period_name[jj]+'_pha']
            for dataset in self.dataset_list:
                for jj in range(0,self.scv.n_periods):
                    ind_list.extend(self.variable_list[dataset.name_ref][self.scv.period_name[jj]+'_amp'])
                    if (dataset.periods_flag[jj]):
                            ind_list.extend(self.variable_list[dataset.name_ref][self.scv.period_name[jj]+'_off'])

        if np.size(ind_list)==0: return population

        for var_ind in ind_list:
            self.bounds[var_ind,:] = pop_mean[var_ind] + self.bounds[var_ind,:]
            fix_sel = (population[:,var_ind]<=self.bounds[var_ind,0]) | (population[:,var_ind]>=self.bounds[var_ind,1])
            population[fix_sel,var_ind] = pop_mean[var_ind]

        return population





def Get_PyOrbit_Input(file_conf,mc):
    file = open(file_conf, 'r')
    row = file.readlines()

    pl_row = 0
    pr_row = 0
    in_row = 0

    for line in row:
        info = line.split()
        if line.find("Input_"+ `in_row` + " ") > -1:
            input_name = info[0]
            input_file = info[1]
            input_kind = info[2]
            input_model = info[3:]
            mc.dataset_list.append(Dataset(in_row,input_kind,input_file,input_model))
            if (in_row>0):
                mc.dataset_list[in_row].common_Tref( mc.dataset_list[0].Tref)
            in_row+=1

        if line.find("Output") > -1:
            mc.planet_name = info[1]
        #if line.find("Nplanets") > -1:
        #    n_planets = np.asarray(info[1], dtype=np.int64)
        if line.find("Planet_" + `pl_row` + " ") > -1:
            mc.pcv.add_planet(info[0],np.asarray(info[1], dtype=np.int64))

            Pmin = np.asarray(info[2], dtype=np.double)
            Pmax = np.asarray(info[3], dtype=np.double)
            mc.pcv.logP_bounds[info[0]]=np.asarray([np.log2(Pmin), np.log2(Pmax)])

            Kmin = np.asarray(info[4], dtype=np.double)
            Kmax = np.asarray(info[5], dtype=np.double)
            mc.pcv.logK_bounds[info[0]]=np.asarray([np.log2(Kmin), np.log2(Kmax)])

            if len(info) > 6:
                mc.pcv.e_bounds[info[0]]=np.asarray([info[6],info[7]], dtype=np.double)
            else:
                mc.pcv.e_bounds[info[0]]=np.asarray([0., 1.],dtype=np.double)
            pl_row += 1
            # if line.find("LongTerm") > -1:
            ## Two possible drift are possible: one caused by the instrument
            ## (and this one is coded in the data through the rv_l flag)
            ## and one due to the star, and common to all the datasets
            ## Use this option in the second case
            ## info= line.split()
            # n_trends = np.asarray(info[2], dtype=np.int64)
        if line.find("Prot") > -1:
            ##info= line.split()
            Pmin = np.asarray(info[2], dtype=np.double)
            Pmax = np.asarray(info[3], dtype=np.double)
            #kind_sin = np.asarray(info[1], dtype=np.int64)

            mc.scv.Prot_bounds = [Pmin, Pmax]

        ## The number of chuncks for the activity fit
        #if line.find("Nperiods") > -1:
        #    scv.n_periods = np.asarray(info[1], dtype=np.int64)

        # Each chunk has a different rotational period, still within the boundaries
        # defined in  Prot, and it is defined as an interval in BJD
        if line.find("Period_" + `pr_row` + " ") > -1:
            P0_tmp = np.asarray(info[1], dtype=np.double)
            P1_tmp = np.asarray(info[2], dtype=np.double)
            mc.scv.add_period_range([P0_tmp,P1_tmp],info[3:])
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
            #scv.n_spec = np.asarray(info[1], dtype=np.int64)
            mc.scv.phase = np.asarray(info[1:], dtype=np.int64)

        if line.find("Recenter_Bounds") > -1:
            # required to avoid a small bug in the code
            # if the dispersion of PyDE walkers around the median value is too broad,
            # then emcee walkers will start outside the bounds, causing an error
            if np.asarray(info[1], dtype=np.int64) < 1:
                mc.recenter_bounds_flag = False


    if (mc.scv.n_periods != pr_row):
        print 'PROBLEM'

