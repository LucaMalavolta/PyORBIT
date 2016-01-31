import numpy as np
import emcee
from pyde.de import DiffEvol
import h5py
import kepler_exo as kp

# The structure of the THETA variable:
# if flag_rv:
#    n_orbpams[pp]  -> 3 or 5 parameters for each planet
#    n_rvj
#    n_rvo
#    n_rvl
# if flag_p1:
#    n_p1j
#    n_p1o
#    n_p1l
# if flag_p2:
#    n_p2j
#    n_p2o
#    n_p2l
# if flag_fw:
#    n_fwj
#    n_fwo
#    n_fwl
# if flag_bs:
#    n_bsj
#    n_bso
#    n_bsl
# Now with the kind_soin flag:
# if kind_sin==1:
#    Prot (1)
# if kind_sin==2 :
#    Prot, (n_pof=0), n_pha (tot:2)
#    if flag_rv:
#       n_amp
#    if flag_p1:
#       n_amp
#    if flag_p2:
#       n_amp
#    if flag_fw:
#       n_amp
#    if flag_bs:
#       n_amp
# if kind_sin==3:
#    Prot, n_pof, n_pha=1
#    if flag_rv:
#       n_amp
#    if flag_ph:
#       n_amp
#    if flag_fw:
#       n_amp
#    if flag_bs:
#       n_amp
# if kind_sin==4:
#    Prot, (n_pof=0), n_pha=1,
#    if flag_rv:
#       n_amp
#    if flag_ph:
#       n_amp
#    if flag_fw:
#       n_amp
#    if flag_bs:
#       n_amp
# if kind_sin==5:
#    Prot, n_pof, n_pha
#    if flag_rv:
#       n_amp
#    if flag_ph:
#       n_amp
#    if flag_fw:
#       n_amp
#    if flag_bs:
#       n_amp

def lnlike(theta):
    ## derived_orb = np.empty(0, dtype=np.double)
    derived_orb = []
    ip = 0
    if flag_rv:
        rv_model = np.zeros(rv_n, dtype=np.double)
        rv_j_add = np.zeros(rv_n, dtype=np.double)

        ## Going trough the orbital parameters of each individual planet
        ## to compute the non-interacting Keplerian RV
        for pp in xrange(0, n_planets):
            e     = np.asarray(0., dtype=np.double)
            omega = np.asarray(0., dtype=np.double)

            ## P and K are converted from their logarithmic value
            P = 10.0 ** (theta[ip + 0])
            k = 10.0 ** (theta[ip + 1])
            ## Their values are stored in output to have the chain and correlation plot
            ## in a readable form
            derived_orb.append(P)
            derived_orb.append(k)

            if (n_orbpams[pp] == 5):
                ## for non-null eccentricity, sqrt(e)*sin(omega) and sqrt(e)*cos(omega)
                ## are converted to e and omega, and their value are given in output
                ## for readable display of chains and triangle diagram
                # theta[ip + 3] = sqrt(e)*sin(omega)
                # theta[ip + 4] = sqrt(e)*cos(omega)
                e = theta[ip + 3] ** 2 + theta[ip + 4] ** 2
                omega = np.arctan2(theta[ip + 3], theta[ip + 4])
                derived_orb.append(e)
                derived_orb.append(omega)
            rv_model += kp.kepler_RV_T0P(rv_x0, theta[ip + 2], P, k, e, omega)
            ip += n_orbpams[pp]

        # for ii in xrange(0, n_trends):
        #   rv_model += rv_x0 ** ii * theta[ip + ii]
        # ip += n_trends

        for ii in xrange(0, n_rvj): rv_j_add += theta[ip + ii] * rvj_mask[:, ii]
        ip += n_rvj

        for ii in xrange(0, n_rvo): rv_model += theta[ip + ii] * rvo_mask[:, ii]
        ip += n_rvo

        for ii in xrange(0, n_rvl): rv_model += theta[ip + ii] * rv_x0 * rvl_mask[:, ii]
        ip += n_rvl

    if flag_p1:
        p1_model = np.zeros(p1_n, dtype=np.double)
        p1_j_add = np.zeros(p1_n, dtype=np.double)

        for ii in xrange(0, n_p1j): p1_j_add += theta[ip + ii] * p1j_mask[:, ii]
        ip += n_p1j

        for ii in xrange(0, n_p1o): p1_model += theta[ip + ii] * p1o_mask[:, ii]
        ip += n_p1o

        for ii in xrange(0, n_p1l): p1_model += theta[ip + ii] * p1_x0 * p1l_mask[:, ii]
        ip += n_p1l

    if flag_p2:
        p2_model = np.zeros(p2_n, dtype=np.double)
        p2_j_add = np.zeros(p2_n, dtype=np.double)

        for ii in xrange(0, n_p2j): p2_j_add += theta[ip + ii] * p2j_mask[:, ii]
        ip += n_p2j

        for ii in xrange(0, n_p2o): p2_model += theta[ip + ii] * p2o_mask[:, ii]
        ip += n_p2o

        for ii in xrange(0, n_p2l): p2_model += theta[ip + ii] * p2_x0 * p2l_mask[:, ii]
        ip += n_p2l

    if flag_fw:
        fw_model = np.zeros(fw_n, dtype=np.double)
        fw_j_add = np.zeros(fw_n, dtype=np.double)

        for ii in xrange(0, n_fwj): fw_j_add += theta[ip + ii] * fwj_mask[:, ii]
        ip += n_fwj

        for ii in xrange(0, n_fwo): fw_model += theta[ip + ii] * fwo_mask[:, ii]
        ip += n_fwo

        for ii in xrange(0, n_fwl): fw_model += theta[ip + ii] * fw_x0 * fwl_mask[:, ii]
        ip += n_fwl

    if flag_bs:
        bs_model = np.zeros(bs_n, dtype=np.double)
        bs_j_add = np.zeros(bs_n, dtype=np.double)

        for ii in xrange(0, n_bsj): bs_j_add += theta[ip + ii] * bsj_mask[:, ii]
        ip += n_bsj

        for ii in xrange(0, n_bso): bs_model += theta[ip + ii] * bso_mask[:, ii]
        ip += n_bso

        for ii in xrange(0, n_bsl): bs_model += theta[ip + ii] * bs_x0 * bsl_mask[:, ii]
        ip += n_bsl

    ## kind_sin option 1 (SPLINE fit) has been removed from this code
    if kind_sin > 1:
        Prot_run = theta[ip]
        ip += 1

        ph_off = np.zeros(n_pof+1, dtype=np.double)

        if (phase_coherence):
            ph_off[0] = theta[ip + n_pof]
            ph_off[1:n_pof+1] = theta[ip:ip + n_pof]

        if ((not phase_synchro) and (not phase_coherence)):
            ph_off[1:n_pof+1] = theta[ip:ip + n_pof]

        ip += n_pof

        for nk in xrange(0, n_periods):

            ph_sin = np.zeros(n_pha, dtype=np.double)

            if (phase_synchro):
                ph_sin[:] = theta[ip]

            if ((not phase_synchro) and (not phase_coherence)):
                ph_sin[:] = theta[ip:ip + n_pha]


            ip += n_pha
            ip_pof = 0

            if rvp_flag[nk]:
                rv_xph = (rv_x0 / Prot_run) % 1
                for jj in xrange(0, n_rva):
                    for ii in xrange(0, n_amp[jj]):
                        rv_model += rvp_mask[:,nk] * rva_mask[:, jj] * theta[ip + ii] * np.sin(
                            ((ii + 1.) * rv_xph + ph_sin[ii] + ph_off[ip_pof])* 2. * np.pi)
                    ip += n_amp[jj]
                ip_pof += 1
            if p1p_flag[nk]:
                p1_xph = (p1_x0 / Prot_run) % 1
                for ii in range(0,n_pho):
                    p1_model += p1p_mask[:,nk] * theta[ip + ii] * np.sin(
                        ((ii + 1.) * p1_xph + ph_sin[ii] + ph_off[ip_pof]) * 2. * np.pi)
                ip += n_pho
                ip_pof += 1
            if p2p_flag[nk]:
                p2_xph = (p2_x0 / Prot_run) % 1
                tp_pof = ip_pof - 1
                for ii in range(0, n_pho):
                    p2_model += p2p_mask[:, nk] * theta[ip + ii] * np.sin(
                        ((ii + 1.) * p2_xph + ph_sin[ii] + ph_off[tp_pof]) * 2. * np.pi)
                ip += n_pho
                ## ip_pof += 1
            if fwp_flag[nk]:
                fw_xph = (fw_x0 / Prot_run) % 1
                for jj in xrange(0, n_fwa):
                    for ii in range(0, n_amp):
                        ## fw_model += fwp_mask[:, nk] * (
                        fw_model += fwp_mask[:, nk] * fwa_mask[:, jj] * theta[ip + ii] * np.sin(
                            ((ii + 1.) * fw_xph + ph_sin[ii] + ph_off[ip_pof]) * 2. * np.pi)
                    ip += n_amp
                ip_pof += 1
            if bsp_flag[nk]:
                bs_xph = (bs_x0 / Prot_run) % 1
                for jj in xrange(0, n_bsa):
                    for ii in range(0, n_amp[jj]):
                        ## bs_model += bsp_mask[:, nk] * (
                        bs_model += bsp_mask[:, nk] * bsa_mask[:, jj] * theta[ip + ii] * np.sin(
                            ((ii + 1.) * bs_xph + ph_sin[ii] + ph_off[ip_pof]) * 2. * np.pi)
                    ip += n_amp[jj]
                ip_pof += 1

    chi2_OUT = 0.
    if flag_rv:
        rve_env = 1.0 / (rv_e ** 2.0 + rv_j_add ** 2.0)
        chi2_OUT += np.sum((rv_y - rv_model) ** 2 * rve_env - np.log(rve_env))
    if flag_p1:
        p1e_env = 1.0 / (p1_e ** 2.0 + p1_j_add ** 2)
        chi2_OUT += np.sum((p1_y - p1_model) ** 2 * p1e_env - np.log(p1e_env))
    if flag_p2:
        p2e_env = 1.0 / (p2_e ** 2.0 + p2_j_add ** 2)
        chi2_OUT += np.sum((p2_y - p2_model) ** 2 * p2e_env - np.log(p2e_env))
    if flag_fw:
        fwe_env = 1.0 / (fw_e ** 2.0 + fw_j_add ** 2)
        chi2_OUT += np.sum((fw_y - fw_model) ** 2 * fwe_env - np.log(fwe_env))
    if flag_bs:
        bse_env = 1.0 / (bs_e ** 2.0 + bs_j_add ** 2)
        chi2_OUT += np.sum((bs_y - bs_model) ** 2 * bse_env - np.log(bse_env))
    ## print ii_tot,   chi2_bin, chi2_bin/ii_tot, theta[ip]
    return -0.5 * (chi2_OUT), derived_orb


def lnprior(theta):
    ## The parameters are stored as a vector of values, so unpack them
    ## We're using only uniform priors
    ip = 0
    period_stack = []
    for pp in xrange(0, n_planets):

        if not bounds[ip + 0, 0] < theta[ip + 0] < bounds[ip + 0, 1]:
            return -np.inf
        ## this is to avoid period crossing, which results in anomalous
        ## semiamplitdes compensated by opposite phases
        period_stack.append(theta[ip + 0])
        for pk in xrange(0, pp):
            if np.abs(period_stack[pp] - period_stack[pk]) < 0.06:
                return -np.inf

        if not bounds[ip + 1, 0] < theta[ip + 1] < bounds[ip + 1, 1]:
            return -np.inf

        if not bounds[ip + 2, 0] <= theta[ip + 2] < bounds[ip + 2, 1]:
            return -np.inf

        if (n_orbpams[pp] > 3):
            if not bounds[ip + 3, 0] < theta[ip + 3] < bounds[ip + 3, 1]:
                return -np.inf

            if not bounds[ip + 4, 0] < theta[ip + 4] < bounds[ip + 4, 1]:
                return -np.inf

            e = theta[ip + 3] ** 2 + theta[ip + 4] ** 2
            if not e_bounds[pp][0] <= e < e_bounds[pp][1]:
                return -np.inf

            omega = np.arctan2(theta[ip + 3], theta[ip + 4])
            if not -np.pi < omega < np.pi:
                return -np.inf

        ip += n_orbpams[pp]

    for ie in xrange(ip, ip + n_trends + n_syspams + n_extpams):
        if not bounds[ie, 0] < theta[ie] < bounds[ie, 1]:
            return -np.inf

    ## check if the periods of the planets cross the
    ## rotational period
    ip += n_trends + n_syspams
    if (kind_sin > 0):
        logP_rot = np.log10(theta[ip])
        for pp in xrange(0, n_planets):
            if np.abs(period_stack[pp] - logP_rot) < 0.06:
                return -np.inf

    return 0.0


def lnprob(theta):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf, [45]
    lk, der_orb = lnlike(theta)
    if not np.isfinite(lk):
        return -np.inf, der_orb
    return lp + lk, der_orb


def lnprob_PyDE(theta):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    lk, der_orb = lnlike(theta)
    if not np.isfinite(lk):
        return -np.inf
    return lp + lk


n_orbpams = []
## User-specific bound for each planet
logP_bounds = []
logK_bounds = []
e_bounds = []  ## this is used to set a constraint on ecccentricity


print 'Waiting for .config file name...'

file_conf = raw_input()
## file_conf="OC102_2pE_spl.config"
file = open(file_conf, 'r')
row = file.readlines()

kind_sin = 0
n_planets = 0
n_trends = 0
n_sin = 0

rv_n = 0
p1_n = 0
p2_n = 0
fw_n = 0
bs_n = 0

n_rvo = 0
n_rvj = 0
n_rvl = 0
n_rva = 0

n_p1j = 0
n_p1o = 0
n_p1l = 0
n_p2j = 0
n_p2o = 0
n_p2l = 0

n_fwj = 0
n_fwo = 0
n_fwl = 0
n_fwa = 0

n_bsj = 0
n_bso = 0
n_bsl = 0
n_bsa = 0

n_pin = 0
n_amp = 0
n_pha = 0
n_pof = 0

spline_ord = 3

flag_p1 = False
flag_p2 = False
flag_rv = False
flag_fw = False
flag_bs = False
pl_row = 1

period_range = []
pr_row = 0
period_sel   = False
sa_row = 1

thin = 1

n_spec = 1
phase_coherence = False
phase_synchro   = False
recenter_bounds = True

pof_bounds = np.asarray([0., 1.0], dtype=np.double)
pha_bounds = np.asarray([0., 1.0], dtype=np.double)
rva_bounds = np.asarray([0., 100.], dtype=np.double)
p1a_bounds = np.asarray([0., 0.5], dtype=np.double)
p2a_bounds = np.asarray([0., 0.5], dtype=np.double)
fwa_bounds = np.asarray([0., 2000.], dtype=np.double)
bsa_bounds = np.asarray([0., 60.], dtype=np.double)

for line in row:
    info = line.split()
    if line.find("Input_RV") > -1:
        input_rv = info[1]
        flag_rv = True
    if line.find("Input_P1") > -1:
        input_p1 = info[1]
        flag_p1 = True
    if line.find("Input_P2") > -1:
        input_p2 = info[1]
        flag_p2 = True
    if line.find("Input_BS") > -1:
        input_bs = info[1]
        flag_bs = True
    if line.find("Input_FW") > -1:
        input_fw = info[1]
        flag_fw = True
    if line.find("Output") > -1:
        planet_name = info[1]
    if line.find("Nplanets") > -1:
        n_planets = np.asarray(info[1], dtype=np.int64)
    if line.find("Planet" + `pl_row` + " ") > -1:
        n_orbpams.append(np.asarray(info[1], dtype=np.int64))

        Pmin = np.asarray(info[2], dtype=np.double)
        Pmax = np.asarray(info[3], dtype=np.double)
        logP_bounds.append([np.log10(Pmin), np.log10(Pmax)])

        Kmin = np.asarray(info[4], dtype=np.double)
        Kmax = np.asarray(info[5], dtype=np.double)
        logK_bounds.append([np.log10(Kmin), np.log10(Kmax)])

        if len(info) > 6:
            e_bounds.append([np.asarray(info[6], dtype=np.double), np.asarray(info[7], dtype=np.double)])
        else:
            e_bounds.append([0., 1.])
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
        kind_sin = np.asarray(info[1], dtype=np.int64)
        Prot_bounds = [Pmin, Pmax]
        n_periods = 1
    ## The number of chuncks for the activity fit
    if line.find("Nperiods") > -1:
        n_periods = np.asarray(info[1], dtype=np.int64)
        period_sel = True
        # Each chunk has a different rotational period, still within the boundaries
        # defined in  Prot, and it is defined as an interval in BJD

    if line.find("Period" + `pr_row` + " ") > -1:
        P0_tmp = np.asarray(info[1], dtype=np.double)
        P1_tmp = np.asarray(info[2], dtype=np.double)

        period_range.append(np.asarray([P0_tmp, P1_tmp], dtype=np.double))
        # print period_range[pr_row,1]
        # print period_range[pr_row ][0]
        # print period_range[pr_row ][1]

        pr_row += 1

    if line.find("Ngen") > -1:
        ngen = np.asarray(info[1], dtype=np.int64)
    if line.find("Nsteps") > -1:
        nsteps = np.asarray(info[1], dtype=np.int64)
    if line.find("Nburn") > -1:
        nburn = np.asarray(info[1], dtype=np.int64)
    if line.find("Npop_mult") > -1:
        npop_mult = np.asarray(info[1], dtype=np.int64)
    if line.find("Thin") > -1:
        thin = np.asarray(info[1], dtype=np.int64)

    if line.find("Spec_dataset") > -1:
        # this is the number of spectral dataset, followd by the number of
        # sinusoids that must be used for each dataset for the activity fit
        # So, if we are using two dataset where the first one has greater precision
        # we can specifiy to use 3 sinusoids for the first dataset and only two
        # for the second one
        # Additional activity indicators associated to RVs must the same order and number of sinusoids,
        #
        n_spec = np.asarray(info[1], dtype=np.int64)
        n_amp = np.asarray(info[2:n_spec + 2], dtype=np.int64)

    if line.find("Phase_coherence") > -1:
        # choose if the sinusoids at Prot harmonics must have all the same phase
        # (Phase_coherence = 1) or not (default)
        if np.asarray(info[1], dtype=np.int64) >= 1:
            phase_coherence = True

    if line.find("Phase_synchro") > -1:
        # choose if Phot, Rv and other indicators must share the same phases
        # (Phase_synchro = 1) or not (default)
        if np.asarray(info[1], dtype=np.int64) >= 1:
            phase_synchro = True

    if line.find("Recenter_Bounds") > -1:
        # required to avoid a small bug in the code
        # if the dispersion of PyDE walkers around the median value is too broad,
        # then emcee walkers will start outside the bounds, causing an error
        if np.asarray(info[1], dtype=np.int64) < 1:
            recenter_bounds = False

if (period_sel == False):
    n_periods = 1
    period_range.append([np.asarray(0.0, dtype=np.double), np.asarray(5000000.0, dtype=np.double)])
if (n_periods != pr_row):
    print 'PROBLEM'


print 'KIND: ', kind_sin
if (kind_sin == 0):
    flag_p1 = False
    flag_p2 = False
    flag_fw = False
    flag_bs = False
    n_periods = 0
    ## It doesn't make sense to include photometric data if
    ## we don't use them to model the RV jitter

if (not flag_rv) and n_planets > 0:
    print 'No RVs but n_planets>0'
    print 'The code cannot handle this'

rvp_flag = np.zeros(n_periods,dtype=bool)
p1p_flag = np.zeros(n_periods,dtype=bool)
p2p_flag = np.zeros(n_periods,dtype=bool)
fwp_flag = np.zeros(n_periods,dtype=bool)
bsp_flag = np.zeros(n_periods,dtype=bool)


## To disable one of the additional offset or rescaling factor
## They must have -1 as value (Python conv. starting from 0)
if flag_rv:
    print 'Opening: ', input_rv
    data_rv = np.loadtxt(input_rv)
    rv_x = np.asarray(data_rv[:, 0], dtype=np.double)
    rv_y = np.asarray(data_rv[:, 1], dtype=np.double)
    rv_e = np.asarray(data_rv[:, 2], dtype=np.double)
    rv_j = np.asarray(data_rv[:, 3], dtype=np.double)
    ## fit for different RV jitters
    rv_o = np.asarray(data_rv[:, 4], dtype=np.double)
    rv_l = np.asarray(data_rv[:, 5], dtype=np.double)
    ## use different offsets for the data
    ## off must start from zero
    ## -1 values for rv_j and rv_l mean that these

    rv_a = np.asarray(data_rv[:, 6], dtype=np.double)
    # Flag for activity fit: we can choose to fit the same RV amplitude
    # of the activity signal for all the datasets, use dfferent values
    #for each dataset or exclude some of the datasets

    ## Model for RV systematics
    n_rvo = np.max(rv_o.astype(np.int64)) + 1
    n_rvj = np.max(rv_j.astype(np.int64)) + 1
    n_rvl = np.max(rv_l.astype(np.int64)) + 1
    n_rva = np.max(rv_a.astype(np.int64)) + 1

    rv_n = np.size(rv_x)
    print 'RV_N = ', rv_n
    print 'RV_N jitter = ', n_rvj
    print 'RV_N offset = ', n_rvo
    print 'RV_N linear = ', n_rvl
    print 'RV_N activ. = ', n_rva
    print
    ## logK_bounds = np.asarray([np.log10(0.1),np.log10(np.max(rv_y)-np.min(rv_y))+1.],dtype=np.double)

## Reading photometric data in the first filter, if present:
if flag_p1:
    print 'Opening: ', input_p1
    data_p1 = np.loadtxt(input_p1)
    p1_x = np.asarray(data_p1[:, 0], dtype=np.double)
    p1_y = np.asarray(data_p1[:, 1], dtype=np.double)
    p1_e = np.asarray(data_p1[:, 2], dtype=np.double)
    p1_j = np.asarray(data_p1[:, 3], dtype=np.double)
    p1_o = np.asarray(data_p1[:, 4], dtype=np.double)
    p1_l = np.asarray(data_p1[:, 5], dtype=np.double)

    ## Model for PH systematics
    n_p1o = np.max(p1_o.astype(np.int64)) + 1
    n_p1j = np.max(p1_j.astype(np.int64)) + 1
    n_p1l = np.max(p1_l.astype(np.int64)) + 1

    p1_n = np.size(p1_x)
    print 'P1_N = ', p1_n
    print 'P1_N jitter = ', n_p1j
    print 'P1_N offset = ', n_p1o
    print 'P1_N linear = ', n_p1l
    print

# photometric data in a different filter
if flag_p2:
    print 'Opening: ', input_p2
    data_p2 = np.loadtxt(input_p2)
    p2_x = np.asarray(data_p2[:, 0], dtype=np.double)
    p2_y = np.asarray(data_p2[:, 1], dtype=np.double)
    p2_e = np.asarray(data_p2[:, 2], dtype=np.double)
    p2_j = np.asarray(data_p2[:, 3], dtype=np.double)
    p2_o = np.asarray(data_p2[:, 4], dtype=np.double)
    p2_l = np.asarray(data_p2[:, 5], dtype=np.double)

    ## Model for PH systematics
    n_p2o = np.max(p2_o.astype(np.int64)) + 1
    n_p2j = np.max(p2_j.astype(np.int64)) + 1
    n_p2l = np.max(p2_l.astype(np.int64)) + 1

    p2_n = np.size(p2_x)
    print 'P2_N = ', p2_n
    print 'P2_N jitter = ', n_p2j
    print 'P2_N offset = ', n_p2o
    print 'P2_N linear = ', n_p2l
    print


## Reading FWHM data, if present:
if flag_fw:
    print 'Opening: ', input_fw
    data_fw = np.loadtxt(input_fw)
    fw_x = np.asarray(data_fw[:, 0], dtype=np.double)
    fw_y = np.asarray(data_fw[:, 1], dtype=np.double)
    fw_e = np.asarray(data_fw[:, 2], dtype=np.double)
    fw_j = np.asarray(data_fw[:, 3], dtype=np.double)
    fw_o = np.asarray(data_fw[:, 4], dtype=np.double)
    fw_l = np.asarray(data_fw[:, 5], dtype=np.int64)
    ## A linear trend in the HARPS-N data has been detected in the period
    ## --- to --- BJD
    ## since the program doesn't know in principle to which instrument the
    ## data are coming from, the use of the FWHM linear term ust be specified in the
    ## data itself as an additional keyword (fw_l)

    fw_a = np.asarray(data_fw[:, 6], dtype=np.int64)

    # Model for FWHM systematics
    n_fwo = np.max(fw_o.astype(np.int64)) + 1
    n_fwj = np.max(fw_j.astype(np.int64)) + 1
    n_fwl = np.max(fw_l.astype(np.int64)) + 1
    n_fwa = np.max(fw_a.astype(np.int64)) + 1
    ## if np.sum(fw_l) > 0: n_fwl = 1

    fw_n = np.size(fw_x)
    print 'FW_N = ', fw_n
    print 'FW_N jitter = ', n_fwj
    print 'FW_N offset = ', n_fwo
    print 'FW_N activ. = ', n_fwa
    # print 'FW_N linear = ', n_fwl
    print

## Reading BIS data, if present:
if flag_bs:
    print 'Opening: ', input_bs
    data_bs = np.loadtxt(input_bs)
    bs_x = np.asarray(data_bs[:, 0], dtype=np.double)
    bs_y = np.asarray(data_bs[:, 1], dtype=np.double)
    bs_e = np.asarray(data_bs[:, 2], dtype=np.double)
    bs_j = np.asarray(data_bs[:, 3], dtype=np.double)
    bs_o = np.asarray(data_bs[:, 4], dtype=np.double)
    bs_l = np.asarray(data_bs[:, 5], dtype=np.double)
    bs_a = np.asarray(data_bs[:, 6], dtype=np.double)

    ## Model for BIS systematics
    n_bso = np.max(bs_o.astype(np.int64)) + 1
    n_bsj = np.max(bs_j.astype(np.int64)) + 1
    n_bsl = np.max(bs_l.astype(np.int64)) + 1
    n_bsa = np.max(bs_a.astype(np.int64)) + 1

    bs_n = np.size(bs_x)
    print 'BS_N = ', bs_n
    print 'BS_N jitter = ', n_bsj
    print 'BS_N offset = ', n_bso
    print 'BS_N linear = ', n_bsl
    print 'BS_N activ. = ', n_bsa
    print



# BOUND definition
if flag_rv:

    rvj_bounds = np.asarray([0.,10*np.max(rv_e)],dtype=np.double)
    rvj_mask = np.zeros([rv_n, n_rvj], dtype=np.double)
    for ii in xrange(0, n_rvj):
        rvj_mask[(abs(rv_j - ii) < 0.1), ii] = 1

    rvo_bounds = np.asarray([np.min(rv_y), np.max(rv_y)], dtype=np.double)
    rvo_mask = np.zeros([rv_n, n_rvo], dtype=np.double)
    for ii in xrange(0, n_rvo):
        rvo_mask[(abs(rv_o - ii) < 0.1), ii] = 1

    # extremely conservative values for the coefficient of the drift

    rvl_bounds = np.asarray([-1., 1.], dtype=np.double)
    rvl_mask = np.zeros([rv_n, n_rvl], dtype=np.double)
    for ii in xrange(0, n_rvl):
        rvl_mask[(abs(rv_l - ii) < 0.1), ii] = 1

    rva_mask = np.zeros([rv_n, n_rva], dtype=np.double)
    for ii in xrange(0, n_rva):
        rva_mask[(abs(rv_a - ii) < 0.1), ii] = 1

    ## The rotational period is a general property of the system, not
    ## of each dataset, so the boundaries are specified at a higher level
    ## However a mask must be created for each parameter, since they may have
    ## a different set of observations
    rvp_mask = np.zeros([rv_n, n_periods], dtype=np.double)
    for ii in xrange(0, n_periods):
        rvp_sel = (period_range[ii][0] < rv_x) & (rv_x < period_range[ii][1])
        rvp_mask[rvp_sel, ii] = 1
        if (np.sum(rvp_sel)>0): rvp_flag[ii]= True

# While the RV are "just physics", the amplitude and zero point
# of photometry, bis span and FWHM may change for systematics effects
# of the instrument or for just a different way to measure them

if flag_p1:
    p1j_bounds = np.asarray([0., 50 * np.max(p1_e)], dtype=np.double)
    p1j_mask = np.zeros([p1_n, n_p1j], dtype=np.double)
    for ii in xrange(0, n_p1j):
        p1j_mask[(abs(p1_j - ii) < 0.1), ii] = 1

    p1o_bounds = np.asarray([np.min(p1_y), np.max(p1_y)], dtype=np.double)
    p1o_mask = np.zeros([p1_n, n_p1o], dtype=np.double)
    for ii in xrange(0, n_p1o):
        p1o_mask[(abs(p1_o - ii) < 0.1), ii] = 1

    p1l_bounds = np.asarray([-1., 1.], dtype=np.double)
    p1l_mask = np.zeros([p1_n, n_p1l], dtype=np.double)
    for ii in xrange(0, n_p1l):
        p1l_mask[(abs(p1_l - ii) < 0.1), ii] = 1

    p1p_mask = np.zeros([p1_n, n_periods], dtype=np.double)
    for ii in xrange(0, n_periods):
        p1p_sel = (period_range[ii][0] < p1_x) & (p1_x < period_range[ii][1])
        p1p_mask[p1p_sel, ii] = 1
        if (np.sum(p1p_sel)>0): p1p_flag[ii]= True

if flag_p2:
    p2j_bounds = np.asarray([0., 50 * np.max(p2_e)], dtype=np.double)
    p2j_mask = np.zeros([p2_n, n_p2j], dtype=np.double)
    for ii in xrange(0, n_p2j):
        p2j_mask[(abs(p2_j - ii) < 0.1), ii] = 1

    p2o_bounds = np.asarray([np.min(p2_y), np.max(p2_y)], dtype=np.double)
    p2o_mask = np.zeros([p2_n, n_p2o], dtype=np.double)
    for ii in xrange(0, n_p2o):
        p2o_mask[(abs(p2_o - ii) < 0.1), ii] = 1

    p2l_bounds = np.asarray([-1., 1.], dtype=np.double)
    p2l_mask = np.zeros([p2_n, n_p2l], dtype=np.double)
    for ii in xrange(0, n_p2l):
        p2l_mask[(abs(p2_l - ii) < 0.1), ii] = 1

    p2p_mask = np.zeros([p2_n, n_periods], dtype=np.double)
    for ii in xrange(0, n_periods):
        p2p_sel = (period_range[ii][0] < p2_x) & (p2_x < period_range[ii][1])
        p2p_mask[p2p_sel, ii] = 1
        if (np.sum(p2p_sel)>0): p2p_flag[ii]= True

if flag_fw:
    fwj_bounds = np.asarray([0., 50 * np.max(fw_e)], dtype=np.double)
    fwj_mask = np.zeros([fw_n, n_fwj], dtype=np.double)
    for ii in xrange(0, n_fwj):
        fwj_mask[(abs(fw_j - ii) < 0.1), ii] = 1

    fwo_bounds = np.asarray([np.min(fw_y), np.max(fw_y)], dtype=np.double)
    fwo_mask = np.zeros([fw_n, n_fwo], dtype=np.double)
    for ii in xrange(0, n_fwo):
        fwo_mask[(abs(fw_o - ii) < 0.1), ii] = 1

    fwl_bounds = np.asarray([-1., 1.], dtype=np.double)
    fwl_mask = np.zeros([fw_n, n_fwl], dtype=np.double)
    for ii in xrange(0, n_fwl):
        fwl_mask[(abs(fw_l - ii) < 0.1), ii] = 1

    fwa_bounds = np.asarray([0., 5.*(np.max(fw_y) - np.min(fw_y))], dtype=np.double)
    fwa_mask = np.zeros([fw_n, n_fwa], dtype=np.double)
    for ii in xrange(0, n_fwa):
        fwa_mask[(abs(fw_a - ii) < 0.1), ii] = 1  #fwp_mask = np.zeros([fw_n, n_periods], dtype=np.double)

    fwp_mask = np.zeros([fw_n, n_periods], dtype=np.double)
    for ii in xrange(0, n_periods):
        fwp_sel = (period_range[ii][0] < fw_x) & (fw_x < period_range[ii][1])
        fwp_mask[fwp_sel, ii] = 1
        if (np.sum(fwp_sel)>0): fwp_flag[ii]= True
        # print 'FW N jitter = ',n_fwj
        # print 'FW N offset = ',n_fwo

if flag_bs:
    bsj_bounds = np.asarray([0., 50 * np.max(bs_e)], dtype=np.double)
    bsj_mask = np.zeros([bs_n, n_bsj], dtype=np.double)
    for ii in xrange(0, n_bsj):
        bsj_mask[(abs(bs_j - ii) < 0.1), ii] = 1

    bso_bounds = np.asarray([np.min(bs_y), np.max(bs_y)], dtype=np.double)
    bso_mask = np.zeros([bs_n, n_bso], dtype=np.double)
    for ii in xrange(0, n_bso):
        bso_mask[(abs(bs_o - ii) < 0.1), ii] = 1

    bsl_bounds = np.asarray([-1., 1.], dtype=np.double)
    bsl_mask = np.zeros([bs_n, n_bsl], dtype=np.double)
    for ii in xrange(0, n_bsl):
        bsl_mask[(abs(bs_l - ii) < 0.1), ii] = 1

    bsa_bounds = np.asarray([0., 5.*(np.max(bs_y) - np.min(bs_y))], dtype=np.double)
    bsa_mask = np.zeros([bs_n, n_bsa], dtype=np.double)
    for ii in xrange(0, n_bsa):
        bsa_mask[(abs(bs_a - ii) < 0.1), ii] = 1

    bsp_mask = np.zeros([bs_n, n_periods], dtype=np.double)
    for ii in xrange(0, n_periods):
        bsp_sel = (period_range[ii][0] < bs_x) & (bs_x < period_range[ii][1])
        bsp_mask[bsp_sel, ii] = 1
        if (np.sum(bsp_sel)>0): bsp_flag[ii]= True

# Arbitrary T0 point
# This is a simple check to be sure that the Tref is
# determined from one dataset at least
if flag_bs: Tref = np.mean(bs_x, dtype=np.double)
if flag_fw: Tref = np.mean(fw_x, dtype=np.double)
if flag_p2: Tref = np.mean(p2_x, dtype=np.double)
if flag_p1: Tref = np.mean(p1_x, dtype=np.double)
if flag_rv: Tref = np.mean(rv_x, dtype=np.double)

if flag_rv: rv_x0 = rv_x - Tref
if flag_p1: p1_x0 = p1_x - Tref
if flag_p2: p2_x0 = p2_x - Tref
if flag_fw: fw_x0 = fw_x - Tref
if flag_bs: bs_x0 = bs_x - Tref

n_dataset = 0
n_amplitude = 0
if kind_sin > 0:
    n_pin = 1
    # RV and photometry share the same period (the rotational one of the star)
    # awlays otherwise all these calculation are meaningless
    if flag_rv: n_dataset += 1
    if flag_p1: n_dataset += 1
    if flag_p2: n_dataset += 1
    if flag_fw: n_dataset += 1
    if flag_bs: n_dataset += 1

    n_pho = np.amax(n_amp)

    if flag_rv: n_amplitude += np.sum(n_amp) * np.sum(rvp_flag)
    if flag_p1: n_amplitude += n_pho * np.sum(p1p_flag)
    if flag_p2: n_amplitude += n_pho * np.sum(p2p_flag)
    if flag_fw: n_amplitude += np.sum(n_amp[0:n_fwa]) * np.sum(fwp_flag)
    if flag_bs: n_amplitude += np.sum(n_amp[0:n_bsa]) * np.sum(bsp_flag)


    if (phase_coherence):
        n_pha = 1
    else:
        n_pha = np.amax(n_amp)

    if (phase_synchro):
        n_pof = 0
    else:
        n_pof = 1 * (n_dataset - 1)

    #the two photometric curves should have the same phase
    if flag_p2: n_pof -= 1

n_syspams = n_rvo + n_rvj + n_rvl + \
            n_p1j + n_p1o + n_p1l + \
            n_p2j + n_p2o + n_p2l + \
            n_fwj + n_fwo + n_fwl + \
            n_bsj + n_bso + n_bsl
n_extpams = n_pin + n_pof + n_amplitude + n_pha*n_periods
ndim = np.sum(n_orbpams, dtype=np.int64) + n_trends + n_syspams + n_extpams

# Parameters por PyDE and emcee
npop = ndim * npop_mult
# ngen   =  1000
# nsteps = 30000
# nburn  =  5000
if (npop%2==1): npop+=1
print 'Dimensions = ', ndim
print '   N orb = ', np.sum(n_orbpams, dtype=np.int64)
print '   N sys = ', n_syspams
print '   N ext = ', n_extpams
print 'Npop = ', npop
print npop


bounds = np.zeros([ndim, 2], dtype=np.double)
ip = 0
if flag_rv:
    for pp in xrange(0, n_planets):
        # fix here if n_orbpams is changed
        bounds[ip + 0, :] = logP_bounds[pp][:]
        bounds[ip + 1, :] = logK_bounds[pp][:]
        #bounds[ip + 1, :] = logK_bounds[:]
        bounds[ip + 2, :] = [0, 2*np.pi]
        if (n_orbpams[pp] == 5):
            bounds[ip + 3, :] = [-1.0, 1.0]
            bounds[ip + 4, :] = [-1.0, 1.0]
        ip += n_orbpams[pp]

    for ii in xrange(0, n_trends): bounds[ip + ii, :] = [-1., 1.0]
    ip += n_trends

    for ii in xrange(0, n_rvj): bounds[ip + ii, :] = rvj_bounds[:]
    ip += n_rvj

    for ii in xrange(0, n_rvo): bounds[ip + ii, :] = rvo_bounds[:]
    ip += n_rvo

    for ii in xrange(0, n_rvl): bounds[ip + ii, :] = rvl_bounds[:]
    ip += n_rvl

if flag_p1:
    for ii in xrange(0, n_p1j): bounds[ip + ii, :] = p1j_bounds[:]
    ip += n_p1j

    for ii in xrange(0, n_p1o): bounds[ip + ii, :] = p1o_bounds[:]
    ip += n_p1o

    for ii in xrange(0, n_p1l): bounds[ip + ii, :] = p1l_bounds[:]
    ip += n_p1l

if flag_p2:
    for ii in xrange(0, n_p2j): bounds[ip + ii, :] = p2j_bounds[:]
    ip += n_p2j

    for ii in xrange(0, n_p2o): bounds[ip + ii, :] = p2o_bounds[:]
    ip += n_p2o

    for ii in xrange(0, n_p2l): bounds[ip + ii, :] = p2l_bounds[:]
    ip += n_p2l

if flag_fw:
    for ii in xrange(0, n_fwj): bounds[ip + ii, :] = fwj_bounds[:]
    ip += n_fwj

    for ii in xrange(0, n_fwo): bounds[ip + ii, :] = fwo_bounds[:]
    ip += n_fwo

    for ii in xrange(0, n_fwl): bounds[ip + ii, :] = fwl_bounds[:]
    ip += n_fwl

if flag_bs:
    for ii in xrange(0, n_bsj): bounds[ip + ii, :] = bsj_bounds[:]
    ip += n_bsj

    for ii in xrange(0, n_bso): bounds[ip + ii, :] = bso_bounds[:]
    ip += n_bso

    for ii in xrange(0, n_bsl): bounds[ip + ii, :] = bsl_bounds[:]
    ip += n_bsl

if (kind_sin > 0):
    bounds[ip, :] = Prot_bounds[:]
    ip += 1
    for ii in range(0, n_pof):
        bounds[ip + ii, :] = pof_bounds[:]
    ip += n_pof

    for nk in xrange(0, n_periods):

        for ii in range(0, n_pha):
            bounds[ip + ii, :] = pha_bounds[:]
        ip += n_pha

        if rvp_flag[nk]:
            for ii in xrange(0, np.sum(n_amp)):
                bounds[ip + ii, :] = rva_bounds[:]
            ip += np.sum(n_amp)

        if p1p_flag[nk]:
            for ii in xrange(0, n_pho):
                bounds[ip + ii, :] = p1a_bounds[:]
            ip += n_pho

        if p2p_flag[nk]:
            for ii in xrange(0, n_pho):
                bounds[ip + ii, :] = p2a_bounds[:]
            ip += n_pho

        if fwp_flag[nk]:
            for ii in xrange(0, np.sum(n_amp[0:n_fwa])):
                bounds[ip + ii, :] = fwa_bounds[:]
            ip += np.sum(n_amp[0:n_fwa])

        if bsp_flag[nk]:
            for ii in xrange(0, np.sum(n_amp[0:n_bsa])):
                bounds[ip + ii, :] = bsa_bounds[:]
            ip +=  np.sum(n_amp[0:n_bsa])

print 'BOUNDS'
print bounds

print 'PyDE'
de = DiffEvol(lnprob_PyDE, bounds, npop, maximize=True)
de.optimize(ngen)
print 'PyDE completed'

np.savetxt('output/' + planet_name + '_pyDEout_output_bounds.dat',bounds)
np.savetxt('output/' + planet_name + '_pyDEout_output_pops.dat',de.population)


pyde_mean = np.mean(de.population, axis=0)
print pyde_mean
np.savetxt('output/' + planet_name + '_pyDEout_mean.dat',pyde_mean)

# fix for PyDE anomalous results
for ii in xrange(0,ndim):
    if np.amax(de.population[:,ii])-np.amin(de.population[:,ii]) < 10e-7 :
        range_restricted = (bounds[ii,1]-bounds[ii,0])/1000.
        min_bound = np.maximum((pyde_mean[ii]-range_restricted/2.0),bounds[ii,0])
        max_bound = np.minimum((pyde_mean[ii]+range_restricted/2.0),bounds[ii,1])
        de.population[:,ii] =  np.random.uniform(min_bound,max_bound,npop)

# Centering the phase variables around the derived value

if recenter_bounds:
    ip = 0

    for pp in xrange(0, n_planets):
        ii = ip + 2
        bounds[ii,0] = pyde_mean[ii]-np.pi
        bounds[ii,1] = pyde_mean[ii]+np.pi
        fix_sel = (de.population[:,ii]<=bounds[ii,0]) | (de.population[:,ii]>=bounds[ii,1])
        de.population[fix_sel,ii]=pyde_mean[ii]
        #for jk in xrange(0,npop):
        #    if de.population[jk,ii]<bounds[ii,0] or de.population[jk,ii]>bounds[ii,1]: de.population[jk,ii]=pyde_mean[ii]
        ip += n_orbpams[pp]

    ip += n_syspams
    if (kind_sin > 0):
        ip += 1
        for ii in range(ip, ip+n_pof):
            bounds[ii,0] = pyde_mean[ii]-0.5
            bounds[ii,1] = pyde_mean[ii]+0.5
            fix_sel = (de.population[:,ii]<=bounds[ii,0]) | (de.population[:,ii]>=bounds[ii,1])
            de.population[fix_sel,ii]=pyde_mean[ii]
            #for jk in xrange(0,npop):
            #    if de.population[jk,ii]<bounds[ii,0] or de.population[jk,ii]>bounds[ii,1]: de.population[jk,ii]=pyde_mean[ii]

        ip += n_pof
        for nk in xrange(0, n_periods):
            for ii in range(ip, ip+n_pha):
                bounds[ii,0] = pyde_mean[ii]-0.5
                bounds[ii,1] = pyde_mean[ii]+0.5
                fix_sel = (de.population[:,ii]<=bounds[ii,0]) | (de.population[:,ii]>=bounds[ii,1])
                de.population[fix_sel,ii]=pyde_mean[ii]
                #for jk in xrange(0,npop):
                #    if de.population[jk,ii]<bounds[ii,0] or de.population[jk,ii]>bounds[ii,1]: de.population[jk,ii]=pyde_mean[ii]
            ip += n_pha
            if rvp_flag[nk]: ip += np.sum(n_amp)
            if p1p_flag[nk]: ip += n_pho
            if p2p_flag[nk]: ip += n_pho
            if fwp_flag[nk]: ip += np.sum(n_amp[0:n_fwa])
            if bsp_flag[nk]: ip += np.sum(n_amp[0:n_bsa])

    print 'REDEFINED BOUNDS'
    print bounds

    np.savetxt('output/' + planet_name + '_pyDEout_redefined_bounds.dat',bounds)
    np.savetxt('output/' + planet_name + '_pyDEout_redefined_pops.dat',de.population)

print 'emcee'
nwalkers = npop
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=24)
sampler.run_mcmc(de.population, nsteps, thin=thin)

print 'emcee completed'

h5f = h5py.File('output/' + planet_name + '.hdf5', "w")

data_grp = h5f.create_group("data")
data_grp.attrs.create("flag_rv", data=flag_rv)
data_grp.attrs.create("flag_p1", data=flag_p1)
data_grp.attrs.create("flag_p2", data=flag_p2)
data_grp.attrs.create("flag_fw", data=flag_fw)
data_grp.attrs.create("flag_bs", data=flag_bs)
data_grp.attrs.create("Tref", data=Tref)

if flag_rv: data_grp.create_dataset("data_rv", data=data_rv, compression="gzip")
if flag_p1: data_grp.create_dataset("data_p1", data=data_p1, compression="gzip")
if flag_p2: data_grp.create_dataset("data_p2", data=data_p2, compression="gzip")
if flag_fw: data_grp.create_dataset("data_fw", data=data_fw, compression="gzip")
if flag_bs: data_grp.create_dataset("data_bs", data=data_bs, compression="gzip")

if flag_rv:
    data_grp.create_dataset("rvj_mask", data=rvj_mask, compression="gzip")
    data_grp.create_dataset("rvo_mask", data=rvo_mask, compression="gzip")
    data_grp.create_dataset("rvl_mask", data=rvl_mask, compression="gzip")
    data_grp.create_dataset("rva_mask", data=rva_mask, compression="gzip")
    data_grp.create_dataset("rvp_mask", data=rvp_mask, compression="gzip")
if flag_p1:
    data_grp.create_dataset("p1j_mask", data=p1j_mask, compression="gzip")
    data_grp.create_dataset("p1o_mask", data=p1o_mask, compression="gzip")
    data_grp.create_dataset("p1l_mask", data=p1l_mask, compression="gzip")
    data_grp.create_dataset("p1p_mask", data=p1p_mask, compression="gzip")
if flag_p2:
    data_grp.create_dataset("p2j_mask", data=p2j_mask, compression="gzip")
    data_grp.create_dataset("p2o_mask", data=p2o_mask, compression="gzip")
    data_grp.create_dataset("p2l_mask", data=p2l_mask, compression="gzip")
    data_grp.create_dataset("p2p_mask", data=p2p_mask, compression="gzip")
if flag_fw:
    data_grp.create_dataset("fwj_mask", data=fwj_mask, compression="gzip")
    data_grp.create_dataset("fwo_mask", data=fwo_mask, compression="gzip")
    data_grp.create_dataset("fwl_mask", data=fwl_mask, compression="gzip")
    data_grp.create_dataset("fwa_mask", data=fwa_mask, compression="gzip")
    data_grp.create_dataset("fwp_mask", data=fwp_mask, compression="gzip")
if flag_bs:
    data_grp.create_dataset("bsj_mask", data=bsj_mask, compression="gzip")
    data_grp.create_dataset("bso_mask", data=bso_mask, compression="gzip")
    data_grp.create_dataset("bsl_mask", data=bsl_mask, compression="gzip")
    data_grp.create_dataset("bsa_mask", data=bsa_mask, compression="gzip")
    data_grp.create_dataset("bsp_mask", data=bsp_mask, compression="gzip")

data_grp.create_dataset("rvp_flag", data=rvp_flag, compression="gzip")
data_grp.create_dataset("p1p_flag", data=p1p_flag, compression="gzip")
data_grp.create_dataset("p2p_flag", data=p2p_flag, compression="gzip")
data_grp.create_dataset("fwp_flag", data=fwp_flag, compression="gzip")
data_grp.create_dataset("bsp_flag", data=bsp_flag, compression="gzip")

data_grp.attrs.create("rv_n", data=rv_n)
data_grp.attrs.create("p1_n", data=p1_n)
data_grp.attrs.create("p2_n", data=p2_n)
data_grp.attrs.create("fw_n", data=fw_n)
data_grp.attrs.create("bs_n", data=bs_n)

data_grp.attrs.create("n_rvj", data=n_rvj)
data_grp.attrs.create("n_rvo", data=n_rvo)
data_grp.attrs.create("n_rvl", data=n_rvl)
data_grp.attrs.create("n_rva", data=n_rva)
data_grp.attrs.create("n_p1j", data=n_p1j)
data_grp.attrs.create("n_p1o", data=n_p1o)
data_grp.attrs.create("n_p1l", data=n_p1l)
data_grp.attrs.create("n_p2j", data=n_p2j)
data_grp.attrs.create("n_p2o", data=n_p2o)
data_grp.attrs.create("n_p2l", data=n_p2l)
data_grp.attrs.create("n_fwj", data=n_fwj)
data_grp.attrs.create("n_fwo", data=n_fwo)
data_grp.attrs.create("n_fwl", data=n_fwl)
data_grp.attrs.create("n_fwa", data=n_fwl)
data_grp.attrs.create("n_bsj", data=n_bsa)
data_grp.attrs.create("n_bso", data=n_bso)
data_grp.attrs.create("n_bsl", data=n_bsl)
data_grp.attrs.create("n_bsa", data=n_bsa)

emcee_grp = h5f.create_group("emcee")

emcee_grp.attrs.create("kind_sin", data=kind_sin)
emcee_grp.attrs.create("n_planets", data=n_planets)
emcee_grp.attrs.create("n_periods", data=n_periods)
emcee_grp.attrs.create("n_trends", data=n_trends)
emcee_grp.attrs.create("n_orbpams", data=n_orbpams)
emcee_grp.attrs.create("n_syspams", data=n_syspams)
emcee_grp.attrs.create("n_extpams", data=n_extpams)
emcee_grp.attrs.create("n_dataset", data=n_dataset)
emcee_grp.attrs.create("n_pin", data=n_pin)
emcee_grp.attrs.create("n_amp", data=n_amp)
emcee_grp.attrs.create("n_pha", data=n_pha)
emcee_grp.attrs.create("n_pho", data=n_pha)
emcee_grp.attrs.create("n_pof", data=n_pof)

emcee_grp.attrs.create("ndim", data=ndim)
emcee_grp.attrs.create("nsteps", data=nsteps)
emcee_grp.attrs.create("ngen", data=ngen)
emcee_grp.attrs.create("npop", data=npop)
emcee_grp.attrs.create("nburn", data=nburn)
emcee_grp.attrs.create("thin", data=thin)

emcee_grp.create_dataset("bound", data=bounds, compression="gzip")
emcee_grp.create_dataset("chain", data=sampler.chain, compression="gzip")

blobs_out = np.asarray(sampler.blobs, dtype=np.double)

emcee_grp.create_dataset("blobs", data=blobs_out, compression="gzip")

## FLATchains are not saved because the BURNin values are mixed inside, making them
# basically useless


emcee_grp.create_dataset("lnprobability", data=sampler.lnprobability, compression="gzip")
emcee_grp.create_dataset("acceptance_fraction", data=sampler.acceptance_fraction, compression="gzip")
emcee_grp.create_dataset("acor", data=sampler.acor, compression="gzip")
# emcee_grp.create_dataset("get_autorr_time", data=sampler.get_autocorr_time(), compression="gzip")


h5f.close()
