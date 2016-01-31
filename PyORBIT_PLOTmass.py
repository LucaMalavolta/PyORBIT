import h5py
import numpy as np
import kepler_exo as kp
from matplotlib import pyplot as plt
import triangle
#from scipy.interpolate import *

print 'Waiting for .config file name...'

file_conf = raw_input()
#file_conf = 'OC101_fast.config'
file=open(file_conf,'r')
row = file.readlines()


for line in row:
  if line.find("Output") > -1:
    info = line.split()
    planet_name = info[1]

h5f = h5py.File('output/'+planet_name+'.hdf5', "r")

for item in h5f.attrs.keys():
  print item + ":", h5f.attrs[item]
  print item + ":", h5f.attrs[item]
h5f_data  = h5f['/data']
h5f_emcee = h5f['/emcee']
for item in h5f_data.attrs.keys():
  print item + ":", h5f_data.attrs[item]

  if item == 'flag_rv': flag_rv = h5f_data.attrs[item]
  if item == 'flag_p1': flag_p1 = h5f_data.attrs[item]
  if item == 'flag_p2': flag_p2 = h5f_data.attrs[item]
  if item == 'flag_fw': flag_fw = h5f_data.attrs[item]
  if item == 'flag_bs': flag_bs = h5f_data.attrs[item]
  if item == 'Tref'   : Tref    = h5f_data.attrs[item]
  if item == 'rv_n' : rv_n = h5f_data.attrs[item]
  if item == 'p1_n' : p1_n = h5f_data.attrs[item]
  if item == 'p2_n' : p2_n = h5f_data.attrs[item]
  if item == 'fw_n' : fw_n = h5f_data.attrs[item]
  if item == 'bs_n' : bs_n = h5f_data.attrs[item]

#  if item == 'n_rvj' : n_rvj = h5f_data.attrs[item]
#  if item == 'n_rvo' : n_rvo = h5f_data.attrs[item]
#  if item == 'n_rvl' : n_rvl = h5f_data.attrs[item]
#  if item == 'n_rva' : n_rva = h5f_data.attrs[item]
#  if item == 'n_p1j' : n_p1j = h5f_data.attrs[item]
#  if item == 'n_p1o' : n_p1o = h5f_data.attrs[item]
#  if item == 'n_p1l' : n_p1l = h5f_data.attrs[item]
#  if item == 'n_p2j' : n_p2j = h5f_data.attrs[item]
#  if item == 'n_p2o' : n_p2o = h5f_data.attrs[item]
#  if item == 'n_p2l' : n_p2l = h5f_data.attrs[item]
#  if item == 'n_fwj' : n_fwj = h5f_data.attrs[item]
#  if item == 'n_fwo' : n_fwo = h5f_data.attrs[item]
#  if item == 'n_fwl' : n_fwl = h5f_data.attrs[item]
#  if item == 'n_fwa' : n_fwa = h5f_data.attrs[item]
#  if item == 'n_bsj' : n_bsj = h5f_data.attrs[item]
#  if item == 'n_bso' : n_bso = h5f_data.attrs[item]
#  if item == 'n_bsl' : n_bsl = h5f_data.attrs[item]
#  if item == 'n_bsa' : n_bsa = h5f_data.attrs[item]

for item in h5f_emcee.attrs.keys():
  print item + ":", h5f_emcee.attrs[item]
#  if item == 'kind_sin'          : kind_sin  = h5f_emcee.attrs[item]
  if item == 'n_planets'         : n_planets = h5f_emcee.attrs[item]
#  if item == 'n_periods'         : n_periods = h5f_emcee.attrs[item]
  if item == 'n_trends'          : n_trends = h5f_emcee.attrs[item]
  if item == 'n_orbpams'         : n_orbpams = h5f_emcee.attrs[item]
#  if item == 'n_syspams'         : n_syspams = h5f_emcee.attrs[item]
#  if item == 'n_sinpams'         : n_sinpams = h5f_emcee.attrs[item]
#  if item == 'n_extpams'         : n_extpams = h5f_emcee.attrs[item]
#  if item == 'n_dataset'          : n_dataset  = h5f_emcee.attrs[item]
#  if item == 'n_pin'             : n_pin     = h5f_emcee.attrs[item]
#  if item == 'n_amp'             : n_amp     = h5f_emcee.attrs[item]
#  if item == 'n_pha'             : n_pha     = h5f_emcee.attrs[item]
#  if item == 'n_pof'             : n_pof     = h5f_emcee.attrs[item]
  if item == 'ndim'              : ndim      = h5f_emcee.attrs[item]
  if item == 'nsteps'            : nsteps    = h5f_emcee.attrs[item]
  if item == 'ngen'              : ngen      = h5f_emcee.attrs[item]
  if item == 'npop'              : npop      = h5f_emcee.attrs[item]
  if item == 'nburn'             : nburn     = h5f_emcee.attrs[item]
  if item == 'thin'              : thin      = h5f_emcee.attrs[item]


nburn = nburn/thin
nsteps = nsteps/thin

#if flag_rv:
#  data_rv  = h5f['/data/data_rv']
#  rv_x = np.asarray(data_rv[:,0], dtype=np.double)
#  rv_y = np.asarray(data_rv[:,1], dtype=np.double)
#  rv_e = np.asarray(data_rv[:,2], dtype=np.double)
#  rv_j = np.asarray(data_rv[:,3], dtype=np.double)
#  rv_o = np.asarray(data_rv[:,4], dtype=np.double)
#  rv_l = np.asarray(data_rv[:,5], dtype=np.double)
#
#  rvj_mask = np.asarray(h5f['/data/rvj_mask'], dtype=np.double)
#  rvo_mask = np.asarray(h5f['/data/rvo_mask'], dtype=np.double)
#  rvl_mask = np.asarray(h5f['/data/rvl_mask'], dtype=np.double)
#  rva_mask = np.asarray(h5f['/data/rva_mask'], dtype=np.double)
#
### Reading photometric data, if present:
#if flag_p1:
#  data_p1  = h5f['/data/data_p1']
#  p1_x = np.asarray(data_p1[:,0], dtype=np.double)
#  p1_y = np.asarray(data_p1[:,1], dtype=np.double)
#  p1_e = np.asarray(data_p1[:,2], dtype=np.double)
#  p1_j = np.asarray(data_p1[:,3], dtype=np.double)
#  p1_o = np.asarray(data_p1[:,4], dtype=np.double)
#  p1_l = np.asarray(data_p1[:,5], dtype=np.double)
#
#  p1j_mask = np.asarray(h5f['/data/p1j_mask'], dtype=np.double)
#  p1o_mask = np.asarray(h5f['/data/p1o_mask'], dtype=np.double)
#  p1l_mask = np.asarray(h5f['/data/p1l_mask'], dtype=np.double)
#
## photometric data in a different filter
#if flag_p2:
#  data_p2  = h5f['/data/data_p2']
#  p2_x = np.asarray(data_p2[:,0], dtype=np.double)
#  p2_y = np.asarray(data_p2[:,1], dtype=np.double)
#  p2_e = np.asarray(data_p2[:,2], dtype=np.double)
#  p2_j = np.asarray(data_p2[:,3], dtype=np.double)
#  p2_o = np.asarray(data_p2[:,4], dtype=np.double)
#  p2_l = np.asarray(data_p2[:,5], dtype=np.double)
#
#  p2j_mask = np.asarray(h5f['/data/p2j_mask'], dtype=np.double)
#  p2o_mask = np.asarray(h5f['/data/p2o_mask'], dtype=np.double)
#  p2l_mask = np.asarray(h5f['/data/p2l_mask'], dtype=np.double)
### Reading FWHM data, if present:
#if flag_fw:
#  data_fw  = h5f['/data/data_fw']
#  fw_x = np.asarray(data_fw[:,0], dtype=np.double)
#  fw_y = np.asarray(data_fw[:,1], dtype=np.double)
#  fw_e = np.asarray(data_fw[:,2], dtype=np.double)
#  fw_j = np.asarray(data_fw[:,3], dtype=np.double)
#  fw_o = np.asarray(data_fw[:,4], dtype=np.double)
#  fw_l = np.asarray(data_fw[:,5], dtype=np.double)
#
#  fwj_mask = np.asarray(h5f['/data/fwj_mask'], dtype=np.double)
#  fwo_mask = np.asarray(h5f['/data/fwo_mask'], dtype=np.double)
#  fwl_mask = np.asarray(h5f['/data/fwl_mask'], dtype=np.double)
#  fwa_mask = np.asarray(h5f['/data/fwa_mask'], dtype=np.double)
#
### Reading BIS data, if present:
#if flag_bs:
#  data_bs  = h5f['/data/data_bs']
#  bs_x = np.asarray(data_bs[:,0], dtype=np.double)
#  bs_y = np.asarray(data_bs[:,1], dtype=np.double)
#  bs_e = np.asarray(data_bs[:,2], dtype=np.double)
#  bs_j = np.asarray(data_bs[:,3], dtype=np.double)
#  bs_o = np.asarray(data_bs[:,4], dtype=np.double)
#  bs_l = np.asarray(data_bs[:,5], dtype=np.double)
#
#  bsj_mask = np.asarray(h5f['/data/bsj_mask'], dtype=np.double)
#  bso_mask = np.asarray(h5f['/data/bso_mask'], dtype=np.double)
#  bsl_mask = np.asarray(h5f['/data/bsl_mask'], dtype=np.double)
#  bsa_mask = np.asarray(h5f['/data/bsa_mask'], dtype=np.double)

bounds = h5f['/emcee/bound']
chain  = h5f['/emcee/chain']
blobs  = h5f['/emcee/blobs']
lnprobability  = h5f['/emcee/lnprobability']
acceptance_fraction  = h5f['/emcee/acceptance_fraction']
acor                 = h5f['/emcee/acor']

print 'Acceptance Fraction for all walkers:'
print acceptance_fraction[:]



## We start by plotting the chains
#out_med = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(sampler.flatchain[nburn:,:], #[16, 50, 84],axis=0)))
#for ii in xrange(0,ndim):
#  print out_med[ii]


# chain is transposed
chain_T =  np.ndarray([nsteps,npop,ndim],dtype=np.double)
for ii in xrange(0,ndim):
  chain_T[:,:,ii] = chain[:,:,ii].T



lnprb_T =  lnprobability[:][:].T


chain_burnt = chain_T[nburn:,:,:]
blobs_burnt = blobs[nburn:,:,:]
lnprb_burnt = lnprb_T[nburn:,:]

#INSERIRE ELEMENTO DI THINNING!!!!

#flattening the blobs
s = chain_burnt.shape
flatchain = chain_burnt.reshape(s[1] * s[0], s[2])
s = blobs_burnt[:].shape
flatblobs = blobs_burnt.reshape(s[1] * s[0], s[2])
s = lnprb_burnt[:].shape
flatlnprb = lnprb_burnt.reshape(s[1] * s[0])

n_kept = s[1] * s[0]

#index of maximum probability values
maxlnprob_ind = flatlnprb.argmax()


#Obtaining the 1-sigma errors on values
chain_med = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), \
              zip(*np.percentile(flatchain[:,:], [15.865, 50, 84.135],axis=0)))
blobs_med = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), \
              zip(*np.percentile(flatblobs[:,:], [15.865, 50, 84.135],axis=0)))


ip = 0
id = 0

id = 0
ip = 0

phase_factor = 180.0/np.pi
omega_factor = 180.0/np.pi
orbit_pam = np.zeros((n_planets,5),dtype=np.double)

if flag_rv:
  for pp in xrange(0,n_planets):
    # Period
    orbit_pam[pp,0] = blobs_med[id][0]
    id += 1


    orbit_pam[pp,1] = blobs_med[id][0]

    id += 1

    # Phase
    orbit_pam[pp,2] = chain_med[ip+2][0] * phase_factor 

    if (n_orbpams[pp] == 5):
      # eccentricity

      orbit_pam[pp,3] = blobs_med[id][0]
      id += 1

      # Omega

      orbit_pam[pp,4] = blobs_med[id][0]*omega_factor
      id += 1
    ip += n_orbpams[pp]
    print orbit_pam[pp,:]






G_grav=6.67398e-11
M_sun =1.98892e30
M_jup = 1.89813e27
M_ratio = M_sun/M_jup


M_star1 = 0.940
M_star1_err = 0.02
import scipy.optimize

def get_mass(M_star2,M_star1,Period,K1,e0):
  # M_star1, M_star2 in solar masses
  # P in days -> Period is converted in seconds in the routine
  # i in degrees
  # Gravitational constant is given in m^3 kg^-1 s^-2
  # output in m/s
  output =  K1 - (2.*np.pi*G_grav*M_sun/86400.0)**(1.0/3.0) * (1.000 / np.sqrt(1.0-e0**2.0)) * (Period) ** (-1.0/3.0) * ( M_star2*(M_star1+M_star2)**(-2.0/3.0))
  return output

Mu_sun = 132712440018.9
seconds_in_day = 86400
AU_km = 1.4960*10**8
#Triangle plot for each planet
id = 0
ip = 0
x0 = 1./150

include_semiaxis = False
## RVs plot
labels_all = []
labels_plt = []

if include_semiaxis:
  labels_orb5 = ["P", "K", "Phase", "e", "$\omega$", "M$_{Jup} \sin i$","a"]
  labels_orb3 = ["P", "K", "Phase", "M$_{Jup} \sin i","a" ]
  additional_pams = 2
else:
  labels_orb5 = ["P", "K", "Phase", "e", "$\omega$", "M$_{Jup} \sin i$"]
  labels_orb3 = ["P", "K", "Phase", "M$_{Jup} \sin i" ]
  additional_pams = 1

for pp in xrange(0,n_planets):
  labels_plt = []
  plot_truths = np.zeros(n_orbpams[pp]+additional_pams)

  M_star1_rand = np.random.normal(M_star1,M_star1_err,n_kept)
  keep_a = np.zeros(n_kept)

  sample_plan = np.zeros([n_kept,n_orbpams[pp]+additional_pams])
  sample_plan[:,0] = flatblobs[:n_kept,id+0]
  sample_plan[:,1] = flatblobs[:n_kept,id+1]
  sample_plan[:,2] = flatchain[:n_kept,ip+2]*phase_factor

  plot_truths[0:3]=orbit_pam[pp,0:3]
  print orbit_pam[pp,2], orbit_pam[pp,2]<0.00
  if (orbit_pam[pp,2]<0.00):
      plot_truths[2] += 1.0
      sample_plan[:,2]  += 1.0


  id += 2
  if (n_orbpams[pp]==5) :
    sample_plan[:,3] = flatblobs[:n_kept,id+0]
    sample_plan[:,4] = flatblobs[:n_kept,id+1]*omega_factor
    plot_truths[3] = orbit_pam[pp,3]
    plot_truths[4] = orbit_pam[pp,4]

    peri_time = Tref + (-sample_plan[:,2]+sample_plan[:,4])/360.00 *sample_plan[:,0]
    peri_med = np.percentile(peri_time[:], [15.865, 50, 84.135],axis=0)
    print ' Time of periastron: ', peri_med[1], peri_med[2]-peri_med[1],peri_med[1]-peri_med[0]
    

    id += 2
    labels_plt +=  labels_orb5
    for ii in xrange(0,n_kept):
      sample_plan[ii,5] =M_ratio* scipy.optimize.fsolve(get_mass, x0, args=(M_star1_rand[ii],sample_plan[ii,0],sample_plan[ii,1],sample_plan[ii,3]))
      keep_a[ii] = np.power( (Mu_sun*np.power(sample_plan[ii,0]*seconds_in_day/(2*np.pi),2.00)/(AU_km**3)) * M_star1_rand[ii],1.00/3.00 )

      #if include_semiaxis: sample_plan[ii,6] = np.power( (Mu_sun*np.power(sample_plan[ii,0]*seconds_in_day/(2*np.pi),2.00)/(AU_km**3)) * M_star1_rand[ii],1.00/3.00 )
    if include_semiaxis: sample_plan[:,6] = keep_a

    pams_med = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), \
              zip(*np.percentile(sample_plan[:,:], [15.865, 50, 84.135],axis=0)))

    plot_truths[5] = pams_med[5][0]
    if include_semiaxis: plot_truths[6] = pams_med[6][0]

  else:
    labels_plt +=  labels_orb3
    for ii in xrange(0,n_kept):
      sample_plan[ii,3] = M_ratio*scipy.optimize.fsolve(get_mass, x0, args=(M_star1_rand[ii],sample_plan[ii,0],sample_plan[ii,1],0.00000000))
      keep_a[ii] = np.power( (Mu_sun*np.power(sample_plan[ii,0]*seconds_in_day/(2*np.pi),2.00)/(AU_km**3)) * M_star1_rand[ii],1.00/3.00 )

      #if include_semiaxis: sample_plan[ii,4] = np.power( (Mu_sun*np.power(sample_plan[ii,0]*seconds_in_day/(2*np.pi),2.00)/(AU_km**3)) * M_star1_rand[ii],1.00/3.00 )
    if include_semiaxis: sample_plan[:,4] = keep_a
    pams_med = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), \
              zip(*np.percentile(sample_plan[:,:], [15.865, 50, 84.135],axis=0)))
    plot_truths[3] = pams_med[3][0]
    if include_semiaxis: plot_truths[4] = pams_med[4][0]

  print 'Planet ',pp, labels_plt
  print pams_med
  print 'A -> precentiles', np.percentile(keep_a, [15.865, 50, 84.135])

  fig = triangle.corner(sample_plan, labels=labels_plt, truths=plot_truths)
  fig.savefig("output/"+planet_name+"_planet"+`pp`+"_RVtriangle_Mass_wp_wt.pdf", bbox_inches='tight')
  plt.close()
  #
  fig = triangle.corner(sample_plan, labels=labels_plt, truths=plot_truths, plot_contours=True, plot_datapoints=False)
  fig.savefig("output/"+planet_name+"_planet"+`pp`+"_RVtriangle_Mass.pdf", bbox_inches='tight')
  plt.close()

  ip += n_orbpams[pp]

#

#OC101_fast.config
