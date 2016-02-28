import h5py
import numpy as np
import kepler_exo as kp
from matplotlib import pyplot as plt
import triangle
#from scipy.interpolate import *

from getdist import plots, MCSamples
import getdist, IPython
print'Version: ',getdist.__version__


print 'Waiting for .config file name...'

file_conf = raw_input()
#file_conf = 'OC101_fast.config'
file=open(file_conf,'r')
row = file.readlines()

bin_num = 20
print 'bin_num: ',bin_num

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

  if item == 'n_rvj' : n_rvj = h5f_data.attrs[item]
  if item == 'n_rvo' : n_rvo = h5f_data.attrs[item]
  if item == 'n_rvl' : n_rvl = h5f_data.attrs[item]
  if item == 'n_rva' : n_rva = h5f_data.attrs[item]
  if item == 'n_p1j' : n_p1j = h5f_data.attrs[item]
  if item == 'n_p1o' : n_p1o = h5f_data.attrs[item]
  if item == 'n_p1l' : n_p1l = h5f_data.attrs[item]
  if item == 'n_p2j' : n_p2j = h5f_data.attrs[item]
  if item == 'n_p2o' : n_p2o = h5f_data.attrs[item]
  if item == 'n_p2l' : n_p2l = h5f_data.attrs[item]
  if item == 'n_fwj' : n_fwj = h5f_data.attrs[item]
  if item == 'n_fwo' : n_fwo = h5f_data.attrs[item]
  if item == 'n_fwl' : n_fwl = h5f_data.attrs[item]
  if item == 'n_fwa' : n_fwa = h5f_data.attrs[item]
  if item == 'n_bsj' : n_bsj = h5f_data.attrs[item]
  if item == 'n_bso' : n_bso = h5f_data.attrs[item]
  if item == 'n_bsl' : n_bsl = h5f_data.attrs[item]
  if item == 'n_bsa' : n_bsa = h5f_data.attrs[item]

for item in h5f_emcee.attrs.keys():
  print item + ":", h5f_emcee.attrs[item]
  if item == 'kind_sin'          : kind_sin  = h5f_emcee.attrs[item]
  if item == 'n_planets'         : n_planets = h5f_emcee.attrs[item]
  if item == 'n_periods'         : n_periods = h5f_emcee.attrs[item]
  if item == 'n_trends'          : n_trends = h5f_emcee.attrs[item]
  if item == 'n_orbpams'         : n_orbpams = h5f_emcee.attrs[item]
  if item == 'n_syspams'         : n_syspams = h5f_emcee.attrs[item]
  if item == 'n_sinpams'         : n_sinpams = h5f_emcee.attrs[item]
  if item == 'n_rotpams'         : n_rotpams = h5f_emcee.attrs[item]
  if item == 'n_extpams'         : n_extpams = h5f_emcee.attrs[item]
  if item == 'n_dataset'          : n_dataset  = h5f_emcee.attrs[item]
  if item == 'n_pin'             : n_pin     = h5f_emcee.attrs[item]
  if item == 'n_amp'             : n_amp     = h5f_emcee.attrs[item]
  if item == 'n_pha'             : n_pha     = h5f_emcee.attrs[item]
  if item == 'n_pof'             : n_pof     = h5f_emcee.attrs[item]
  if item == 'ndim'              : ndim      = h5f_emcee.attrs[item]
  if item == 'nsteps'            : nsteps    = h5f_emcee.attrs[item]
  if item == 'ngen'              : ngen      = h5f_emcee.attrs[item]
  if item == 'npop'              : npop      = h5f_emcee.attrs[item]
  if item == 'nburn'             : nburn     = h5f_emcee.attrs[item]
  if item == 'thin'              : thin      = h5f_emcee.attrs[item]


nburn = nburn/thin
nsteps = nsteps/thin

n_pho = np.amax(n_amp)


bounds = h5f['/emcee/bound']
chain  = h5f['/emcee/chain']
lnprobability  = h5f['/emcee/lnprobability']
acceptance_fraction  = h5f['/emcee/acceptance_fraction']
acor                 = h5f['/emcee/acor']

print 'Acceptance Fraction for all walkers:'
print acceptance_fraction[:]


## We start by plotting the chains
#out_med = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(sampler.flatchain[nburn:,:], #[16, 50, 84],axis=0)))
#for ii in xrange(0,ndim):
#  print out_med[ii]


### Start plotting the chains
print 'Tref = ',Tref


# chain is transposed
chain_T =  np.ndarray([nsteps,npop,ndim],dtype=np.double)
for ii in xrange(0,ndim):
  chain_T[:,:,ii] = chain[:,:,ii].T

lnprb_T =  lnprobability[:][:].T

lnprob_median = np.median(lnprb_T[:,:])
fig=plt.plot(lnprb_T[:200,:], '-', color='k', alpha=0.3)
plt.axhline(lnprob_median)
plt.savefig('output/'+planet_name+'_LNprob_chain.png', bbox_inches='tight')
plt.close()
print ' LNprob median = ', lnprob_median


chain_burnt = chain_T[:nburn,:,:]
burnt_size = np.size(chain_burnt, axis=0)

print 'Burnt size: ', burnt_size
#lnprb_burnt = lnprb_T[nburn:,:]

#INSERIRE ELEMENTO DI THINNING!!!!

#flattening the chain
#s = chain_burnt.shape
#flatchain = chain_burnt.reshape(s[1] * s[0], s[2])

#s = lnprb_burnt[:].shape
#flatlnprb = lnprb_burnt.reshape(s[1] * s[0])

#n_kept = s[1] * s[0]

orbit_pam = np.zeros((n_planets,5),dtype=np.double)
system_pam = np.zeros(n_rvj+n_rvo+n_rvl)


def GelmanRubin(chains_input):
    n_iters, n_chain = np.shape(chains_input)
    W = np.asarray(0.,dtype=np.double)
    z_pc =  np.sum(chains_input,axis=0)/n_iters  #eq 20
    for nn in xrange(0,n_chain):
        W += (np.sum(np.power(chains_input[:,nn] - z_pc[nn],2))) / ((n_iters-1) * n_chain)  #eq 21
    z_pp = np.sum(chains_input) / (n_chain*n_iters)
    B = np.sum(np.power(z_pc-z_pp,2))*(n_chain/(n_iters-1))
    var = W*(n_chain-1)/n_chain + B/n_chain
    return np.sqrt(var/W)

for nd in xrange(0,ndim): #(0,ndim):
    out_lines = np.zeros(nburn)
    for ii in xrange(20,nburn):
        out_lines[ii]=GelmanRubin(chain_T[:ii,:,nd])

    plt.ylim(0.95,2.3 )
    plt.plot(out_lines[20:], '-', color='k')
    plt.axhline(1.01)
    plt.savefig('output/'+planet_name+'_GRtrace_pam_'+`nd`+'.png', bbox_inches='tight')
    plt.close()


###   print 'Starting the GetDist section'
###
###   #Get the getdist MCSamples objects for the samples, specifying same parameter
###   #na###   mes and labels; if not specified weights are assumed to all be unity
###   names = ["x%s"%i for i in range(ndim)]
###   labels =  ["x_%s"%i for i in range(ndim)]
###
###   samp_out = np.zeros([nsteps-nburn,ndim+2])
###   samp_out[:,0]  = 1
###
###   for nk in xrange(0,npop):
###       samp_out[:,2:] = chain_burnt[:,nk,:]
###       if nk==0:
###           np.savetxt('./output/tmp/samps.txt',samp_out,fmt='%.18e')
###       else:
###           np.savetxt('./output/tmp/samps_'+`nk`+'.txt',samp_out,fmt='%.18e')
###
###   chan = getdist.mcsamples.loadMCSamples('./output/tmp/samps')
###   chan.getSeparateChains()
###   print chan.getConvergeTests()
###
###   ## >>> aa=np.ones([3,5])
###   ## >>> print aa
###   ## [[ 1.  1.  1.  1.  1.]
###   ##  [ 1.  1.  1.  1.  1.]
###   ##  [ 1.  1.  1.  1.  1.]]
###   ## >>> np.savetxt('aaa.out',aa,fmt='%.18e')
###   ## >>> quit()
###   ## ~/CODE/PyORBIT $ more aaa.out
###   ## 1.000000000000000000e+00 1.000000000000000000e+00 1.000000000000000000e+00 1.000000000000000000e+00 1.000000000000000000e+00
###   ## 1.000000000000000000e+00 1.000000000000000000e+00 1.000000000000000000e+00 1.000000000000000000e+00 1.000000000000000000e+00
###   ## 1.000000000000000000e+00 1.000000000000000000e+00 1.000000000000000000e+00 1.000000000000000000e+00 1.000000000000000000e+00
###   ##