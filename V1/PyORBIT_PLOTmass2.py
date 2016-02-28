import h5py
import numpy as np
import kepler_exo as kp
from matplotlib import pyplot as plt
import corner
import csv
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

h5f_data  = h5f['/data']
h5f_emcee = h5f['/emcee']
for item in h5f_data.attrs.keys():

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


for item in h5f_emcee.attrs.keys():
    if item == 'n_planets'         : n_planets = h5f_emcee.attrs[item]
    if item == 'n_trends'          : n_trends = h5f_emcee.attrs[item]
    if item == 'n_orbpams'         : n_orbpams = h5f_emcee.attrs[item]
    if item == 'ndim'              : ndim      = h5f_emcee.attrs[item]
    if item == 'nsteps'            : nsteps    = h5f_emcee.attrs[item]
    if item == 'ngen'              : ngen      = h5f_emcee.attrs[item]
    if item == 'npop'              : npop      = h5f_emcee.attrs[item]
    if item == 'nburn'             : nburn     = h5f_emcee.attrs[item]
    if item == 'thin'              : thin      = h5f_emcee.attrs[item]


nburn = nburn/thin
nsteps = nsteps/thin

bounds = h5f['/emcee/bound']
chain  = h5f['/emcee/chain']
blobs  = h5f['/emcee/blobs']
lnprobability  = h5f['/emcee/lnprobability']
acceptance_fraction  = h5f['/emcee/acceptance_fraction']
acor                 = h5f['/emcee/acor']

#Chains are transposed
chain_T =  np.ndarray([nsteps,npop,ndim],dtype=np.double)
for ii in xrange(0,ndim):
    chain_T[:,:,ii] = chain[:,:,ii].T


lnprb_T =  lnprobability[:][:].T

chain_burnt = chain_T[nburn:,:,:]
blobs_burnt = blobs[nburn:,:,:]
lnprb_burnt = lnprb_T[nburn:,:]

#flattening the blobs
s = chain_burnt.shape
flatchain = chain_burnt.reshape(s[1] * s[0], s[2])
s = blobs_burnt[:].shape
flatblobs = blobs_burnt.reshape(s[1] * s[0], s[2])
s = lnprb_burnt[:].shape
flatlnprb = lnprb_burnt.reshape(s[1] * s[0])

n_kept = s[1] * s[0]

# index of maximum probability values
maxlnprob_ind = flatlnprb.argmax()


# Obtaining the 1-sigma errors on values
chain_med = np.percentile(flatchain[:,:], 50,axis=0)
blobs_med = np.percentile(flatblobs[:,:], 50,axis=0)


id_list = 0
ip_list = 0

phase_factor = 180.0/np.pi
omega_factor = 180.0/np.pi
orbit_pam = np.zeros((n_planets,5),dtype=np.double)

if flag_rv:
    for pp in xrange(0,n_planets):
        # Period
        orbit_pam[pp,0] = blobs_med[id_list]
        id_list += 1


        orbit_pam[pp,1] = blobs_med[id_list]

        id_list += 1

        # Phase
        orbit_pam[pp,2] = chain_med[ip_list + 2]

        if (n_orbpams[pp] == 5):
            # eccentricity

            orbit_pam[pp,3] = blobs_med[id_list]
            id_list += 1

            # Omega

            orbit_pam[pp,4] = blobs_med[id_list]
            id_list += 1
        ip_list += n_orbpams[pp]
        print orbit_pam[pp,:]


G_grav=6.67398e-11
M_sun =1.98892e30
M_jup = 1.89813e27
M_ratio = M_sun/M_jup


M_star1 = 0.935
M_star1_err = 0.013
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
id_list = 0
ip_list = 0
x0 = 1./150


include_semiaxis = False
## RVs plot
labels_all = []
labels_plt = []


labels_orb5 = ["P", "K", "e", "M$_{Jup} \sin i$"]
labels_orb3 = ["P", "K", "M$_{Jup} \sin i" ]
additional_pams = 2

for pp in xrange(0,n_planets):
    labels_plt = []

    if (n_orbpams[pp]==5) :
        output_n = 4
    else:
        output_n = 4

    plot_truths = np.zeros(output_n)
    output_plan = np.zeros([n_kept, output_n])


    M_star1_rand = np.random.normal(M_star1,M_star1_err,n_kept)
    keep_a = np.zeros(n_kept)

    sample_plan = np.zeros([n_kept,5+additional_pams])
    sample_plan[:,0] = flatblobs[:n_kept, id_list + 0]
    sample_plan[:,1] = flatblobs[:n_kept, id_list + 1]
    sample_plan[:,2] = flatchain[:n_kept, ip_list + 2]

    plot_truths[0]=orbit_pam[pp,0]
    plot_truths[1]=orbit_pam[pp,1]
    print orbit_pam[pp,2], orbit_pam[pp,2]<0.00
    if (orbit_pam[pp,2]<0.00):
        sample_plan[:,2]  += 1.0

    ip_list += n_orbpams[pp]
    id_list += 2

    if (n_orbpams[pp]==5) :

        sample_plan[:,3] = flatblobs[:n_kept, id_list + 0]
        sample_plan[:,4] = flatblobs[:n_kept, id_list + 1]
        id_list += 2

        peri_time = Tref + (-sample_plan[:,2]+sample_plan[:,4])/360.00 *sample_plan[:,0]
        peri_med = np.percentile(peri_time[:], [15.865, 50, 84.135],axis=0)
        print ' Time of periastron: ', peri_med[1], peri_med[2]-peri_med[1],peri_med[1]-peri_med[0]


    for ii in xrange(0,n_kept):
        sample_plan[ii,5] =M_ratio* scipy.optimize.fsolve(get_mass, x0, args=(M_star1_rand[ii],sample_plan[ii,0],sample_plan[ii,1],sample_plan[ii,3]))
        sample_plan[ii,6] = np.power( (Mu_sun*np.power(sample_plan[ii,0]*seconds_in_day/(2*np.pi),2.00)/(AU_km**3)) * M_star1_rand[ii],1.00/3.00 )


    if (n_orbpams[pp]==5) :
        plot_truths[2] = orbit_pam[pp,3]
        plot_truths[3] = np.percentile(sample_plan[:,5],50)

        labels_plt +=  labels_orb5

        output_plan[:,0] = sample_plan[:,0]
        output_plan[:,1] = sample_plan[:,1]
        output_plan[:,2] = sample_plan[:,3]
        output_plan[:,3] = sample_plan[:,5]

    else:
        plot_truths[2] = np.percentile(sample_plan[:,3],50)
        labels_plt +=  labels_orb3

        #pams_med = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), \
        #               zip(*np.percentile(sample_plan[:,:], [15.865, 50, 84.135],axis=0)))
        output_plan[:,0] = sample_plan[:,0]
        output_plan[:,1] = sample_plan[:,1]
        output_plan[:,2] = sample_plan[:,5]


    print 'Planet ',pp, labels_plt
    print 'M -> precentiles', np.percentile(sample_plan[:,-2], [0.135, 2.275, 15.865, 50, 84.135, 97.725, 99.865])
    print 'A -> precentiles', np.percentile(sample_plan[:,-1], [0.135, 2.275, 15.865, 50, 84.135, 97.725, 99.865])

    pams_limits = np.zeros([n_orbpams[pp],2])
    for ii in xrange(0,n_orbpams[pp]):
        pams_limits[ii,:] = np.percentile(sample_plan[:,ii], [0.135, 99.865])
        print ii, pams_limits[ii,:]


    #fig = corner.corner(output_plan, labels=labels_plt, truths=plot_truths)
    #fig.savefig("output/"+planet_name+"_planet"+`pp`+"_RVtriangle_Mass2_wp_wt.pdf", bbox_inches='tight')
    #plt.close()
    #
    #fig = corner.corner(output_plan, labels=labels_plt, truths=plot_truths, plot_contours=True, plot_datapoints=False)
    #fig.savefig("output/"+planet_name+"_planet"+`pp`+"_RVtriangle_Mass2.pdf", bbox_inches='tight')
    #plt.close()

    ## Creating data faile for MCMC distributions
    list_labels=['P', 'K', 'e', 'm']

    n_int = 4
    n_bins = 60 + 1

    h5f = h5py.File('output/'+planet_name+'_planet'+`pp`+'_hist1d.hdf5', "w")
    data_grp = h5f.create_group("hist1d")

    data_lim = np.zeros([n_kept,2],dtype=np.double)
    data_edg = np.zeros([n_kept, n_bins],dtype=np.double)
    for ii in xrange(0,n_int):
        data_lim[ii,:] = [np.amin(output_plan[:,ii]),np.amax(output_plan[:,ii])]
        data_edg[ii,:] = np.linspace(data_lim[ii,0],data_lim[ii,1],n_bins)

    for ii in xrange(0,n_int):
        for jj in xrange(ii,n_int):
            x_data = output_plan[:,ii]
            y_data = output_plan[:,jj]
            x_edges=data_edg[ii,:]
            y_edges=data_edg[jj,:]

            if (ii!=jj):
                #hist2d = np.histogram2d(x_data, y_data, bins=[x_edges,y_edges], normed=True)
                hist2d = np.histogram2d(x_data, y_data, bins=[x_edges,y_edges], normed=True)
                hist1d_y = np.histogram(y_data, bins=y_edges, normed=True)

                Hflat = hist2d[0].flatten()
                inds = np.argsort(Hflat)[::-1]
                Hflat = Hflat[inds]
                sm = np.cumsum(Hflat)
                sm /= sm[-1]

                #levels = 1.0 - np.exp(-0.5 * np.arange(1.0, 3.1, 1.0) ** 2)

                #V = np.empty(len(levels))
                #for i, v0 in enumerate(levels):
                #    try:
                #        V[i] = Hflat[sm <= v0][-1]
                #    except:
                #        V[i] = Hflat[0]
                #print list_labels[ii], '_',list_labels[jj]
                #print '   ',np.amax(hist2d[0])
                #print '   ',V
                #print '   ',V/np.amax(hist2d[0])

                x_edges_1d = (x_edges[1:]+ x_edges[:-1])/2
                y_edges_1d = (y_edges[1:]+ y_edges[:-1])/2
                h2d_out = np.zeros([n_bins,n_bins])
                h2d_out[0,1:] = x_edges_1d
                h2d_out[1:,0] = y_edges_1d
                h2d_out[1:,1:] = hist2d[0].T *1. / np.amax(hist2d[0])

                h2d_list =  h2d_out.tolist()
                h2d_list[0][0] = ''
                csvfile = 'output/'+planet_name+'_planet'+`pp`+'_hist2d_'+list_labels[ii]+'_'+list_labels[jj]+'.csv'
                with open(csvfile, "w") as output:
                    writer = csv.writer(output, lineterminator='\n')
                    writer.writerows(h2d_list)
            else:
                hist1d = np.histogram(x_data, bins=x_edges)
                hist1d_norm = hist1d[0]*1. / n_kept
                x_edges_1d = (x_edges[1:]+ x_edges[:-1])/2
                data_grp.create_dataset(list_labels[ii]+'_x', data=x_edges_1d, compression="gzip")
                data_grp.create_dataset(list_labels[ii]+'_y', data=hist1d_norm, compression="gzip")

    #OC101_fast.config

    random_n = 20000
    ii_kep = np.random.randint(low=0, high=n_kept, size=random_n)

    x_kep = np.arange(5900, 7400, 2)
    y_kep = np.zeros([np.size(x_kep), 2])


    x_pha = np.arange(-1.00,2.00,0.005,dtype=np.double)
    y_pha = np.zeros([np.size(x_pha), 2])

    y_flg = np.ones(random_n, dtype=bool)

    y_kep_tmp = np.zeros(np.size(x_kep))
    y_pha_tmp = np.zeros(np.size(x_pha))

    if pams_limits[3,0] < 0.02: pams_limits[3,0]=0.00
    for ii in xrange(0,4):
        y_flg = y_flg & (sample_plan[ii_kep, ii]>=pams_limits[ii,0]) & ( sample_plan[ii_kep, ii]<=pams_limits[ii,1])

    y_kep[:,0] = kp.kepler_RV_T0P(x_kep - Tref, orbit_pam[pp,2], orbit_pam[pp,0], orbit_pam[pp,1], orbit_pam[pp,3], orbit_pam[pp,4])
    y_kep[:,1] = y_kep[:,0]
        # value initialization
    y_pha[:,0] = kp.kepler_RV_T0P(x_pha*orbit_pam[pp,0], orbit_pam[pp,2], orbit_pam[pp,0], orbit_pam[pp,1], orbit_pam[pp,3], orbit_pam[pp,4])
    y_pha[:,1] = y_pha[:,0]

    fig1 =plt.plot(x_kep,y_kep[:,0],c='g')
    fig2 =plt.plot(x_pha,y_pha[:,0],c='g')

    print 'Accepted randomizations: ', np.sum(y_flg)
    for ii in xrange(0, random_n):
        ir = ii_kep[ii]
        if (y_flg[ii]):
            y_kep_tmp = kp.kepler_RV_T0P(x_kep - Tref, sample_plan[ir, 2], sample_plan[ir, 0], sample_plan[ir, 1], sample_plan[ir, 3], sample_plan[ir, 4])
            y_pha_tmp = kp.kepler_RV_T0P(x_pha*sample_plan[ir, 0], sample_plan[ir, 2], sample_plan[ir, 0], sample_plan[ir, 1], sample_plan[ir, 3], sample_plan[ir, 4])

            y_kep[:,0] = np.minimum(y_kep[:,0],y_kep_tmp)
            y_kep[:,1] = np.maximum(y_kep[:,1],y_kep_tmp)

            y_pha[:,0] = np.minimum(y_pha[:,0],y_pha_tmp)
            y_pha[:,1] = np.maximum(y_pha[:,1],y_pha_tmp)

    fig1 =plt.plot(x_kep,y_kep[:,0],c='r')
    fig1 =plt.plot(x_kep,y_kep[:,1],c='b')
    fig2 =plt.plot(x_pha,y_pha[:,0],c='r')
    fig2 =plt.plot(x_pha,y_pha[:,1],c='b')

    #h5f = h5py.File('output/'+planet_name+"planet"+`pp`+'_kep.hdf5', "w")
    #data_grp = h5f.create_group("data")
    #data_grp.create_dataset("x_kep", data=x_kep, compression="gzip")
    #data_grp.create_dataset("y_kep", data=y_kep, compression="gzip")
    #h5f.close()
    #h5f = h5py.File('output/'+planet_name+"planet"+`pp`+'_pha.hdf5', "w")
    #data_grp = h5f.create_group("data")
    #data_grp.create_dataset("x_pha", data=x_kep, compression="gzip")
    #data_grp.create_dataset("y_pha", data=y_kep, compression="gzip")
    #h5f.close()

    fileout = open('output/'+planet_name+"planet"+`pp`+'_kep.dat','w')
    fileout.write('descriptor x_kep y_kep,+- \n')
    for ii in xrange(0,np.size(x_kep)):
        fileout.write('{0:14f} {1:14f} {2:14f}  \n'.format(\
        x_kep[ii], (y_kep[ii,1]+y_kep[ii,0])/2,(y_kep[ii,1]-y_kep[ii,0])/2))
    fileout.close()

    fileout = open('output/'+planet_name+"planet"+`pp`+'_pha.dat','w')
    fileout.write('descriptor x_pha y_pha,+- \n')
    for ii in xrange(0,np.size(x_pha)):
        fileout.write('{0:14f} {1:14f} {2:14f}  \n'.format(\
        x_pha[ii], (y_pha[ii,1]+y_pha[ii,0])/2,(y_pha[ii,1]-y_pha[ii,0])/2))
    fileout.close()

    ##fig = plt.plot(x_list,y_list[:,ii])
    #fig1 = plt.savefig("output/"+planet_name+"_planet"+`pp`+"_sample.pdf", bbox_inches='tight')
    #fig1 = plt.close()
    #
    #fig2 = plt.savefig("output/"+planet_name+"_planet"+`pp`+"_sample.pdf", bbox_inches='tight')
    #fig2 = plt.close()

