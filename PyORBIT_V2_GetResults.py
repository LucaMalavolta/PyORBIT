from PyORBIT_V2_Classes import *
import numpy as np
import h5py
import cPickle as pickle
import scipy.optimize
import csv
import os
import corner
import argparse
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

G_grav = 6.67398e-11
M_sun = 1.98892e30
M_jup = 1.89813e27
M_ratio = M_sun / M_jup
EJ_ratio = 317.83
Mu_sun = 132712440018.9
seconds_in_day = 86400
AU_km = 1.4960 * 10 ** 8


def GelmanRubin(chains_input):
    n_iters, n_chain = np.shape(chains_input)
    W = np.asarray(0., dtype=np.double)
    z_pc = np.sum(chains_input, axis=0) / n_iters  # eq 20
    for nn in xrange(0, n_chain):
        W += (np.sum(np.power(chains_input[:, nn] - z_pc[nn], 2))) / ((n_iters - 1) * n_chain)  # eq 21
    z_pp = np.sum(chains_input) / (n_chain * n_iters)
    B = np.sum(np.power(z_pc - z_pp, 2)) * (n_chain / (n_iters - 1))
    var = W * (n_chain - 1) / n_chain + B / n_chain
    return np.sqrt(var / W)


def get_mass(M_star2, M_star1, Period, K1, e0):
    # M_star1, M_star2 in solar masses
    # P in days -> Period is converted in seconds in the routine
    # i in degrees
    # Gravitational constant is given in m^3 kg^-1 s^-2
    # output in m/s
    output = K1 - (2. * np.pi * G_grav * M_sun / 86400.0) ** (1.0 / 3.0) * (1.000 / np.sqrt(1.0 - e0 ** 2.0)) * (
                                                                                                                    Period) ** (
                                                                                                                    -1.0 / 3.0) * (
                      M_star2 * (M_star1 + M_star2) ** (-2.0 / 3.0))
    return output


parser = argparse.ArgumentParser(prog='PyORBIT_V2_GetResults.py', description='Extract results from output MCMC')
# parser.add_argument('-l', type=str, nargs='+', help='line identificator')
parser.add_argument('-i', type=str, nargs='+', required=True, help='config file')
parser.add_argument('-v', type=str, nargs='?', default='False', help='Create Veusz ancillary files')
parser.add_argument('-t', type=str, nargs='?', default='False', help='Create GR traces')
parser.add_argument('-nburn', type=int, nargs='?', default=0, help='emcee burn-ins')

args = parser.parse_args()

file_conf = args.i[0]

# file_conf = raw_input()

mc = ModelContainer()
get_pyorbit_input(file_conf, mc)


M_star1 = mc.star_mass_val
M_star1_err = mc.star_mass_err


dir_output = './' + mc.planet_name + '/'

mc.variable_list = pickle.load(open(dir_output + 'vlist.pick', 'rb'))
mc.scv.use_offset = pickle.load(open(dir_output + 'scv_offset.pick', 'rb'))

print mc.variable_list

h5f = h5py.File(dir_output + mc.planet_name + '.hdf5', "r")

h5f_data = h5f['/data']
h5f_emcee = h5f['/emcee']

for item in h5f_emcee.attrs.keys():
    print item + ":", h5f_emcee.attrs[item]

    if item == 'nwalkers': mc.nwalkers = h5f_emcee.attrs[item]
    if item == 'ndim': mc.ndim = h5f_emcee.attrs[item]

mc.bounds = h5f['/emcee/bound']
chain = h5f['/emcee/chain']
# blobs  = h5f['/emcee/blobs']
lnprobability = h5f['/emcee/lnprobability']
acceptance_fraction = h5f['/emcee/acceptance_fraction']
acor = h5f['/emcee/acor']

print
print '*************************************************************'
print
print 'Acceptance Fraction for all walkers:'
print acceptance_fraction[:]
print
print '*************************************************************'
print

if args.nburn > 0:
    mc.nburn = args.nburn

mc.nsample = mc.nsteps / mc.thin
mc.nburnin = mc.nburn / mc.thin

lnprb_T = lnprobability[:][:].T

chain_T = np.ndarray([mc.nsample, mc.nwalkers, mc.ndim], dtype=np.double)
for ii in xrange(0, mc.ndim):
    chain_T[:, :, ii] = chain[:, :, ii].T

chain_burnt = chain_T[mc.nburnin:, :, :]
s = chain_burnt.shape
flatchain = chain_burnt.reshape(s[1] * s[0], s[2])
n_kept = s[1] * s[0]

chain_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                           zip(*np.percentile(flatchain[:, :], [15.865, 50, 84.135], axis=0))))

mc.results_resumen(chain_med[:, 0])

print
print '*************************************************************'
print

lnprob_median = np.median(lnprb_T[:, :])
fig = plt.plot(lnprb_T[:, :], '-', alpha=0.5)
plt.axhline(lnprob_median)
plt.savefig(dir_output + 'LNprob_chain.png', bbox_inches='tight')
plt.close()
print 'LNprob median = ', lnprob_median

print
print 'Tref: ', mc.Tref
print
print '*************************************************************'
print

# Plotting out the chains
for ii in xrange(0, mc.ndim):
    print mc.pam_names[ii], chain_med[ii, 0], ' +\sigma ', chain_med[ii, 1], ' -\sigma ', chain_med[ii, 2]
    file_name = dir_output + 'chain_' + repr(ii) + '_' + mc.pam_names[ii] + '.png'
    plt.plot(chain_T[:, :, ii], '-', alpha=0.5)
    plt.axvline(mc.nburnin, c='r')
    plt.axhline(chain_med[ii, 0], c='k')
    plt.savefig(file_name, bbox_inches='tight')
    plt.close()

print
print '*************************************************************'
print

if args.t != 'False':
    for nd in xrange(0, mc.ndim):  # (0,ndim):
        out_absc = np.arange(0, mc.nburnin, 1)
        out_lines = np.zeros(mc.nburnin)
        for ii in xrange(20, mc.nburnin):
            out_lines[ii] = GelmanRubin(chain_T[:ii, :, nd])

        plt.ylim(0.95, 2.3)
        plt.plot(out_absc[20:], out_lines[20:], '-', color='k')
        plt.axhline(1.01)
        plt.savefig(dir_output + 'GRtrace_pam_' + repr(nd) + '.png', bbox_inches='tight')
        plt.close()

    print
    print '*************************************************************'
    print

x0 = 1. / 150

M_star1_rand = np.random.normal(M_star1, M_star1_err, n_kept)

if 'kepler' in mc.model_list:
    for planet_name in mc.pcv.name_ref:

        print planet_name

        if mc.pcv.n_orbpams[planet_name] == 5:

            sample_plan = np.zeros([n_kept, 8])
            sample_plan[:, 0], sample_plan[:, 1], sample_plan[:, 2], sample_plan[:, 3], sample_plan[:, 4] = \
                mc.pcv.convert_params(flatchain[:, mc.variable_list[planet_name]['logP']],
                                      flatchain[:, mc.variable_list[planet_name]['logK']],
                                      flatchain[:, mc.variable_list[planet_name]['phase']],
                                      flatchain[:, mc.variable_list[planet_name]['esino']],
                                      flatchain[:, mc.variable_list[planet_name]['ecoso']])

            # Time of periastron
            sample_plan[:, 5] = mc.Tref + (-sample_plan[:, 2] + sample_plan[:, 4]) / 360.00 * sample_plan[:, 0]

            for ii in xrange(0, n_kept):
                # Planet mass
                sample_plan[ii, 6] = M_ratio * scipy.optimize.fsolve(get_mass, x0, args=(
                    M_star1_rand[ii], sample_plan[ii, 0], sample_plan[ii, 1], sample_plan[ii, 3]))
                # semi-major axis
                sample_plan[ii, 7] = np.power(
                    (Mu_sun * np.power(sample_plan[ii, 0] * seconds_in_day / (2 * np.pi), 2) / (AU_km ** 3.0)) *
                    M_star1_rand[ii], 1.00 / 3.00)

            sample_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                        zip(*np.percentile(sample_plan[:, :], [15.865, 50, 84.135], axis=0))))

            print 'Period = ', sample_med[0, 0], ' +\sigma ', sample_med[0, 1], ' -\sigma ', sample_med[0, 2]
            print 'K      = ', sample_med[1, 0], ' +\sigma ', sample_med[1, 1], ' -\sigma ', sample_med[1, 2]
            print 'phase  = ', sample_med[2, 0], ' +\sigma ', sample_med[2, 1], ' -\sigma ', sample_med[2, 2]
            print 'e      = ', sample_med[3, 0], ' +\sigma ', sample_med[3, 1], ' -\sigma ', sample_med[3, 2]
            print 'o      = ', sample_med[4, 0], ' +\sigma ', sample_med[4, 1], ' -\sigma ', sample_med[4, 2]
            print 'Tperi  = ', sample_med[5, 0], ' +\sigma ', sample_med[5, 1], ' -\sigma ', sample_med[5, 2]
            print 'Mass_J = ', sample_med[6, 0], ' +\sigma ', sample_med[6, 1], ' -\sigma ', sample_med[6, 2]
            print 'Mass_E = ', sample_med[6, 0]*EJ_ratio, \
                ' +\sigma ', sample_med[6, 1]*EJ_ratio, ' -\sigma ', sample_med[6, 2]*EJ_ratio
            print 'a      = ', sample_med[7, 0], ' +\sigma ', sample_med[7, 1], ' -\sigma ', sample_med[7, 2]

            sel_list = [0, 1, 2, 3, 4, 6]
            sel_label = ['P', 'K', 'phase', 'e', 'omega', 'Mass_J']

        else:

            sample_plan = np.zeros([n_kept, 8])
            sample_plan[:, 0], sample_plan[:, 1], sample_plan[:, 2], sample_plan[:, 3], sample_plan[:, 4] = \
                mc.pcv.convert_params(flatchain[:, mc.variable_list[planet_name]['logP']],
                                      flatchain[:, mc.variable_list[planet_name]['logK']],
                                      flatchain[:, mc.variable_list[planet_name]['phase']],
                                      np.zeros(n_kept), np.zeros(n_kept))

            for ii in xrange(0, n_kept):
                # Planet mass
                sample_plan[ii, 6] = M_ratio * \
                                     scipy.optimize.fsolve(get_mass, x0,
                                                           args=(
                                                                 M_star1_rand[ii], sample_plan[ii, 0],
                                                                 sample_plan[ii, 1], sample_plan[ii, 3]))
                # semi-major axis
                sample_plan[ii, 7] = np.power(
                    (Mu_sun * np.power(sample_plan[ii, 0] * seconds_in_day / (2 * np.pi), 2) / (AU_km ** 3.0)) *
                    M_star1_rand[ii], 1.00 / 3.00)

            sample_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                        zip(*np.percentile(sample_plan[:, :], [15.865, 50, 84.135], axis=0))))

            print 'Period = ', sample_med[0, 0], ' +\sigma ', sample_med[0, 1], ' -\sigma ', sample_med[0, 2]
            print 'K      = ', sample_med[1, 0], ' +\sigma ', sample_med[1, 1], ' -\sigma ', sample_med[1, 2]
            print 'phase  = ', sample_med[2, 0], ' +\sigma ', sample_med[2, 1], ' -\sigma ', sample_med[2, 2]
            print 'Mass_J = ', sample_med[6, 0], ' +\sigma ', sample_med[6, 1], ' -\sigma ', sample_med[6, 2]
            print 'Mass_E = ', sample_med[6, 0]*EJ_ratio, \
                ' +\sigma ', sample_med[6, 1]*EJ_ratio, ' -\sigma ', sample_med[6, 2]*EJ_ratio
            print 'a      = ', sample_med[7, 0], ' +\sigma ', sample_med[7, 1], ' -\sigma ', sample_med[7, 2]

            sel_list = [0, 1, 2, 6]
            sel_label = ['P', 'K', 'phase', 'Mass_J']

        fig = corner.corner(sample_plan[:, sel_list], labels=sel_label, truths=sample_med[sel_list, 0])
        fig.savefig(dir_output + planet_name + "_corners.pdf", bbox_inches='tight')
        plt.close()

        print
        print '-----------------------------'
        print

        #print ' Makes the phase plot'
        #mc(chain_med[0,:])
        #for dataset in mc.dataset_list:

        if args.v != 'False':

            list_labels = ['P', 'K', 'e', 'm']

            n_int = 4

            output_plan = np.zeros([n_kept, n_int], dtype=np.double)
            output_plan[:, 0] = sample_plan[:, 0]
            output_plan[:, 1] = sample_plan[:, 1]
            output_plan[:, 2] = sample_plan[:, 4]
            output_plan[:, 3] = sample_plan[:, 6]
            plot_truths = np.percentile(output_plan[:, :], [15.865, 50, 84.135], axis=0)

            veusz_dir = dir_output + '/Veuz_plot/'
            if not os.path.exists(veusz_dir):
                os.makedirs(veusz_dir)

            n_bins = 60 + 1

            h5f = h5py.File(veusz_dir + planet_name + '_hist1d.hdf5', "w")
            data_grp = h5f.create_group("hist1d")

            data_lim = np.zeros([n_int, 2], dtype=np.double)
            data_edg = np.zeros([n_int, n_bins], dtype=np.double)
            for ii in xrange(0, n_int):
                data_lim[ii, :] = [np.amin(output_plan[:, ii]), np.amax(output_plan[:, ii])]
                data_edg[ii, :] = np.linspace(data_lim[ii, 0], data_lim[ii, 1], n_bins)

            for ii in xrange(0, n_int):
                for jj in xrange(ii, n_int):
                    x_data = output_plan[:, ii]
                    y_data = output_plan[:, jj]
                    x_edges = data_edg[ii, :]
                    y_edges = data_edg[jj, :]

                    print 'X_edges', x_edges
                    print 'Y_edges', y_edges
                    if ii != jj:
                        hist2d = np.histogram2d(x_data, y_data, bins=[x_edges, y_edges], normed=True)
                        hist1d_y = np.histogram(y_data, bins=y_edges, normed=True)

                        Hflat = hist2d[0].flatten()
                        inds = np.argsort(Hflat)[::-1]
                        Hflat = Hflat[inds]
                        sm = np.cumsum(Hflat)
                        sm /= sm[-1]

                        x_edges_1d = (x_edges[1:] + x_edges[:-1])/2
                        y_edges_1d = (y_edges[1:] + y_edges[:-1])/2
                        h2d_out = np.zeros([n_bins, n_bins])
                        h2d_out[0, 1:] = x_edges_1d
                        h2d_out[1:, 0] = y_edges_1d
                        h2d_out[1:, 1:] = hist2d[0].T *1. / np.amax(hist2d[0])

                        h2d_list =  h2d_out.tolist()
                        h2d_list[0][0] = ''
                        csvfile = veusz_dir + planet_name + '_hist2d_' + list_labels[ii] + '_' + list_labels[jj] + '.csv'
                        with open(csvfile, "w") as output:
                            writer = csv.writer(output, lineterminator='\n')
                            writer.writerows(h2d_list)
                    else:
                        hist1d = np.histogram(x_data, bins=x_edges)
                        hist1d_norm = hist1d[0]*1. / n_kept
                        x_edges_1d = (x_edges[1:]+ x_edges[:-1])/2
                        data_grp.create_dataset(list_labels[ii]+'_x', data=x_edges_1d, compression="gzip")
                        data_grp.create_dataset(list_labels[ii]+'_y', data=hist1d_norm, compression="gzip")

            #plot lower-upper limits for the mass

            pams_limits = np.zeros([mc.pcv.n_orbpams[planet_name], 2])
            for ii in xrange(0, mc.pcv.n_orbpams[planet_name]):
                pams_limits[ii, :] = np.percentile(sample_plan[:, ii], [0.135, 99.865])
                print ii, pams_limits[ii, :]

            random_n = 2000
            ii_kep = np.random.randint(low=0, high=n_kept, size=random_n)

            x_kep = np.arange(5900, 7900, 2)
            y_kep = np.zeros([np.size(x_kep), 2])

            x_pha = np.arange(-1.00, 2.00, 0.005,dtype=np.double)
            y_pha = np.zeros([np.size(x_pha), 2])

            y_flg = np.ones(random_n, dtype=bool)

            y_kep_tmp = np.zeros(np.size(x_kep))
            y_pha_tmp = np.zeros(np.size(x_pha))

            if pams_limits[3, 0] < 0.02 : pams_limits[3, 0] = 0.00
            for ii in xrange(0, 4):
                y_flg = y_flg & (sample_plan[ii_kep, ii] >= pams_limits[ii,0]) & (sample_plan[ii_kep, ii] <= pams_limits[ii, 1])

            y_kep[:, 0] = kp.kepler_RV_T0P(x_kep - mc.Tref, sample_med[2, 0], sample_med[0, 0], sample_med[1, 0],
                                           sample_med[3, 0], sample_med[4, 0])
            y_kep[:, 1] = y_kep[:,0]
            # value initialization
            y_pha[:, 0] = kp.kepler_RV_T0P(x_pha*sample_med[0, 0], sample_med[2, 0], sample_med[0, 0], sample_med[1, 0],
                                           sample_med[3, 0], sample_med[4, 0])
            y_pha[:, 1] = y_pha[:, 0]

            fig1 = plt.plot(x_kep, y_kep[:, 0], c='g')
            fig2 = plt.plot(x_pha, y_pha[:, 0], c='g')

            print 'Accepted randomizations: ', np.sum(y_flg)
            for ii in xrange(0, random_n):
                ir = ii_kep[ii]
                if y_flg[ii]:
                    y_kep_tmp = kp.kepler_RV_T0P(x_kep - mc.Tref, sample_plan[ir, 2], sample_plan[ir, 0], sample_plan[ir, 1], sample_plan[ir, 3], sample_plan[ir, 4])
                    y_pha_tmp = kp.kepler_RV_T0P(x_pha*sample_plan[ir, 0], sample_plan[ir, 2], sample_plan[ir, 0], sample_plan[ir, 1], sample_plan[ir, 3], sample_plan[ir, 4])

                    y_kep[:, 0] = np.minimum(y_kep[:, 0],y_kep_tmp)
                    y_kep[:, 1] = np.maximum(y_kep[:, 1],y_kep_tmp)

                    y_pha[:, 0] = np.minimum(y_pha[:, 0],y_pha_tmp)
                    y_pha[:, 1] = np.maximum(y_pha[:, 1],y_pha_tmp)

            fig1 = plt.plot(x_kep, y_kep[:, 0], c='r')
            fig1 = plt.plot(x_kep, y_kep[:, 1], c='b')
            fig2 = plt.plot(x_pha, y_pha[:, 0], c='r')
            fig2 = plt.plot(x_pha, y_pha[:, 1], c='b')

            # h5f = h5py.File('output/'+planet_name+"planet"+`pp`+'_kep.hdf5', "w")
            # data_grp = h5f.create_group("data")
            # data_grp.create_dataset("x_kep", data=x_kep, compression="gzip")
            # data_grp.create_dataset("y_kep", data=y_kep, compression="gzip")
            # h5f.close()
            # h5f = h5py.File('output/'+planet_name+"planet"+`pp`+'_pha.hdf5', "w")
            # data_grp = h5f.create_group("data")
            # data_grp.create_dataset("x_pha", data=x_kep, compression="gzip")
            # data_grp.create_dataset("y_pha", data=y_kep, compression="gzip")
            # h5f.close()

            fileout = open(veusz_dir + planet_name + '_kep.dat','w')
            fileout.write('descriptor x_kep y_kep,+- \n')
            for ii in xrange(0,np.size(x_kep)):
                fileout.write('{0:14f} {1:14f} {2:14f}  \n'.format(\
                    x_kep[ii], (y_kep[ii, 1]+y_kep[ii, 0])/2, (y_kep[ii, 1]-y_kep[ii, 0])/2))
            fileout.close()

            fileout = open(veusz_dir + planet_name + '_pha.dat','w')
            fileout.write('descriptor x_pha y_pha,+- \n')
            for ii in xrange(0,np.size(x_pha)):
                fileout.write('{0:14f} {1:14f} {2:14f}  \n'.format(\
                    x_pha[ii], (y_pha[ii, 1]+y_pha[ii, 0])/2, (y_pha[ii, 1]-y_pha[ii, 0])/2))
            fileout.close()

