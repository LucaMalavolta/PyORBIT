from classes.model_container import ModelContainerMultiNest
from classes.input_parser import yaml_parser
import numpy as np
import os
import sys
import argparse
import json
import threading, subprocess
import pymultinest

def show(filepath):
    """ open the output (pdf) file for the user """
    if os.name == 'mac': subprocess.call(('open', filepath))
    elif os.name == 'nt': os.startfile(filepath)

parser = argparse.ArgumentParser(prog='PyORBIT_V3_MultiNest.py', description='PyDE+emcee runner')
parser.add_argument('config_file', type=str, nargs=1, help='config file')


args = parser.parse_args()

file_conf = args.config_file[0]

mc = ModelContainerMultiNest()

yaml_parser(file_conf, mc)


dir_output = './' + mc.planet_name + '/'
if not os.path.exists(dir_output):
    os.makedirs(dir_output)
if not os.path.exists(dir_output+'chains'):
    os.makedirs(dir_output+'chains')

mc.create_bounds()


if bool(mc.pcv.dynamical):
        mc.pcv.prepare_dynamical(mc)

print 'Dimensions = ', mc.ndim
print '   '
print 'Variable list:', mc.variable_list
print
print 'Variable bounds:', mc.bounds
print


'''
    On Linux system (BASH):
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib
    export LD_PRELOAD=/usr/lib/openmpi/lib/libmpi.so:$LD_PRELOAD
    export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgfortran.so.3
    mpirun -np 4 python run_PyPolyChord.py

    on Mac:
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib
    export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgfortran.so.3
    export LD_PRELOAD=/opt/local/lib/openmpi//lib/libmpi.so:$LD_PRELOAD
    mpirun -np 4 python run_PyPolyChord.py

'''


# number of dimensions our problem has
n_params = mc.ndim
print '--->' , mc.pam_names

mc.results_resumen(np.ones(n_params))
parameters = mc.pam_names


progress = pymultinest.ProgressPlotter(n_params=n_params, outputfiles_basename=dir_output + 'chains/2-');
progress.start()
threading.Timer(2, show, [dir_output + "chains/2-phys_live.points.pdf"]).start()  # delayed opening
# run MultiNest
pymultinest.run(mc.multinest_call, mc.multinest_priors, n_params, resume=True,
                sampling_efficiency='parameter',
                evidence_tolerance=0.2, max_iter=0,
                verbose=True, n_live_points=1000, outputfiles_basename=dir_output + 'chains/2-')
progress.stop()

# lets analyse the results
a = pymultinest.Analyzer(n_params=n_params, outputfiles_basename=dir_output + 'chains/2-')
s = a.get_stats()


#print s
print a.get_best_fit()
print s['marginals'][0]['median'], s['marginals'][0]['3sigma']

print()
print("-" * 30, 'ANALYSIS', "-" * 30)
print("Global Evidence:\n\t%.15e +- %.15e" % (
s['nested sampling global log-evidence'], s['nested sampling global log-evidence error']))

import matplotlib.pyplot as plt

plt.clf()

# run MultiNest
# pymultinest.run(mc.pymultinest_call, mc.pymultinest_priors, mc.ndim, outputfiles_basename=dir_output, resume = False, verbose = True)
# json.dump(parameters, open(dir_output+'params.json', 'w')) # save parameter names

# Here we will plot all the marginals and whatnot, just to show off
# You may configure the format of the output here, or in matplotlibrc
# All pymultinest does is filling in the data of the plot.

# Copy and edit this file, and play with it.

p = pymultinest.PlotMarginalModes(a)
plt.figure(figsize=(5 * n_params, 5 * n_params))
# plt.subplots_adjust(wspace=0, hspace=0)
for i in range(n_params):
    plt.subplot(n_params, n_params, n_params * i + i + 1)
    p.plot_marginal(i, with_ellipses=True, with_points=True, grid_points=50)
    plt.ylabel("Probability")
    plt.xlabel(parameters[i])

    for j in range(i):
        plt.subplot(n_params, n_params, n_params * j + i + 1)
        # plt.subplots_adjust(left=0, bottom=0, right=0, top=0, wspace=0, hspace=0)
        p.plot_conditional(i, j, with_ellipses=False, with_points=True, grid_points=30)
        plt.xlabel(parameters[i])
        plt.ylabel(parameters[j])

plt.savefig(dir_output + "chains/marginals_multinest.pdf")  # , bbox_inches='tight')
show(dir_output + "chains/marginals_multinest.pdf")

for i in range(n_params):
    outfile = '%s-mode-marginal-%d.pdf' % (a.outputfiles_basename, i)
    p.plot_modes_marginal(i, with_ellipses=True, with_points=True)
    plt.ylabel("Probability")
    plt.xlabel(parameters[i])
    plt.savefig(outfile, format='pdf', bbox_inches='tight')
    plt.close()

    outfile = '%s-mode-marginal-cumulative-%d.pdf' % (a.outputfiles_basename, i)
    p.plot_modes_marginal(i, cumulative=True, with_ellipses=True, with_points=True)
    plt.ylabel("Cumulative probability")
    plt.xlabel(parameters[i])
    plt.savefig(outfile, format='pdf', bbox_inches='tight')
    plt.close()

print("Take a look at the pdf files in chains/")
