from __future__ import print_function, division
import numpy as np
#%matplotlib notebook
import matplotlib.pyplot as plt
import collections
import matplotlib.gridspec as gridspec # GRIDSPEC !

dir_base = './'
dir_mods = 'K2-141_example01_LC_1p'
dir_plot = '/emcee_plot/MAP_model_files/'

filename = 'K2-141_example01_LC_1p_transit'

plot_properties = {
    'x_limits': [-0.165, 0.165],
    'y_limits': None,
    'y_residual_limits': None
}

planet_list = ['b']

def plots_in_grid():
    # Partially taken from here:
    # http://www.sc.eso.org/~bdias/pycoffee/codes/20160407/gridspec_demo.html

    gs = gridspec.GridSpec(2,1, height_ratios=[3.0,1.0])
    # Also make sure the margins and spacing are apropriate
    gs.update(left=0.2, right=0.95, bottom=0.08, top=0.93, wspace=0.02, hspace=0.03)

    ax_0 = plt.subplot(gs[0])
    ax_1 = plt.subplot(gs[1])

    # Adding minor ticks only to x axis
    from matplotlib.ticker import AutoMinorLocator
    minorLocator = AutoMinorLocator()
    ax_0.xaxis.set_minor_locator(minorLocator)
    ax_1.xaxis.set_minor_locator(minorLocator)

    # Disabling the offset on top of the plot
    ax_0.ticklabel_format(useOffset=False)
    ax_1.ticklabel_format(useOffset=False)
    return ax_0, ax_1


LC_full = np.genfromtxt(dir_base + dir_mods + dir_plot + 'LCdata_full.dat', skip_header=1)
LC_modb = np.genfromtxt(dir_base + dir_mods + dir_plot + 'LCdata_lc_model_b.dat', skip_header=1)
LC_modb_full = np.genfromtxt(dir_base + dir_mods + dir_plot + 'LCdata_lc_model_b_full.dat', skip_header=1)

t1_hours = LC_modb[:,1]*24.
transit_phase = LC_modb[:,1]


diff_array = LC_modb_full[1:,1]-LC_modb_full[:-1,1]
ind_start = np.where(diff_array<0)[0][0]+1
ind_stop = np.where(diff_array<0)[0][1]+1
#index_sort = np.argsort(LC_modb_full[:,1])

fig = plt.figure(1, figsize=(8,8))

ax_0, ax_1 = plots_in_grid()

ax_0.errorbar(transit_phase, LC_modb[:,8]+1., yerr=LC_modb[:,9], color='black', markersize=0, alpha=0.25, fmt='o', zorder=0)
ax_0.scatter(transit_phase, LC_modb[:,8]+1. , c='C0', s=4, zorder=1, label = 'First transit')

ax_0.plot(LC_modb_full[ind_start:ind_stop,1],LC_modb_full[ind_start:ind_stop,3]+1., color='k', linestyle='-', zorder=2)

ax_1.errorbar(transit_phase, LC_modb[:,10], yerr=LC_modb[:,11], color='black', markersize=0, alpha=0.25, fmt='o', zorder=1)

ax_1.scatter(transit_phase, LC_modb[:,10], c='C0', s=4, zorder=2)
ax_1.axhline(0.000, c='k', zorder=3)

#ax_0.legend( loc='lower left', fontsize='medium')

if plot_properties.get('x_limits', False):
    ax_0.set_xlim(plot_properties['x_limits'][0], plot_properties['x_limits'][1])
    ax_1.set_xlim(plot_properties['x_limits'][0], plot_properties['x_limits'][1])
if plot_properties.get('y_limits', False):
    ax_0.set_ylim(plot_properties['y_limits'][0], plot_properties['y_limits'][1])
if plot_properties.get('y_residual_limits', False):
    ax_1.set_ylim(plot_properties['y_residual_limits'][0], plot_properties['y_residual_limits'][1])

#ax_0.tick_params(axis='x', which='minor')
#ax_0.yaxis.set_tick_params(which='minor', right = 'off')
#ax_0.minorticks_on()

#ax_0.set_ylim(0.9898,1.0022)
ax_0.set_ylabel('Normalized flux')
ax_1.set_xlabel('Time from mid-transit [h]')
ax_1.set_ylabel('Residuals')
plt.savefig(filename+'.pdf', dpi=300, bbox_inches='tight')
