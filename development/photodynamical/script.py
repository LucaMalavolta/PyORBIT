import numpy as np
import os
import sys
# import glob
# import time as timer
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table
from pytransit import RoadRunnerModel
import matplotlib.pyplot as plt

from pytrades import constants as cst
from pytrades import ancillary as anc
from pytrades import pytrades

set_unit_base = anc.set_unit_base
anc.set_rcParams()

here_folder = os.path.abspath(".")

photometry_folder = os.path.join(here_folder, "photometry")
# radial_velocities_folder = os.path.join(here_folder, "radial_velocities")
# transit_times_folder = os.path.join(here_folder, "transit_times")

out_folder = os.path.join(here_folder, "output")
os.makedirs(out_folder, exist_ok=True)

map_out_folder = os.path.join(out_folder, "01_map")
os.makedirs(map_out_folder, exist_ok=True)

# tess_s13_file = os.path.join(
#     photometry_folder, 
#     "hlsp_qlp_tess_ffi_s0013-0000000254113311_tess_v01_llc.txt"
# )
# time_s13, flux_sap_s13, flux_sap_err_s13 = np.genfromtxt(tess_s13_file, delimiter=",", unpack=True)

tess_s13_file = os.path.join(
    photometry_folder, 
    "hlsp_qlp_tess_ffi_s0013-0000000254113311_tess_v01_llc.fits"
)
with fits.open(tess_s13_file) as hdul:
    print(hdul.info())
    s13_header = hdul[1].header
    tess_s13 = Table(hdul[1].data)

# print(s13_header)
texp = s13_header["TIMEDEL"]
print("exposure time = {:.6f}d == {:.0f}".format(texp, texp*cst.day2sec))
n_over = int(texp*cst.day2sec / 120.0)+1

# We are associating the following variables to columns within the tess_s13 dataset
# As always, we are only interested in the Time, Flux and Fluxx error columns

time_s13, flux_sap_s13, flux_sap_err_s13 = tess_s13["TIME"], tess_s13["KSPSAP_FLUX"], tess_s13["KSPSAP_FLUX_ERR"]
quality = tess_s13["QUALITY"]
ok = np.logical_and(
    np.isnan(flux_sap_s13) == False,
    quality == 0
)
time_s13, flux_sap_s13, flux_sap_err_s13 = time_s13[ok], flux_sap_s13[ok], flux_sap_err_s13[ok]

fileout = open('photometry/tess_v01_llc_PyORBIT.dat', 'w')
for b, v, e in zip(time_s13, flux_sap_s13, flux_sap_err_s13 ):
    fileout.write('{0:16.8f} {1:12.8f} {2:12.8f} 0 -1 -1 \n'.format(b,v,e))
fileout.close()

# a transit of planet b
file_cheops = os.path.join(
    photometry_folder,
    "CHEOPS-PIPE_TOI-1130b_v04_CH_PR100031_TG042201_V0200_detrended.dat"
)
time_cheops_b, flux_cheops_b, flux_cheops_err_b = np.genfromtxt(file_cheops, usecols=(0,1,2), unpack=True)

# a transit of planet c
file_cheops = os.path.join(
    photometry_folder,
    "CHEOPS-PIPE_TOI-1130c_v03_CH_PR100015_TG018201_V0200_detrended.dat"
)
time_cheops_c, flux_cheops_c, flux_cheops_err_c = np.genfromtxt(file_cheops, usecols=(0,1,2), unpack=True)

# the transit of planet b and c
file_cheops = os.path.join(
    photometry_folder,
    "CHEOPS-PIPE_TOI-1130c_v06_CH_PR120053_TG004701_V0200_normalized.dat"
)
time_cheops_bc, flux_cheops_bc, flux_cheops_err_bc = np.genfromtxt(file_cheops, usecols=(0,1,2), unpack=True)

time_all = np.concatenate((time_s13, time_cheops_b, time_cheops_c, time_cheops_bc))
#time_all is just an array made up by "joinig" (concatenating) the "Time" arrays of each of the four data sets.

t_epoch = 1657.0 # I BELIEVE the t_epoch is the time at which we define our orbital parameters
t_start = np.min(time_all) - 10.0 # The start time of the simulation is just a little bit before the first (smallest) timestamp within the time_all array
# t_start = t_epoch
# t_end = Time("2023-08-30T00:00:00", format="isot", scale="tdb").jd - cst.btjd
t_end = np.max(time_all) + 1.0 # Finishing time is a little bit after the last (highest) timestamp in time_all
t_int = t_end - t_start

print("t_epoch = {}".format(t_epoch))
print("t_start = {}".format(t_start))
print("t_end   = {}".format(t_end))
print("total integration time: {} days".format(t_int))

body_names = ["star", "b", "c"]
n_body = len(
    body_names
)  # number of bodies (NOT PLANETS) in the system, that is star + planets

sim = pytrades.PhotoTRADES(
    n_body,
    t_epoch,
    t_start,
    t_int,
    duration_check=1,
    encounter_check=True,
    do_hill_check=False,
    amd_hill_check=False,
    rv_res_gls=False,
)



cheops = {}
anc_cheops = {}
cheops[0] =  pytrades.set_photometry_portion(
    time_cheops_b, flux_cheops_b, flux_cheops_err_b,
            n_oversample=1,
            t_exp_d=60.0 / 86400.)
cheops[1] =  pytrades.set_photometry_portion(
    time_cheops_c, flux_cheops_c, flux_cheops_err_c,
            n_oversample=1,
            t_exp_d=60.0 / 86400.)
cheops[2] =  pytrades.set_photometry_portion(
    time_cheops_bc, flux_cheops_bc, flux_cheops_err_bc,
            n_oversample=1,
            t_exp_d=60.0 / 86400.)
anc_cheops[0] = None
anc_cheops[1] = None
anc_cheops[2] = None

tess = {}
anc_tess = {}

tess[0] = pytrades.set_photometry_portion(
            time_s13, flux_sap_s13, flux_sap_err_s13,
            n_oversample=1,
            t_exp_d=60.0 / 86400.)
anc_tess[0] = None
#print(tess[0])


chiron_file = os.path.join('photoTRADES_TOI-1130/radial_velocities/', "Huang2020_obsRV.dat")
t_rv_chiron, rv_chiron, erv_chiron = np.genfromtxt(chiron_file, unpack=True)
chiron = {"time": t_rv_chiron, "rv": rv_chiron, "rv_err": erv_chiron}


sim.add_photometry("cheops", cheops, anc_cheops)
#sim.add_photometry("tess", tess, anc_tess)
sim.add_radial_velocity("CHIRON", chiron, 1)  # chiron will have the source id 1


sim.set_radial_velocity_sorting()  # sort the radial velocities\

sim.update_n_data()

# These orbital parameters are defined at one specific instant. Is it at t_epoch?

M_msun = np.array(
    [0.745059, 19.833346 * cst.Mears, 335.603435 * cst.Mears]
)  # Masses in Solar unit
R_rsun = np.array(
    [0.697470, 3.657000 * cst.Rears, 12.983016 * cst.Rears]
)  # Radii in Solar unit
P_day = np.array([0.0, 4.074554, 8.350190])  # Periods in days
ecc_val = np.array([0.0, 0.0521624, 0.039773])  # eccentricities
argp_deg = np.array([0.0, 141.11112, 182.502357])  # argument of pericenters in degrees
mA_deg = np.array([0.0, 159.696701, 233.068994])  # mean anonalies in degrees
inc_deg = np.array([0.0, 87.494901, 87.613475])  # inclinations in degrees
lN_deg = np.array([0.0, 180.0, 179.993043])  # longitude of ascending nodes in degrees


sim.add_default_physical_parameters(
    M_msun, R_rsun, P_day, ecc_val, argp_deg, mA_deg, inc_deg, lN_deg
)

(
    time_steps,
    orbits,
    transits,
    durations,
    lambda_rm,
    kep_elem,
    body_flag,
    rv_sim,
    stable,
) = sim.orbital_parameters_to_transits(
    M_msun, R_rsun, P_day, ecc_val, argp_deg, mA_deg, inc_deg, lN_deg
)
phot_sim_v0 = sim.get_simulate_flux(
    R_rsun, ld_quads, transits, durations, body_flag, kep_elem,
    time_key="time"
)

phot_sim_v0, rv_sim, sim_transits = sim.full_photodyn(
    M_msun, R_rsun, P_day, ecc_val, argp_deg, mA_deg, inc_deg, lN_deg, ld_quads
)