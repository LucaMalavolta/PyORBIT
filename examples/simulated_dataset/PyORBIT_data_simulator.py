import numpy as np

import sys
import os
import matplotlib.pyplot as plt

sys.path.insert(0, '../pyorbit/')
import pyorbit.classes.kepler_exo as kp
import pyorbit.classes.constants as constants


"""
Initial realization of the parameters

np.random.seed(12345678)

planet_b = {
    'P': np.random.uniform(6,9),
    'K': np.random.uniform(10,20),
    'f': np.random.uniform(0,np.pi),
    'e': np.random.uniform(0,0.2),
    'o': np.random.uniform(0,np.pi),
    'i': np.random.uniform(85,90),
}

planet_c = {
    'P': np.random.uniform(12,18),
    'K': np.random.uniform(5,15),
    'f': np.random.uniform(0,np.pi),
    'e': np.random.uniform(0,0.3),
    'o': np.random.uniform(0,np.pi),
    'i': np.random.uniform(85,90),
}

instrument = {
    'RV_precision': 1.000,
    'RV_offset1': np.random.uniform(4000.0, 5000.0),
    'RV_offset2': np.random.uniform(4000.0, 5000.0),
    'T0_precision': 0.001
}

"""

bjd_syn = np.arange(6000, 6100, 1, dtype=np.double)
bjd_syn -= (bjd_syn-6025.0)*(4./1440)

bjd_obs = np.random.normal(np.arange(6000, 6100, 1, dtype=np.double), 0.2)
Tref = 6000.0000

planet_b  = {
    'P': 6.737412701474708,
    'K': 15.964286053214767,
    'Tc': 6002.473649024,
    'e': 0.07578202101359831,
    'i': 86.19154785304445,
    'o': 0.07681624494196555}

planet_c  = {
    'P': 14.327606002800504,
    'K': 11.805410430293563,
    'Tc': 6010.0937947242,
    'e': 0.22822106564776393,
    'i': 85.317208074894,
    'o': 0.6966217861505585}

instrument  = {
    'RV_precision': 1.0,
    'RV_offset1': 4779.443746523752,
    'RV_offset2': 4721.741371149122,
    'T0_precision': 0.001
}


planet_b['f'] = kp.kepler_Tc2phase_Tref(planet_b['P'],
                                        planet_b['Tc']-Tref,
                                        planet_b['e'],
                                        planet_b['o'])
# 5.396785202049886

planet_c['f'] = kp.kepler_Tc2phase_Tref(planet_c['P'],
                                        planet_c['Tc']-Tref,
                                        planet_c['e'],
                                        planet_c['o'])
# 3.114055661251522

# print(kp.get_planet_mass(planet_b['P'], planet_b['K'], planet_b['e'], 1.0000, approximation_limit=1.) * constants.Msear)
# print(kp.get_planet_mass(planet_c['P'], planet_c['K'], planet_c['e'], 1.0000, approximation_limit=1.) * constants.Msear)
# [47.02042486]
# [43.65936615]

bjd0 = bjd_obs - Tref

y_pla = instrument['RV_offset1'] + \
        kp.kepler_RV_T0P(bjd0,
                         planet_b['f'],
                         planet_b['P'],
                         planet_b['K'],
                         planet_b['e'],
                         planet_b['o']) + \
        kp.kepler_RV_T0P(bjd0,
                         planet_c['f'],
                         planet_c['P'],
                         planet_c['K'],
                         planet_c['e'],
                         planet_c['o'])

mod_pl = np.random.normal(y_pla, instrument['RV_precision'])

Tcent_b =  np.random.normal(
    np.arange(0,5)*planet_b['P'] + kp.kepler_phase2Tc_Tref(planet_b['P'], planet_b['f'], planet_b['e'], planet_b['o']) + Tref,
    instrument['T0_precision'])

Tcent_c =  np.random.normal(
    np.arange(0,10)*planet_c['P'] + kp.kepler_phase2Tc_Tref(planet_c['P'], planet_c['f'], planet_c['e'], planet_c['o']) + Tref,
    instrument['T0_precision'])


fileout = open('simulated_2p_RV.dat', 'w')
for b, r in zip(bjd_obs, mod_pl):
    fileout.write('{0:f}  {1:.2f}  {2:.2f}  {3:d}  {4:d}  {5:d}\n'.format(b, r, instrument['RV_precision'], 0, 0, -1))
fileout.close()

fileout = open('simulated_2p_Tc_b.dat', 'w')
for i_Tc, v_Tc in enumerate(Tcent_b):
    fileout.write('{0:d}  {1:.4f}  {2:.4f}  {3:d}\n'.format(i_Tc, v_Tc, instrument['T0_precision'], 0))
fileout.close()

fileout = open('simulated_2p_Tc_c.dat', 'w')
for i_Tc, v_Tc in enumerate(Tcent_c):
    fileout.write('{0:d}  {1:.4f}  {2:.4f}  {3:d}\n'.format(i_Tc, v_Tc, instrument['T0_precision'], 0))
fileout.close()




###################################################################
