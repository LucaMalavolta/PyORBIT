import numpy as np

import sys
import os
import matplotlib.pyplot as plt

sys.path.insert(0, '../pyorbit/classes/')
import kepler_exo as kp


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

planet_b  = {
    'P': 6.737412701474708,
    'K': 15.964286053214767,
    'f': 1.1271771440198723,
    'e': 0.07578202101359831,
    'i': 86.19154785304445,
    'o': 0.07681624494196555}

planet_c  = {
    'P': 14.327606002800504,
    'K': 11.805410430293563,
    'f': 2.63686698774923,
    'e': 0.22822106564776393,
    'i': 85.317208074894,
    'o': 0.6966217861505585}

instrument  = {
    'RV_precision': 1.0,
    'RV_offset1': 4779.443746523752,
    'RV_offset2': 4721.741371149122,
    'T0_precision': 0.001
}



bjd_syn = np.arange(6000, 6050, 1, dtype=np.double)
bjd_syn -= (bjd_syn-6025.0)*(4./1440)

bjd_obs = np.random.normal(np.arange(6000, 6050, 1, dtype=np.double), 0.2)
Tref = 6025.0000

photometry = {
    'phot_precision': 0.0001,
    'phot_bjd': np.arange(6000.0, 6070.0, (29.45/60./24.)*10.0) #pseudo-K2 range
}


activity = {
    'Prot': 10.34,
    'Pdec': 55.23,
    'Oamp': 0.34,
    'Hamp_RV1': 9.24,
    'Hamp_RV2': 15.80,
    'Hamp_PH': 0.023,
}

polynomial_trend = {'c1': 0.5454, 'c2': 0.08934}



"""
TestCase01: single planet, no external contraints
TestCase02: single planet, transit times
TestCase03: two planets, transit times of the first one
TestCase04: two planets, transit times of the first one, dynamical integration
TestCase05: photometric curve analyzed with Gaussian Process
TestCase06: RV + photometric light curve analyzed with Gaussian Process
TestCase07: two RV datasets analyzed with Gaussian Process, with priors
"""

def create_testcase01():

    bjd0 = bjd_obs - Tref

    y_pla = kp.kepler_RV_T0P(bjd0,
                             planet_b['f'],
                             planet_b['P'],
                             planet_b['K'],
                             planet_b['e'],
                             planet_b['o']) + instrument['RV_offset1']

    mod_pl = np.random.normal(y_pla, instrument['RV_precision'])

    fileout = open('TestCase01_RV.dat', 'w')
    for b, r in zip(bjd_obs, mod_pl):
        fileout.write('{0:f}  {1:.2f}  {2:.2f}  {3:d}  {4:d}  {5:d}\n'.format(b, r, instrument['RV_precision'], 0, 0, -1))
    fileout.close()


def create_testcase02():

    bjd0 = bjd_obs - Tref

    y_pla = kp.kepler_RV_T0P(bjd0,
                             planet_b['f'],
                             planet_b['P'],
                             planet_b['K'],
                             planet_b['e'],
                             planet_b['o']) + instrument['RV_offset1']

    mod_pl = np.random.normal(y_pla, instrument['RV_precision'])

    Tcent_b =  np.random.normal(
        np.arange(0, 5)*planet_b['P'] + kp.kepler_phase2Tc_Tref(planet_b['P'], planet_b['f'], planet_b['e'], planet_b['o']) + Tref,
        instrument['T0_precision'])

    fileout = open('TestCase02_RV.dat', 'w')
    for b, r in zip(bjd_obs, mod_pl):
        fileout.write('{0:f}  {1:.2f}  {2:.2f}  {3:d}  {4:d}  {5:d}\n'.format(b, r, instrument['RV_precision'], 0, 0, -1))
    fileout.close()

    fileout = open('TestCase02_Tcent_b.dat', 'w')
    for i_Tc, v_Tc in enumerate(Tcent_b):
        fileout.write('{0:d}  {1:.4f}  {2:.4f}  {3:d}\n'.format(i_Tc, v_Tc, instrument['T0_precision'], 0))
    fileout.close()


def create_testcase03():

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

    #Tcent_b = kp.kepler_phase2Tc_Tref(planet_b['P'], planet_b['f'], planet_b['e'], planet_b['o']) + Tref
    Tcent_b =  np.random.normal(
        np.arange(0,5)*planet_b['P'] + kp.kepler_phase2Tc_Tref(planet_b['P'], planet_b['f'], planet_b['e'], planet_b['o']) + Tref,
        instrument['T0_precision'])

    fileout = open('TestCase03_RV.dat', 'w')
    for b, r in zip(bjd_obs, mod_pl):
        fileout.write('{0:f}  {1:.2f}  {2:.2f}  {3:d}  {4:d}  {5:d}\n'.format(b, r, instrument['RV_precision'], 0, 0, -1))
    fileout.close()

    #fileout = open('TestCase01_Tcent_b.dat', 'w')
    #for ii in xrange(0, np.size(Transit_Time)):
    #    fileout.write('{0:d} {1:f} {2:f} {3:d} \n'.format(ii, Tcent_b, 0.01, 0))
    #fileout.close()

    fileout = open('TestCase03_Tcent_b.dat', 'w')
    for i_Tc, v_Tc in enumerate(Tcent_b):
        fileout.write('{0:d}  {1:.4f}  {2:.4f}  {3:d}\n'.format(i_Tc, v_Tc, instrument['T0_precision'], 0))
    fileout.close()


def create_testcase04():

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

    #Tcent_b = kp.kepler_phase2Tc_Tref(planet_b['P'], planet_b['f'], planet_b['e'], planet_b['o']) + Tref
    Tcent_b =  np.random.normal(
        np.arange(0,10)*planet_b['P'] + kp.kepler_phase2Tc_Tref(planet_b['P'], planet_b['f'], planet_b['e'], planet_b['o']) + Tref,
        instrument['T0_precision'])

    Tcent_c =  np.random.normal(
        np.arange(0,10)*planet_c['P'] + kp.kepler_phase2Tc_Tref(planet_c['P'], planet_c['f'], planet_c['e'], planet_c['o']) + Tref,
        instrument['T0_precision'])

    fileout = open('TestCase04_RV.dat', 'w')
    for b, r in zip(bjd_obs, mod_pl):
        fileout.write('{0:f}  {1:.2f}  {2:.2f}  {3:d}  {4:d}  {5:d}\n'.format(b, r, instrument['RV_precision'], 0, 0, -1))
    fileout.close()

    fileout = open('TestCase04_Tcent_b.dat', 'w')
    for i_Tc, v_Tc in enumerate(Tcent_b):
        fileout.write('{0:d}  {1:.4f}  {2:.4f}  {3:d}\n'.format(i_Tc, v_Tc, instrument['T0_precision'], 0))
    fileout.close()

    fileout = open('TestCase04_Tcent_c.dat', 'w')
    for i_Tc, v_Tc in enumerate(Tcent_c):
        fileout.write('{0:d}  {1:.4f}  {2:.4f}  {3:d}\n'.format(i_Tc, v_Tc, instrument['T0_precision'], 0))
    fileout.close()


def create_testcase05():

    import george

    bjd0 = photometry['phot_bjd'] - Tref
    err = bjd0*0 + photometry['phot_precision']

    """ Conversion of the physical parameter to the internally defined parameter
    to be passed to george
    """
    gp_pams = np.zeros(4)
    gp_pams[0] = np.log(activity['Hamp_PH'])*2
    gp_pams[1] = np.log(activity['Pdec'])*2
    gp_pams[2] = 1. / (2*activity['Oamp'] ** 2)
    gp_pams[3] = np.log(activity['Prot'])

    kernel = np.exp(gp_pams[0]) * \
                      george.kernels.ExpSquaredKernel(metric=np.exp(gp_pams[1])) * \
                      george.kernels.ExpSine2Kernel(gamma=gp_pams[2], log_period=gp_pams[3])

    gp = george.GP(kernel)

    gp.compute(bjd0, err)
    prediction = gp.sample(bjd0)
    obs_photometry = np.random.normal(prediction, photometry['phot_precision'])

    fileout = open('TestCase05_photometry.dat', 'w')
    for b, p in zip(photometry['phot_bjd'], obs_photometry):
        fileout.write('{0:14f} {1:14f} {2:14f} {3:5d} {4:5d} {5:5d} \n'.format(
            b, p, photometry['phot_precision'], 0, 0, -1))
    fileout.close()


def create_testcase06():
    import george

    bjd0 = photometry['phot_bjd'] - Tref
    err = bjd0 * 0 + photometry['phot_precision']

    """ Conversion of the physical parameter to the internally defined parameter
    to be passed to george
    """
    gp_pams = np.zeros(4)
    gp_pams[0] = np.log(activity['Hamp_PH']) * 2
    gp_pams[1] = np.log(activity['Pdec']) * 2
    gp_pams[2] = 1. / (2 * activity['Oamp'] ** 2)
    gp_pams[3] = np.log(activity['Prot'])

    kernel = np.exp(gp_pams[0]) * \
             george.kernels.ExpSquaredKernel(metric=np.exp(gp_pams[1])) * \
             george.kernels.ExpSine2Kernel(gamma=gp_pams[2], log_period=gp_pams[3])

    gp = george.GP(kernel)

    gp.compute(bjd0, err)
    prediction = gp.sample(bjd0)
    obs_photometry = np.random.normal(prediction, photometry['phot_precision'])

    fileout = open('TestCase06_photometry.dat', 'w')
    for b, p in zip(photometry['phot_bjd'], obs_photometry):
        fileout.write('{0:14f} {1:14f} {2:14f} {3:5d} {4:5d} {5:5d} \n'.format(
            b, p, photometry['phot_precision'], 0, 0, -1))
    fileout.close()

    bjd0 = bjd_obs - Tref

    err = bjd0 * 0 + instrument['RV_precision']

    gp_pams[0] = np.log(activity['Hamp_RV1']) * 2
    kernel = np.exp(gp_pams[0]) * \
             george.kernels.ExpSquaredKernel(metric=np.exp(gp_pams[1])) * \
             george.kernels.ExpSine2Kernel(gamma=gp_pams[2], log_period=gp_pams[3])

    gp = george.GP(kernel)

    gp.compute(bjd0, err)
    prediction = gp.sample(bjd0)

    y_pla = kp.kepler_RV_T0P(bjd0,
                             planet_b['f'],
                             planet_b['P'],
                             planet_b['K'],
                             planet_b['e'],
                             planet_b['o']) + instrument['RV_offset1']

    mod_pl = np.random.normal(y_pla + prediction, instrument['RV_precision'])

    Tcent_b =  np.random.normal(
        np.arange(0,1)*planet_b['P'] + kp.kepler_phase2Tc_Tref(planet_b['P'], planet_b['f'], planet_b['e'], planet_b['o']) + Tref,
        instrument['T0_precision'])

    fileout = open('TestCase06_Tcent_b.dat', 'w')
    for i_Tc, v_Tc in enumerate(Tcent_b):
        fileout.write('{0:d}  {1:.4f}  {2:.4f}  {3:d}\n'.format(i_Tc, v_Tc, instrument['T0_precision'], 0))
    fileout.close()

    fileout = open('TestCase06_RV.dat', 'w')
    for b, r in zip(bjd_obs, mod_pl):
        fileout.write('{0:f}  {1:.2f}  {2:.2f}  {3:d}  {4:d}  {5:d}\n'.format(b, r, instrument['RV_precision'], 0, 0, -1))
    fileout.close()


def create_testcase07():
    import george

    bjd0 = bjd_obs - Tref

    err = bjd0 * 0 + instrument['RV_precision']


    gp_pams = np.zeros(4)
    gp_pams[0] = np.log(activity['Hamp_RV1']) * 2
    gp_pams[1] = np.log(activity['Pdec'])*2
    gp_pams[2] = 1. / (2*activity['Oamp'] ** 2)
    gp_pams[3] = np.log(activity['Prot'])

    kernel = np.exp(gp_pams[0]) * \
             george.kernels.ExpSquaredKernel(metric=np.exp(gp_pams[1])) * \
             george.kernels.ExpSine2Kernel(gamma=gp_pams[2], log_period=gp_pams[3])

    gp = george.GP(kernel)

    gp.compute(bjd0, err)
    prediction1 = gp.sample(bjd0)

    gp_pams[0] = np.log(activity['Hamp_RV2']) * 2
    kernel = np.exp(gp_pams[0]) * \
             george.kernels.ExpSquaredKernel(metric=np.exp(gp_pams[1])) * \
             george.kernels.ExpSine2Kernel(gamma=gp_pams[2], log_period=gp_pams[3])

    gp = george.GP(kernel)

    gp.compute(bjd0, err)
    prediction2 = gp.sample(bjd0)

    y_pla = kp.kepler_RV_T0P(bjd0,
                             planet_b['f'],
                             planet_b['P'],
                             planet_b['K'],
                             planet_b['e'],
                             planet_b['o'])

    mod_pl1 = np.random.normal(y_pla + prediction1 + instrument['RV_offset1'], instrument['RV_precision'])
    mod_pl2 = np.random.normal(y_pla + prediction2 + instrument['RV_offset2'], instrument['RV_precision'])

    Tcent_b =  np.random.normal(
        np.arange(0,1)*planet_b['P'] + kp.kepler_phase2Tc_Tref(planet_b['P'], planet_b['f'], planet_b['e'], planet_b['o']) + Tref,
        instrument['T0_precision'])

    fileout = open('TestCase07_RV_dataset1.dat', 'w')
    for b, r in zip(bjd_obs, mod_pl1):
        fileout.write('{0:f}  {1:.2f}  {2:.2f}  {3:d}  {4:d}  {5:d}\n'.format(b, r, instrument['RV_precision'], 0, 0, -1))
    fileout.close()

    fileout = open('TestCase07_RV_dataset2.dat', 'w')
    for b, r in zip(bjd_obs, mod_pl2):
        fileout.write('{0:f}  {1:.2f}  {2:.2f}  {3:d}  {4:d}  {5:d}\n'.format(b, r, instrument['RV_precision'], 0, 0, -1))
    fileout.close()

    fileout = open('TestCase07_Tcent_b.dat', 'w')
    for i_Tc, v_Tc in enumerate(Tcent_b):
        fileout.write('{0:d}  {1:.4f}  {2:.4f}  {3:d}\n'.format(i_Tc, v_Tc, instrument['T0_precision'], 0))
    fileout.close()


def create_testcase08():

    bjd0 = bjd_obs - Tref

    y_pla = kp.kepler_RV_T0P(bjd0,
                             planet_b['f'],
                             planet_b['P'],
                             planet_b['K'],
                             planet_b['e'],
                             planet_b['o']) \
            + instrument['RV_offset1'] \
            + polynomial_trend['c1']*bjd0 + polynomial_trend['c2']*bjd0*bjd0

    mod_pl = np.random.normal(y_pla, instrument['RV_precision'])

    Tcent_b =  np.random.normal(
        np.arange(0, 5)*planet_b['P'] + kp.kepler_phase2Tc_Tref(planet_b['P'], planet_b['f'], planet_b['e'], planet_b['o']) + Tref,
        instrument['T0_precision'])

    fileout = open('TestCase08_RV.dat', 'w')
    for b, r in zip(bjd_obs, mod_pl):
        fileout.write('{0:f}  {1:.2f}  {2:.2f}  {3:d}  {4:d}  {5:d}\n'.format(b, r, instrument['RV_precision'], 0, 0, -1))
    fileout.close()

    fileout = open('TestCase08_Tcent_b.dat', 'w')
    for i_Tc, v_Tc in enumerate(Tcent_b):
        fileout.write('{0:d}  {1:.4f}  {2:.4f}  {3:d}\n'.format(i_Tc, v_Tc, instrument['T0_precision'], 0))
    fileout.close()



create_testcase01()
create_testcase02()
create_testcase03()
create_testcase04()
create_testcase05()
create_testcase06()
create_testcase07()
create_testcase08()

###################################################################
