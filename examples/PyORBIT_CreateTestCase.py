import numpy as np

import sys
import os
import matplotlib.pyplot as plt

sys.path.insert(0, '../pyorbit/classes/')
import kepler_exo as kp


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
    'RV_offset1': np.random.uniform(4000.0,5000.0),
    'RV_offset2': np.random.uniform(4000.0,5000.0),
    'T0_precision': 0.001
}

bjd_obs = np.random.normal(np.arange(6000, 6050, 1, dtype=np.double), 0.2)
Tref = 6025.0000

"""

# planet_b {
    'P': 6.737412701474708,
    'K': 15.964286053214767,
    'f': 1.1271771440198723,
    'e': 0.07578202101359831,
    'i': 86.19154785304445,
    'o': 0.07681624494196555}

# planet_c {
    'P': 14.327606002800504,
    'K': 11.805410430293563,
    'f': 2.63686698774923,
    'e': 0.22822106564776393,
    'i': 85.317208074894,
    'o': 0.6966217861505585}

# instrument {
    'RV_precision': 1.0,
    'RV_offset1': 4779.443746523752,
    'RV_offset2': 4721.741371149122}

"""

""" 
TestCase01: single planet, no external contraints
TestCase02: single planet, transit times
TestCase03: two planets, transit times of the first one

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

    #Tcent_b = kp.kepler_Tcent_T0P(planet_b['P'], planet_b['f'], planet_b['e'], planet_b['o']) + Tref
    Tcent_b =  np.random.normal(
        np.arange(0,5)*kp.kepler_Tcent_T0P(planet_b['P'], planet_b['f'], planet_b['e'], planet_b['o']) + Tref,
        instrument['T0_precision'])

    fileout = open('TestCase02_RV.dat', 'w')
    for b, r in zip(bjd_obs, mod_pl):
        fileout.write('{0:f}  {1:.2f}  {2:.2f}  {3:d}  {4:d}  {5:d}\n'.format(b, r, instrument['RV_precision'], 0, 0, -1))
    fileout.close()

    #fileout = open('TestCase01_Tcent_b.dat', 'w')
    #for ii in xrange(0, np.size(Transit_Time)):
    #    fileout.write('{0:d} {1:f} {2:f} {3:d} \n'.format(ii, Tcent_b, 0.01, 0))
    #fileout.close()

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

    #Tcent_b = kp.kepler_Tcent_T0P(planet_b['P'], planet_b['f'], planet_b['e'], planet_b['o']) + Tref
    Tcent_b =  np.random.normal(
        np.arange(0,5)*kp.kepler_Tcent_T0P(planet_b['P'], planet_b['f'], planet_b['e'], planet_b['o']) + Tref,
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




create_testcase01()
create_testcase02()
create_testcase03()

###################################################################


























def create_test_1planet_circular():
    x = np.arange(6000, 6100, 1, dtype=np.double)
    x = np.random.normal(x, 0.4)
    Tref = np.mean(x, dtype=np.double)
    x0 = x - Tref

    print Tref
    P = 23.4237346
    K = 43.47672
    phase = 0.34658203
    e = 0.000
    omega = np.pi/2.
    offset = 45605

    y_pla = kp.kepler_RV_T0P(x0, phase, P, K, e, omega) + offset

    mod_pl = np.random.normal(y_pla, 2)

    Transit_Time = kp.kepler_Tcent_T0P(P, phase, e, omega) + Tref
    print 'Transit Time:', Transit_Time
    print 'transit Time - Tref', Transit_Time - Tref

    print Transit_Time, np.floor((Transit_Time-Tref+0.01) / P) * P + Tref + \
            kp.kepler_Tcent_T0P(P, phase, e, omega)
    quit()

    plt.scatter(x, mod_pl)
    plt.axvline(Transit_Time)
    plt.show()

    fileout = open('test_1planet_circular_RV.dat', 'w')
    for ii in xrange(0, np.size(x)):
        fileout.write('{0:14f} {1:14f} {2:14f} {3:5d} {4:5d} {5:5d} \n'.format(x[ii], mod_pl[ii], 2., 0, 0, -1))
    fileout.close()

    fileout = open('test_1planet_circular_Tcent_0.dat', 'w')
    for ii in xrange(0, np.size(Transit_Time)):
        fileout.write('{0:6d} {1:14f} {2:14f} \n'.format(ii, Transit_Time, 0.01))
    fileout.close()


def create_test_1planet_GP():

    import george
    x = np.arange(6000, 6200, 1., dtype=np.double)
    x = np.random.normal(x, 0.4)
    Tref = np.mean(x, dtype=np.double)
    x0 = x - Tref

    err = x*0.0 + 2.0

    print Tref
    P = 23.4237346
    K = 43.47672
    phase = 0.34658203
    e = 0.13
    omega = 0.673434
    offset = 0 #45605

    Transit_Time = kp.kepler_Tcent_T0P(P, phase, e, omega) + Tref + P*trip
    print Transit_Time

    y_pla = kp.kepler_RV_T0P(x0, phase, P, K, e, omega) + offset

    input_pams = {'Prot': 17.37462, 'Pdec': 50.93840, 'Oamp': 0.35, 'Hamp': 50.03 }

    gp_pams = np.zeros(4)
    gp_pams[0] = np.log(input_pams['Hamp'])*2
    gp_pams[1] = np.log(input_pams['Pdec'])*2
    gp_pams[2] = 1. / (2*input_pams['Oamp'] ** 2)
    gp_pams[3] = np.log(input_pams['Prot'])


    kernel = np.exp(gp_pams[0]) * \
                      george.kernels.ExpSquaredKernel(metric=np.exp(gp_pams[1])) * \
                      george.kernels.ExpSine2Kernel(gamma=gp_pams[1], log_period=gp_pams[2])
    gp = george.GP(kernel, solver=george.HODLRSolver, seed=42)
    #gp = george.GP(kernel)
    gp.compute(x0, err)
    pred = gp.sample(x)

    mod_pl = np.random.normal(y_pla+pred, 2)

    plt.plot(x, y_pla, c='r')
    plt.plot(x, pred, c='b')
    plt.scatter(x, mod_pl, c='y')

    vert_lines = np.arange(-50.0, 50.0, input_pams['Prot'])+Tref
    for v in vert_lines:
        plt.axvline(v, c='b')

    for v in Transit_Time:
        plt.axvline(v, c='r')
    plt.show()

    fileout = open('test_1planet_GP_RV.dat', 'w')
    for ii in xrange(0, np.size(x)):
        fileout.write('{0:14f} {1:14f} {2:14f} {3:5d} {4:5d} {5:5d} \n'.format(x[ii], mod_pl[ii], 2., 0, 0, -1))
    fileout.close()

    fileout = open('test_1planet_GP_Tcent_0.dat', 'w')
    for ii in xrange(0, np.size(Transit_Time)):
        fileout.write('{0:6d} {1:14f} {2:14f} \n'.format(ii, Transit_Time[ii], 0.01))
    fileout.close()


def create_test_1planet_multijit_multioff():
    x = np.arange(6000, 6100, 1, dtype=np.double)
    x = np.random.normal(x, 0.4)
    Tref = np.mean(x, dtype=np.double)
    x0 = x - Tref

    print Tref
    P = 23.4237346
    K = 43.47672
    phase = 0.34658203
    e = 0.13
    omega = 0.673434
    offset = 45605

    y_pla = kp.kepler_RV_T0P(x0, phase, P, K, e, omega) + offset

    mod_pl = np.random.normal(y_pla, 2)

    Transit_Time = kp.kepler_Tcent_T0P(P, phase, e, omega) + Tref
    print 'Transit Time:', Transit_Time
    print 'transit Time - Tref', Transit_Time - Tref

    jit = np.zeros(np.size(x))
    jit[x>6030.0]+=1
    jit[x>6070.0]+=1

    off = np.zeros(np.size(x))
    off[x>6030.0]+=1
    off[x>6070.0]+=1

    plt.scatter(x, mod_pl)
    plt.axvline(Transit_Time)
    plt.show()

    fileout = open('test_1planet_RV_multijit_multioff.dat', 'w')
    for ii in xrange(0, np.size(x)):
        fileout.write('{0:f} {1:f} {2:f} {3:f} {4:f} {5:f} \n'.format(x[ii], mod_pl[ii], 2., jit[ii], off[ii], -1))
    fileout.close()

    fileout = open('test_1planet_Tcent_multijit_multioff.dat', 'w')
    for ii in xrange(0, np.size(Transit_Time)):
        fileout.write('{0:6d} {1:14f} {2:14f} {3:d} \n'.format(ii, Transit_Time, 0.01, 0))
    fileout.close()


def create_test_2planets_multijit_multioff():
    x = np.arange(6000, 6100, 1, dtype=np.double)
    x = np.random.normal(x, 0.4)
    Tref = np.mean(x, dtype=np.double)
    x0 = x - Tref
    offset = 45605

    print Tref
    P = 23.4237346
    K = 43.47672
    phase = 0.34658203
    e = 0.13
    omega = 0.673434

    y_pla = kp.kepler_RV_T0P(x0, phase, P, K, e, omega)
    Transit_Time1 = kp.kepler_Tcent_T0P(P, phase, e, omega) + Tref
    print 'Transit Time:', Transit_Time1
    print 'transit Time - Tref', Transit_Time1 - Tref

    P = 3.4237346
    K = 1.47672
    phase = 2.0658203
    e = 0.01
    omega = 5.573434

    y_pla += kp.kepler_RV_T0P(x0, phase, P, K, e, omega)
    Transit_Time2 = kp.kepler_Tcent_T0P(P, phase, e, omega) + Tref
    print 'Transit Time:', Transit_Time2
    print 'transit Time - Tref', Transit_Time2 - Tref

    y_pla += offset

    mod_pl = np.random.normal(y_pla, 2)


    jit = np.zeros(np.size(x))
    jit[x>6030.0]+=1
    jit[x>6070.0]+=1

    off = np.zeros(np.size(x))
    off[x>6030.0]+=1
    off[x>6070.0]+=1

    fileout = open('test_2planets_RV_multijit_multioff.dat', 'w')
    for ii in xrange(0, np.size(x)):
        fileout.write('{0:f} {1:f} {2:f} {3:f} {4:f} {5:f} \n'.format(x[ii], mod_pl[ii], 2., jit[ii], off[ii], -1))
    fileout.close()

    fileout = open('test_2planets_Tcent_b_multijit_multioff.dat', 'w')
    for ii in xrange(0, np.size(Transit_Time1)):
        fileout.write('{0:6d} {1:14f} {2:14f} {3:d} \n'.format(ii, Transit_Time1, 0.01, 0))
    fileout.close()

    fileout = open('test_2planets_Tcent_c_multijit_multioff.dat', 'w')
    for ii in xrange(0, np.size(Transit_Time2)):
        fileout.write('{0:6d} {1:14f} {2:14f} {3:d} \n'.format(ii, Transit_Time2, 0.01, 0))
    fileout.close()


def create_test_1planet_polynomial_trend():
    x = np.arange(6000, 6100, 1, dtype=np.double)
    x = np.random.normal(x, 0.4)
    Tref = np.mean(x, dtype=np.double)
    x0 = x - Tref

    print Tref
    P = 23.4237346
    K = 43.47672
    phase = 0.34658203
    e = 0.13
    omega = 0.673434
    offset = 45605

    order = 2
    import numpy.polynomial.polynomial

    coeff = np.zeros(order+1)
    coeff[1] = 0.256
    coeff[2] = 0.0120

    rv_trend = numpy.polynomial.polynomial.polyval(x0, coeff)

    y_pla = kp.kepler_RV_T0P(x0, phase, P, K, e, omega) + offset + rv_trend

    mod_pl = np.random.normal(y_pla, 2)

    Transit_Time = kp.kepler_Tcent_T0P(P, phase, e, omega) + Tref
    print 'Transit Time:', Transit_Time
    print 'transit Time - Tref', Transit_Time - Tref

    plt.scatter(x, mod_pl)
    plt.axvline(Transit_Time)
    plt.show()

    fileout = open('test_1planet_polynomial_trend_RV.dat', 'w')
    for ii in xrange(0, np.size(x)):
        fileout.write('{0:14f} {1:14f} {2:14f} {3:5d} {4:5d} {5:5d} \n'.format(x[ii], mod_pl[ii], 2., 0, 0, -1))
    fileout.close()

    fileout = open('test_1planet_polynomial_trend_Tcent_0.dat', 'w')
    for ii in xrange(0, np.size(Transit_Time)):
        fileout.write('{0:6d} {1:14f} {2:14f} \n'.format(ii, Transit_Time, 0.01))
    fileout.close()



