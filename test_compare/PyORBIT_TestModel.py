import numpy as np

import sys
import os
import matplotlib.pyplot as plt

sys.path.insert(0, '../pyorbit/classes/')
import kepler_exo as kp


def create_test_sinusoids():

    x = np.arange(6000, 6500, 1, dtype=np.double)
    Tref = np.mean(x, dtype=np.double)
    x0 = x - Tref

    P = 23.4237346
    K = 43.47672
    phase = 0.34658203
    e = 0.13
    omega = 0.673434
    offset = 45605

    y_pla = kp.kepler_RV_T0P(x0, phase, P, K, e, omega) + offset

    mod_pl = np.random.normal(y_pla, 2)

    n_periods = 2
    n_pha = 4
    # Period1

    Prot = 9.234
    pha = np.zeros([n_periods, n_pha], dtype=np.double)
    pha[0, :] = [0.452, 0.894, 0.622, 0.344]
    pha[1, :] = [0.236, 0.522, 0.974, 0.334]

    off = np.zeros(n_periods, dtype=np.double)
    off[:] = 0.56424

    amp = np.zeros([n_periods, n_pha], dtype=np.double)
    amp[0, :] = [35.0, 24.0, 17.0, 7.0]
    amp[1, :] = [30.0, 19.0, 12.0, 8.0]

    p_mask = np.zeros([n_periods, np.size(x)], dtype=np.bool)
    p_mask[0, :] = (x > 5000) & (x < 6250)
    p_mask[1, :] = (x > 6250) & (x < 7000)

    har = np.arange(1, np.size(amp, axis=1) + 1, 1., dtype=np.double)

    mod_rv = np.zeros(np.size(x), dtype=np.double)
    mod_ph1 = np.zeros(np.size(x), dtype=np.double)
    mod_ind = np.zeros([n_periods, np.size(x)], dtype=np.double)
    xph = (x0 / Prot) % 1

    for ii in xrange(0, n_periods):
        for jj in xrange(0, n_pha):
            mod_ind[ii, :] = p_mask[ii, :] * amp[ii, jj] * np.sin((har[jj] * xph + pha[ii, jj] + off[ii]) * 2. * np.pi)
            mod_rv += mod_ind[ii, :]
            # plt.plot(x[p_mask[ii,:]],mod_ind[ii,p_mask[ii,:]])

    amp = np.zeros([n_periods, n_pha], dtype=np.double)
    amp[0, :] = [0.0350, 0.0240, 0.0170, 0.0070]
    amp[1, :] = [0.0300, 0.0190, 0.0120, 0.0080]

    for ii in xrange(0, n_periods):
        for jj in xrange(0, n_pha):
            mod_ind[ii, :] = p_mask[ii, :] * amp[ii, jj] * np.sin((har[jj] * xph + pha[ii, jj]) * 2. * np.pi)
            mod_ph1 += mod_ind[ii, :]
            # plt.plot(x[p_mask[ii,:]],mod_ind[ii,p_mask[ii,:]])

    mod_ph = np.random.normal(mod_ph1, 0.002)

    fileout = open('example/test_sinusoids_RV.dat', 'w')
    for ii in xrange(0, np.size(x)):
        fileout.write(
            '{0:14f} {1:14f} {2:14f} {3:5d} {4:5d} {5:5d} \n'.format(x[ii], mod_rv[ii] + mod_pl[ii], 1., 0, 0, -1))
    fileout.close()

    fileout = open('test1/test_sinusoid_PH.dat', 'w')
    for ii in xrange(0, np.size(x)):
        fileout.write('{0:14f} {1:14f} {2:14f} {3:5d} {4:5d} {5:5d} \n'.format(x[ii], mod_ph[ii], 0.002, 0, 0, -1))
    fileout.close()


def create_test_1planet():
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

    trip = np.arange(10, 3, -1)
    Transit_Time = kp.kepler_Tcent_T0P(P, phase, e, omega) + Tref
    print 'Transit Time:', Transit_Time
    print 'transit Time - Tref', Transit_Time - Tref

    plt.scatter(x, mod_pl)
    plt.axvline(Transit_Time)
    plt.show()

    fileout = open('test_1planet_RV.dat', 'w')
    for ii in xrange(0, np.size(x)):
        fileout.write('{0:14f} {1:14f} {2:14f} {3:5d} {4:5d} {5:5d} \n'.format(x[ii], mod_pl[ii], 2., 0, 0, -1))
    fileout.close()

    fileout = open('test_1planet_Tcent_0.dat', 'w')
    for ii in xrange(0, np.size(Transit_Time)):
        fileout.write('{0:6d} {1:14f} {2:14f} \n'.format(ii, Transit_Time, 0.01))
    fileout.close()


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

    trip = np.arange(10, 3, -1)
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

    trip = np.arange(-5, 6, 1)
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

    trip = np.arange(10, 3, -1)
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

    rv_trend = numpy.polynomial.polynomial.polyval(x-Tref-22.5, coeff)

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


create_test_1planet_polynomial_trend()

#create_test_1planet_GP()
#create_test_1planet_circular()
#create_test_1planet_multijit_multioff()
#create_test_2planets_multijit_multioff()
