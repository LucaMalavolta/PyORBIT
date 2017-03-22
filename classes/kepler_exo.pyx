import numpy as np

# +
# NAME:
#    exofast_keplereq
# PURPOSE: 
#    Solve Kepler's Equation
# DESCRIPTION:
#    Solve Kepler's Equation. Method by S. Mikkola (1987) Celestial
#       Mechanics, 40 , 329-334. 
#    result from Mikkola then used as starting value for
#       Newton-Raphson iteration to extend the applicability of this
#       function to higher eccentricities


G_grav = 6.67428e-11 # Gravitational Constants in SI system [m^3/kg/s^2]
M_sun = 1.9884e30 # Value from TRADES


def kepler_E(M_in, ec):
    E = 0.0
    E0 = 0.0
    M = np.atleast_1d(M_in)
    ecc = np.asarray(ec, dtype=np.double)
    eccanom = np.zeros(np.size(M), dtype=np.double)

    for ii in xrange(0, np.size(M)):
        # -np.pi < M < np.pi
        mx = M[ii]
        if mx > np.pi:
            mx = mx % (2. * np.pi)
            if mx > np.pi:
                mx = mx - (2. * np.pi)
        if mx <= -np.pi:
            mx = mx % (2. * np.pi)
            if mx < -np.pi:
                mx += (2. * np.pi)

        if ecc < 1e-10:
            eccanom[ii] = mx
        else:

            # equation 9a
            aux = 4.0 * ecc + 0.50
            alpha = (1.0 - ecc) / aux
            beta = mx / (2.0 * aux)

            # equation 9b
            ## the actual equation 9b is much much slower, but gives the same
            ## answer (probably because more refinement necessary)
            aux = np.sqrt(beta * beta + alpha * alpha * alpha)
            z = beta + aux
            if z < 0.:
                z = beta - aux
            z = z ** (1. / 3.)

            if abs(z) < 1e-8:
                s0 = 0.
            else:
                s0 = z - alpha / z

            s1 = s0 - (0.078 * s0 ** 5) / ((1.) + ecc)
            e0 = mx + ecc * (3. * s1 - 4. * s1 ** 3.)

            se0 = np.sin(e0)
            ce0 = np.cos(e0)

            f = e0 - ecc * se0 - mx
            f1 = (1.0) - ecc * ce0
            f2 = ecc * se0
            f3 = ecc * ce0
            u1 = -f / f1
            u2 = -f / (f1 + 0.5 * f2 * u1)
            u3 = -f / (f1 + 0.5 * f2 * u2 + (1. / 6.) * f3 * u2 ** 2.)
            u4 = -f / (f1 + 0.5 * f2 * u3 + (1. / 6.) * f3 * u3 ** 2 - (1. / 24.) * f2 * u3 ** 3)

            ecan_tmp = e0 + u4

            if ecan_tmp >= 2. * np.pi:
                ecan_tmp = ecan_tmp - 2. * np.pi

            if ecan_tmp < 0.:
                ecan_tmp = ecan_tmp + 2. * np.pi

            ## Now get more precise solution using Newton Raphson method
            ## for those times when the Kepler equation is not yet solved
            ## to better than 1e-10
            ## (modification J. Wilms)

            if mx < 0.:
                mx = mx + 2. * np.pi

            ## calculate the differences
            diff = abs(ecan_tmp - ecc * np.sin(ecan_tmp) - mx)
            if diff > abs(diff - 2 * np.pi):
                diff = abs(diff - 2 * np.pi)

            thresh1 = 1e-8
            thresh2 = 10000
            countt = 0

            while (diff > thresh1 and countt < thresh2):
                ## E-e sinE-M
                fe = (ecan_tmp - ecc * np.sin(ecan_tmp) - mx) % (2 * np.pi)
                ## f' = 1-e*cosE
                fs = (1. - ecc * np.cos(ecan_tmp)) % (2 * np.pi)
                oldval = ecan_tmp
                ecan_tmp = (oldval - fe / fs)
                diff = abs(oldval - ecan_tmp)
                countt += 1
            ## range reduction
            if ecan_tmp >= 2. * np.pi:
                ecan_tmp = ecan_tmp % 2. * np.pi
            if ecan_tmp < 0.:
                ecan_tmp = ecan_tmp % 2. * np.pi + 2. * np.pi

            eccanom[ii] = ecan_tmp

    return eccanom


def kepler_K1(M_star1, M_star2, Period, i, e0):
    # M_star1, M_star2 in solar masses
    # P in days -> Period is converted in seconds in the routine
    # i in degrees
    # Gravitational constant is given in m^3 kg^-1 s^-2
    # output in m/s
    K1 = (2. * np.pi * G_grav * M_sun / 86400.0) ** (1.0 / 3.0) * (
    np.sin(i * np.pi / 180.0) / np.sqrt(1.0 - e0 ** 2.0)) * (Period) ** (-1.0 / 3.0) * (
         M_star2 * (M_star1 + M_star2) ** (-2.0 / 3.0))
    return K1


def kepler_RV(BJD, TPeri, Period, gamma, K, e0, omega0):
    # omega = argument of pericenter
    # Mean Anomaly
    #
    MeAn = 2. * np.pi * (1. + (BJD - TPeri) / Period % 1.)

    if abs(e0) < 1e-3:
        TrAn = np.asarray(MeAn, dtype=np.double)
        e = np.asarray(0., dtype=np.double)
        omega = np.asarray(0., dtype=np.double)
    else:
        if e0 < 0.:
            e = np.asarray(-e0, dtype=np.double)
            omega = np.asarray(omega0, dtype=np.double) + np.pi
        else:
            e = np.asarray(e0, dtype=np.double)
            omega = np.asarray(omega0, dtype=np.double)

        # Eccentric Anomaly
        EccAn = kepler_E(MeAn, e)
        TrAn = 2. * np.arctan(np.sqrt((1.0 + e) / (1.0 - e)) * np.tan(EccAn / 2.0))

    rv = K * (np.cos(TrAn + omega) + e * np.cos(omega)) + gamma

    return rv


def kepler_RV_T0P(BJD0, phase, Period, K, e0, omega0):
    # BJD0 is given as BJD-T0, where T0 is arbitrarily defined by the user
    # Tperi_ is substituted by _phase_, which is the phase of the orbit where
    #        BJD0+T0+phase*Period = Tperi
    # omega = argument of pericenter
    # Mean Anomaly
    #
    # MeAn = 2.*np.pi*(1. + (BJD - TPeri)/Period % 1.)
    # MeAn = 2.*np.pi*(1. + (BJD0 - phase*Period)/Period % 1.)

    # MeAn = 2.*np.pi*(1. + ((BJD0/Period)-phase) % 1.)
    omega = np.asarray(omega0, dtype=np.double)
    e = np.asarray(e0, dtype=np.double)
    MeAn = 2. * np.pi * (1. + ((BJD0 / Period) + (phase - omega0) / (2 * np.pi)) % 1.)

    if abs(e0) < 1e-3:
        TrAn = np.asarray(MeAn, dtype=np.double)
        e = np.asarray(0., dtype=np.double)
    else:
        if e0 < 0.:
            e = -1 * e
            omega += np.pi

        # Eccentric Anomaly
        EccAn = kepler_E(MeAn, e)
        TrAn = 2. * np.arctan(np.sqrt((1.0 + e) / (1.0 - e)) * np.tan(EccAn / 2.0))

    rv = K * (np.cos(TrAn + omega) + e * np.cos(omega))

    return rv


def kepler_Tcent_T0P(Period, phase, e0, omega0):
    # The closest Tcent after Tref is given back
    TrAn = np.pi / 2 - omega0
    EccAn = 2. * np.arctan(np.sqrt((1.0 - e0) / (1.0 + e0)) * np.tan(TrAn / 2.0))
    MeAn = EccAn - e0 * np.sin(EccAn)
    return (MeAn - phase + omega0) / (2 * np.pi) * Period % Period
