from pyorbit.classes.common import np
from scipy.optimize import fsolve
import pyorbit.classes.constants as constants

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

__all__ = ["kepler_K1", "kepler_RV", "kepler_Tc2phase_Tref", "kepler_phase2Tc_Tref", "get_planet_mass"]


def kepler_E(M_in, ec):
    E = 0.0
    E0 = 0.0
    M = np.atleast_1d(M_in)
    ecc = np.asarray(ec, dtype=np.double)
    eccanom = np.zeros(np.size(M), dtype=np.double)

    for ii in range(0, np.size(M)):
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


def kepler_K1(m_star1, m_star2, period, i, e0):
    """ Computes the radial velocity semi-amplitude of the primary star

    :param m_star1: mass of the primary, in Solar mass units
    :param m_star2: mass of the secondary/planet, in Solar mass units
    :param period: orbital period of star2, in [d]
    :param i: orbital inclination of star2 wrt the observer (0=face on), in [deg]
    :param e0: orbital eccentricity of star2
    :return: k1, the observed radial velocity semi-amplitude of the primary, in [m s^-1]
    """
    # period must be given in days, conversion factor to seconds are included in the routine
    # constants.Gsi: Gravitational constant in SI system [m^3 kg^-1 s^-2]
    # constants.Msun: Sun mass in SI system [kg]
    # 86400. / constants.d2s: seconds in a day

    return (2. * np.pi * constants.Gsi * constants.Msun / 86400.) ** (1. / 3.) \
           * (np.sin(i * np.pi / 180.0) / np.sqrt(1. - e0 ** 2.)) * period ** (-1. / 3.) \
           * (m_star2 * (m_star1 + m_star2) ** (-2. / 3.))


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
    #

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


def kepler_phase2Tc_Tref(Period, phase, e0, omega0):
    # The closest Tcent after Tref is given back
    TrAn = np.pi / 2 - omega0
    EccAn = 2. * np.arctan(np.sqrt((1.0 - e0) / (1.0 + e0)) * np.tan(TrAn / 2.0))
    MeAn = EccAn - e0 * np.sin(EccAn)
    return (MeAn - phase + omega0) / (2 * np.pi) * Period % Period


def kepler_Tc2phase_Tref(Period, Tcent , e0, omega0):
    # The closest Tcent after Tref is given back
    TrAn = np.pi / 2 - omega0
    EccAn = 2. * np.arctan(np.sqrt((1.0 - e0) / (1.0 + e0)) * np.tan(TrAn / 2.0))
    MeAn = EccAn - e0 * np.sin(EccAn)
    return (omega0 + MeAn - Tcent / Period * 2 * np.pi) % (2 * np.pi)


def f_get_mass(m_star2, m_star1, period, e0, k1):
    """ Computes the difference between the input radial velocity semi-amplitude
    of the primary star and the value corresponding to the provided orbital parameters.
    Supporting function to get_planet_mass subroutine

    :param m_star2: mass of the secondary/planet, in Solar mass units
    :param m_star1: mass of the primary, in Solar mass units
    :param period: orbital period of star2, in [d]
    :param e0: orbital eccentricity of star2
    :param k1: observed RV semi-amplitude of the primary
    :return: the difference between the observed and theoretical RV semi-amplitude of the primary, in [m s^-1]
    """
    # period must be given in days, conversion factor to seconds are included in the routine
    # constants.Gsi: Gravitational constant in SI system [m^3 kg^-1 s^-2]
    # constants.Msun: Sun mass in SI system [kg]
    # 86400. / constants.d2s: seconds in a day

    # M_star1, M_star2 in solar masses
    # P in days -> Period is converted in seconds in the routine
    # inclination assumed to be 90 degrees
    # Gravitational constant in SI system [in m^3 kg^-1 s^-2]
    # output in m/s

    return k1 \
           - ((2. * np.pi * constants.Gsi * constants.Msun / 86400.0) ** (1. / 3.)
              * (1. / np.sqrt(1. - e0 ** 2.))
              * period ** (-1. / 3.)
              * (m_star2 * (m_star1 + m_star2) ** (-2. / 3.)))

def get_approximate_mass(period, k1, e0, m_star1):
    """ Return the approximate mass of the planet in Solar mass units, in the assumption that M_planet << M_star

    :param period: orbital period of star2, in [d]
    :param k1: observed RV semi-amplitude of the primary
    :param e0: orbital eccentricity of star2
    :param m_star1: mass of the primary, in Solar mass units
    :return: mass of the planet, in Solar mass units
    """
    return k1 / ((2. * np.pi * constants.Gsi * constants.Msun / 86400.0) ** (1. / 3.)
         * (1. / np.sqrt(1. - e0 ** 2.))
         * period ** (-1. / 3.)
         * (m_star1 ** (-2. / 3.)))


def get_planet_mass(P, K, e, Mstar, approximation_limit=30.):

    n = np.size(K)
    if n == 1:
        M_approx = min(get_approximate_mass(P, K, e, Mstar), 2*constants.Msear)
        return fsolve(f_get_mass, M_approx, args=(Mstar, P, e, K))

    M_approx = get_approximate_mass(P, K, e, Mstar)

    if np.average(M_approx) > approximation_limit/constants.Msear:
        print('Computing exact mass of the planet (average approximate mass larger than {0:3.1f} Me)'.format(approximation_limit))
        M_init = np.average(M_approx)
        for i in range(0, n):
            M_approx[i] = fsolve(f_get_mass, np.average(M_init), args=(Mstar[i], P[i], e[i], K[i]))
    else:
        print('Computing planetary mass under the approximation M_planet << M_star (threshold at {0:3.1f} Me)'.format(approximation_limit))

    return M_approx
