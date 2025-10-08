import sys
#sys.path.insert(0, '/Users/malavolta/Astro/CODE/PyORBIT/')

from pyorbit.subroutines.common import np
from scipy.optimize import fsolve, newton
import pyorbit.subroutines.constants as constants


__all__ = ["kepler_K1", "kepler_RV", "kepler_RV_T0P", "kepler_Tc2phase_Tref", "kepler_phase2Tc_Tref", "get_planet_mass", "kepler_true_anomaly_orbital_distance", 
            "kepler_compute_rv_semiamplitude", "kepler_compute_rv", "kepler_compute_rv_deltabjd", "kepler_compute_deltaTc_from_meanlong",
            "kepler_compute_meanlong_from_deltaTc", "kepler_compute_deltaTperi_from_deltaTc", "kepler_compute_deltaTperi_from_meanlong",
            "kepler_get_planet_mass", "kepler_compute_trueanomaly_orbitaldistance"]

#TODO: user commented definition below in PyORBIT version 12

#__all__ = ["kepler_compute_rv_semiamplitude", "kepler_compute_rv", "kepler_compute_rv_deltabjd", "kepler_compute_deltaTc_from_meanlong",
#           "kepler_compute_meanlong_from_deltaTc", "kepler_compute_deltaTperi_from_deltaTc", "kepler_compute_deltaTperi_from_meanlong",
#           "kepler_get_planet_mass", "kepler_compute_trueanomaly_orbitaldistance"]

def kepler_K1(mass_primary, mass_secondary, period, inclination, eccentricity):
    """Alias for kepler_compute_rv_semiamplitude, to ensure back-compatibility
    #TODO: remove in PyORBIT version 12
    """
    return kepler_compute_rv_semiamplitude(mass_primary, mass_secondary, period, inclination, eccentricity)

def kepler_RV(bjd, Tperi, period, rv_semiamplitude, ecc, omega_deg):
    """Alias for kepler_compute_rv, to ensure back-compatibility
    #TODO:  remove in PyORBIT version 12
    """
    return kepler_compute_rv(bjd, Tperi, period, rv_semiamplitude, ecc, omega_deg)

def kepler_RV_T0P(bjd_tref, mean_long, period, rv_semiamplitude, ecc=0., omega_deg=90., Omega_deg=0.0):
    """Alias for kepler_compute_rv_deltabjd, to ensure back-compatibility
    #TODO: remove in PyORBIT version 12 
    """
    return kepler_compute_rv_deltabjd(bjd_tref, rv_semiamplitude, period, mean_long, ecc, omega_deg, Omega_deg)

def kepler_phase2Tc_Tref(period, mean_long, ecc=0., omega_deg=90., Omega_deg=0.0):
    """Alias for kepler_compute_deltaTc_from_meanlong, to ensure back-compatibility
    #TODO: remove in PyORBIT version 12 
    """
    return kepler_compute_deltaTc_from_meanlong(period, mean_long, ecc, omega_deg, Omega_deg)

def kepler_Tc2phase_Tref(period, delta_Tc, ecc, omega_deg, Omega_deg=0.0):
    """Alias for kepler_compute_meanlong_from_deltaTc, to ensure back-compatibility
    #TODO: remove in PyORBIT version 12 
    """
    return kepler_compute_meanlong_from_deltaTc(period, delta_Tc, ecc, omega_deg, Omega_deg)

def kepler_Tc2Tperi_Tref(period, delta_Tc, ecc, omega_deg):
    """Alias for kepler_compute_meanlong_from_deltaTc, to ensure back-compatibility
    #TODO: remove in PyORBIT version 12 
    """
    return kepler_compute_deltaTperi_from_deltaTc(period, delta_Tc , ecc, omega_deg)

def kepler_phase2Tperi_Tref(period, mean_long, ecc, omega_deg, Omega_deg):
    """Alias for kepler_compute_meanlong_from_deltaTc, to ensure back-compatibility
    #TODO: remove in PyORBIT version 12 
    """
    return kepler_compute_deltaTperi_from_meanlong(period, mean_long, ecc, omega_deg, Omega_deg)

def get_planet_mass(period, rv_semiamplitude, ecc, mass_star, approximation_limit=30.):
    """Alias for kepler_get_planet_mass, to ensure back-compatibility
    #TODO: remove in PyORBIT version 12 
    """
    return kepler_get_planet_mass(period, rv_semiamplitude, ecc, mass_star, approximation_limit)

def kepler_true_anomaly_orbital_distance(bjd_tref, delta_Tc, period, ecc, omega_deg, semimajor_axis, Omega_deg=0.0):
    """Alias for kepler_compute_trueanomaly_orbitaldistance, to ensure back-compatibility
    #TODO: remove in PyORBIT version 12 
    """
    return kepler_compute_trueanomaly_orbitaldistance(bjd_tref, semimajor_axis, delta_Tc, period, ecc, omega_deg, Omega_deg)


def f0_keplerE(ecan_tmp, ecc, mx):
    """ Support function for kepler_E """
    return (ecan_tmp - ecc * np.sin(ecan_tmp) - mx) % (2 * np.pi)

def f1_keplerE(ecan_tmp, ecc, mx):
    """ Support function for kepler_E """
    ## f' = 1-e*cosE
    return (1. - ecc * np.cos(ecan_tmp)) % (2 * np.pi)

def kepler_E(M_in, ecc):
    """ Solves Kepler's equation to find the eccentric anomaly, E, given the mean anomaly, M, and eccentricity, e.
    """
    mx = np.atleast_1d(M_in) % (2. * np.pi)

    """
    #if ecc < 1e-10:
    #    return mx

    # equation 9a
    aux = 4.0 * ecc + 0.50
    alpha = (1.0 - ecc) / aux
    beta = mx / (2.0 * aux)

    # equation 9b
    ## the actual equation 9b is much much slower, but gives the same
    ## answer (probably because more refinement necessary)
    aux = np.sqrt(beta * beta + alpha * alpha * alpha)
    z = beta + aux
    sel_z = (z < 0.)

    z[sel_z] = beta[sel_z] + aux[sel_z]

    z = z ** (1. / 3.)

    s0 = mx * 0.
    sel_z = (abs(z) > 1e-8)
    s0[sel_z] = z - alpha / z

    s1 = s0 - (0.078 * s0 ** 5) / ((1.) + ecc)

    # first guess of e0
    e0 = mx + ecc * (3. * s1 - 4. * s1 ** 3.)

    # difference from first estimate of E
    e0_diff = e0 - ecc * np.sin(e0) - mx

    """

    # Using Fulton scheme
    k = 0.85

    e0 = mx + np.sign(np.sin(mx)) * k * ecc  # first guess at E
    # fiarr should go to zero when converges
    e0_diff = ( e0 - ecc * np.sin(e0) - mx)


    conv = 1.0e-12  # convergence criterion
    convd = np.where(np.abs(e0_diff) > conv)[0]  # which indices have not converged
    nd = len(convd)  # number of unconverged elements
    count = 0

    while nd > 0:  # while unconverged elements exist
        count += 1

        e0_conv = e0[convd]# just the unconverged elements ...

        se0 = np.sin(e0_conv)
        ce0 = np.cos(e0_conv)

        f = e0_conv - ecc * se0 - mx[convd]
        f1 = 1. - ecc * ce0
        f2 = ecc * se0
        f3 = ecc * ce0
        u1 = -f / f1
        u2 = -f / (f1 + 0.5 * f2 * u1)
        u3 = -f / (f1 + 0.5 * f2 * u2 + (1. / 6.) * f3 * u2 ** 2.)
        u4 = -f / (f1 + 0.5 * f2 * u3 + (1. / 6.) * f3 * u3 ** 2 - (1. / 24.) * f2 * u3 ** 3)

        e0_conv += u4

        e0[convd] = e0_conv
        e0_diff = e0 - ecc * np.sin( e0 ) - mx
        convd = np.abs(e0_diff) > conv  # test for convergence
        nd = np.sum(convd is True)


        #ecan_tmp = newton(f0_keplerE, ecan_tmp, fprime=f1_keplerE, args=(ecc, mx))
        #ecan_tmp = ecan_tmp % (2. * np.pi)
    ecan_tmp = e0 % (2. * np.pi)

    return ecan_tmp


def kepler_compute_rv_semiamplitude(mass_primary, mass_secondary, period, inclination, eccentricity):
    """ Computes the radial velocity semi-amplitude of the primary star
    Period must be given in days, conversion factor to seconds are included in the routine
    constants.Gsi: Gravitational constant in SI system [m^3 kg^-1 s^-2]
    constants.Msun: Sun mass in SI system [kg]

    :param mass_primary: mass of the primary (in Solar mass units)
    :param mass_secondary: mass of the secondary/planet (in Solar mass units)
    :param period: orbital period of star2, (in days)
    :param inclination: orbital inclination of star2 wrt the observer (0=face on) (in degrees)
    :param eccentricity: orbital eccentricity of star2
    :return: k1, the observed radial velocity semi-amplitude of the primary (in m/s)
    """

    # 86400. / constants.d2s: seconds in a day

    return (2. * np.pi * constants.Gsi * constants.Msun / 86400.) ** (1. / 3.) \
           * (np.sin(inclination * constants.deg2rad) / np.sqrt(1. - eccentricity ** 2.)) * period ** (-1. / 3.) \
           * (mass_secondary * (mass_primary + mass_secondary) ** (-2. / 3.))


def kepler_compute_rv(bjd, Tperi, period, rv_semiamplitude, ecc, omega_deg):
    """ Compute the radial velocity at a given time (bjd) if the time of periastron passage (Tperi)
    is known, assuming a Keplerian orbit with given period,
    eccentricity (ecc), and argument of pericenter (omega_deg).
    :param bjd: Barycentric Julian Date at which to compute the radial velocity (in days)
    :param Tperi: Time of periastron passage (in days)
    :param period: Orbital period of the planet (in days)
    :param rv_semiamplitude: Radial velocity semi-amplitude (in m/s)
    :param ecc: Orbital eccentricity (default is 0)
    :param omega_deg: Argument of pericenter (in degrees, default is 90)
    :return: Radial velocity at time bjd (in m/s)
    """

    MeAn = 2. * np.pi * (1. + (bjd - Tperi) / period % 1.)

    if abs(ecc) < 1e-3:
        TrAn = np.asarray(MeAn, dtype=np.double)
        ecc_eq = np.asarray(0., dtype=np.double)
        omega = np.asarray(0., dtype=np.double)
    else:
        if ecc < 0.:
            ecc_eq = np.asarray(-ecc, dtype=np.double)
            omega = np.asarray(omega_deg, dtype=np.double)*constants.deg2rad + np.pi
        else:
            ecc_eq = np.asarray(ecc, dtype=np.double)
            omega = np.asarray(omega_deg, dtype=np.double)*constants.deg2rad

        # Eccentric Anomaly
        EccAn = kepler_E(MeAn, ecc_eq)
        TrAn = 2. * np.arctan(np.sqrt((1.0 + ecc_eq) / (1.0 - ecc_eq)) * np.tan(EccAn / 2.0))

    rv = rv_semiamplitude * (np.cos(TrAn + omega) + ecc_eq * np.cos(omega))

    return rv


def kepler_compute_rv_deltabjd(bjd_tref, rv_semiamplitude, period, mean_long, ecc, omega_deg, Omega_deg=0.0):
    """
    Compute the radial velocity at a given time difference from a reference epoch (bjd_tref=BJD-Tref)
    from the mean longitude at reference epoch (mean_long) assuming a Keplerian orbit with given period, 
    eccentricity (ecc), argument of pericenter (omega_deg), and longitude of the ascending node (Omega_deg).
    :param bjd_tref: Time difference from reference epoch (in days)
    :param rv_semiamplitude: Radial velocity semi-amplitude (in m/s)
    :param period: Orbital period of the planet (in days)
    :param mean_long: Mean longitude at reference epoch (in degrees)
    :param ecc: Orbital eccentricity (default is 0)
    :param omega: Argument of pericenter (in degrees, default is 90)
    :param Omega: Longitude of the ascending node (in degrees, default is 0)
    :return: Radial velocity at time bjd_tref (in m/s)
    """

    omega = np.asarray(omega_deg, dtype=np.double)* constants.deg2rad
    Omega = np.asarray(Omega_deg, dtype=np.double)* constants.deg2rad
    ecc_eq = np.asarray(ecc, dtype=np.double)
    MeAn = (2*np.pi *(bjd_tref / period) + (mean_long*constants.deg2rad - omega - Omega)) % (2 * np.pi)

    if abs(ecc) < 1e-3:
        TrAn = np.asarray(MeAn, dtype=np.double)
        ecc_eq = np.asarray(0., dtype=np.double)
    else:
        if ecc < 0.:
            ecc_eq = -1 * ecc_eq
            omega += np.pi

        # Eccentric Anomaly
        EccAn = kepler_E(MeAn, ecc_eq)
        TrAn = 2. * np.arctan(np.sqrt((1.0 + ecc_eq) / (1.0 - ecc_eq)) * np.tan(EccAn / 2.0))

    rv = rv_semiamplitude * (np.cos(TrAn + omega) + ecc_eq * np.cos(omega))

    return rv


def kepler_compute_deltaTc_from_meanlong(period, mean_long, ecc=0., omega_deg=90., Omega_deg=0.0):
    """
    Compute the difference between the first *time of inferior conjunction* (Tc) and the reference time Tref
    from the mean longitude at reference epoch (mean_long) assuming a Keplerian orbit with given period, 
    eccentricity (ecc), argument of pericenter (omega_deg), and longitude of the ascending node (Omega_deg).
    The *time of inferior conjunction* is defined as the time when the planet is closest to the line of sight
    between the observer and the star, i.e. when the true anomaly (TrAn) is equal to pi/2 - omega.
    For circular orbits, the time of inferior conjunction corresponds to the central time of transit.

    :param period: Orbital period of the planet (in days)
    :param mean_long: Mean longitude at reference epoch (in degrees)
    :param ecc: Orbital eccentricity (default is 0)
    :param omega: Argument of pericenter (in degrees, default is 90)
    :param Omega: Longitude of the ascending node (in degrees, default is 0)
    :return: minimum positive difference between the time of inferior conjunction (Tv) and reference time (Tref)
    """
    omega = omega_deg * constants.deg2rad
    Omega = Omega_deg * constants.deg2rad
    TrAn = np.pi / 2 - omega
    EccAn = 2. * np.arctan(np.sqrt((1.0 - ecc) / (1.0 + ecc)) * np.tan(TrAn / 2.0))
    MeAn = EccAn - ecc * np.sin(EccAn)
    MeAn_Tref = mean_long*constants.deg2rad - omega - Omega
    return (MeAn - MeAn_Tref) / (2 * np.pi) * period % period

def kepler_compute_meanlong_from_deltaTc(period, delta_Tc, ecc, omega_deg, Omega_deg=0.0):
    """
    Compute the mean longitude at reference epoch (mean_long) from delta_Tc, i.e., the difference between the first
    *time of inferior conjunction* (Tc) and the reference time Tref, assuming a Keplerian orbit with given period, 
    eccentricity (ecc), argument of pericenter (omega_deg), and longitude of the ascending node (Omega_deg).
    The *time of inferior conjunction* is defined as the time when the planet is closest to the line of sight
    between the observer and the star, i.e. when the true anomaly (TrAn) is equal to pi/2 - omega.  

    :param period: Orbital period of the planet (in days)
    :param delta_Tc: Difference between the time of inferior conjunction and reference time (in days)
    :param ecc: Orbital eccentricity
    :param omega: Argument of pericenter (in degrees)
    :param Omega: Longitude of the ascending node (in degrees, default is 0)
    :return: Mean longitude at reference epoch (in degrees, in the range 0-360)
    """
    omega = omega_deg * constants.deg2rad
    Omega = Omega_deg * constants.deg2rad
    TrAn = np.pi / 2 - omega
    EccAn = 2. * np.arctan(np.sqrt((1.0 - ecc) / (1.0 + ecc)) * np.tan(TrAn / 2.0))
    MeAn = EccAn - ecc * np.sin(EccAn)
    return ((omega + Omega + MeAn - delta_Tc / period * 2 * np.pi) * constants.rad2deg) % 360.

def kepler_compute_deltaTperi_from_deltaTc(period, delta_Tc , ecc, omega_deg):
    """
    Compute difference between the *time of periastron passage* (Tperi) and the reference time (Tref)
    starting from delta_Tc, i.e., the difference between the first *time of inferior conjunction* (Tc) 
    and the reference time (Tref), assuming a Keplerian orbit with given period, 
    eccentricity (ecc), and argument of pericenter (omega_deg).
    :param period: Orbital period of the planet (in days)
    :param delta_Tc: Difference between the time of inferior conjunction and reference time (in days)
    :param ecc: Orbital eccentricity
    :param omega: Argument of pericenter (in degrees)
    :return: minimum positive difference between the time of periastron passage (Tperi) and reference time (Tref)
    """
    TrAn = np.pi / 2 - omega_deg * constants.deg2rad
    EccAn = 2. * np.arctan(np.sqrt((1.0 - ecc) / (1.0 + ecc)) * np.tan(TrAn / 2.0))
    MeAn = EccAn - ecc * np.sin(EccAn)
    return delta_Tc - MeAn/(2*np.pi)*period

def kepler_compute_deltaTperi_from_meanlong(period, mean_long, ecc, omega_deg, Omega_deg=0.0):
    """
    Compute difference between the *time of periastron passage* (Tperi) and the reference time (Tref)
    starting from the mean longitude at reference epoch (mean_long) assuming a Keplerian orbit with given period, 
    eccentricity (ecc), argument of pericenter (omega_deg), and longitude of the ascending node (Omega_deg).
    :param period: Orbital period of the planet (in days)
    :param mean_long: Mean longitude at reference epoch (in degrees)
    :param ecc: Orbital eccentricity
    :param omega: Argument of pericenter (in degrees)
    :param Omega: Longitude of the ascending node (in degrees, default is 0)
    :return: minimum positive difference between the time of periastron passage (Tperi) and reference time (Tref)
    """
    omega = omega_deg * constants.deg2rad
    Omega = Omega_deg * constants.deg2rad
    TrAn = np.pi / 2 - omega
    EccAn = 2. * np.arctan(np.sqrt((1.0 - ecc) / (1.0 + ecc)) * np.tan(TrAn / 2.0))
    MeAn = EccAn - ecc * np.sin(EccAn)
    MeAn_Tref = mean_long*constants.deg2rad - omega - Omega
    Tcent = (MeAn - MeAn_Tref) / (2 * np.pi) * period % period
    return Tcent - MeAn/(2*np.pi)*period

def f_get_mass(mass_secondary, mass_primary, period, ecc, rv_semiamplitude):
    """ Compute the difference between the input radial velocity semi-amplitude
    of the primary star and the value corresponding to the provided orbital parameters.
    Supporting function to *kepler_get_planet_mass* subroutine.

    constants.Gsi: Gravitational constant in SI system [m^3 kg^-1 s^-2]
    constants.Msun: Sun mass in SI system [kg]

    :param mass_secondary: mass of the secondary/planet (in Solar mass units)
    :param mass_primary: mass of the primary star (in Solar mass units)
    :param period: orbital period of the secondary (in days)
    :param ecc: orbital eccentricity of the planet
    :param rv_semiamplitude: observed RV semi-amplitude of the primary (in m/s)
    :return: the difference between the observed and theoretical RV semi-amplitude of the primary (in m/s)
    """

    return rv_semiamplitude - ((2. * np.pi * constants.Gsi * constants.Msun / 86400.0) ** (1. / 3.)
            * (1. / np.sqrt(1. - ecc ** 2.))
            * period ** (-1. / 3.)
            * (mass_secondary * (mass_primary + mass_secondary) ** (-2. / 3.)))

def get_approximate_mass(period, rv_semiamplitude, ecc, mass_primary):
    """ Return the approximate mass of the planet in Solar mass units, under the assumption that
    the mass of the planet is negligible compared to the mass of the star.
    Supporting function to *kepler_get_planet_mass* subroutine.
    For a more precise calculation, decrease the *approximation_limit* value of the function
    *kepler_get_planet_mass* to enable the Newton-Raphson method. 

    :param period: orbital period of the planet (in days)
    :param rv_semiamplitude: observed RV semi-amplitude of the primary (in m/s)
    :param ecc: orbital eccentricity of star2
    :param mass_primary: mass of the primary star (in Solar mass units)
    :return: mass of the planet (in Solar mass units)
    """
    return rv_semiamplitude / ((2. * np.pi * constants.Gsi * constants.Msun / 86400.0) ** (1. / 3.)
         * (1. / np.sqrt(1. - ecc ** 2.))
         * period ** (-1. / 3.)
         * (mass_primary ** (-2. / 3.)))


def kepler_get_planet_mass(period, rv_semiamplitude, ecc, mass_star, approximation_limit=30.):
    """ Compute the mass of the planet in Solar mass units, given the orbital period (in days),
    the observed RV semi-amplitude (in m/s), the orbital eccentricity, and the mass of the primary star (in Solar mass units).
    If the approximate mass of the planet (computed under the assumption that the mass of the planet is negligible compared to the mass of the star)
    is larger than the specified approximation limit, the function will use the more accurate method to compute the planet's mass. 
    The default value for the approximation limit is 30 Earth masses, which is a good compromise between speed and accuracy.
    To force the more accurate method, decrease the *approximation_limit* value.
    To convert the planetary mass,
    ```python
    import pyorbit.subroutines.constants as constants
    mass_planet_earth = mass_planet_solar * constants.Msear
    mass_planet_jupiter = mass_planet_solar * constants.Msjup
    ```
    :param period: orbital period of the planet (in days)
    :param rv_semiamplitude: observed RV semi-amplitude of the primary (in m/s)
    :param ecc: orbital eccentricity of the planet
    :param mass_star: mass of the primary star (in Solar mass units)
    :param approximation_limit: threshold in Earth mass units to switch between the approximate and the more accurate method (default is 30 Earth masses)
    :return: mass of the planet (in Solar mass units)
    """

    n = np.size(rv_semiamplitude)
    if n == 1:
        M_approx = min(get_approximate_mass(period, rv_semiamplitude, ecc, mass_star), 2*constants.Msear)
        return fsolve(f_get_mass, M_approx, args=(mass_star, period, ecc, rv_semiamplitude))

    M_approx = get_approximate_mass(period, rv_semiamplitude, ecc, mass_star)

    if np.average(M_approx) > approximation_limit/constants.Msear:
        print('Computing exact mass of the planet (mean of approximate mass distribution larger than {0:3.1f} Me)'.format(approximation_limit))
        M_init = np.average(M_approx)
        for i in range(0, n):
            M_approx[i] = fsolve(f_get_mass, np.average(M_init), args=(mass_star[i], period[i], ecc[i], rv_semiamplitude[i]))
    else:
        print('Computing planetary mass under the approximation M_planet << M_star (threshold at {0:3.1f} Me)'.format(approximation_limit))

    return M_approx

def kepler_compute_trueanomaly_orbitaldistance(bjd_tref, semimajor_axis, delta_Tc, period, ecc, omega_deg, Omega_deg=0.0):
    """Compute the true anomaly and orbital distance of at the time bjd_tref (BJD-Tref) given its orbital parameters.
    :param bjd_tref: Time difference from reference epoch (in days)
    :param semimajor_axis: Semi-major axis of the planet (either stellar units or AU)
    :param delta_Tc: Difference between the time of inferior conjunction and reference time (in days)
    :param period: Orbital period of the planet (in days)
    :param ecc: Orbital eccentricity
    :param omega_deg: Argument of pericenter (in degrees)
    :param Omega_deg: Longitude of the ascending node (in degrees, default is 0)
    :return: True anomaly (in radians) and orbital distance (same input as semi-major axis) at time bjd_tref
    """

    mean_long = kepler_compute_meanlong_from_deltaTc(period, delta_Tc, ecc, omega_deg, Omega_deg) * constants.deg2rad

    omega = np.asarray(omega_deg, dtype=np.double) * constants.deg2rad
    Omega = np.asarray(Omega_deg, dtype=np.double) * constants.deg2rad
    ecc_array = np.asarray(ecc, dtype=np.double)
    MeAn = (2*np.pi *(bjd_tref / period) + (mean_long - omega - Omega)) % (2 * np.pi)

    if abs(ecc) < 1e-3:
        TrAn = np.asarray(MeAn, dtype=np.double)
        ecc_array = np.asarray(0., dtype=np.double)
        r_orb = semimajor_axis
    else:
        if ecc < 0.:
            ecc_array = -1 * ecc_array
            omega += np.pi

        # Eccentric Anomaly
        EccAn = kepler_E(MeAn, ecc_array)
        TrAn = 2. * np.arctan(np.sqrt((1.0 + ecc_array) / (1.0 - ecc_array)) * np.tan(EccAn / 2.0))
        r_orb = semimajor_axis * (1. - ecc_array ** 2) / (1. + ecc_array * np.cos(TrAn))
    return TrAn, r_orb
