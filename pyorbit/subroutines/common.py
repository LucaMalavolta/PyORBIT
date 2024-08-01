from __future__ import print_function
import os
import sys
from scipy import stats
from scipy.interpolate import splrep, splev
import pyorbit.subroutines.constants as constants

import numpy as np
from ordered_set import OrderedSet
np.seterr(invalid='ignore')


# old base 2 logarithm


def get_var_log(var, fix, i):
    if len(np.shape(var)) == 1:
        return np.log2(var[i], dtype=np.double)
    else:
        return np.log2(var[:, i], dtype=np.double)


def get_var_exp(var, fix, i):
    if len(np.shape(var)) == 1:
        return np.exp2(var[i], dtype=np.double)
    else:
        return np.exp2(var[:, i], dtype=np.double)


def get_var_log_base2(var, fix, i):
    if len(np.shape(var)) == 1:
        return np.log2(var[i], dtype=np.double)
    else:
        return np.log2(var[:, i], dtype=np.double)


def get_var_exp_base2(var, fix, i):
    if len(np.shape(var)) == 1:
        return np.exp2(var[i], dtype=np.double)
    else:
        return np.exp2(var[:, i], dtype=np.double)


def get_var_log_base10(var, fix, i):
    if len(np.shape(var)) == 1:
        return np.log10(var[i], dtype=np.double)
    else:
        return np.log10(var[:, i], dtype=np.double)


def get_var_exp_base10(var, fix, i):
    if len(np.shape(var)) == 1:
        return 10**(var[i])
    else:
        return 10**(var[:, i])


def get_var_log_natural(var, fix, i):
    if len(np.shape(var)) == 1:
        return np.log(var[i], dtype=np.double)
    else:
        return np.log(var[:, i], dtype=np.double)


def get_var_exp_natural(var, fix, i):
    if len(np.shape(var)) == 1:
        return np.exp(var[i], dtype=np.double)
    else:
        return np.exp(var[:, i], dtype=np.double)

def get_var_sin(var, fix, i):
    if len(np.shape(var)) == 1:
        return np.sin(var[i] * constants.deg2rad, dtype=np.double)
    else:
        return np.sin(var[:, i]* constants.deg2rad, dtype=np.double)


def get_var_arcsine(var, fix, i):
    if len(np.shape(var)) == 1:
        return np.arcsin(var[i], dtype=np.double)*constants.rad2deg
    else:
        return np.arcsin(var[:, i], dtype=np.double)*constants.rad2deg


def get_var_sine(var, fix, i):
    if len(np.shape(var)) == 1:
        return np.sin(var[i] * constants.deg2rad, dtype=np.double)
    else:
        return np.sin(var[:, i]* constants.deg2rad, dtype=np.double)


def get_var_arccosine(var, fix, i):
    if len(np.shape(var)) == 1:
        return np.arccos(var[i], dtype=np.double)*constants.rad2deg
    else:
        return np.arccos(var[:, i], dtype=np.double)*constants.rad2deg

def get_var_cosine(var, fix, i):
    if len(np.shape(var)) == 1:
        return np.cos(var[i] * constants.deg2rad, dtype=np.double)
    else:
        return np.cos(var[:, i]* constants.deg2rad, dtype=np.double)


def get_var_val(var, fix, i):
    if len(np.shape(var)) == 1:
        return var[i]
    else:
        return var[:, i]


def get_fix_val(var, fix, i):
    if len(np.shape(fix)) == 1:
        return fix[i]
    else:
        return fix[:, i]


def get_2var_sre(var, fix, i):
    if len(np.shape(var)) == 1:
        ecoso = var[i[0]]
        esino = var[i[1]]
    else:
        ecoso = var[:, i[0]]
        esino = var[:, i[1]]
    return np.square(ecoso, dtype=np.double) + np.square(esino, dtype=np.double)


def get_2var_e(var, fix, i):
    if len(np.shape(var)) == 1:
        ecoso = var[i[0]]
        esino = var[i[1]]
    else:
        ecoso = var[:, i[0]]
        esino = var[:, i[1]]
    return np.sqrt(np.square(ecoso, dtype=np.double) + np.square(esino, dtype=np.double))


def get_2var_o(var, fix, i):
    if len(np.shape(var)) == 1:
        ecoso = var[i[0]]
        esino = var[i[1]]
    else:
        ecoso = var[:, i[0]]
        esino = var[:, i[1]]
    return np.arctan2(esino, ecoso, dtype=np.double) * constants.rad2deg


# Kipping transformation for quadratic limb darkening coefficients
def get_2var_c1(var, fix, i):
    if len(np.shape(var)) == 1:
        q1 = var[i[0]]
        q2 = var[i[1]]
    else:
        q1 = var[:, i[0]]
        q2 = var[:, i[1]]
    return 2.0*np.sqrt(q1, dtype=np.double) * q2


def get_2var_c2(var, fix, i):
    if len(np.shape(var)) == 1:
        q1 = var[i[0]]
        q2 = var[i[1]]
    else:
        q1 = var[:, i[0]]
        q2 = var[:, i[1]]
    return np.sqrt(q1, dtype=np.double) * (1.0 - 2.0*q2)


def get_2darray_from_val(val):
    out = np.zeros(2, dtype=np.double)
    try:
        if np.shape(val) == 1:
            out[0] = val[0]
        else:
            out[:] = val[0:2]
    except:
        try:
            out[0] = val
        except:
            out = val
    return out


def get_2var_rho(var, fix, i):
    if len(np.shape(var)) == 1:
        mass = var[i[0]]
        radius = var[i[1]]
    else:
        mass = var[:, i[0]]
        radius = var[:, i[1]]
    return mass/radius**3

def get_2var_radius(var, fix, i):
    if len(np.shape(var)) == 1:
        mass = var[i[0]]
        density = var[i[1]]
    else:
        mass = var[:, i[0]]
        density = var[:, i[1]]
    return np.power(mass/density, 1./3.)


def get_2var_mass(var, fix, i):
    if len(np.shape(var)) == 1:
        density = var[i[0]]
        radius = var[i[1]]
    else:
        density = var[:, i[0]]
        radius = var[:, i[1]]
    return density*radius**3


def get_2var_product(var, fix, i):
    if len(np.shape(var)) == 1:
        out = var[i[0]] * var[i[1]]
    else:
        out = var[:, i[0]] * var[:, i[1]]
    return out


def get_2var_veq_istar_vsini(var, fix, i):
    # first parameter: v_eq
    # second parameter: i_star
    # output: vsini
    if len(np.shape(var)) == 1:
        out = var[i[0]] * np.sin(var[i[1]] * constants.deg2rad)
    else:
        out = var[:, i[0]] * np.sin(var[:, i[1]] * constants.deg2rad)
    return out

def get_2var_veq_cosistar_vsini(var, fix, i):
    # first parameter: v_eq
    # second parameter: cosi_star
    # output: vsini
    if len(np.shape(var)) == 1:
        istar = np.arccos(var[i[1]])
        out = var[i[0]] * np.sin(istar)
    else:
        istar = np.arccos(var[:, i[1]])
        out = var[:, i[0]] * np.sin(istar)
    return out

def get_2var_veq_rot_radius(var, fix, i):
    # first parameter: v_eq
    # second parameter: p_rot
    # R = v*P / 2pi
    # R[R_sun] = v_eq[km/s] * (P[d]* 86400[s/d] ) / (2*np.pi) / R_sun[km]
    if len(np.shape(var)) == 1:
        out = var[i[0]] * (var[i[1]] * constants.d2s / (2*np.pi) / constants.Rsun)
    else:
        out = var[:, i[0]] * (var[:, i[1]] * constants.d2s / (2*np.pi) / constants.Rsun)
    return out

def get_2var_veq_radius_rot(var, fix, i):
    # first parameter: v_eq
    # second parameter: radius
    # P = 2pi R / v
    # P[d]  = (2*np.pi) R[R_sun]*R_sun[km] /  v_eq[km/s] /  86400[s/d]
    if len(np.shape(var)) == 1:

        out = (2*np.pi) * var[i[1]]  / var[i[0]] * (constants.Rsun /constants.d2s)
        out = var[i[0]] * (var[i[1]] * constants.d2s / (2*np.pi) / constants.Rsun)
    else:
        out = (2*np.pi) * var[:, i[1]]  / var[:, i[0]] * (constants.Rsun /constants.d2s)
    return out




def get_2var_prot_rstar_veq(var, fix, i):
    #first parameter: p_rot
    #second parameter: r_star
    #  v = 2piR / P
    #  v[km/s] = 2pi R[km] / P[d]
    if len(np.shape(var)) == 1:
        out =  2* np.pi * var[i[1]] * constants.Rsun / (var[i[0]] * constants.d2s )
    else:
        out =  2* np.pi * var[:,i[1]] * constants.Rsun / (var[:,i[0]] * constants.d2s )
    return out


def get_3var_prot_rstar_istar_veq(var, fix, i):
    #first parameter: p_rot
    #second parameter: r_star
    #third parameter: i_star
    #  v = 2piR / P
    #  v[km/s] = 2pi R[km] / P[d]
    if len(np.shape(var)) == 1:
        out =  2* np.pi * var[i[1]] * constants.Rsun / (var[i[0]] * constants.d2s ) * np.sin(var[i[2]] * constants.deg2rad)
    else:
        out =  2* np.pi * var[:,i[1]] * constants.Rsun / (var[:,i[0]] * constants.d2s ) *np.sin(var[:, i[2]] * constants.deg2rad)
    return out

def get_3var_prot_rstar_cosistar_veq(var, fix, i):
    #first parameter: p_rot
    #second parameter: r_star
    #third parameter: i_star
    #  v = 2piR / P
    #  v[km/s] = 2pi R[km] / P[d]
    if len(np.shape(var)) == 1:
        i_star = np.arccos(var[i[2]])
        out =  2* np.pi * var[i[1]] * constants.Rsun / (var[i[0]] * constants.d2s ) * np.sin(i_star)
    else:
        i_star = np.arccos(var[:, i[2]])
        out =  2* np.pi * var[:,i[1]] * constants.Rsun / (var[:,i[0]] * constants.d2s ) *np.sin(i_star)
    return out




def giveback_priors(kind, bounds, pams, val):
    """ The code is supposed to give -np.inf log-likelihood when the
    parameters are outside the boundaries,
    so this case is not encompassed in the definition of the priors """

    if kind == 'File':
        """ Special case where 'pams' is actually a KDE """
        return np.log(pams.pdf(val)[0])

    if kind == 'None':
        return 0.00

    if kind == 'Gaussian':
        # return np.log(stats.norm.pdf(val, loc=pams[0], scale=pams[1]))
        return -(val - pams[0]) ** 2 / (2 * pams[1] ** 2) - 0.5 * np.log(2*np.pi) - np.log(pams[1])

    if kind == 'HalfGaussian' or kind == 'PositiveHalfGaussian':
        # return -(val - pams[0]) ** 2 / (2 * pams[1] ** 2)
        if val >= pams[0]:
            return np.log(2*stats.norm.pdf(val, loc=pams[0], scale=pams[1]))
        else:
            return -np.inf

    if kind == 'NegativeHalfGaussian':
        # return -(val - pams[0]) ** 2 / (2 * pams[1] ** 2)
        if val <= pams[0]:
            return np.log(2*stats.norm.pdf(val, loc=pams[0], scale=pams[1]))
        else:
            return -np.inf

    if kind == 'Uniform':
        return np.log(1./(bounds[1]-bounds[0]))

    if kind == 'Jeffreys' or kind == 'TruncatedJeffreys':
        return np.log(1./(val*np.log(bounds[1]/bounds[0])))

    if kind == 'ModifiedJeffreys' or kind == 'TruncatedModifiedJeffreys':
        """ Used for the semi-amplitude of the RV curve
            bounds[1] = Kmax (suggested 999 m/s)
            pams[0] = K_0 (suggested 1 m/s)
        """
        return np.log(1./(pams[0]*(1. + val/pams[0])) * 1./np.log(1.+bounds[1]/pams[0]))

    if kind == 'TruncatedRayleigh':
        """ bounds[1] = e_max
            pams[0] = sigma_2 (suggested 0.2)
        """
        return np.log((val/pams[0]**2 * np.exp(-0.5 * (val/pams[0])**2))/(1. - np.exp(-0.5*(bounds[1]/pams[0])**2)))

    if kind == "WhiteNoisePrior":
        """ bounds[1] = noise_max (99 m/s)
            pams[0] = noise_0 (suggested 1 m/s)
        """
        return np.log(1./(pams[0]*(1.0 + val/pams[0])) * 1.0/np.log(1.0+bounds[1]/pams[0]))

    if kind == 'BetaDistribution' or kind == 'Beta':
        return np.log(stats.beta.pdf((val-bounds[0])/(bounds[1]-bounds[0]), pams[0], pams[1]))


"""
DEPRECATED
Special subroutine to transform MultiNest/PolyChord priors, i.e., trasnform the datacube from [0:1] to physical
values while taking into account the priors

def nested_sampling_prior_transformation(kind, bounds, pams):

    #x = np.linspace(0.000000, 1.000000, num=1001, endpoint=True)
    x_var = np.linspace(0.000000, 1.000000, num=10001, endpoint=True, dtype=np.double)*(bounds[1]-bounds[0]) + bounds[0]
    area = np.zeros(len(x_var), dtype=np.double)

    for x_num, x_val in enumerate(x_var):
        area[x_num:] += np.exp(giveback_priors(kind, bounds, pams, x_val)) * (1. / 10000.) + 0.000000000001
    area[0] = 0
    area /= area[-1]

    return interp1d(area, x_var, kind='cubic')
"""


def nested_sampling_prior_prepare(kind, bounds, pams, space):
    """
    This subroutine computes the coefficient of the spline interpolation of
    the inverse cumulative function
    In some special cases, it returns the parameters required by the intrinsic
    function, e.g. scipi.stats.norm.icf according to their implementation in
    nested_sampling_prior_compute()

    :param kind: type of prior
    :param bounds: list/array with lower and upper limits for parameter exploration
    :param pams: parameters relative to the prior probability function
    :return:
    """

    if kind == 'Uniform':
        return bounds

    if kind in ['Gaussian', 'HalfGaussian',  'PositiveHalfGaussian', 'NegativeHalfGaussian','BetaDistribution', 'Beta']:
        return pams

    """ All the following priors are defined only if the variable is sampled in the Natural space"""
    if space != 'Linear':
        print()
        print(' *** ERROR in the YAML file ***')
        print(' You are using a prior that is not supported in a non-Linear sampling space')
        print(' add this keyword in the YAML file for each parameter not sampled in the ')
        print(' Linear space and with a prior other than Uniform or Gaussian')
        print('   spaces: ')
        print('       pam: Linear')
        print()
        quit()

    x_var = np.linspace(0.000000, 1.000000, num=10001, endpoint=True, dtype=np.double)*(bounds[1]-bounds[0]) + bounds[0]
    area = np.zeros(len(x_var), dtype=np.double)

    for x_num, x_val in enumerate(x_var):
        area[x_num:] += np.exp(giveback_priors(kind, bounds, pams, x_val)) * (1. / 10000.) + 0.000000000001
    area[0] = 0
    area /= area[-1]

    return splrep(area, x_var)


def nested_sampling_prior_compute(val, kind, coeff, space):
    """
    In same cases ()

    :param val:
    :param kind:
    :param coeff:
    :param space:
    :return:
    """

    if kind == 'Uniform':
        return val * (coeff[1] - coeff[0]) + coeff[0]

    if kind == 'Gaussian':
        x_new = stats.norm.ppf(val, coeff[0], coeff[1])

        if space == 'Linear':
            return x_new
        elif space in ['Log_Base2', 'Logarithmic']:
            return np.log2(x_new)
        elif space == 'Log_Base10':
            return np.log10(x_new)
        elif space == 'Log_Natural':
            return np.log(x_new)
        elif space == 'Sine_Angle':
            return np.sin(x_new*constants.rad2deg)
        elif space == 'Cosine_Angle':
            return np.cos(x_new*constants.rad2deg)

    if kind in ['HalfGaussian', 'PositiveHalfGaussian']:

        x_new = stats.halfnorm.ppf(val, coeff[0], coeff[1])

        if space == 'Linear':
            return x_new
        elif space in ['Log_Base2', 'Logarithmic']:
            return np.log2(x_new)
        elif space == 'Log_Base10':
            return np.log10(x_new)
        elif space == 'Log_Natural':
            return np.log(x_new)
        elif space == 'Sine_Angle':
            return np.sin(x_new*constants.rad2deg)
        elif space == 'Cosine_Angle':
            return np.cos(x_new*constants.rad2deg)

    if kind in ['NegativeHalfGaussian']:

        x_new = coeff[0] - stats.halfnorm.ppf(val, 0, coeff[1])

        if space == 'Linear':
            return x_new
        elif space in ['Log_Base2', 'Logarithmic']:
            return np.log2(x_new)
        elif space == 'Log_Base10':
            return np.log10(x_new)
        elif space == 'Log_Natural':
            return np.log(x_new)
        elif space == 'Sine_Angle':
            return np.sin(x_new*constants.rad2deg)
        elif space == 'Cosine_Angle':
            return np.cos(x_new*constants.rad2deg)

    if kind in ['BetaDistribution', 'Beta']:
        x_new = stats.beta.ppf(val, coeff[0], coeff[1])

        if space == 'Linear':
            return x_new
        elif space in ['Log_Base2', 'Logarithmic']:
            return np.log2(x_new)
        elif space == 'Log_Base10':
            return np.log10(x_new)
        elif space == 'Log_Natural':
            return np.log(x_new)
        elif space == 'Sine_Angle':
            return np.sin(x_new*constants.rad2deg)
        elif space == 'Cosine_Angle':
            return np.cos(x_new*constants.rad2deg)

    return splev(val, coeff)


def compute_value_sigma(samples):
    if np.size(np.shape(samples)) == 1:
        sample_med = np.zeros(3)
        sample_tmp = np.percentile(samples, [15.865, 50, 84.135], axis=0)
        sample_med[0] = sample_tmp[1]
        sample_med[1] = sample_tmp[2] - sample_tmp[1]
        sample_med[2] = sample_tmp[0] - sample_tmp[1]

    elif np.size(np.shape(samples)) == 2:
        sample_med = np.asarray(
            list(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]), zip(*np.percentile(samples, [15.865, 50, 84.135], axis=0)))))

    else:
        print('ERROR!!! ')
        return None
    return sample_med


def pick_MAP_parameters(samples, lnprob):

    indmax = np.argmax(lnprob)
    if np.size(np.shape(samples)) == 1:
        return samples[indmax], lnprob[indmax]
    elif np.size(np.shape(samples)) == 2:
        return samples[indmax, :], lnprob[indmax]
    else:
        print('ERROR!!! ')
        return None


def pick_sampleMED_parameters(samples, lnprob):
    pams_perc = np.percentile(samples, [15.865, 50, 84.135], axis=0)
    dist_array = 0.0
    for i_pam in range(0, len(samples[0, :])):
        # dist_array += ((samples[:, i_pam]-pams_perc[1, i_pam])/(pams_perc[2,i_pam]-pams_perc[0,i_pam]))**2
        dist_array += ((samples[:, i_pam]-pams_perc[1, i_pam])/(np.ptp(samples[:, i_pam])))**2

    id_val = np.argmin(dist_array)
    return samples[id_val, :], lnprob[id_val]


def convert_rho_to_ars(P, rho):

    return np.power(constants.Gsi * (constants.d2s * constants.d2s) * (P**2)
                    * rho * constants.rho_Sun / (3. * np.pi), 1./3.)


def convert_ars_to_a(a, Rs):
    return a * Rs * constants.RsunAU


def convert_PMsMp_to_a(P, Ms, Mp):
    # planet mass in solar masses
    # output in Astronomical Units
    return np.power(P**2 * constants.Giau * (Ms + Mp/constants.Msear) / (2 * np.pi)**2., 1./3.)


def convert_b_to_i(b, e, o, a):
    o_rad  = o * constants.deg2rad

    rho_e = (1. - e ** 2) / (1. + e * np.sin(o_rad))
    arccos_argument = b / a / rho_e
    if np.size(arccos_argument) <= 1:

        if arccos_argument > 1.:
            arccos_argument = 1
        if arccos_argument < -1.:
            arccos_argument = -1
    else:
        arccos_argument = [1. if b > 1. else b for b in arccos_argument]
        arccos_argument = [-1. if b < -1. else b for b in arccos_argument]

    return np.arccos(arccos_argument)*180./np.pi


def convert_RTaAU_to_insol(R, Ts, a_abs):
    return R**2 * (Ts/constants.Sun_temperature)**4 / a_abs**2 * constants.Sun_constant
