from pyorbit.classes.common import \
    get_var_val,\
    get_2var_c1,\
    get_2var_c2,\
    nested_sampling_prior_prepare

from pyorbit.models.abstract_common import AbstractCommon

from pyorbit.models.limb_darkening import LimbDarkening_1Pam, LimbDarkening_2Pam, LimbDarkening_4Pam


class Batman_LimbDarkening_Linear(LimbDarkening_1Pam):
    model_class = 'batman_ld_linear'
    ld_type = 'linear'


class Batman_LimbDarkening_Quadratic(LimbDarkening_2Pam):
    model_class = 'batman_ld_quadratic'
    ld_type = 'quadratic'


class Batman_LimbDarkening_SquareRoot(LimbDarkening_2Pam):
    model_class = 'batman_ld_square-root'
    ld_type = 'square-root'


class Batman_LimbDarkening_Logarithmic(LimbDarkening_2Pam):
    model_class = 'batman_ld_logarithmic'
    ld_type = 'logarithmic'


class Batman_LimbDarkening_Exponential(LimbDarkening_2Pam):
    model_class = 'batman_ld_exponential'
    ld_type = 'exponential'


class Batman_LimbDarkening_Power2(LimbDarkening_1Pam):
    model_class = 'batman_ld_power2'
    ld_type = 'power2'


class Batman_LimbDarkening_NonLinear(LimbDarkening_4Pam):
    model_class = 'batman_ld_nonlinear'
    ld_type = 'nonlinear'
