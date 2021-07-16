
from pyorbit.models.limb_darkening import LimbDarkening_1Pam, LimbDarkening_2Pam, LimbDarkening_4Pam


#class PyTransit_LimbDarkening_Linear(LimbDarkening_1Pam):
#    model_class = 'pytransit_ld_linear'
#    ld_type = 'linear'


class PyTransit_LimbDarkening_Quadratic(LimbDarkening_2Pam):
    model_class = 'pytransit_ld_quadratic'
    ld_type = 'quadratic'


#class PyTransit_LimbDarkening_SquareRoot(LimbDarkening_2Pam):
#    model_class = 'pytransit_ld_square-root'
#    ld_type = 'square-root'


#class PyTransit_LimbDarkening_Logarithmic(LimbDarkening_2Pam):
#    model_class = 'pytransit_ld_logarithmic'
#    ld_type = 'logarithmic'


#class PyTransit_LimbDarkening_Exponential(LimbDarkening_2Pam):
#    model_class = 'pytransit_ld_exponential'
#    ld_type = 'exponential'


class PyTransit_LimbDarkening_Power2(LimbDarkening_1Pam):
    model_class = 'pytransit_ld_power2'
    ld_type = 'power2'


#class PyTransit_LimbDarkening_NonLinear(LimbDarkening_4Pam):
#    model_class = 'pytransit_ld_nonlinear'
#    ld_type = 'nonlinear'
