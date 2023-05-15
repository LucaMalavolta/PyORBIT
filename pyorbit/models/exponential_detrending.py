from pyorbit.subroutines.common import np
from pyorbit.models.abstract_model import AbstractModel
from pyorbit.models.polynomial_detrending import PolynomialDetrending

class ExponentialDetrending(PolynomialDetrending):

    def __init__(self, *args, **kwargs):
        AbstractModel.__init__(self, *args, **kwargs)
        PolynomialDetrending.__init__(self, *args, **kwargs)
        self.exponential_detrending = True