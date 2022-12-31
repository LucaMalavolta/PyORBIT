from pyorbit.subroutines.common import *
from pyorbit.classes.model_container_abstract import ModelContainer


class ModelContainerOptimize(ModelContainer):

    def __init__(self):
        super(self.__class__, self).__init__()

        self.optimize_dir_output = None
        self.optimize_parameters = {'method': 'Nelder-Mead',
                                    'maxiter': 1000000.0,
                                    'adaptive': True,
                                    'disp': True}

        self.current_neg_loglikelihood = 10000000.0
        """
        Check SciPy documentation for available methods:
        https://docs.scipy.org/doc/scipy/reference/optimize.html
        """

    def negative_log_priors_likelihood(self, theta):

        neg_loglikelihood = -self.__call__(theta)

        # workaround to avoid RuntimeWarning when the one of the variable is outside the boundaries
        if np.isfinite(neg_loglikelihood):
            self.current_neg_loglikelihood = neg_loglikelihood
            return neg_loglikelihood
        else:
            return self.current_neg_loglikelihood*10.

    """ 
    def check_bounds(self, theta):
        print(type(theta))
        for ii in range(0, self.ndim):

            print(theta[ii])
            #if not (self.bounds[ii, 0] < theta[ii] < self.bounds[ii, 1]):
            #    return False
    """
