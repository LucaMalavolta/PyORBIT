from pyorbit.subroutines.common import np, OrderedSet

try:
    import celerite
    from celerite.terms import Term



    class SHOTerm(Term):
        r"""
        A term representing a stochastically-driven, damped harmonic oscillator
        As in celertie, but accpeting the physical parameters instead of their
        logarithm

        The PSD of this term is
        .. math::
            S(\omega) = \sqrt{\frac{2}{\pi}} \frac{S_0\,\omega_0^4}
            {(\omega^2-{\omega_0}^2)^2 + {\omega_0}^2\,\omega^2/Q^2}
        with the parameters ``log_S0``, ``log_Q``, and ``log_omega0``.
        Args:
            S0 (float): parameter :math:`S_0`.
            Q (float): parameter :math:`Q`.
            omega0 (float): parameter :math:`\omega_0`.
        """

        parameter_names = ("S0", "Q", "w0")

        def __repr__(self):
            return "SHOTerm({0.S0}, {0.Q}, {0.w0})".format(self)

        def get_real_coefficients(self, params):
            S0, Q, w0 = params
            if Q >= 0.5:
                return np.empty(0), np.empty(0)

            f = np.sqrt(1.0 - 4.0 * Q**2)
            return (
                0.5*S0*w0*Q*np.array([1.0+1.0/f, 1.0-1.0/f]),
                0.5*w0/Q*np.array([1.0-f, 1.0+f])
            )

        def get_complex_coefficients(self, params):
            S0, Q, w0 = params
            if Q < 0.5:
                return np.empty(0), np.empty(0), np.empty(0), np.empty(0)

            f = np.sqrt(4.0 * Q**2-1)
            return (
                S0 * w0 * Q,
                S0 * w0 * Q / f,
                0.5 * w0 / Q,
                0.5 * w0 / Q * f,
            )

except (ModuleNotFoundError,ImportError):
    pass