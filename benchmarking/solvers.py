# Analytical Peters solution
import numpy as np
from scipy.constants import G, c
from scipy.integrate import quad, solve_ivp

M_SUN = 1.989e30
AU    = 1.496e11
SEC_PER_YEAR = 365.25 * 24 * 3600

def peters_merger_time(m1, m2, a0, e0):
    """Merger time from Peters (1964) Eq. 5.14.

    ----------
    m1, m2 : float
        Component masses in solar masses.
    a0 : float
        Initial semi-major axis in AU.
    e0 : float
        Initial eccentricity  (0 <= e0 < 1).

    Returns
    -------
    T_merge_yr : float
        Merger time in years.
    """
    m1_kg = m1 * M_SUN
    m2_kg = m2 * M_SUN
    M_kg  = m1_kg + m2_kg
    a0_m  = a0 * AU

    # Peters (1964) Eq. 5.6 prefactor
    beta = (64.0 / 5.0) * G**3 * m1_kg * m2_kg * M_kg / c**5

    if e0 == 0:
        # Circular orbit
        T_s = a0_m**4 / (4.0 * beta)
    else:
        # General case — Peters (1964) Eq. 5.14
        # c0 = a0 * (1 - e0^2) / e0^(12/19) * (1 + 121/304 * e0^2)^(-870/2299)
        c0 = (a0_m * (1.0 - e0**2) / e0**(12.0 / 19.0)
              * (1.0 + (121.0 / 304.0) * e0**2)**(-870.0 / 2299.0))

        def integrand(e):
            return (e**(29.0 / 19.0)
                    * (1.0 + (121.0 / 304.0) * e**2)**(1181.0 / 2299.0)
                    / (1.0 - e**2)**1.5)

        integral, _ = quad(integrand, 0.0, e0)
        T_s = (12.0 / 19.0) * c0**4 / beta * integral

    return T_s / SEC_PER_YEAR


# non-transformed ODE solver

class PetersGW():


    def __init__(self, m1, m2, a0, e0):
        self.m1 = m1 * M_SUN
        self.m2 = m2 * M_SUN
        self.M = self.m1 + self.m2
        self.a0 = a0 * AU
        self.e0 = e0

        # Pre-compute the constant prefactor beta (Peters 1964, Eq. 5.6)
        # beta = (64/5) * G^3 * m1 * m2 * (m1+m2) / c^5
        self.beta = (64.0 / 5.0) * G**3 * self.m1 * self.m2 * self.M / c**5

        self._solution = None


    def _derivs(self, t, y):
        """Right-hand side of Peters' ODEs.

        Parameters
        ----------
        t : float
            Time (seconds).
        y : array_like
            State vector [a, e] where a is in metres and e is dimensionless.

        Returns
        -------
        dydt : list
            [da/dt, de/dt]
        """
        a, e = y

        # Guard against unphysical values during integration
        if a <= 0 or e < 0 or e >= 1:
            return [0.0, 0.0]

        e2 = e * e
        one_minus_e2 = 1.0 - e2

        # da/dt  (Peters 1964, Eq. 5.6)
        dadt = (-self.beta / (a**3 * one_minus_e2**(3.5))
                * (1.0 + (73.0 / 24.0) * e2 + (37.0 / 96.0) * e2**2))

        # de/dt  (Peters 1964, Eq. 5.7)
        # coefficient: (304/15) / (64/5) = 19/12
        dedt = (-self.beta * (19.0 / 12.0) * e
                / (a**4 * one_minus_e2**(2.5))
                * (1.0 + (121.0 / 304.0) * e2))

        return [dadt, dedt]


    @staticmethod
    def _merger_event(t, y):
        """Event function: triggers when the semi-major axis reaches a
        small fraction of its initial value (effectively merger)."""
        a = y[0]
        return a  # Crossing zero from above → merger

    _merger_event.terminal = True
    _merger_event.direction = -1

    def integrate(self, t_max=14e9, rtol=1e-10, atol=1e-12,
                  max_step=np.inf, dense_output=True):
        """Integrate the orbital evolution until merger or *t_max*.

        Parameters
        ----------
        t_max : float, optional
            Maximum integration time in **years** (default: 14 Gyr).
        rtol, atol : float
            Relative / absolute tolerances passed to ``solve_ivp``.
        max_step : float
            Maximum internal step size in seconds.
        dense_output : bool
            Whether to request dense (interpolated) output.

        Returns
        -------
        sol : OdeResult
            The ``scipy.integrate.solve_ivp`` result object.
        """

        y0 = [self.a0, self.e0]
        t_max_s = t_max * SEC_PER_YEAR

        sol = solve_ivp(
            self._derivs,
            (0.0, t_max_s),
            y0,
            method="RK45",
            events=self._merger_event,
            rtol=rtol,
            atol=atol,
            max_step=max_step,
            dense_output=dense_output,
        )
        self._solution = sol
        return sol

    @property
    def merger_time_yr(self):
        return self._solution.t[-1] / SEC_PER_YEAR
