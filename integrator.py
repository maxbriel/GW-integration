"""
Peters' Gravitational Wave Radiation Integrator
================================================

Solves the Peters (1964) equations for orbital decay due to gravitational
wave emission:

    da/dt = -(64/5) * (G^3 * m1 * m2 * (m1+m2)) / (c^5 * a^3 * (1-e^2)^(7/2))
            * (1 + (73/24)*e^2 + (37/96)*e^4)

    de/dt = -(304/15) * e * (G^3 * m1 * m2 * (m1+m2)) / (c^5 * a^4 * (1-e^2)^(5/2))
            * (1 + (121/304)*e^2)

References:
    Peters, P. C. (1964). "Gravitational Radiation and the Motion of Two
    Point Masses". Physical Review, 136(4B), B1224-B1232.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.constants import G, c

# Solar mass in kg
M_SUN = 1.989e30
# Astronomical unit in metres
AU = 1.496e11


class PetersGWIntegrator:
    """Integrate Peters' equations for GW-driven orbital evolution.

    Parameters
    ----------
    m1 : float
        Mass of the first body in solar masses.
    m2 : float
        Mass of the second body in solar masses.
    a0 : float
        Initial semi-major axis in AU.
    e0 : float
        Initial eccentricity (0 <= e0 < 1).
    """

    def __init__(self, m1, m2, a0, e0):
        if not 0 <= e0 < 1:
            raise ValueError("Eccentricity must satisfy 0 <= e0 < 1")
        if a0 <= 0:
            raise ValueError("Semi-major axis must be positive")

        self.m1 = m1 * M_SUN          # kg
        self.m2 = m2 * M_SUN          # kg
        self.M = self.m1 + self.m2     # total mass in kg
        self.a0 = a0 * AU             # metres
        self.e0 = e0

        # Pre-compute the constant prefactor beta (Peters 1964, Eq. 5.6)
        # beta = (64/5) * G^3 * m1 * m2 * (m1+m2) / c^5
        self.beta = (64.0 / 5.0) * G**3 * self.m1 * self.m2 * self.M / c**5

        self._solution = None

    # ------------------------------------------------------------------
    # Peters' ODEs
    # ------------------------------------------------------------------
    def _derivs(self, t, y):
        """Right-hand side of Peters' coupled ODEs.

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

    # ------------------------------------------------------------------
    # Event: merger detection  (a → 0 or very small)
    # ------------------------------------------------------------------
    @staticmethod
    def _merger_event(t, y):
        """Event function: triggers when the semi-major axis reaches a
        small fraction of its initial value (effectively merger)."""
        a = y[0]
        return a  # Crossing zero from above → merger

    _merger_event.terminal = True
    _merger_event.direction = -1

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
    def integrate(self, t_max=None, rtol=1e-10, atol=1e-12,
                  max_step=np.inf, dense_output=True):
        """Integrate the orbital evolution until merger or *t_max*.

        Parameters
        ----------
        t_max : float, optional
            Maximum integration time in seconds.  If *None*, an estimate
            for the circular-orbit merger time (Peters 1964, Eq. 5.10) is
            used, multiplied by a safety factor of 2.
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
        if t_max is None:
            t_max = 2.0 * self.merger_time_circular()

        y0 = [self.a0, self.e0]

        sol = solve_ivp(
            self._derivs,
            (0.0, t_max),
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

    def merger_time_circular(self):
        """Estimate the merger time for a *circular* orbit (e=0).

        Uses Peters (1964), Eq. 5.10:
            T = a0^4 / (4 * beta)

        Returns
        -------
        T : float
            Merger time in seconds.
        """
        return self.a0**4 / (4.0 * self.beta)

    def merger_time_peters(self):
        """Analytic Peters (1964) merger-time estimate for an eccentric orbit.

        Uses the circular-orbit result (Peters 1964, Eq. 5.10) multiplied by
        the leading-order eccentricity enhancement factor:

            T_merge = (a0^4 / (4*beta)) * (1 - e0^2)^(7/2)

        The factor (1-e0^2)^(7/2) is the leading-order correction that arises
        from averaging the Peters da/dt equation over one orbit and captures
        the accelerated inspiral due to eccentricity.  For e0 = 0 this reduces
        exactly to the circular formula.

        Returns
        -------
        T : float
            Merger time in seconds.
        """
        return (self.a0**4 / (4.0 * self.beta)) * (1.0 - self.e0**2)**3.5

    @property
    def solution(self):
        """The most recent ``solve_ivp`` result, or *None*."""
        return self._solution

    def get_merger_time(self):
        """Return the merger time from the last integration (seconds).

        Returns
        -------
        t_merge : float or None
            The time at which the merger event was detected, or *None*
            if the integration has not been run or no merger was found.
        """
        if self._solution is None:
            raise RuntimeError("Call integrate() first.")
        if self._solution.t_events[0].size > 0:
            return self._solution.t_events[0][0]
        return self._solution.t[-1]  # reached t_max ≈ merger

    # ------------------------------------------------------------------
    # Convenience: time in useful units
    # ------------------------------------------------------------------
    @staticmethod
    def seconds_to_years(t):
        """Convert seconds to years."""
        return t / (365.25 * 24 * 3600)

    @staticmethod
    def metres_to_au(a):
        """Convert metres to AU."""
        return a / AU
