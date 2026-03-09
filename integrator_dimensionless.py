"""
Dimensionless Peters (1964) GW Integrator — log-separation variable
====================================================================

Re-formulates the Peters orbital-decay equations in dimensionless
variables with the log-orbital-separation  s = -ln(alpha)  as the
independent variable, where  alpha = a / a0.

Definitions
-----------
    alpha = a / a0           dimensionless semi-major axis
    τ = t / t0               dimensionless time
    s = -ln(alpha)           log-separation (independent variable)
    t0 = a0^4 / (4 beta)     natural timescale  (circular-orbit merger time)
    beta  = (64/5) G^3 m1 m2 (m1+m2) / c^5

(Following Andrews+19)

ODE system (state vector y = [tau, e], independent variable s)
-------------------------------------------------------------
    dtau/ds =  4 exp(-4s) (1 - e^2)^(7/2) / g(e)
    de/ds = -(19/12) e (1 - e^2) h(e) / g(e)

where
    g(e) = 1 + (73/24) e^2 + (37/96) e^4
    h(e) = 1 + (121/304) e^2

Note that de/ds depends *only* on e — the eccentricity evolution as a
function of log-separation is universal (mass- and scale-independent).

Integration starts at  s = 0  (alpha = 1, i.e. a = a0)  and proceeds to
large s  (alpha -> 0, merger).

References
----------
    Peters, P. C. (1964). "Gravitational Radiation and the Motion of Two
    Point Masses". Physical Review, 136(4B), B1224-B1232.


"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.constants import G, c

# ---------- physical constants in SI ----------
M_SUN = 1.989e30   # solar mass [kg]
AU = 1.496e11       # astronomical unit [m]

# ----------- eccentricity factors -------------
def _g(e2):
    """g(e) = 1 + (73/24) e^2 + (37/96) e^4"""
    return 1.0 + (73.0 / 24.0) * e2 + (37.0 / 96.0) * e2 * e2


def _h(e2):
    """h(e) = 1 + (121/304) e^2"""
    return 1.0 + (121.0 / 304.0) * e2


# ==================================================================
class PetersGWDimensionless:
    """Dimensionless Peters (1964) integrator using s = -ln(a/a0).

    Parameters
    ----------
    m1 : float
        Mass of the first body in solar masses.
    m2 : float
        Mass of the second body in solar masses.
    a0 : float
        Initial semi-major axis in AU.
    e0 : float
        Initial eccentricity (0 ≤ e0 < 1).
    """

    def __init__(self, m1, m2, a0, e0):
        if not 0 <= e0 < 1:
            raise ValueError("Eccentricity must satisfy 0 <= e0 < 1")
        if a0 <= 0:
            raise ValueError("Semi-major axis must be positive")

        # Store inputs in SI
        self.m1_solar = m1
        self.m2_solar = m2
        self.a0_au = a0
        self.e0 = e0

        self.m1 = m1 * M_SUN          # [kg]
        self.m2 = m2 * M_SUN          # [kg]
        self.M = self.m1 + self.m2     # total mass [kg]
        self.a0 = a0 * AU             # [m]

        # beta = (64/5) G^3 m1 m2 (m1+m2) / c^5
        self.beta = (64.0 / 5.0) * G**3 * self.m1 * self.m2 * self.M / c**5

        # Characteristic timescale t0 = a0^4 / (4 beta)
        # This is the circular-orbit merger time.
        self.t0 = self.a0**4 / (4.0 * self.beta)

        self._solution = None

    # ------------------------------------------------------------------
    # ODE right-hand side
    # ------------------------------------------------------------------
    @staticmethod
    def _derivs(s, y):
        """RHS of the dimensionless Peters ODEs.

        Parameters
        ----------
        s : float
            Independent variable  s = -ln(alpha).
        y : array_like
            State vector [tau, e].

        Returns
        -------
        dyds : list
            [dtau/ds, de/ds]
        """
        tau, e = y

        # Guard against unphysical eccentricities
        if e < 0 or e >= 1:
            return [0.0, 0.0]

        e2 = e * e
        one_minus_e2 = 1.0 - e2
        g = _g(e2)
        h = _h(e2)

        # dtau/ds = 4 exp(-4s) (1-e^2)^(7/2) / g(e)
        dtau_ds = 4.0 * np.exp(-4.0 * s) * one_minus_e2**3.5 / g

        # de/ds = -(19/12) e (1-e^2) h(e) / g(e)
        de_ds = -(19.0 / 12.0) * e * one_minus_e2 * h / g

        return [dtau_ds, de_ds]

    # ------------------------------------------------------------------
    # Integration
    # ------------------------------------------------------------------
    def integrate(self, s_max=1000.0, rtol=1e-12, atol=1e-14,
                  max_step=np.inf, dense_output=True):
        """Integrate from s = 0 (a = a0) to s = s_max.

        Parameters
        ----------
        s_max : float, optional
            Maximum value of s = -ln(alpha).  Default 100, giving
            alpha = exp(-30) ≈ 9.4 \times 10^{-14}.
        rtol, atol : float
            Tolerances for ``solve_ivp``.
        max_step : float
            Maximum step in s.
        dense_output : bool
            Request interpolating solution.

        Returns
        -------
        sol : OdeResult
            ``scipy.integrate.solve_ivp`` result.  ``sol.t`` contains
            the s-values; ``sol.y[0]`` is tau(s) and ``sol.y[1]`` is e(s).
        """
        y0 = [0.0, self.e0]        # tau(0) = 0, e(0) = e0

        sol = solve_ivp(
            self._derivs,
            (0.0, s_max),
            y0,
            method="RK45",
            rtol=rtol,
            atol=atol,
            max_step=max_step,
            dense_output=dense_output,
        )

        self._solution = sol
        return sol

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------
    @property
    def solution(self):
        """Most recent ``solve_ivp`` result, or *None*."""
        return self._solution

    def get_s(self):
        """Return the array of s values from the last integration."""
        self._check_solved()
        return self._solution.t

    def get_alpha(self):
        """Return α(s) = exp(-s), the dimensionless semi-major axis."""
        return np.exp(-self.get_s())

    def get_tau(self):
        """Return τ(s), the dimensionless time."""
        self._check_solved()
        return self._solution.y[0]

    def get_eccentricity(self):
        """Return e(s)."""
        self._check_solved()
        return self._solution.y[1]

    # ------------------------------------------------------------------
    # Physical-unit convenience methods
    # ------------------------------------------------------------------
    def get_time_seconds(self):
        """Return t(s) in seconds:  t = tau \times t0."""
        return self.get_tau() * self.t0

    def get_time_years(self):
        """Return t(s) in years."""
        return self.get_time_seconds() / (365.25 * 24.0 * 3600.0)

    def get_semi_major_axis_au(self):
        """Return a(s) in AU:  a = alpha \times a0."""
        return self.get_alpha() * self.a0_au

    def get_semi_major_axis_metres(self):
        """Return a(s) in metres:  a = alpha \times a0."""
        return self.get_alpha() * self.a0

    def get_merger_time_seconds(self):
        """Return the merger time in seconds (tau_final \times t0)."""
        self._check_solved()
        return self._solution.y[0, -1] * self.t0

    def get_merger_time_years(self):
        """Return the merger time in years."""
        return self.get_merger_time_seconds() / (365.25 * 24.0 * 3600.0)

    @property
    def t0_years(self):
        """Characteristic timescale t₀ in years (= circular merger time)."""
        return self.t0 / (365.25 * 24.0 * 3600.0)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _check_solved(self):
        if self._solution is None:
            raise RuntimeError("Call integrate() first.")

    def __repr__(self):
        return (f"PetersGWDimensionless(m1={self.m1_solar} M☉, "
                f"m2={self.m2_solar} M☉, a0={self.a0_au} AU, "
                f"e0={self.e0})")
