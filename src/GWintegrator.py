import numpy as np
from scipy.constants import G, c
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# ---------- physical constants in SI ----------
M_SUN = 1.989e30   # solar mass [kg]
AU = 1.496e11       # astronomical unit [m]
SEC_PER_YEAR = 365.25 * 24 * 3600  # seconds in a year


# ----------- eccentricity factors -------------
def _g(e2):
    """g(e) = 1 + (73/24) e^2 + (37/96) e^4"""
    return 1.0 + (73.0 / 24.0) * e2 + (37.0 / 96.0) * e2 * e2


def _f(e2):
    """f(e) = 1 + (121/304) e^2"""
    return 1.0 + (121.0 / 304.0) * e2


class GWIntegrator:
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

        # Characteristic timescale
        # circular merger time
        self.t0 = self.a0**4 / (self.beta)
        self._solution = None

    @staticmethod
    def _rhs(s,y):
        """RHS of Peters' equations in ln-space.

        Parameters
        ----------
        s : float
            Independent variable  s = -ln(alpha).
        y : array_like
            State vector [tau, l],
            where tau = t / t0 (dimensionless time) and l = ln(e)

        Returns
        -------
        dyds : list
            [dtau/ds, dl/ds]
        """

        tau, l = y

        e = np.exp(l)
        e2 = e * e

        if e2 >= 1.0:   # guard against e >= 1
            return [0.0, 0.0]

        one_minus_e2 = 1.0 - e2
        G = _g(e2)
        F = _f(e2)

        dtau_ds = np.exp(-4.0 * s) * one_minus_e2**3.5 / G
        dl_ds = -(19.0 / 12.0) * one_minus_e2 * F / G

        return [dtau_ds, dl_ds]


    def __call__(self, t, s_max=1e3, rtol=1e-12, atol=1e-12,
                  max_step=np.inf):
        """Return the orbital configuration at one or more times.

        Integrates the systems of ODEs for the maximum in the given t array.
        It then uses the `dense_output` from `solve_ivp` to request the other t values.
        To this ends, it inverts the s->tau relationship using scipy.brentq.

        For time values beyond the merger time of the system, np.nan is returned.

        Parameters
        ----------
        t : float or array_like of float
            Times (in years) at which to evaluate the orbit.
        s_max : float, optional
            Maximum s = -ln(alpha) to integrate to (default: 1000).
            This is a safeguard against infinite integration if t_max is None.
        rtol, atol : float
            Relative / absolute tolerances passed to ``solve_ivp``.
        max_step : float
            Maximum internal step size in seconds.

        Returns
        -------
        a : np.ndarray
            Semi-major axis in AU at each requested time.
        e : np.ndarray
            Eccentricity at each requested time.
        """
        t = np.atleast_1d(np.asarray(t, dtype=float))

        sol = self.integrate(t_max=float(np.max(t)), s_max=s_max,
                             rtol=rtol, atol=atol, max_step=max_step,
                             dense_output=True)

        # check if s_max reached instead of t_max
        if sol.status == 0:
            print("Merger conditions have been reached")
            print("Some t values might not return a value")

        # Compare merger_time against t array
        mask = t <= self.merger_time_yr
        masked_t = t[mask]

        s_lo, s_hi = sol.t[0], sol.t[-1]

        out_a = np.empty(len(t))
        out_e = np.empty(len(t))
        for i, ti in enumerate(masked_t):
            tau_target = ti * SEC_PER_YEAR / self.t0
            # tau(s) is monotone -> brentq converges
            # inverting s -> tau relationship
            s_star = brentq(lambda s: sol.sol(s)[0] - tau_target, s_lo, s_hi,
                            xtol=1e-14, rtol=1e-14)
            out_a[i] = np.exp(-s_star) * self.a0 / AU
            out_e[i] = np.exp(sol.sol(s_star)[1])

        out_a[~mask] = np.nan
        out_e[~mask] = np.nan

        return out_a, out_e

    def integrate(self, t_max=None, s_max=1e3, rtol=1e-12, atol=1e-12,
                  max_step=np.inf, dense_output=True):
        """Integrate Peters' equations until merger or t_max.

        Parameters
        ----------
        t_max : float, optional
            Maximum integration time in yrs.
            If None, integrate until merger.
        s_max : float, optional
            Maximum s = -ln(alpha) to integrate to (default: 1000).
            This is a safeguard against infinite integration if t_max is None.
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
        self.t_max = t_max * SEC_PER_YEAR if t_max is not None else None

        # Stop evolution based on max time
        def final_time_event(s, y):
            """Event function: triggers when tau reaches t_max / t0.
            Calculated in seconds

            Parameters
            ----------
            s : float
                Independent variable  s = -ln(alpha).
            y : array_like
                [tau, l]
            """
            tau = y[0]
            return tau - (self.t_max / self.t0) if self.t_max is not None else -1.0

        final_time_event.terminal = True
        final_time_event.direction = 0

        # Set l to a machine small value if given e0==0
        l0 = np.log(self.e0) if self.e0 > 0 else np.log(np.finfo(float).tiny)
        y0 = [0.0, l0]

        sol = solve_ivp(
            self._rhs,
            (0.0, s_max),
            events=final_time_event,
            y0=y0,
            method="RK45",
            rtol=rtol,
            atol=atol,
            max_step=max_step,
            dense_output=dense_output,
        )

        self._solution = sol

        return sol

    ### Main output variables

    @property
    def solution(self):
        """The most recent ``solve_ivp`` result"""
        return self._solution

    @property
    def time_array_yr(self):
        """Get the time array in years."""
        self._check_solved()
        return self._solution.y[0] * self.t0 / SEC_PER_YEAR

    @property
    def separation_array_AU(self):
        """Return a(s) = alpha * a0 the semi-major axis in AU."""
        self._check_solved()
        return self.get_alpha() * self.a0 / AU

    @property
    def eccentricity_array(self):
        """Return e(s) = exp(l(s))."""
        self._check_solved()
        return np.exp(self._solution.y[1])


    # Internal output variables

    def get_alpha(self):
        """Return alpha(s) = exp(-s), the dimensionless semi-major axis."""
        self._check_solved()
        return np.exp(-self.get_s())

    def get_s(self):
        """Return the array of s values"""
        self._check_solved()
        return self._solution.t

    @property
    def merger_time_yr(self):
        """Get the merger time in years."""
        self._check_solved()
        self._check_merger()
        return self._solution.y[0, -1] * self.t0 / SEC_PER_YEAR

    # Helper functions
    def _check_solved(self):
        if self._solution is None:
            raise RuntimeError("Call integrate() first.")

    # Check merger
    def _check_merger(self):
        # sol status = 0 means the end of the integrtaion interval was reached.
        if self._solution.status != 0:
            raise RuntimeError('Merger not found.')


    def __repr__(self):
        return (f"GWIntegrator"
                f"    m1={self.m1_solar} M☉"
                f"    m2={self.m2_solar} M☉"
                f"    a0={self.a0_au} AU"
                f"    e0={self.e0}")
