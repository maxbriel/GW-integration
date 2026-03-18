import numpy as np
from scipy.integrate import solve_ivp
from scipy.constants import G, c

# ---------- physical constants in SI ----------
M_SUN = 1.989e30   # solar mass [kg]
AU = 1.496e11       # astronomical unit [m]
SEC_PER_YEAR = 365.25 * 24 * 3600  # seconds in a year


# ----------- eccentricity factors -------------
def _g(e2):
    """g(e) = 1 + (73/24) e^2 + (37/96) e^4"""
    return 1.0 + (73.0 / 24.0) * e2 + (37.0 / 96.0) * e2 * e2


def _h(e2):
    """h(e) = 1 + (121/304) e^2"""
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
        self.t0 = self.a0**4 / (4.0 * self.beta)  
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
        
        if e2 >= 1.0:   # guard against e >= 1 (shouldn't occur for e0 < 1)
            return [0.0, 0.0]
        one_minus_e2 = 1.0 - e2
        
        g = _g(e2)
        h = _h(e2)
        
        dtau_ds = 4.0 * np.exp(-4.0 * s) * one_minus_e2**3.5 / g
        dl_ds = -(19.0 / 12.0) * one_minus_e2 * h / g
        
        return [dtau_ds, dl_ds]
    
    
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
        t_max = t_max * SEC_PER_YEAR if t_max is not None else None
        
        def final_time_event(s, y):
            """Event function: triggers when tau reaches t_max / t0."""
            tau = y[0]
            return tau - (t_max / self.t0) if t_max is not None else -1.0
        
        final_time_event.terminal = True
        final_time_event.direction = 0
        
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
    
    @property
    def time_years(self):
        """Get the time array in years."""
        self._check_solved()
        return self._solution.y[0] * self.t0 / SEC_PER_YEAR    
    
    # Helper functions
    
    def _check_solved(self):
        if self._solution is None:
            raise RuntimeError("Call integrate() first.")