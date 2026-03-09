#!/usr/bin/env python3
"""
Compare the two Peters (1964) integrators.

Produces a two-panel figure (a vs t and e vs t) overlaying:
    • Original integrator  (physical variables, solid lines)
    • Dimensionless integrator  (s-variable, dashed lines)

A residual panel is included for each quantity.

Usage
-----
    python plot_evolution.py

Adjust *m1*, *m2*, *a0*, and *e0* below to explore different binaries.
"""

import numpy as np
import matplotlib.pyplot as plt
from integrator import PetersGWIntegrator
from integrator_dimensionless import PetersGWDimensionless

# ── Binary parameters ─────────────────────────────────────────────────
m1 = 1.4          # solar masses  (e.g. neutron star)
m2 = 1.4          # solar masses
a0 = 0.01         # initial semi-major axis in AU  (~1.5e9 m)
e0 = 0.6          # initial eccentricity
# ──────────────────────────────────────────────────────────────────────

def main():
    # ── Peters analytic merger time (no integration, Peters 1964 Eq. 5.10) ──
    _tmp = PetersGWIntegrator(m1, m2, a0, e0)
    t_merge_analytic_yr = _tmp.seconds_to_years(_tmp.merger_time_peters())
    t_max = _tmp.merger_time_peters()   # use analytic estimate as upper bound

    # ── Run both integrators ──────────────────────────────────────────
    # 1) Original (physical variables)
    gw_phys = PetersGWIntegrator(m1, m2, a0, e0)
    sol_phys = gw_phys.integrate(t_max=t_max)

    t_yr_phys = gw_phys.seconds_to_years(sol_phys.t)
    a_au_phys = gw_phys.metres_to_au(sol_phys.y[0])
    e_phys = sol_phys.y[1]
    
    t_merge_phys = gw_phys.seconds_to_years(gw_phys.get_merger_time())

    
    print(f"Physical integrator: {len(t_yr_phys)} time steps")
    print(f"Physical integrator: status {sol_phys.message}")
    print(f" -> a_final = {a_au_phys[-1]} AU")
    print(f" -> a0 = {a0} AU")
    print(f" -> alpha_final = {a_au_phys[-1]/a0}")
    print(f" -> e_final = {e_phys[-1]:.2e}")
    print(f" -> merger time = {t_merge_phys} yr (physical)")

    # 2) Dimensionless (s-variable)
    gw_dim = PetersGWDimensionless(m1, m2, a0, e0)
    gw_dim.integrate()
    print(f"Dimensionless integrator: {len(gw_dim.get_time_years())} time steps\n"
          f" (s_max = {gw_dim._solution.t[-1]:.2f})\n"
            f" -> alpha_final = {gw_dim.get_alpha()[-1]:.2e})\n"
            f" -> e_final = {gw_dim.get_eccentricity()[-1]:.2e})\n"
            f" -> merger time = {gw_dim.get_merger_time_years():.6e} yr\n"
            f" (status {gw_dim._solution.message})")

    t_yr_dim = gw_dim.get_time_years()
    a_au_dim = gw_dim.get_semi_major_axis_au()
    e_dim = gw_dim.get_eccentricity()

    t_merge_dim = gw_dim.get_merger_time_years()
    
    print(f"Dimensionless merger time: {t_merge_dim} yr")
    
    print(t_merge_dim, t_merge_phys)

    # ── Print summary ─────────────────────────────────────────────────
    print(f"Binary: {m1} M_sun + {m2} M_sun")
    print(f"Initial orbit:  a0 = {a0} AU,  e0 = {e0}")
    print(f"Merger time (Peters analytic): {t_merge_analytic_yr:.6e} yr")
    print(f"Merger time (physical):        {t_merge_phys:.6e} yr")
    print(f"Merger time (dimensionless):   {t_merge_dim:.6e} yr")
    print(f"Relative difference:           "
          f"{abs(t_merge_phys - t_merge_dim) / t_merge_dim:.2e}")

    # ── Interpolate dimensionless solution onto physical time grid ────
    # for a fair residual comparison
    from scipy.interpolate import interp1d

    # Remove duplicate time entries (tau saturates near merger)
    _, unique_idx = np.unique(t_yr_dim, return_index=True)
    t_yr_dim_u = t_yr_dim[unique_idx]
    a_au_dim_u = a_au_dim[unique_idx]
    e_dim_u = e_dim[unique_idx]

    interp_a = interp1d(t_yr_dim_u, a_au_dim_u, kind="cubic",
                        bounds_error=False, fill_value=np.nan)
    interp_e = interp1d(t_yr_dim_u, e_dim_u, kind="cubic",
                        bounds_error=False, fill_value=np.nan)

    a_dim_on_phys = interp_a(t_yr_phys)
    e_dim_on_phys = interp_e(t_yr_phys)

    # ── Figure ────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(12, 7),
                             gridspec_kw={"width_ratios": [3, 1]})

    # ---- Semi-major axis (left) ----
    ax_a = axes[0, 0]
    ax_a.plot(t_yr_phys, a_au_phys,
              color="steelblue", lw=1.8, label="Physical")
    ax_a.plot(t_yr_dim, a_au_dim,
              color="orange", lw=1.8, ls="--", label="Dimensionless")
    ax_a.set_ylabel("Semi-major axis  $a$  [AU]")
    ax_a.set_title(
        f"Peters (1964) — "
        f"$m_1={m1}\\,M_\\odot$, $m_2={m2}\\,M_\\odot$, "
        f"$a_0={a0}$ AU, $e_0={e0}$"
    )
    ax_a.legend(fontsize=9)
    ax_a.grid(True, ls="--", alpha=0.4)

    # ---- Semi-major axis residual (right) ----
    ax_a_res = axes[0, 1]
    mask_a = np.isfinite(a_dim_on_phys) & (a_dim_on_phys > 0)
    rel_a = (a_au_phys[mask_a] - a_dim_on_phys[mask_a]) / a_dim_on_phys[mask_a]
    ax_a_res.plot(t_yr_phys[mask_a], rel_a, color="steelblue", lw=0.8)
    ax_a_res.set_ylabel("Relative diff  $\\Delta a / a$")
    ax_a_res.axhline(0, color="k", lw=0.5)
    ax_a_res.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    ax_a_res.grid(True, ls="--", alpha=0.4)
    ax_a_res.set_title("Residual")

    # ---- Eccentricity (left) ----
    ax_e = axes[1, 0]
    ax_e.plot(t_yr_phys, e_phys,
              color="firebrick", lw=1.8, label="Physical")
    ax_e.plot(t_yr_dim, e_dim,
              color="orange", lw=1.8, ls="--", label="Dimensionless")
    ax_e.set_ylabel("Eccentricity  $e$")
    ax_e.set_xlabel("Time  [yr]")
    ax_e.set_ylim(-0.02, min(1.0, e0 * 1.15))
    ax_e.legend(fontsize=9)
    ax_e.grid(True, ls="--", alpha=0.4)

    # ---- Eccentricity residual (right) ----
    ax_e_res = axes[1, 1]
    mask_e = np.isfinite(e_dim_on_phys) & (np.abs(e_dim_on_phys) > 1e-10)
    rel_e = (e_phys[mask_e] - e_dim_on_phys[mask_e]) / e_dim_on_phys[mask_e]
    ax_e_res.plot(t_yr_phys[mask_e], rel_e, color="firebrick", lw=0.8)
    ax_e_res.set_ylabel("Relative diff  $\\Delta e / e$")
    ax_e_res.set_xlabel("Time  [yr]")
    ax_e_res.axhline(0, color="k", lw=0.5)
    ax_e_res.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    ax_e_res.grid(True, ls="--", alpha=0.4)
    ax_e_res.set_title("Residual")

    for ax in axes.flatten():
        ax.set_xscale('log')
        #ax.set_xlim(0, t_merge_phys * 1.05)

    plt.tight_layout()
    plt.savefig("orbital_evolution.png", dpi=150)
    print("Saved figure -> orbital_evolution.png")
    plt.show()


if __name__ == "__main__":
    main()
