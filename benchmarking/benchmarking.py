"""Script to benchmark the GWintegrator against the analytical solution for merger time
and against the non-transformed differential equations."""
#
import matplotlib.pyplot as plt
import numpy as np
from solvers import PetersGW, peters_merger_time

from GWintegrator import GWIntegrator

# ----------
# Example 1
# ----------
m10 = 30  # Msun
m20 = 10  # Msun
a0 = 0.05 # AU
e = 0
print("m1, m2, a0, e0")
print(f'{m10:.2e}, {m20:.2e}, {a0:.2e}, {e:.2f}')

t_merger = peters_merger_time(m10, m20, a0, e)

integrator = GWIntegrator(m10, m20, a0, e)
sol_new = integrator.integrate()
nfev_new = sol_new.nfev

integrator2 = PetersGW(m10, m20, a0, e)
sol_old = integrator2.integrate()
nfev_old = sol_old.nfev

print("Merger time")
print('analytical', t_merger)
print('old       ', integrator2.merger_time_yr)
print('new       ', integrator.merger_time_yr)

print('Solver messages:')
print('old:', integrator2._solution.message)
print('new:', integrator._solution.message)

print("Function evaluations")
print('old       ', nfev_old)
print('new       ', nfev_new)
print('ratio     ', nfev_new/nfev_old)
print("#"*80)

# ----------
# Example 2
# ----------
m10 = 30  # Msun
m20 = 10  # Msun
a0 = 0.01 # AU
e = 0.5
print("m1, m2, a0, e0")
print(f'{m10:.2e}, {m20:.2e}, {a0:.2e}, {e:.2f}')

t_merger = peters_merger_time(m10, m20, a0, e)

integrator = GWIntegrator(m10, m20, a0, e)
sol_new = integrator.integrate()
nfev_new = sol_new.nfev

integrator2 = PetersGW(m10, m20, a0, e)
sol_old = integrator2.integrate()
nfev_old = sol_old.nfev

print("Merger time")
print('analytical', t_merger)
print('old       ', integrator2.merger_time_yr)
print('new       ', integrator.merger_time_yr)


print('Solver messages:')
print('old:', integrator2._solution.message)
print('new:', integrator._solution.message)

print("Function evaluations")
print('old       ', nfev_old)
print('new       ', nfev_new)
print('ratio     ', nfev_new/nfev_old)
print("#"*80)




# ----------
# Example 3
# ----------
m10 = 30  # Msun
m20 = 1e4  # Msun
a0 = 0.01 # AU
e = 0.9
print("m1, m2, a0, e0")
print(f'{m10:.2e}, {m20:.2e}, {a0:.2e}, {e:.2f}')

t_merger = peters_merger_time(m10, m20, a0, e)

integrator = GWIntegrator(m10, m20, a0, e)
sol_new = integrator.integrate()
nfev_new = sol_new.nfev

integrator2 = PetersGW(m10, m20, a0, e)
sol_old = integrator2.integrate()
nfev_old = sol_old.nfev


print("Merger time")
print('analytical', t_merger)
print('old       ', integrator2.merger_time_yr)
print('new       ', integrator.merger_time_yr)

print('Solver messages:')
print('old:', integrator2._solution.message)
print('new:', integrator._solution.message)

print("Function evaluations")
print('old       ', nfev_old)
print('new       ', nfev_new)
print('ratio     ', nfev_new/nfev_old)
print("#"*80)




# ----------
# Example 4
# ----------
m10 = 3  # Msun
m20 = 1e7  # Msun
a0 = 0.01 # AU
e = 0.1
print("m1, m2, a0, e0")
print(f'{m10:.2e}, {m20:.2e}, {a0:.2e}, {e:.2f}')

t_merger = peters_merger_time(m10, m20, a0, e)

integrator = GWIntegrator(m10, m20, a0, e)
sol_new = integrator.integrate()
nfev_new = sol_new.nfev

integrator2 = PetersGW(m10, m20, a0, e)
sol_old = integrator2.integrate()
nfev_old = sol_old.nfev

print("Merger time")
print('analytical', t_merger)
print('old       ', integrator2.merger_time_yr)
print('new       ', integrator.merger_time_yr)

print('Solver messages:')
print('old:', integrator2._solution.message)
print('new:', integrator._solution.message)

print("Function evaluations")
print('old       ', nfev_old)
print('new       ', nfev_new)
print('ratio     ', nfev_new/nfev_old)
print("#"*80)


# ----------
# Example 5
# ----------
m10 = 1e4  # Msun
m20 = 1e7  # Msun
a0 = 1 # AU
e = 0.9
print("m1, m2, a0, e0")
print(f'{m10:.2e}, {m20:.2e}, {a0:.2e}, {e:.2f}')

t_merger = peters_merger_time(m10, m20, a0, e)

integrator = GWIntegrator(m10, m20, a0, e)
sol_new = integrator.integrate()
nfev_new = sol_new.nfev

integrator2 = PetersGW(m10, m20, a0, e)
sol_old = integrator2.integrate()
nfev_old = sol_old.nfev

print("Merger time")
print('analytical', t_merger)
print('old       ', integrator2.merger_time_yr)
print('new       ', integrator.merger_time_yr)

print('Solver messages:')
print('old:', integrator2._solution.message)
print('new:', integrator._solution.message)

print("Function evaluations")
print('old       ', nfev_old)
print('new       ', nfev_new)
print('ratio     ', nfev_new/nfev_old)
print("#"*80)



# ----------
# Example 6
# ----------
m10 = 0.3  # Msun
m20 = 0.5  # Msun
a0 = 1e-4 # AU
e = 0.9
print("m1, m2, a0, e0")
print(f'{m10:.2e}, {m20:.2e}, {a0:.2e}, {e:.2f}')

t_merger = peters_merger_time(m10, m20, a0, e)

integrator = GWIntegrator(m10, m20, a0, e)
sol_new = integrator.integrate()
nfev_new = sol_new.nfev

integrator2 = PetersGW(m10, m20, a0, e)
sol_old = integrator2.integrate()
nfev_old = sol_old.nfev

print("Merger time")
print('analytical', t_merger)
print('old       ', integrator2.merger_time_yr)
print('new       ', integrator.merger_time_yr)

print('Solver messages:')
print('old:', integrator2._solution.message)
print('new:', integrator._solution.message)

print("Function evaluations")
print('old       ', nfev_old)
print('new       ', nfev_new)
print('ratio     ', nfev_new/nfev_old)
print("#"*80)



# ══════════════════════════════════════════════════════════════════════
# Figure: orbital evolution comparison  (separation + eccentricity)
# x-axis = t [yr],  y-axes on log scale for separation
# No dependence on analytical merger time.
# ══════════════════════════════════════════════════════════════════════
m10 = 1e4  # Msun
m20 = 1e7  # Msun
a0 = 1     # AU
e = 0.9
print("Figure case — m1, m2, a0, e0")
print(f'{m10:.2e}, {m20:.2e}, {a0:.2e}, {e:.2f}')

integrator = GWIntegrator(m10, m20, a0, e)
integrator.integrate()

integrator2 = PetersGW(m10, m20, a0, e)
integrator2.integrate()

t_old = integrator2.time_array_yr
t_new = integrator.time_array_yr
a_old = integrator2.separation_array_AU
a_new = integrator.separation_array_AU
e_old = integrator2.eccentricity_array
e_new = integrator.eccentricity_array

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# ── Left panel: semi-major axis ───────────────────────────────────────
ax1.plot(t_old, a_old, color='#d62728', lw=1.5,
         label='Original Peters ODE')
ax1.plot(t_new, a_new, color='#1f77b4', lw=1.5, ls='--',
         label='Transformed (ln-space) ODE')

# Mark where old solver stopped
ax1.plot(t_old[-1], a_old[-1], 'x', color='#d62728', ms=4, mew=2.5,
         zorder=5)
ax1.annotate(f'Old solver fails\n$a$ = {a_old[-1]:.2e} AU',
             xy=(t_old[-1], a_old[-1]),
             xytext=(t_old[-1] * 0.6, a_old[-1] * 0.1),
             arrowprops=dict(arrowstyle='->', color='#d62728', lw=1),
             fontsize=9, color='#d62728')

ax1.set_yscale('log')
ax1.set_xlabel('Time [yr]')
ax1.set_ylabel('Semi-major axis [AU]')
ax1.legend(loc='lower left', fontsize=9)
ax1.set_title(f'$m_1=10^4\\,M_\\odot$, $m_2=10^7\\,M_\\odot$, '
              f'$a_0={a0}$ AU, $e_0={e}$', fontsize=10)

# ── Right panel: eccentricity ─────────────────────────────────────────
ax2.plot(t_old, e_old, color='#d62728', lw=1.5,
         label='Original Peters ODE')
ax2.plot(t_new, e_new, color='#1f77b4', lw=1.5, ls='--',
         label='Transformed (ln-space) ODE')

ax2.plot(t_old[-1], e_old[-1], 'x', color='#d62728', ms=4, mew=2.5,
         zorder=5)

ax2.set_xlabel('Time [yr]')
ax2.set_ylabel('Eccentricity')
ax2.set_yscale('log')
ax2.set_ylim(1e-10, 1)
ax2.legend(loc='lower left', fontsize=9)
ax2.set_title('Eccentricity evolution', fontsize=10)

# nfev text box
nfev_txt = (f"nfev original:     {integrator2._solution.nfev}\n"
            f"nfev transformed: {integrator._solution.nfev}\n"
            f"ratio: {integrator2._solution.nfev / integrator._solution.nfev:.1f}×")
ax2.text(0.03, 0.50, nfev_txt, transform=ax2.transAxes,
         fontsize=8, verticalalignment='top', horizontalalignment='left',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='wheat', alpha=0.8))

fig.tight_layout()
fig.savefig('orbital_evolution_comparison.png', dpi=150)
print("\nSaved: orbital_evolution_comparison.png")
plt.show()
