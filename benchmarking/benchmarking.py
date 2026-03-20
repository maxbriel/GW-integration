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



# Make a figure of an example systems approaching the merger time

t = np.linspace(0,1e9,1000000000)
m10 = 1e4  # Msun
m20 = 1e7  # Msun
a0 = 1 # AU
e = 0.9
print("m1, m2, a0, e0")
print(f'{m10:.2e}, {m20:.2e}, {a0:.2e}, {e:.2f}')

t_merger = peters_merger_time(m10, m20, a0, e)

integrator = GWIntegrator(m10, m20, a0, e)
integrator.integrate()


integrator2 = PetersGW(m10, m20, a0, e)
integrator2.integrate()


plt.figure()
#plt.axvline(t_merger, color='black')
plt.plot(integrator2.time_array_yr[-1] - integrator2.time_array_yr, integrator2.separation_array_AU, label='old', ls='-',)
plt.plot(integrator.time_array_yr[-1] - integrator.time_array_yr, integrator.separation_array_AU, label='new', ls='--')
#plt.xlim(0.001, None)
plt.xscale('log')
plt.ylim(1e-6, 1e1)
plt.xlim(1e-18, 1e1)
plt.yscale('log')
plt.xlabel('time till final value [yr]')
plt.ylabel('separation [AU]')
plt.legend()
plt.show()
