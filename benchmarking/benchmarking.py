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

t_merger = peters_merger_time(m10, m20, a0, e)

integrator = GWIntegrator(m10, m20, a0, e)
integrator.integrate()

integrator2 = PetersGW(m10, m20, a0, e)
integrator2.integrate()

print("Merger time")
print('analytical', t_merger)
print('old       ', integrator2.merger_time_yr)
print('new       ', integrator.merger_time_yr)


print('Solver messages:')
print('old:', integrator2._solution.message)
print('new:', integrator._solution.message)

print("#"*80)

# ----------
# Example 2
# ----------
m10 = 30  # Msun
m20 = 10  # Msun
a0 = 0.01 # AU
e = 0.5

t_merger = peters_merger_time(m10, m20, a0, e)

integrator = GWIntegrator(m10, m20, a0, e)
integrator.integrate()


integrator2 = PetersGW(m10, m20, a0, e)
integrator2.integrate()
print("Merger time")
print('analytical', t_merger)
print('old       ', integrator2.merger_time_yr)
print('new       ', integrator.merger_time_yr)


print('Solver messages:')
print('old:', integrator2._solution.message)
print('new:', integrator._solution.message)
print("#"*80)




# ----------
# Example 3
# ----------
m10 = 30  # Msun
m20 = 1e4  # Msun
a0 = 0.01 # AU
e = 0.9

t_merger = peters_merger_time(m10, m20, a0, e)

integrator = GWIntegrator(m10, m20, a0, e)
integrator.integrate()


integrator2 = PetersGW(m10, m20, a0, e)
integrator2.integrate()
print("Merger time")
print('analytical', t_merger)
print('old       ', integrator2.merger_time_yr)
print('new       ', integrator.merger_time_yr)

print('Solver messages:')
print('old:', integrator2._solution.message)
print('new:', integrator._solution.message)
print("#"*80)




# ----------
# Example 3
# ----------
m10 = 3  # Msun
m20 = 1e7  # Msun
a0 = 0.01 # AU
e = 0.1
print("m1, m2, a0, e0")
print(f'{m10:.2e}, {m20:.2e}, {a0:.2e}, {e:.2f}')

t_merger = peters_merger_time(m10, m20, a0, e)

integrator = GWIntegrator(m10, m20, a0, e)
integrator.integrate()


integrator2 = PetersGW(m10, m20, a0, e)
integrator2.integrate()
print("Merger time")
print('analytical', t_merger)
print('old       ', integrator2.merger_time_yr)
print('new       ', integrator.merger_time_yr)

print('Solver messages:')
print('old:', integrator2._solution.message)
print('new:', integrator._solution.message)
print("#"*80)



# ----------
# Example 3
# ----------
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
print("Merger time")
print('analytical', t_merger)
print('old       ', integrator2.merger_time_yr)
print('new       ', integrator.merger_time_yr)

print('Solver messages:')
print('old:', integrator2._solution.message)
print('new:', integrator._solution.message)
print("#"*80)
