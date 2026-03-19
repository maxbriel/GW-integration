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
print(t_merger)

integrator = GWIntegrator(m10, m20, a0, e)
integrator.integrate()

print(integrator.merger_time_yr)


integrator2 = PetersGW(m10, m20, a0, e)
integrator2.integrate()

print(integrator2.merger_time_yr)
