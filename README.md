# GWintegrator

## Introduction

GWintegrator is a Python package for computing the orbital evolution of gravitationally bound binary systems approaching merger. It transforms the orbital equations based on [Peters (1964)](https://link.aps.org/doi/10.1103/PhysRev.136.B1224) to achieve more stable numerical integration. Unlike other methods that focus solely on merger time, GWintegrator provides the orbital evolution—allowing you to determine the orbital configuration (semi-major axis and eccentricity) at any time before merger.

## Installation instructions

- via pip
- via github

## How to use

## general case
import
```from GWintegrator import GWintegrator```

setup
```
m10 = 100  # Msun
m20 = 100  # Msun
a0 = 0.01  # AU
e0 = 0.43
integrator = GWintegrator(m1, m2, a, e)
integrator.integrate()

print(integrator.merger_time)
```


getting values at different t:
```
t = [1e3, 1e4, 1e5, 1e7]
a, e = integrator(t)
```

available commands:
- integrate()
- solution
- time_array_yr
- separation_array_AU
- eccentricity_array
- merger_time_yr (if reached)




# Benchmarks

- Peters merger time calculation

## Circular orbits

Analytical solution to integration of Peters' equations:

$$
T_{m} = \frac{}
$$

We can use this analytical solution to benchmark the ODE integrator.

Example setups:
```
m10 = 30  # Msun
m20 = 10  # Msun
a0 = 0.05 # AU
e = 0
```

## Eccentric orbits

For eccentric orbits we need to solve the integral given by Peters:


Example setups
```
m10 = 100  # Msun
m20 = 100  # Msun
a0 = 0.01  # AU
e0 = 0.43
```

```
m10 = 30  # Msun
m20 = 10  # Msun
a0 = 0.05 # AU
e = 0.9999
```




# Caveats

- Cannot handle e=0, instead a small offset (machine precision offset) is added.
