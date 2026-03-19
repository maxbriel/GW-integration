# GWintegrator

## Introduction

GWintegrator is a Python package for computing the orbital evolution of gravitationally bound binary systems approaching merger.
It transforms the orbital equations based on [Peters (1964)](https://link.aps.org/doi/10.1103/PhysRev.136.B1224) to achieve more stable numerical integration.
Unlike other methods that focus solely on merger time, GWintegrator provides the orbital evolution—allowing you to determine the orbital configuration (semi-major axis and eccentricity) at any time before merger.

## Installation instructions

- via pip
- via github

## How to use

### General case

**Step 1: Import the integrator**

```python
from GWintegrator import GWintegrator
```

**Step 2: Define your binary system parameters**

Specify the masses of both objects (in solar masses), the initial semi-major
axis (in AU), and the initial eccentricity:

```python
m1 = 100  # Primary mass in Msun
m2 = 100  # Secondary mass in Msun
a0 = 0.01  # Initial semi-major axis in AU
e0 = 0.43  # Initial eccentricity (between 0 and 1)
```

**Step 3: Create an integrator instance and compute the orbital evolution**

Initialize the `GWintegrator` object with your parameters and call
the `integrate()` method to compute the complete orbital evolution until merger:

```python
integrator = GWintegrator(m1, m2, a0, e0)
integrator.integrate()

# Get the merger time
print(integrator.merger_time)
```

**Step 4: Query orbital parameters at specific times**

You can retrieve the semi-major axis and eccentricity at any time before merger
by calling the integrator with an array of times (in years):

```python
t = [1e3, 1e4, 1e5, 1e7]  # Times in years
a, e = integrator(t)
```
This will run the integration using the maximum time (`t_max`) in the array (1e7 years in the example)
and return the orbital parameters at the specified times based on the `dense_output` from
`solve_ivp`. If the system merges before reaching `t_max`, orbital values with times
after the merger time will return `np.nan`.

### Available Methods and Properties

Once you've created and integrated a `GWintegrator` object, you can access the following methods and properties:

**`integrate()`** — Computes the orbital evolution of the binary system from the initial conditions until merger.

**`solution`** — The raw solution object returned by `scipy.integrate.solve_ivp`. Contains all interpolation data and can be called like a function to retrieve `(tau, l)` at specified $s$ values.

**`time_array_yr`** — Returns a NumPy array of time values (in years) where the solution was evaluated.

**`separation_array_AU`** — Returns a NumPy array of semi-major axis values (in AU) corresponding to the time array.

**`eccentricity_array`** — Returns a NumPy array of eccentricity values corresponding to the time array.

**`merger_time_yr`** — Returns the time (in years) at which the merger occurs. Only available after calling `integrate()` or having called the integrator. If the system has not yet merged within the integration time, this will raise an error.

Additional methods are available to retrieve the solutions in other parameter spaces, i.e. $\alpha$ space.

## Benchmarks

- Peters merger time calculation

### Circular orbits

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

### Eccentric orbits

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




## Caveats

- Cannot handle e=0, instead a small offset (machine precision offset) is added.
