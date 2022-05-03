import time

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.integrate import fixed_quad, quad

# Copied from NumII Exercise 11



# TODO: Currently hat functions. Change to polynomials of degree 3.
def get_phi(grid, i):  
  assert i < grid.shape[0] - 1
  return lambda y: (y - grid[i-1])/(grid[i] - grid[i-1]) * (y > grid[i-1])*(y < grid[i]) + (1 - (y - grid[i])/(grid[i+1] - grid[i])) * (y > grid[i])*(y < grid[i+1])

n = 10
x_plot = np.linspace(-1, 1, 1000)

grid = np.linspace(-1, 1, n+2)
fig, ax = plt.subplots()
for i in range(1, n+1):
  phi = get_phi(grid,i)
  ax.plot(x_plot , phi(x_plot))


def get_v(grid, coeff):
  return lambda y: sum([coeff[i-1]*(get_phi(grid, i))(y) for i in range(1, len(coeff)+1)])

def f(x) :
  return 2. * ( 2. * np.pi  ** 2 *  x **  2 - 1. ) * np.sin( 2 * np.pi * x )  \
		- 8 * np.pi * x * np.cos( 2 * np.pi *  x )

def u(x):
  return x ** 2 * np.sin(2 * np.pi * x)

def get_A(grid):
  invEleSize = 1/(grid[1:] - grid[:-1])
  return diags([-invEleSize[1:-1], invEleSize[:-1] + invEleSize[1:], -invEleSize[1:-1]], [-1, 0, 1], format = "csc")

def get_rhs(grid , f):
  return np.array([fixed_quad(lambda x: f(x)*(get_phi(grid, i))(x), grid[i-1], grid[i+1], n = 2500)[0] for i in range(1, len(grid[1:-1])+1)])

n = 20
fig,ax = plt.subplots()
x_plot = np.linspace(-1, 1, 100)
ax.plot(x_plot , u(x_plot), label = "u")

name = "equidistant grid"
grid = np.linspace(-1,1,n + 2)
uh = get_v(grid, spsolve(get_A(grid), get_rhs(grid, f)))
ax.plot(x_plot, uh(x_plot), label = name)

ax.legend()

plt.show()
print('Press enter to continue')
input()
