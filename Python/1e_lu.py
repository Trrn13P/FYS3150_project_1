# Project 1e) LU-decomp. Solving Ax=b.

import sys
import time
import numpy as np
import matplotlib.pyplot as plt

# The non-zero elements of A is given through terminal command:
a = float(sys.argv[1])
b = float(sys.argv[2])
c = float(sys.argv[3])
n = int(sys.argv[4])  # large n causes problems (memory)! 

x = np.linspace(0,1,n+2)
h = 1/(n+1)
f = 100*np.exp(-10*x)

B = h**2*f[1:-1]

A = np.zeros((n,n))        # nxn matrix.

for i in range(n):
    if not i==0:
        A[i][i-1] = a
    A[i][i] = b
    if not i==n-1:
        A[i][i+1] = c

t0 = time.perf_counter_ns()  # clock starts.

V = np.linalg.solve(A,B)   # solves Av = B except endpoints.

t1 = time.perf_counter_ns() - t0  # clock stops.

v = np.zeros(n+2)          # numerical solution.
v[0]   = 0
v[n+1] = 0

for i in range(n):
    v[i+1] = V[i]

u = 1 - (1 - np.exp(-10))*x - np.exp(-10*x)   # analytic solution.

print("Time = {:} ns".format(t1))

plt.plot(x,u)
plt.plot(x,v)
plt.show()
