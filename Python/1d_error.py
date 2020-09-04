# Project 1, d) Computing relative error.
# Av = B.

from numba import jit
import sys
import numpy as np
import matplotlib.pyplot as plt

# The non-zero elements of A is given through terminal command:
a = float(sys.argv[1])
b = float(sys.argv[2])
c = float(sys.argv[3])

@jit

def func(n,C,D,uC,uD,a,b,c):
    # Solves the matrix equation Av=B, including boundaries (u(C) and u(D)).
    # A is a special case nxn tridiagonal matrix.

    x = np.linspace(C,D,n+2)
    h = (D-C)/(n+2)
    f = 100*np.exp(-10*x)

    B = h**2*f

    v = np.zeros(n+2)         # numerical solution.
    v[0]   = uC                # initial condition at x = C.
    v[n+1] = uD                # initial condition at x = D.

    # Row reduction. First step:
    y = np.zeros(n) ; y[0] = b
    z = np.zeros(n) ; z[0] = c
    g = np.zeros(n) ; g[0] = B[1]

    # Forward substitution:
    for i in range(1,n):
        y[i] = y[i-1]*b/a - z[i-1]
        g[i] = y[i-1]*B[i+1]/a - g[i-1]
        z[i] = y[i-1]

    # Solving final set of equations:
    v[n] = g[n-1]/y[n-1]
    # Backward substitution:
    for i in range(1,n):
        j = n - i
        v[j] = (g[j-1] - z[j-1]*v[j+1])/y[j-1]

    u = 1 - (1 - np.exp(-10))*x - np.exp(-10*x)   # analytic solution.

    return h,v,u,n

n = np.asarray([10,100,1000,10**4,10**5,10**6,10**7])

# Relative error for different step sizes:
eps  = np.zeros(len(n))
step = np.zeros(len(n))

for i in range(len(n)):
    h,v,u,m = func(n=n[i],C=0,D=1,uC=0,uD=0,a=a,b=b,c=c)
    rel_error = (v[1:m-1] - u[1:m-1])/u[1:m-1]
    epsi    = np.log10(abs(rel_error))
    eps[i]  = np.max(epsi)   # relative error, log10. Max value.
    step[i] = np.log10(h)    # step size, log10.

plt.plot(step,eps)
plt.xlabel(r"$log_{10}(h)$ - step size")
plt.ylabel(r"$log_{10}$ - relative error")
plt.title("Project 1d) Relative error as a function of step size")
plt.grid()
plt.savefig("1d_error.png")
plt.show()


"""
terminal> python3 1d_error.py -1 2 -1
"""
