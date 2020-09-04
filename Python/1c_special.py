# Project 1, c) Special algorithm.
# Av = B.
"""
The matrix A has identical matrix elements along the diagonal and identical
(but different) values for the non- diagonal elements.
Specialize your algorithm to the special case and find the number of floating
point operations for this specific tri-diagonal matrix.
Compare the CPU time with the general algorithm (1b_general.py).
"""

import sys
import time
import numpy as np
import matplotlib.pyplot as plt

# The non-zero elements of A is given through terminal command:
a = float(sys.argv[1])
b = float(sys.argv[2])
c = float(sys.argv[3])

def func(n,C,D,uC,uD,a,b,c):
    # Solves the matrix equation Av=B, including boundaries (u(C) and u(D)).
    # A is a special case nxn tridiagonal matrix.

    x = np.linspace(C,D,n+2)
    h = (D-C)/(n+1)
    f = 100*np.exp(-10*x)

    B = h**2*f

    v = np.zeros(n+2)         # numerical solution.
    v[0]   = uC                # initial condition at x = C.
    v[n+1] = uD                # initial condition at x = D.

    t0 = time.perf_counter_ns()   # start timer.

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

    t1 = time.perf_counter_ns() - t0   # end timer.

    FLOPS = 7*n   # approximatly the number of floating point operations.

    return x,v,t1,FLOPS

n = np.asarray([10,100,1000,10**4,10**5,10**6])

p = np.linspace(0,1,1000)

u = 1 - (1 - np.exp(-10))*p - np.exp(-10*p)   # analytic solution.

for i in range(len(n)):
    x,v,t1,FLOPS = func(n=n[i],C=0,D=1,uC=0,uD=0,a=a,b=b,c=c)
    plt.plot(x,v,label=r"Numerical solution, $n={:}$ steps".format(n[i]))
    plt.plot(p,u,label="Analytic solution")
    plt.xlabel("x") ; plt.ylabel("u")
    plt.title(r"Project 1c | $FLOPS={:}$ | $CPU-time={:}\:ns$".format(FLOPS,t1))
    plt.legend() ; plt.grid()
    #plt.savefig("1c_{:}.png".format(n[i]))
    plt.show()

"""
terminal> python3 1c_special.py -1 2 -1
"""
