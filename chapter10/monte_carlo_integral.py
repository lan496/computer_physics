import numpy as np
import random

L = 3.0
N = 1e7

def monte_carlo_integral(f):
    V = 8.0 * L**3
    
    inc = sum([f([random.uniform(-L, L), random.uniform(-L,L), random.uniform(-L,L)]) for i in np.arange(N)])
    return V * inc / N

def f(x):
    if 0 < x[0] < 1 and 0 < x[1] < 1 and x[0] + x[1] - 1 < x[2] < x[0] + x[1]:
        return x[2]
    else:
        return 0

print monte_carlo_integral(f)
