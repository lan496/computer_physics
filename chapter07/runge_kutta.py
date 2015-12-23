import numpy as np
import matplotlib.pyplot as plt

def euler(f, t, x, h):
    return x + h * f(t, x)

def runge_kutta(f, t, x, h):
    k1 = f(t, x)
    k2 = f(t + h/2, x + k1/2)
    k3 = f(t + h/2, x + k2/2)
    k4 = f(t + h, x + k3)

    return x + h/6*(k1 + 2*k2 + 2*k3 + k4)

def differential_equation(f, t0, x0, t_end, h):
    # f (t, x[]) -> xx[]
    step = (int)((t_end - t0) / h)
    q = [[0, 0] for i in np.arange(step)]
    q[0] = [t0, x0]
    for i in np.arange(1,step):
        #q[i] = [t0 + i*h, runge_kutta(f, t0 + i*h, q[i-1][1], h)]
        q[i] = [t0 + i*h, euler(f, t0 + i*h, q[i-1][1], h)]
        if i % (step/10) == 0:
            print (i*t_end)/step
    return q

def logistic_equation():
    h = 0.0000001
    q = differential_equation(lambda t,x: x*(1-x), 0, 0.5, 6, h)
    plt.plot([e[0] for e in q], [e[1] for e in q])
    plt.show()

logistic_equation()
