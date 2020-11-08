"""

dy/dx = ax - bxy  ---> f_x(t,x,y)

dy/dt = -cy + dxy  ---> f_y(t,x,y)

x(0) = 10
x(0) = 10
a = 1.1
b = 0.4
c = 0.4
d = 0.1

Calculando inclinações:
K1x = f_x(x_n,y_n,t_n)
K1y = f_y(x_n,y_n,t_n)

K2x = f_x(x_n + h*K1x , y_n + h*K1y , t_n + h)
K2y = f_y(x_n + h*K1x , y_n + h*K1y , t_n + h)

"""

"""

    euler e euler melhorado (heun) para sistemas.

"""

import numpy as np
import matplotlib.pylab as plt


def fx(t, r):  # FUNÇÃO QUE DEFINE A EDO (derivada)
    x, y  = r
    a = 1.1
    b = 0.4
    c = 0.4
    d = 0.1
    return np.array([a*x - b*x*y, -c*y + d*x*y])



def EDO_euler_sistemas(f, r0, t0, NUMBER_OF_STEPS=100, h=0.01):

    NUMBER_OF_EQUATIONS = len(r0)

    r = np.zeros([NUMBER_OF_STEPS, NUMBER_OF_EQUATIONS], dtype=np.float32)
    t = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)

    r[0] = r0
    t[0] = t0

    for n in range(0, NUMBER_OF_STEPS - 1):
        t[n+1] = t[n]+h
        K1= fx(t[n],r[n])

        r[n+1] = r[n] + K1*h

    return (t, r)

def EDO_heun_sistemas(f, r0, t0, NUMBER_OF_STEPS=100, h=0.01):

    NUMBER_OF_EQUATIONS = len(r0)

    r = np.zeros([NUMBER_OF_STEPS, NUMBER_OF_EQUATIONS], dtype=np.float32)
    t = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)

    r[0] = r0
    t[0] = t0

    for n in range(0, NUMBER_OF_STEPS - 1):
        t[n+1] = t[n]+h
        K1= fx(t[n],r[n])
        K2= fx(t[n+1],r[n] + K1*h )

        r[n+1] = r[n] + 0.5*(K1+K2)*h

    return (t, r)

t,r = EDO_heun_sistemas(fx,(10,10), 0 , NUMBER_OF_STEPS = 1000, h = 0.05)

x=r[:,0]
y=r[:,1]

plt.plot(t, x)
plt.plot(t, y)

plt.show()
