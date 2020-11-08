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

import numpy as np
import matplotlib.pylab as plt


def fx(t, x, y):  # FUNÇÃO QUE DEFINE A EDO (derivada)
    a = 1.1
    b = 0.4
    c = 0.4
    d = 0.1
    return np.array([a*x - b*x*y, -c*y + d*x*y])



def EDO_euler(f, y0, t0, NUMBER_OF_STEPS=100, h=0.01):
    # recebe a função derivada, y(0) e t(0),
    # os demais paramentros são definidos, não precisa passar

    y = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)
    t = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)

    y[0] = 1
    t[0] = 0

    for n in range(0, NUMBER_OF_STEPS - 1):
        K1 = f(t[n], y[n])
        y[n+1] = y[n] + K1*h  # fórmula iterativa do método de euler
        t[n+1] = t[n]+h  # próximo passo de tempo

    return (t, y)


def EDO_heun(f, y0, t0, NUMBER_OF_STEPS=100, h=0.01):
    # recebe a função derivada, y(0) e t(0),
    # os demais paramentros são definidos, não precisa passar

    y = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)
    t = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)

    y[0] = 1
    t[0] = 0

    for n in range(0, NUMBER_OF_STEPS - 1):
        t[n+1] = t[n]+h
        K1 = f(t[n], y[n])
        K2 = f(t[n+1], y[n] + K1*h)
        y[n+1] = y[n] + 0.5*(K1 + K2)*h  # fórmula iterativa do método de euler

    return (t, y)


def EDO_heun_sistemas(f, r0, t0, NUMBER_OF_STEPS=100, h=0.01):

    NUMBER_OF_EQUATIONS = len(r0)

    x = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)
    y = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)
    t = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)

    x[0], y[0] = r0
    t[0] = t0

    for n in range(0, NUMBER_OF_STEPS - 1):
        t[n+1] = t[n]+h
        K1x, K1y = fx(t[n],x[n] ,y[n])
        K2x, K2y = fx(t[n+1],x[n] + K1x*h ,y[n] + K1y*h )

        x[n+1] = x[n] + 0.5*(K1x + K2x)*h
        y[n+1] = y[n] + 0.5*(K1y + K2y)*h

    return (t, x, y)

t,x,y = EDO_heun_sistemas(fx,(10,10), 0 , NUMBER_OF_STEPS = 1000, h = 0.05)

plt.plot(t, x)
plt.plot(t, y)

#te, ye = EDO_euler(fx, 1, 0)
#th, yh = EDO_heun(fx, 1, 0)

#plt.plot(te, ye, "ro")
#plt.plot(th, yh, "bs")

plt.show()
