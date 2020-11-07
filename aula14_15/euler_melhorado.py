"""
K1 = f(t_n,y_n)

K2 = f(t_{n+1}, y_n + h*K1)

y_{n+1} = y_n + h*((K1 + K2)/2)

"""

import numpy as np
import matplotlib.pylab as plt


def f(t, y):  # FUNÇÃO QUE DEFINE A EDO (derivada)
    return 1 - t + 4*y


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
        y[n+1] =  y[n] + 0.5*(K1 + K2)*h # fórmula iterativa do método de euler


    return (t, y)

te, ye = EDO_euler(f, 1, 0)
th, yh = EDO_heun(f, 1, 0)

plt.plot(te, ye, "ro")
plt.plot(th, yh, "bs")
plt.show()
