

# método de euler
# y_{n+1} = y_n + f(t_n,y_n)h

# EDO a ser resolvida:
# y' = 1 - t + 4y
# y(0) = 1

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


t, y = EDO_euler(f, 1, 0) 
plt.plot(t, y)
plt.show()
