# método de euler
# y_{n+1} = y_n + f(t_n,y_n)h


import numpy as np
import matplotlib.pylab as plt
import math as mt  # função rt(t) usa exponencial


"""
Exemplo:

# EDO a ser resolvida:
# y' = 1 - t + 4y
# y(0) = 1

K1 = f(t_n,y_n)
y_{n+1} = y_n + K1*h

def f(t, y):  # FUNÇÃO QUE DEFINE A EDO (derivada)
    return 1 - t + 4*y

def EDO_euler(f, y0, t0, NUMBER_OF_STEPS=200, h=0.01):
    # recebe a função derivada, y(0) e t(0),
    # os demais paramentros são definidos por padrão, caso não sejam definidos na chamada da função

    y = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)
    t = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)

    y[0] = y0
    t[0] = t0

    for n in range(0, NUMBER_OF_STEPS - 1):
        K1 = f(t[n], y[n])
        y[n+1] = y[n] + K1*h  # fórmula iterativa do método de euler
        t[n+1] = t[n]+h  # próximo passo de tempo

    return (t, y)
"""


"""
def rt(t):  # O que é essa função?
    return (1/10 * mt.exp(-(t-10)**2/7))


def f(x1, x2, x3, x4, t):  # FUNÇÃO QUE DEFINE A EDO (derivada)
    y1 = x2
    dy1 = -(7.5 + 50) / 25 * x1 - 30/25 * x2 + \
        7.5/25 * x3 + 30/25 * x4 + 50/25 * rt(t)
    y2 = x4
    dy2 = -7.5/300 * x1 + 30/300 * x2 - 7.5/300 * x3 - 30 / 300 * x4

    return(y1, dy1, y2, dy2)
"""

def rt(t):
    return (1/10 * mt.exp(-(t-10)**2/7))

def f(t,r):
    x1,x2,x3,x4 = r

    y1 = x2
    dy1 = -(7.5 + 50) / 25 * x1 - 30/25 * x2 + 7.5/25 * x3 + 30/25 * x4 + 50/25 * rt(t)
    y2 = x4
    dy2 = -7.5/300 * x1 + 30/300 * x2 - 7.5/300 * x3 - 30 / 300 * x4

    return np.array([y1, dy1, y2, dy2])


def EDO_euler_sistemas(f, r0, t0, NUMBER_OF_STEPS=100, h=0.01):
    NUMBER_OF_EQUATIONS = len(r0)

    r = np.zeros( [NUMBER_OF_STEPS,NUMBER_OF_EQUATIONS], dtype=np.float32 )
    t = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)

    r[0] = r0
    t[0] = t0

    for n in range(0, NUMBER_OF_STEPS - 1):
        t[n+1] = t[n]+h
        K1 = f(t[n], r[n])
        r[n+1] = r[n] + K1*h

    return (t, r)

t, r = EDO_euler_sistemas(f,(0,0,0,0),0,NUMBER_OF_STEPS=200, h=0.5)

pos1 = r[:,1] - r[:,3]
vel1 = r[:,0] - r[:,2]

plt.subplot(211)
plt.plot( t, pos1, color = 'red', label = 'Euler')
plt.legend()
plt.show()

"""

def euler(f,h, NUMBER_OF_STEPS):
    y1 = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)
    y2 = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)
    y3 = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)
    y4 = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)
    t = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)
    pos = np.zeros(4)
    vel = np.zeros(4)

    #os valores iniciais de y1,y2,y3,y4 e t, são zero na posição zero.
    #y1[0] = 0
    #y2[0] = 0
    #y3[0] = 0
    #y4[0] = 0
    #t[0] = 0

    for n in range(0, NUMBER_OF_STEPS-1):
        K1 = f(y1[n], y2[n], y3[n], y4[n], t[n])
        y1[n+1] = y1[n] + K1[0] * h
        y2[n+1] = y2[n] + K1[1] * h
        y3[n+1] = y3[n] + K1[2] * h
        y4[n+1] = y4[n] + K1[3] * h
        t[n+1] = t[n] + h

    pos = y2 - y4
    vel = y1 - y3

    return (pos, vel, t)



h = 0.5
NS = 200
#vet1 = euler(h, NS)
pos1, vel1, tempo = euler(h, NS)
#pos1 = vet1[0]
#tempo = vet1[2]
#vel1 = vet1[1]

plt.subplot(211)
plt.plot( tempo, pos1, color = 'red', label = 'Euler')
#plt.show()
plt.legend()

# recebe a função derivada, y(0) e t(0),
# os demais paramentros (NUMBER_OF_STEPS=100, h=0.01) são definidos, não precisa passar
#t, y = EDO_euler(f, 1, 0)
#plt.plot(t, y)
plt.show()

"""