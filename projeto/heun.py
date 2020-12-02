import numpy as np
import matplotlib.pylab as plt
import math as mt


"""

K1 = f(t_n,y_n)

K2 = f(t_{n+1}, y_n + h*K1)

y_{n+1} = y_n + h*((K1 + K2)/2)
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


def EDO_heun_sistemas(f, r0, t0, NUMBER_OF_STEPS=100, h=0.01):
    # recebe a função derivada, y(0) e t(0),
    # os demais paramentros são definidos, não precisa passar
    NUMBER_OF_EQUATIONS = len(r0)

    r = np.zeros( [NUMBER_OF_STEPS,NUMBER_OF_EQUATIONS], dtype=np.float32 )
    t = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)
    #y = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)

    #x=r[:,0]
    #y=r[:,1]

    r[0,:] = r0
    t[0] = t0

    for n in range(0, NUMBER_OF_STEPS - 1):
        t[n+1] = t[n]+h
        K1 = f(t[n], r[n,:])
        K2 = f(t[n+1], r[n,:] + K1*h)
        r[n+1] = r[n,:] + 0.5*(K1 + K2)*h
        #y[n+1] =  y[n] + 0.5*(K1 + K2)*h # fórmula iterativa do método de euler melhorado
        #y[n+1] == r[n+1], calculado acima.

    posicao = r[:,1] - r[:,3]
    velocidade = r[:,0] - r[:,2]
    return (t, posicao, velocidade)

t, r = EDO_heun_sistemas(f,(0,0,0,0),0,NUMBER_OF_STEPS=200, h=0.5)

pos1 = r[:,1] - r[:,3]
vel1 = r[:,0] - r[:,2]

plt.subplot(211)
plt.plot( t, pos1, color = 'red', label = 'Euler Melhorado')
plt.legend()

plt.show()

"""

def rt(t):
    return (1/10 * mt.exp(-(t-10)**2/7))


def f(x1, x2, x3, x4, t):
    y1 = x2
    dy1 = -(7.5 + 50) / 25 * x1 - 30/25 * x2 + \
        7.5/25 * x3 + 30/25 * x4 + 50/25 * rt(t)
    y2 = x4
    dy2 = -7.5/300 * x1 + 30/300 * x2 - 7.5/300 * x3 - 30 / 300 * x4

    return(y1, dy1, y2, dy2)

def eulermelhorado(f,h, NUMBER_OF_STEPS):

    x1 = np.zeros(NUMBER_OF_STEPS)
    x2 = np.zeros(NUMBER_OF_STEPS)
    x3 = np.zeros(NUMBER_OF_STEPS)
    x4 = np.zeros(NUMBER_OF_STEPS)
    t = np.zeros(NUMBER_OF_STEPS)
    vel = np.zeros(4)
    pos = np.zeros(4)

    x1[0] = 0
    x2[0] = 0
    x3[0] = 0
    x4[0] = 0
    t[0] = 0
   
    for n in range(0, NUMBER_OF_STEPS - 1):
        K1 = f(x1[n], x2[n], x3[n], x4[n], t[n])
        K2 = f(x1[n] + K1[0] * h, x2[n] + K1[1] * h, x3[n] + K1[2] * h, x4[n] + K1[3] * h, t[n])
        x1[n+1] = x1[n] + 0.5 * (K1[0] + K2[0]) * h
        x2[n+1] = x2[n] + 0.5 * (K1[1] + K2[1]) * h
        x3[n+1] = x3[n] + 0.5 * (K1[2] + K2[2]) * h
        x4[n+1] = x4[n] + 0.5 * (K1[3] + K2[3]) * h
        t[n+1] = t[n] + h

        pos = x2 - x4
        vel = x1 - x3

    return (pos, vel, t)


h = 0.5
NS = 200

pos1, vel1, tempo = eulermelhorado(f,h, NS)

plt.subplot(211)
plt.plot( tempo, pos1, color = 'red', label = 'Euler Melhorado')
plt.legend()

plt.show()

"""