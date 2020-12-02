import numpy as np
import matplotlib.pylab as plt
import math as mt

"""
for n in range(0, NUMBER_OF_STEPS-1):
t[n+1] = t[n] + h

r[n+1]=r[n]+(1/6)*(K1 + 2*K2 + 2*K3 + K4)


K1 = f(t[n],r[n])
K2 = f(t[n] + (1/2)*h , r[n] + (1/2)*K1*h)
K3 = f(t[n] + (1/2)*h , r[n] + (1/2)*K2*h)
K4 = f(t[n] + h , r[n] + K3*h)


y[n+1]=y[n]+(1/6)*(K1 + 2*K2 + 2*K3 + K4)


K1 = f(t[n],y[n])
K2 = f(t[n] + (1/2)*h , y[n] + (1/2)*K1*h)
K3 = f(t[n] + (1/2)*h , y[n] + (1/2)*K2*h)
K4 = f(t[n] + h , y[n] + K3*h)
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


def EDO_rk_sistemas(f, r0, t0, NUMBER_OF_STEPS=100, h=0.01):
    NUMBER_OF_EQUATIONS = len(r0)

    r = np.zeros( [NUMBER_OF_STEPS,NUMBER_OF_EQUATIONS], dtype=np.float32 )
    t = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)

    r[0,:] = r0
    t[0] = t0

    for n in range(0, NUMBER_OF_STEPS - 1):
        t[n+1] = t[n]+h

        K1 = f(t[n],r[n])
        K2 = f(t[n] + (1/2)*h , r[n] + (1/2)*K1*h)
        K3 = f(t[n] + (1/2)*h , r[n] + (1/2)*K2*h)
        K4 = f(t[n+1] , r[n] + K3*h)
        r[n+1]=r[n]+(1/6)*(K1 + 2*K2 + 2*K3 + K4)

    posicao = r[:,1] - r[:,3]
    velocidade = r[:,0] - r[:,2]
    return (t, posicao, velocidade)

t, r = EDO_rk_sistemas(f,(0,0,0,0),0,NUMBER_OF_STEPS=200, h=0.5)

pos1 = r[:,1] - r[:,3]
vel1 = r[:,0] - r[:,2]

plt.subplot(211)
plt.plot( t, pos1, color = 'g', label = 'Runge-Kutta')
plt.legend()

plt.show()


"""
def rungekutta(f,h, NUMBER_OF_STEPS):

    x1 = np.zeros(NUMBER_OF_STEPS)
    x2 = np.zeros(NUMBER_OF_STEPS)
    x3 = np.zeros(NUMBER_OF_STEPS)
    x4 = np.zeros(NUMBER_OF_STEPS)
    t = np.zeros(NUMBER_OF_STEPS)
    K = np.zeros(4)
    K1 = np.zeros(4)
    K2 = np.zeros(4)
    K3 = np.zeros(4)
    K4 = np.zeros(4)
    vel = np.zeros(4)
    pos = np.zeros(4)

    x1[0] = 0
    x2[0] = 0
    x3[0] = 0
    x4[0] = 0
    t[0] = 0

    for n in range(1, NUMBER_OF_STEPS):

        K1 = f(x1[n-1], x2[n-1], x3[n-1], x4[n-1], t[n-1])
        K2 = f(x1[n-1] + K1[0] * 0.5 * h, x2[n-1] + K1[1] * 0.5 * h, x3[n-1] +
               K1[2] * 0.5 * h, x4[n-1] + K1[3] * 0.5 * h, t[n-1] + 0.5 * h)
        K3 = f(x1[n-1] + K2[0] * 0.5 * h, x2[n-1] + K2[1] * 0.5 * h, x3[n-1] +
               K2[2] * 0.5 * h, x4[n-1] + K2[3] * 0.5 * h, t[n-1] + 0.5 * h)
        K4 = f(x1[n-1] + K3[0] * h, x2[n-1] + K3[1] * h, x3[n-1] +
               K3[2] * h, x4[n-1] + K3[3] * h, t[n-1] + h)

        K[0] = 1/6 * (K1[0] + 2 * K2[0] + 2 * K3[0] + K4[0])
        K[1] = 1/6 * (K1[1] + 2 * K2[1] + 2 * K3[1] + K4[1])
        K[2] = 1/6 * (K1[2] + 2 * K2[2] + 2 * K3[2] + K4[2])
        K[3] = 1/6 * (K1[3] + 2 * K2[3] + 2 * K3[3] + K4[3])

        x1[n] = x1[n-1] + K[0] * h
        x2[n] = x2[n-1] + K[1] * h
        x3[n] = x3[n-1] + K[2] * h
        x4[n] = x4[n-1] + K[3] * h

        t[n] = t[n-1] + h

    pos = x2 - x4
    vel = x1 - x3
    #print("pos rk = ", pos)
    #print("vel rk = ", vel)
    return (pos, vel)
"""