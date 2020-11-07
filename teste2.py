#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''

Métodos numéricos para solução de E.D.O.s: modelagem da suspensão de um
automóvel

'''

import numpy as np
import matplotlib.pylab as plt
import math as mt


def rt(t):
    return (1/10 * mt.exp(-(t-10)**2/7))


def f(x1, x2, x3, x4, t):
    y1 = x2
    dy1 = -(7.5 + 50) / 25 * x1 - 30/25 * x2 + \
        7.5/25 * x3 + 30/25 * x4 + 50/25 * rt(t)
    y2 = x4
    dy2 = -7.5/300 * x1 + 30/300 * x2 - 7.5/300 * x3 - 30 / 300 * x4

    return(y1, dy1, y2, dy2)


def euler(h, NUMBER_OF_STEPS):
    x1 = np.zeros(NUMBER_OF_STEPS)
    x2 = np.zeros(NUMBER_OF_STEPS)
    x3 = np.zeros(NUMBER_OF_STEPS)
    x4 = np.zeros(NUMBER_OF_STEPS)
    t = np.zeros(NUMBER_OF_STEPS)
    pos = np.zeros(4)
    vel = np.zeros(4)

    x1[0] = 0
    x2[0] = 0
    x3[0] = 0
    x4[0] = 0
    t[0] = 0

    for n in range(1, NUMBER_OF_STEPS):
        vet1 = f(x1[n-1], x2[n-1], x3[n-1], x4[n-1], t[n-1])
        x1[n] = x1[n-1] + vet1[0] * h
        x2[n] = x2[n-1] + vet1[1] * h
        x3[n] = x3[n-1] + vet1[2] * h
        x4[n] = x4[n-1] + vet1[3] * h
        t[n] = t[n-1] + h

    pos = x2 - x4
    vel = x1 - x3

    return (pos, vel, t)


def eulermelhorado(h, NUMBER_OF_STEPS):

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

    for n in range(1, NUMBER_OF_STEPS):
        vet1 = f(x1[n-1], x2[n-1], x3[n-1], x4[n-1], t[n-1])
        vet2 = f(x1[n-1] + vet1[0] * h, x2[n-1] + vet1[1] * h,
                 x3[n-1] + vet1[2] * h, x4[n-1] + vet1[3] * h, t[n-1])
        x1[n] = x1[n-1] + 0.5 * (vet1[0] + vet2[0]) * h
        x2[n] = x2[n-1] + 0.5 * (vet1[1] + vet2[1]) * h
        x3[n] = x3[n-1] + 0.5 * (vet1[2] + vet2[2]) * h
        x4[n] = x4[n-1] + 0.5 * (vet1[3] + vet2[3]) * h
        t[n] = t[n-1] + h

        pos = x2-x4
        vel = x1 - x3

    return (pos, vel)


def rungekutta(h, NUMBER_OF_STEPS):

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

        K1 = f(x1[n-1],                     x2[n-1],
               x3[n-1],                     x4[n-1],                    t[n-1])
        K2 = f(x1[n-1] + K1[0] * 0.5 * h, x2[n-1] + K1[1] * 0.5 * h, x3[n-1] +
               K1[2] * 0.5 * h, x4[n-1] + K1[3] * 0.5 * h, t[n-1] + 0.5 * h)
        K3 = f(x1[n-1] + K2[0] * 0.5 * h, x2[n-1] + K2[1] * 0.5 * h, x3[n-1] +
               K2[2] * 0.5 * h, x4[n-1] + K2[3] * 0.5 * h, t[n-1] + 0.5 * h)
        K4 = f(x1[n-1] + K3[0] * h,        x2[n-1] + K3[1] * h,
               x3[n-1] + K3[2] * h,        x4[n-1] + K3[3] * h,        t[n-1] + h)
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
    return (pos, vel)


def rungekuttafehlberg(h, NUMBER_OF_STEPS):
    x1 = np.zeros(NUMBER_OF_STEPS)
    x2 = np.zeros(NUMBER_OF_STEPS)
    x3 = np.zeros(NUMBER_OF_STEPS)
    x4 = np.zeros(NUMBER_OF_STEPS)
    z1 = np.zeros(NUMBER_OF_STEPS)
    z2 = np.zeros(NUMBER_OF_STEPS)
    z3 = np.zeros(NUMBER_OF_STEPS)
    z4 = np.zeros(NUMBER_OF_STEPS)
    t = np.zeros(NUMBER_OF_STEPS)
    vel = np.zeros(4)
    pos = np.zeros(4)

    x1[0] = 0
    x2[0] = 0
    x3[0] = 0
    x4[0] = 0
    t[0] = 0
    tol = 0.00001

    K1 = np.zeros(4)
    K2 = np.zeros(4)
    K3 = np.zeros(4)
    K4 = np.zeros(4)
    K5 = np.zeros(4)
    K6 = np.zeros(4)
    q = np.zeros(4)
    status = 0

    for n in range(1, NUMBER_OF_STEPS):

        while (status == 0):

            K1 = f(x1[n-1], x2[n-1], x3[n-1], x4[n-1], t[n-1])
            K2 = f(x1[n-1] + h/4 * K1[0], x2[n-1] + h/4 * K1[1],
                   x3[n-1] + h/4 * K1[2], x4[n-1] + h/4 * K1[3], t[n-1] + h/4)
            K3 = f(x1[n-1] + h/32 * (3 * K1[0] + 9 * K2[0]), x2[n-1] + h/32 * (3 * K1[1] + 9 * K2[1]), x3[n-1] +
                   h/32 * (3 * K1[2] + 9 * K2[2]), x4[n-1] + h/32 * (3 * K1[3] + 9 * K2[3]), t[n-1] + 3 * h/8)
            K4 = f(x1[n-1] + h/2197 * (1932 * K1[0] - 7200 * K2[0] + 7296 * K3[0]), x2[n-1] + h/2197 * (1932 * K1[1] - 7200 * K2[1] + 7296 * K3[1]), x3[n-1] +
                   h/2197 * (1932 * K1[2] - 7200 * K2[2] + 7296 * K3[2]), x4[n-1] + h/2197 * (1932 * K1[3] - 7200 * K2[3] + 7296 * K3[3]), t[n-1] + 12/13 * h)
            K5 = f(x1[n-1] + h * (439/216 * K1[0] - 8 * K2[0] + 3680/513 * K3[0] - 845/4104 * K4[0]), x2[n-1] + h * (439/216 * K1[1] - 8 * K2[1] + 3680/513 * K3[1] - 845/4104 * K4[1]),
                   x3[n-1] + h * (439/216 * K1[2] - 8 * K2[2] + 3680/513 * K3[2] - 845/4104 * K4[2]), x4[n-1] + h * (439/216 * K1[3] - 8 * K2[3] + 3680/513 * K3[3] - 845/4104 * K4[3]), t[n-1] + h)
            K6 = f(x1[n-1] + h * (-8/27 * K1[0] + 2 * K2[0] - 3544/2565 * K3[0] + 1859/4104 * K4[0] - 11/40 * K5[0]), x2[n-1] + h * (-8/27 * K1[1] + 2 * K2[1] - 3544/2565 * K3[1] + 1859/4104 * K4[1] - 11/40 * K5[1]), x3[n-1] +
                   h * (-8/27 * K1[2] + 2 * K2[2] - 3544/2565 * K3[2] + 1859/4104 * K4[2] - 11/40 * K5[2]), x4[n-1] + h * (-8/27 * K1[3] + 2 * K2[3] - 3544/2565 * K3[3] + 1859/4104 * K4[3] - 11/40 * K5[3]), t[n-1] + 0.5 * h)

            x1[n] = x1[n-1] + (16/135 * K1[0] + 6656/12825 * K3[0] +
                               28561/56430 * K4[0] - 9/50 * K5[0] + 2/55 * K6[0]) * h
            x2[n] = x2[n-1] + (16/135 * K1[1] + 6656/12825 * K3[1] +
                               28561/56430 * K4[1] - 9/50 * K5[1] + 2/55 * K6[1]) * h
            x3[n] = x3[n-1] + (16/135 * K1[2] + 6656/12825 * K3[2] +
                               28561/56430 * K4[2] - 9/50 * K5[2] + 2/55 * K6[2]) * h
            x4[n] = x4[n-1] + (16/135 * K1[3] + 6656/12825 * K3[3] +
                               28561/56430 * K4[3] - 9/50 * K5[3] + 2/55 * K6[3]) * h

            z1[n] = x1[n-1] + (25/216 * K1[0] + 1408/2565 *
                               K3[0] + 2197/4104 * K4[0] - 1/5 * K5[0]) * h
            z2[n] = x2[n-1] + (25/216 * K1[1] + 1408/2565 *
                               K3[1] + 2197/4104 * K4[1] - 1/5 * K5[1]) * h
            z3[n] = x3[n-1] + (25/216 * K1[2] + 1408/2565 *
                               K3[2] + 2197/4104 * K4[2] - 1/5 * K5[2]) * h
            z4[n] = x4[n-1] + (25/216 * K1[3] + 1408/2565 *
                               K3[3] + 2197/4104 * K4[3] - 1/5 * K5[3]) * h

            q[0] = 0.84 * (tol/abs(z1[n] - x1[n]) * h)**(1/4)
            q[1] = 0.84 * (tol/abs(z2[n] - x2[n]) * h)**(1/4)
            q[2] = 0.84 * (tol/abs(z3[n] - x3[n]) * h)**(1/4)
            q[3] = 0.84 * (tol/abs(z4[n] - x4[n]) * h)**(1/4)
            h2 = h

            for i in range(0, 4):
                if q[i] < 1:
                    h1 = h * q[i]
                    if (h1 < h2):
                        h2 = h1
                    status = 0

            if (h == h2):
                h2 = h * q[0]
                for i in range(1, 4):
                    h1 = h * q[i]
                    if (h1 < h2):
                        h2 = h1
                status = 1

            h = h2

    pos = x2 - x4
    vel = x1 - x3
    return (pos, vel)


h = 0.5
NS = 200

vet1 = euler(h, NS)
pos1 = vet1[0]
tempo = vet1[2]
vel1 = vet1[1]
vet2 = eulermelhorado(h, NS)
pos2 = vet2[0]
vel2 = vet2[1]
vet3 = rungekutta(h, NS)
pos3 = vet3[0]
vel3 = vet3[1]
vet4 = rungekuttafehlberg(h, NS)
pos4 = vet4[0]
vel4 = vet4[1]

plt.subplot(211)
#plt.plot( tempo, pos1, color = 'red', label = 'Euler')
#plt.plot (tempo, pos2, color = 'green', label = 'Euler Melhorado')
plt.plot(tempo, pos3, color='blue', label='Runge-Kutta')
plt.plot(tempo, pos3, color='yellow', label='Runge-Kutta-Fehlberg')
plt.legend()
plt.title('POSIÇÃO PERCEBIDA PELA MOLA ks')
plt.xlabel('tempo (s)')
plt.ylabel('posição (m) ')

plt.subplot(212)

plt.plot(tempo, vel1, color='red', label='Euler')
plt.plot(tempo, vel2, color='green', label='Euler Melhorado')
plt.plot(tempo, vel3, color='blue', label='Runge-Kutta')
plt.plot(tempo, vel3, color='yellow', label='Runge-Kutta-Fehlberg')
plt.xlabel('tempo (s)')
plt.ylabel('velocidade (m/s)')
plt.title('VELOCIDADE PERCEBIDA PELA MOLA ks')
plt.legend()

plt.rcParams['figure.figsize'] = (15, 10)
# plt.savefig('grafico.png')

plt.show()
