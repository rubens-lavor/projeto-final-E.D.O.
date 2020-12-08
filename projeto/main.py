import numpy as np
import matplotlib.pylab as plt
import math as mt

import euler
import heun
import runge_kutta as rk
import runge_kutta_fehlberg as rkf


def rt(t):
    return (1/10 * mt.exp(-(t-10)**2/7))


def f(t,r):
    x1,x2,x3,x4 = r

    y1 = x2
    dy1 = -(7.5 + 50) / 25 * x1 - 30/25 * x2 + 7.5/25 * x3 + 30/25 * x4 + 50/25 * rt(t)
    y2 = x4
    dy2 = -7.5/300 * x1 + 30/300 * x2 - 7.5/300 * x3 - 30 / 300 * x4

    return np.array([y1, dy1, y2, dy2])



t, r = euler.EDO_euler_sistemas(f,(0,0,0,0),0,NUMBER_OF_STEPS=200, h=0.5)
tempo = t
pos1 = r[:,1] - r[:,3]
vel1 = r[:,0] - r[:,2]

t, r = heun.EDO_heun_sistemas(f,(0,0,0,0),0,NUMBER_OF_STEPS=200, h=0.5)
pos2 = r[:,1] - r[:,3]
vel2 = r[:,0] - r[:,2]


t, r = rk.EDO_rk_sistemas(f,(0,0,0,0),0,NUMBER_OF_STEPS=200, h=0.5)
pos3 = r[:,1] - r[:,3]
vel3 = r[:,0] - r[:,2]

t, r = rkf.EDO_rkf_sistemas(f,(0,0,0,0),0,NUMBER_OF_STEPS=200, h=0.5)
pos4 = r[:,1] - r[:,3]
vel4 = r[:,0] - r[:,2]


plt.subplot(211)
plt.plot( tempo, pos1, color = 'red', label = 'Euler')
plt.plot (tempo, pos2, color = 'green', label = 'Euler Melhorado')
plt.plot(tempo, pos3, color='blue', label='Runge-Kutta')
plt.plot(tempo, pos4, color='yellow', label='Runge-Kutta-Fehlberg')
plt.legend()
plt.title('POSIÇÃO PERCEBIDA PELA MOLA ks')
plt.xlabel('tempo (s)')
plt.ylabel('posição (m) ')

plt.subplot(212)

plt.plot(tempo, vel1, color='red', label='Euler')
plt.plot(tempo, vel2, color='green', label='Euler Melhorado')
plt.plot(tempo, vel3, color='blue', label='Runge-Kutta')
plt.plot(tempo, vel4, color='yellow', label='Runge-Kutta-Fehlberg')
plt.xlabel('tempo (s)')
plt.ylabel('velocidade (m/s)')
plt.title('VELOCIDADE PERCEBIDA PELA MOLA ks')
plt.legend()

plt.rcParams['figure.figsize'] = (15, 10)
# plt.savefig('grafico.png')

plt.show()
