import matplotlib.pyplot as plt
import numpy as np
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


'''
def feval(funcName, *args):
    return eval(funcName)(*args)


def RKF45(func, yinit, x_range, h):
    m = len(yinit)
    n = int((x_range[-1] - x_range[0])/h)

    x = x_range[0]
    y = yinit

    xsol = np.empty(0)
    xsol = np.append(xsol, x)

    ysol = np.empty(0)
    ysol = np.append(ysol, y)

    for i in range(n):
        k1 = feval(func, x, y)

        yp2 = y + k1*(h/5)

        k2 = feval(func, x+h/5, yp2)

        yp3 = y + k1*(3*h/40) + k2*(9*h/40)

        k3 = feval(func, x+(3*h/10), yp3)

        yp4 = y + k1*(3*h/10) - k2*(9*h/10) + k3*(6*h/5)

        k4 = feval(func, x+(3*h/5), yp4)

        yp5 = y - k1*(11*h/54) + k2*(5*h/2) - k3*(70*h/27) + k4*(35*h/27)

        k5 = feval(func, x+h, yp5)

        yp6 = y + k1*(1631*h/55296) + k2*(175*h/512) + k3 * \
            (575*h/13824) + k4*(44275*h/110592) + k5*(253*h/4096)

        k6 = feval(func, x+(7*h/8), yp6)

        for j in range(m):
            y[j] = y[j] + h*(37*k1[j]/378 + 250*k3[j]/621 +
                             125*k4[j]/594 + 512*k6[j]/1771)

        x = x + h
        xsol = np.append(xsol, x)

        for r in range(len(y)):
            ysol = np.append(ysol, y[r])

    return [xsol, ysol]
'''

def rungekuttafehlberg(f, h, NUMBER_OF_STEPS):

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
                    status = 1

            if (h == h2):
                h2 = h * q[0]
                for i in range(1, 4):
                    h1 = h * q[i]
                    if (h1 < h2):
                        h2 = h1
                status = 0

            h = h2

    pos = x2 - x4
    vel = x1 - x3
    return (pos, vel, t)



h = 0.5
NS = 200

pos1, vel1, tempo = rungekuttafehlberg(f,h, NS)

plt.subplot(211)
plt.plot( pos1, tempo, color = 'red', label = 'pos - rungekuttafehlberg')
plt.plot( vel1, tempo, color = 'green', label = 'vel - rungekuttafehlberg')
plt.legend()

plt.show()
