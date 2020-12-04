import numpy as np

"""
Entendenddo o método:
---------------------

q = alpha*((tolerancia*h)/abs(y4 - y5))**(1/k)

onde:
alpha < 1
k é a ordem da edo


Se q < 1 é preciso reduzir o passo(h) original para ter o erro
abaixo da tolerância, logo o passo y[n+1] precisa ser refeito.

Se q > 1 pode-se aumentar o passo(h) que o erro ainda ficará menor que 
a tolerância

"""


def EDO_rkf_sistemas(f, r0, t0, NUMBER_OF_STEPS=100, h=0.01, alpha = 0.86, tol = 0.00001, k = 4):
    NUMBER_OF_EQUATIONS = len(r0)

    r = np.zeros( [NUMBER_OF_STEPS,NUMBER_OF_EQUATIONS], dtype=np.float32 )
    t = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)
    y4 = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)
    y5 = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)
    q = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)

    q_minimo = 0

    r[0,:] = r0
    t[0] = t0

    for n in range(0, NUMBER_OF_STEPS - 1):

        while(True):
            K1 = f(t[n],r[n])
            K2 = f(t[n] +   (1/4)*h, r[n] + h*(1/4)*K1)
            K3 = f(t[n] +   (3/8)*h, r[n] + h*((3*K1 + 9*K2)/32))
            K4 = f(t[n] + (12/13)*h, r[n] + h*((1932*K1 - 7200*K2 + 7296*K3)/2197))

            K5 = f(t[n] + h, r[n] + h*((439/216)*K1 - 8*K2 + (3680/513)*K3 - (845/4104)*K4))
            K6 = f(t[n] + (1/2)*h, r[n] + h*(-(8/27)*K1 + 2*K2 - (3544/2565)*K3 + (1859/4104)*K4 - (11/40)*K5))

            y5 = r[n] + ((16/135)*K1 + (6656/12825)*K3 + (28561/56430)*K4 - (9/50)*K5 + (2/55)*K6)*h
            y4 = r[n] + ((25/216)*K1 + (1408/2565)*K3 + (2197/4104)*K4 - (1/5)*K5)*h

            q = alpha*((tol*h)/abs(y4 - y5))**(1/k)
            
            q_minimo = q[0]
            for i in range(1, NUMBER_OF_EQUATIONS-1):
                if q_minimo > q[i]:
                    q_minimo = q[i]

            if q_minimo >= 1 : break 
            else: h = q_minimo*h

        t[n+1] = t[n]+h
        r[n+1] = y5
        h = q_minimo*h
        #print ("h =", h)


    return (t, r)
