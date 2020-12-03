import numpy as np

"""

K1 = f(t_n,y_n)

K2 = f(t_{n+1}, y_n + h*K1)

y_{n+1} = y_n + h*((K1 + K2)/2)

"""

def EDO_heun_sistemas(f, r0, t0, NUMBER_OF_STEPS=100, h=0.01):
    NUMBER_OF_EQUATIONS = len(r0)

    r = np.zeros( [NUMBER_OF_STEPS,NUMBER_OF_EQUATIONS], dtype=np.float32 )
    t = np.zeros(NUMBER_OF_STEPS, dtype=np.float32)

    r[0,:] = r0
    t[0] = t0

    for n in range(0, NUMBER_OF_STEPS - 1):
        t[n+1] = t[n]+h
        K1 = f(t[n], r[n,:])
        K2 = f(t[n+1], r[n,:] + K1*h)
        r[n+1] = r[n,:] + 0.5*(K1 + K2)*h

    return (t, r)

