import numpy as np

"""
K1 = f(t_n,y_n)
y_{n+1} = y_n + K1*h

"""

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
        #y[n+1] == r[n+1], calculado acima.
    
    return (t, r)