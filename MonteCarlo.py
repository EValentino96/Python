import numpy as np
import matplotlib.pyplot as plt

##---------- Monte Carlo ---------##
# Here we want to estimate pi

def p(x): # probability function
    if(x < 0 or x > 1):
        return 0
    else:
        return 1

def g(x): # actual function
    return np.sqrt(1 - x**2)

n = 1000

for j in range(0,10):

    xs = []
    
    for i in range(0,n):
        x = np.random.rand()
        xs.append(x)

    integral = 0

    for i in range(0,n):
        integral += g(xs[i])/p(xs[i])

    integral /= n

    print(integral*4) # we're supposed to get pi/4, scale to show pi
