import numpy as np
import matplotlib.pyplot as plt

##---------- Metropolis Hastings ----------##

for i in range(0,5):

    x = 0
    sigma = 2
    mu = 1
    xs = []
    n = 100000
    stepMax = 0.1

    def p(x):
        return 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x - mu)**2/(2*sigma**2)) # gauss curve w/ area of 1

    for i in range(0,n):

        dx = stepMax*2*(np.random.rand() - 0.5) # rand num from -0.1 to 0.1
        testx = x + dx
    
        if(p(testx)/p(x) > 1):
            x = testx
        elif(np.random.rand() < p(testx)/p(x)):
            x = testx
    
        xs.append(x)

    def g(x):
        if(x < 0):
            return 0
        return x**2*np.exp(-x**2)

    integral = 0
    for i in range(0,n):
        integral += g(xs[i])/p(xs[i])
    integral /= n
    print(integral)
