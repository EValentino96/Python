import numpy as np

##---------- Radioactive Decay ----------##

def getNumberPoints(tStart,tEnd,timeStep):
    return np.int((tEnd-tStart)/timeStep + 1)

def radDecay():
    #Speeding through this
    k = 0.1
    r = 10

    tStart = 0
    tEnd = 10
    h = 1e-3
    n = getNumberPoints(tStart,tEnd,h)

    t = np.linspace(tStart,tEnd,n)
    rArray = np.zeros(n)

    rArray[0] = r

    for j in range(0,n):
        mid = r + (h/2)*(-k*r) 
        r += h*(-k*mid)
        rArray[j] = r
    
    print(t,"\n",rArray)

radDecay()
