import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

##---------- Time Independent Schrödinger Equation ----------##

n = 1000
a = 2
hbar = 1
m = 2
V1 = 100
V2 = 200
## we really have a spacing of 2a/(n-3)
hfix = 2*a/(n + 1)
### make x spread out more
x = np.linspace(-4*a,4*a,n)
deltaX = x[1]-x[0]

D2 = np.zeros(shape = (n,n))

for i in range(0,n):
    D2[i,i] = -2
    if(i - 1 > -1):
        D2[i, i - 1]=1
    if(i + 1 < n):
        D2[i, i + 1]=1
D2 /= deltaX**2

## build a V matrix
def Vfunc(x):
    if(x < -a):
        return V1
    elif(x >= -2 and x <= 2):
        return 0
    elif(x > a):
        return V2
    
V = np.zeros(shape = (n,n))
VVec = np.zeros(shape = (n))
for i in range(0,n):
    V[i,i] = Vfunc(x[i])
    VVec[i] = V[i,i]

## add in V to H
H = -hbar**2/2/m * D2 + V
val, vec = la.eig(H)

isort = np.argsort(val)

ishow=0
Area = vec[0,isort[ishow]]**2*deltaX
for i in range(1,n-1):
    Area += vec[i,isort[ishow]]**2*deltaX
Area += vec[n-1,isort[ishow]]**2*deltaX
vec[:,isort[ishow]] /= np.sqrt(Area)

print(np.sort(val)[0:5])
plt.plot(x,100*vec[:,isort[ishow]])
plt.plot(x,VVec)
plt.show()

