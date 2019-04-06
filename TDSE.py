import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import matplotlib.animation as animation

##---------- Time Dependent Schrödinger Equation ----------## 

## Parameters
n = 500
a = 5
Vbarrier = 1000
m = 1
hbar = 1
x = np.linspace(-3*a,3*a,n)
deltaX = x[1]-x[0]
nt = 1000
t = np.linspace(0,10,nt)

# initial conditions
psiR = np.ndarray(n)
psiI = np.ndarray(n)
nLevel = 3
sig = 0.1
mu = 2.5
for i in range(0,n):
    if(np.abs(x[i]) < a):
        psiR[i] = np.sin(nLevel*np.pi/2/a*(x[i]-a))
        #psiR[i] = (1/np.sqrt(2*np.pi*(sig**2))) * np.exp(-((x[i]-mu)**2)/(2*(sig**2)))
    else:
        #psiR[i] = 0.0
        #psiR[i] = np.exp(-20*(x[i])**2)
        psiR[i] = (1/np.sqrt(2*np.pi*(sig**2))) * np.exp(-((x[i] - mu)**2)/(2*(sig**2)))
    psiI[i] = 0.0

# potential function    
V = np.zeros(n)
def Vfunc(x):
    if(x < -a):
        return Vbarrier
    elif(x > a):
        return Vbarrier
    elif(x < -1):
        return 100
    elif(x > 1):
        return 100
    else:
        return 0
    
for i in range(0,n):
    V[i] = Vfunc(x[i])

# matmul(H,psi) as a loop operation
def Hmul(psi):
    result = np.zeros(n)
    #get curvature
    result[0] = -2*psi[0]+psi[1]
    for i in range(1,n-1):
        result[i] = psi[i-1]-2*psi[i]+psi[i+1]
    result[n-1] = psi[n-2]-2*psi[n-1]
    result /= deltaX**2
    # scale D^2
    result *= -hbar**2/2/m
    # add potential energy
    for i in range(0,n):
        result[i] += V[i]*psi[i]
    return result

# timestepping function, rates of change        
def f(psiC,t):
    psiR = psiC[0:n]
    psiI = psiC[n:2*n]
    dpR =  1/hbar*Hmul(psiI)
    dpI =  -1/hbar*Hmul(psiR)
    return np.concatenate((dpR,dpI))

# store real and imag as one long array
psi0 = np.concatenate((psiR,psiI))

# solve problem
sol = odeint(f,psi0,t)

# plot solution

psiR = sol[0][0:n]
psiI = sol[0][n:2*n]
    
fig, ax = plt.subplots()
line0, = ax.plot(x, psiR*psiR+psiI*psiI)
line1, = ax.plot(x, psiR)
line2, = ax.plot(x, psiI)

def animate(i):
    psiR = sol[i][0:n]
    psiI = sol[i][n:2*n]
    line0.set_ydata(psiR*psiR+psiI*psiI)  # update the data
    line1.set_ydata(psiR)  # update the data
    line2.set_ydata(psiI)
    return line0,line1,line2,

ani = animation.FuncAnimation(fig,
                              animate,
                              np.arange(0,nt),
                              interval=1,
                              blit=True)
plt.show()

### suppose you want to determine probability between two values
## normalize last time step
psiRFinal = sol[-1][0:n]
psiIFinal = sol[-1][n:2*n]
## trapezoidal rule integration weights
weight = np.ones(n)
weight[0]=0.5
weight[-1]=0.5
area = np.sum(weight*deltaX*(psiRFinal*psiRFinal+psiIFinal*psiIFinal))
psiRFinal /= np.sqrt(area)
psiIFinal /= np.sqrt(area)

## integrate over range
prob = 0
probA = 1
probB = 5
for i in range(0,n):
    ## integrate with reimann sum
    if(x[i] > probA and x[i] < probB):
        prob += deltaX*(psiRFinal[i]*psiRFinal[i]+psiIFinal[i]*psiIFinal[i])

print(prob)
