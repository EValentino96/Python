import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

actual = 2-(112/np.exp(10))
a = 0
b = np.pi/2
n = 8
x = np.linspace(a,b,n)
h = x[1]-x[0]

##---------- Trapezoidal ----------##
w = np.ones(n)*h
w[0] = w[1] = h/2

def f(x):
    return np.sin(x)*np.sin(x)

sum = 0
for i in range(0,n):
    sum = sum + w[i]*f(x[i])

err = np.abs(sum - np.pi/4)

print(err)

##---------- Quadrature ----------##

w = np.ones(n)*h/3
for i in range(1,n-1):
    if(i%2 == 0):
        w[i] *= 2
    else:
        w[i] *= 4

sum = 0
for i in range(0,n):
    sum = sum + w[i]*f(x[i])

print(actual - sum)

##---------- Monte Carlo ----------##

H = 0.5
bbox = (b -a)*H
nrand = 100
hits = 0

for i in range(0,nrand):
    xr = a+b*np.random.rand()
    yr = np.random.rand()*H
    if(yr < f(xr)):
        hits += 1

integral = bbox*hits/nrand

print(integral - actual)
