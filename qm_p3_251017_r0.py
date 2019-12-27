# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 14:28:48 2017

@author: hyejeonc
"""

import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline

L = 20
hbar = 1
k0 = 20
m = 1
E = hbar**2*k0**2/(2*m)
omega = E/hbar
sigmax = 1
x0= 5
C = 1



Nx = 1000
x=np.linspace(0, L, Nx)
dx = x[1]-x[0]
dt = 0.0001

#Problem1

t = 0
psi = C*np.exp(-(x-x0)**2/(2*sigmax**2))*np.exp(1j*(k0*x-omega*t))
t = dt/2
psi.real = np.real( C*np.exp(-(x-x0)**2/(2*sigmax**2))*np.exp(1j*(k0*x-omega*t)) )


probability=psi*psi.conj()
C = np.sum(probability*dx)
psi = psi / np.sqrt(C)

#Problem3
l = L/50
v0 = E/2

for j in range(Nx):
    if ((dx * j) > (L/2)-(l/2)) and ((dx * j) < (L/2)+(l/2)):
        v[j] =v0

#Problem2
vg = k0*hbar/m
T = (L/2)/vg
Nt = int(T / dt)
#potential
v = np.zeros(Nx) #free particles

#loop of time
for i in range(Nt):
    dt=dt*i
    # loop of position
    for j in range(1, Nx-1):
        psi.imag[j] = psi.imag[j] - dt*(((v[j]/hbar)*psi.real[j])-(hbar/(2*m))*(psi.real[j+1]-2*psi.real[j]+psi.real[j-1])/dx**2)
        
    for j in range(1, Nx-1):
        psi.real[j] = psi.real[j] + dt*(((v[j]/hbar)*psi.imag[j])-(hbar/(2*m))*(psi.imag[j+1]-2*psi.imag[j]+psi.imag[j-1])/dx**2)




probability = psi*psi.conj()
C = np.sum(probability*dx)

plt.plot(x, psi.real)
plt.plot(x, psi.imag)
plt.plot(x, probability)
plt.plot(x, v/200)
plt.plot(x, v)
plt.savefig('problem3.png')

