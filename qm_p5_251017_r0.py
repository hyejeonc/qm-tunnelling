# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 08:46:45 2017

@author: hyejeonc
"""

import numpy as np
import matplotlib.pyplot as plt
#matplotlib inline

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

# Problem1
t = 0
psi = C*np.exp(-(x-x0)**2/(2*sigmax**2))*np.exp(1j*(k0*x-omega*t))
t = dt/2 
psi.real = np.real( C*np.exp(-(x-x0)**2/(2*sigmax**2))*np.exp(1j*(k0*x-omega*t)) )


probability = psi*psi.conj()
C = np.sum(probability*dx)
psi = psi / np.sqrt(C)

probability = psi*psi.conj()
C = np.sum(probability*dx)


#Problem2
vg = k0*hbar/m
T = (L/2)/vg
Nt = int(T / dt)


#Problem3
ll = np.linspace(0,L/20,50)
r_save = np.zeros_like(ll)
t_save = np.zeros_like(ll)

v0 = 9*E/10
v = np.zeros(Nx)
l = L/50


for q in range(50):
    l = ll[q]
    
    for j in range(Nx):
        if ((dx * j) > (L/2)-(l/2)) and ((dx * j) < (L/2)+(l/2)):
            v[j] = v0
    
    t = 0
    psi = C*np.exp(-(x-x0)**2/(2*sigmax**2))*np.exp(1j*(k0*x-omega*t))
    t = dt/2
    psi.real = np.real( C*np.exp(-(x-x0)**2/(2*sigmax**2))*np.exp(1j*(k0*x-omega*t)) )
    
    
    probability=psi*psi.conj()
    C = np.sum(probability*dx)
    psi = psi / np.sqrt(C)
    
    probability = psi*psi.conj()
    C = np.sum(probability*dx)
    
    #loop of time
    for i in range(Nt):
        # loop of position
        for j in range(1, Nx-1):
            psi.imag[j] = psi.imag[j] - dt*(((v[j]/hbar)*psi.real[j])-(hbar/(2*m))*(psi.real[j+1]-2*psi.real[j]+psi.real[j-1])/dx**2)
            
        for j in range(1, Nx-1):
            psi.real[j] = psi.real[j] + dt*(((v[j]/hbar)*psi.imag[j])-(hbar/(2*m))*(psi.imag[j+1]-2*psi.imag[j]+psi.imag[j-1])/dx**2)
    
    
    probability = psi*psi.conj()
    
    

    #problem 5 
    
    #v = np.zeros(Nx)
    #V0= 0   
    
    #probability=psi*psi.conj()
    V0 = np.sum(probability*dx)
    
    reflected = np.sum(probability[0:int(Nx/2)]*dx)
    transmitted = np.sum(probability[int(Nx/2):int(Nx-1)]*dx)

    r_save[q] = reflected
    t_save[q] = transmitted


print(r_save)
print(t_save)

plt.plot(ll, r_save)
plt.plot(ll, t_save)
plt.savefig('problem5.png')


