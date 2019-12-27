# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 03:32:14 2017

@author: HYEJEONG
"""


import matplotlib.pyplot as plt
import numpy as np
#import tensorflow as tf
import math as m
#import sympy as sp

# Constants not changed
hbar = 1
ma = 1
k_0 = 20
L = 20

# Constants for setting
x_s = 5 # Initial position of function x
sig_x = 0.5 # Bandwidth ;0.5 1.0 2.0  & 
# x_s = 15  sig_x = 0.5  1.5  
Nx = 100 # Discrete points
x = np.linspace(x_s, L, Nx)
dx = x[1] - x[0]
dt = 0.001
V = np.zeros(Nx) # Potential, free particles

# Constants needed to calculate
w = m.sqrt(L/ma)
E = (hbar**2) * k_0**2/2 * ma
dx = L/(Nx - 1)

prob_fun = np.exp(-(x - x_s)**2/2 * sig_x**2)**2
psipsi = np.sum(prob_fun*dx)
C = 1/np.sqrt(psipsi)

psi.real = np.zeros(Nx)
psi.imag = np.zeros((Nx), dtype=np.complex)
dr = psi.real[2]-psi.real[1]
di = psi.imag[2]-psi.imag[1]

# Problem 2

t = 0 #stationary state for Psi_imaginary
psi = C * np.exp(-(x - x_s)**2/(2 * sig_x**2)) \
    * np.exp(1j * (k0 * x - omega * t))

psi.imag = np.imag(C*np.exp(-(x-x_s)**2/2*sig_x**2)*np.exp(1j*(k_0*x-w*t)))

t=dt/2 #for Psi_real
psi.real = np.real(C*np.exp(-(x-x_s)**2/2*sig_x**2)*np.exp(1j*(k_0*x-w*t)))
###############################################################################
#New constant for Problem 2

vg = hbar*k_0/ma
T = L/2*vg

V = np.zeros(Nx) #Potential energy, free particles
Nt = T/dt


for l in range(0, int(Nt)):
    l = l + l*dt
    
    for i in range(x_s, Nx):
        di = -dt*(((V/hbar)*psi.real[i])-((hbar/2*ma)*(psi.real[i+1]-2*psi.real[i]+psi.real[i-1])/(dx**2)))
        #psi.imag[i]=psi.imag[i]-dt*((V/hbar)*psi.real[i]-(hbar/2*ma)*(psi.real[i+1]-2*psi.real[i]+psi.real[i-1])/(dx**2))
    psi.imag[i]=psi.imag[i]+di
    
    for i in range(x_s, Nx):
        dr = dt*((V/hbar)*psi.imag[i]-(hbar/2*ma)*(psi.imag[i+1]-2*psi.imag[i]+psi.imag[i-1])/(dx**2))
        #psi.real[i]=psi.real[i]+dt*((V/hbar)*psi.imag[i]-(hbar/2*ma)*(psi.imag[i+1]-2*psi.imag[i]+psi.imag[i-1])/(dx**2))
    psi.real[i] = psi.real[i] + dr


prob = psi*psi.conj()

plt.plot(x, psi.real)
plt.plot(x, psi.imag)
plt.plot(x, prob)
plt.show()

plt.savefig('ploblem2.png')
