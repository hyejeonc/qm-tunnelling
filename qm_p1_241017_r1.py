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
sig_x = 0.5 # Bandwidth
Nx = 100 # Discrete points
x = np.linspace(x_s, L, Nx)
dx = x[1] - x[0]
dt = 0.001

# Constants needed to calculate
w = m.sqrt(L/ma)
E = (hbar**2) * k_0**2 / 2 * ma
dx = L/(Nx - 1)

prob_fun = np.exp(-(x - x_s)**2/2 * sig_x**2)**2
psipsi = np.sum(prob_fun * dx)
C = 1/np.sqrt(psipsi)

# Problem 1

t = 0 
psi = C * np.exp(-(x - x_s)**2/(2 * sig_x**2)) \
    * np.exp(1j * (k0 * x - omega * t))
#for Psi_imag
psi.imag = np.conj(C*np.exp(-(x - x_s)**2/2 * sig_x**2) \
    * np.exp(1j * (k_0 * x - w * t)))

#for Psi_real
t = dt/2 
psi.real = np.real(C * np.exp(-(x - x_s)**2/2 * sig_x**2) * np.exp(1j * (k_0*x-w*t)))

prob=psi*psi.conj()

plt.plot(x,psi.real)
plt.plot(x,psi.imag)
plt.plot(x,prob)
plt.show()
plt.savefig('ploblem1.png')


