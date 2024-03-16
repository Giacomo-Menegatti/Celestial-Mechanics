#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 23:06:43 2023

@author: menegattig
"""

from RungeKutta import RK,RK_Adapt
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

plt.close('all')

#state vector: x1,y1, x2,y2, x3,y3,  vx1,vy1, vx2,vy2, vx3,vy3

m = [3,4,5]
y0 = [1.0,3.0, -2.0,-1.0, 1.0,-1.0,  0.0,0.0, 0.0,0.0,  0.0,0.0]
G = 1

def f(t, y):
    N = 3                   #Number of bodies
    a = np.zeros((N,2))      #Accelerations vector
    for i in range(N):         # i goes form 0 to 2
        for j in range(i+1,N):  # j goes from i+1 to 2
            # rij = ri-rj
            r = np.array([y[2*i]-y[2*j], y[2*i+1]-y[2*j+1]])
            r3 = (r[0]**2+r[1]**2)**1.5
            
            a[i] += -G*m[j]/r3*r
            a[j] += G*m[i]/r3*r
            
    return np.array([y[6],y[7], y[8],y[9], y[10],y[11], a[0][0],a[0][1], a[1][0],a[1][1], a[2][0],a[2][1]])

def hf(y):
    d = 0
    N = 3
    for i in range(N): 
        for j in range(i+1,N):
            r2 = (y[2*i]-y[2*j])**2 + (y[2*i+1]-y[2*j+1])**2
            d += 1/r2
    return d
    
def R(Y):
    # Calculates the radius between the bodies and gives back a list
    # containing r01, r02, r12
    r = np.ndarray((3,len(Y.t)))
    
    r[0] = ((Y.y[0]-Y.y[2])**2 + (Y.y[1]-Y.y[3])**2)**0.5  #r01
    r[1] = ((Y.y[0]-Y.y[4])**2 + (Y.y[1]-Y.y[5])**2)**0.5  #r02
    r[2] = ((Y.y[2]-Y.y[4])**2 + (Y.y[3]-Y.y[5])**2)**0.5  #r12
    
    return r

def L(Y):   #Calculate the angular momentum
    L = m[0]*(Y.y[0]*Y.y[7]-Y.y[1]*Y.y[6])+m[1]*(Y.y[2]*Y.y[9]-Y.y[3]*Y.y[8])+m[2]*(Y.y[4]*Y.y[11]-Y.y[5]*Y.y[10])
    return L
    
def E(Y):   #Calculate the total energy
    r = R(Y)
    E = 0.5*m[0]*(Y.y[6]**2+Y.y[7]**2)+0.5*m[1]*(Y.y[8]**2+Y.y[9]**2)+0.5*m[2]*(Y.y[10]**2+Y.y[11]**2)
    E = -m[0]*m[1]/r[0] -m[0]*m[2]/r[1] -m[1]*m[2]/r[2]
    return E

    
t0, t = 0., 80
h = 1e-3

#Y = RK(f, y0, t0, t, h)
#Y = RK_Adapt(f, y0, t0, t, h, hf)
Y = integrate.solve_ivp(f, [t0,t], y0, method='RK45', rtol=1e-13, at max_step=1e-3)

E = E(Y)
L = L(Y)


fig,ax = plt.subplots()
ax.set(xlabel='x', ylabel='y', aspect = 'equal')
ax.plot(Y.y[0],Y.y[1], label='first body')
ax.plot(Y.y[2],Y.y[3], label='second body')
ax.plot(Y.y[4],Y.y[5], label='third body')
ax.legend()

#%%

#Angular momentum and energy plot
fig, (ax1,ax3) = plt.subplots(1,2,sharex=True)

ax1.set(xlabel = 't')
ax3.set(xlabel = 't')
ax1.plot(Y.t,E,'r-',label='energy')
ax3.plot(Y.t,L,'b-',label='angular momentum')
ax1.legend()
ax3.legend()
plt.show()