#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 23:06:43 2023

@author: menegattig
"""

from RungeKutta import RK, RK_Adapt
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

plt.close('all')

#state vector: x1,y1,vx1,vy1,   x2,y2,vx2,vy2,

x0 = [0.0,0.0, 0.0,0.0,  384400.0,0.0, 0.91647306922544,0.91647306922544]
GM = 398600.4415

def f(t, x):
    # r=x2-x1, points toward the second mass    
    r = np.array([x[4]-x[0],x[5]-x[1]])
    
    r3 = (r[0]**2 + r[1]**2)**1.5
    a = GM/r3*r
    # f(y) = (vx1,vy1, ax1,ay1,  vx2,vy2, ax2,ay2)
    return np.array([x[2],x[3], a[0],a[1],   x[6],x[7], -a[0],-a[1]])

def hf(x):  #Variable time step calculator: f=1/r**2
    r = np.array([x[4]-x[0],x[5]-x[1]])
    r2 = (r[0]**2 + r[1]**2)    
    return 1/r2

def E(Y):   #Calculate the energy
    r = R(Y)
    E = np.ndarray(len(Y.y))
    E = 0.5*(Y.y[2]**2+Y.y[3]**2+Y.y[6]**2+Y.y[7]**2)-GM/r
    return E
        
def R(Y):   #Calculate the radius 
    r = ((Y.y[4]-Y.y[0])**2+(Y.y[5]-Y.y[1])**2)**0.5
    return r
    
    
def L(Y):   #Calculate the specific moment
    L = np.array(len(Y.y))
    L = (Y.y[0]*Y.y[3]-Y.y[1]*Y.y[2])+(Y.y[4]*Y.y[7]-Y.y[5]*Y.y[6])
    return L

def CM(Y):  #Calculate the CM position
    CM = np.zeros((2, len(Y.t)))
    CM[0] = 0.5*(Y.y[0]+Y.y[4]) #x
    CM[1] = 0.5*(Y.y[1]+Y.y[5]) #y
    return CM

t0, t = 0., 3.0e6
h = 1e4

Y = RK(f, x0, t0, t, h)
Y2 = RK_Adapt(f, x0, t0, t, h, hf)
Y1 = integrate.solve_ivp(f, (t0,t), x0, method='RK45', rtol=1e-12)

CM = CM(Y)
E1 = E(Y)
L1 = L(Y)

E2 = E(Y1)
L2 = L(Y1)

E3 = E(Y2)
L3 = L(Y2)

plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.plot(Y.y[0],Y.y[1],label='first body')
plt.plot(Y.y[4],Y.y[5],label='second body')
plt.plot(CM[0],CM[1],label='center of mass')
plt.legend()


#%%
#Angular momentum and energy plot
fig, (ax1,ax3) = plt.subplots(1,2,sharex=True)


ax1.plot(Y.t,E1,'r-',label='energy')
ax3.plot(Y.t,L1,'b-',label='angular momentum')
ax1.legend()
ax3.legend()
plt.show()

#%%
# trajectory plot in the CM reference system

fig,ax = plt.subplots()
ax.set(xlabel='x [km]', ylabel='y [km]', aspect = 'equal')
plt.plot(Y.y[0]-CM[0],Y.y[1]-CM[1],label='first body')
plt.plot(Y.y[4]-CM[0],Y.y[5]-CM[1],label='second body')
plt.legend()

#%%
# trajectory comparison between different methods
fig,ax = plt.subplots()
ax.set(xlabel='x [km]', ylabel='y [km]', aspect = 'equal')
plt.plot(Y.y[0],Y.y[1],label='Runge Kutta')
plt.plot(Y1.y[0],Y1.y[1],label='RK45')
plt.legend()







