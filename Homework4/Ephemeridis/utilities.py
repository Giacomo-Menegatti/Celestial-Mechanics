#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 16:04:21 2023

@author: menegattig
"""

import matplotlib.pyplot as plt
import numpy as np

'''Kepler equation solver using Newton-Raphson'''
def KeplerEq(M,e,rtol=1e-12):
  k = np.floor(M/(2*np.pi))
  M = M-k*2*np.pi           #Subtract all the complete orbits
  if M == 0:  #This case has a known result and it's implemented here to avoid errors later in the divisions
    return 0;
  E = M   #First guess for E
  err = 1 #First guess for the error
  while(err > rtol):
    newE = E - (E-e*np.sin(E)-M)/(1-e*np.cos(E))
    err = np.fabs(newE/E-1)
    E = newE
  return E

'''Rotation Matrix'''
def Rz(theta):
  return np.array([[np.cos(theta),-np.sin(theta),0.0],[np.sin(theta), np.cos(theta), 0.0],[0.0, 0.0, 1.0]])

def Rx(theta):
  return np.array([[1.0, 0.0, 0.0], [0.0, np.cos(theta),-np.sin(theta)],[0.0, np.sin(theta), np.cos(theta)]])

''' Coordinates and velocity conversion from Laplece RF to Ecliptic RF'''

def X(a,e,i,node_lon,peri_lon,mean_lon):
  omega = peri_lon-node_lon #argument of pericenter
  M = mean_lon-peri_lon #mean anomaly
  E = KeplerEq(M, e)
  X = np.array([a*(np.cos(E)-e),a*(1-e**2)**0.5*np.sin(E), 0.0])  #Position in the Laplace RF
  X = np.dot(Rz(omega),X)   #Rotation around the z-axis of the argument of pericenter
  X = np.dot(Rx(i),X)   #Inclination of the orbital plane
  X = np.dot(Rz(node_lon),X)  #Rotation to the nodal line longitude
  return X

def V(n,a,e,i,node_lon,peri_lon,mean_lon):
  omega = peri_lon-node_lon #argument of pericenter
  M = mean_lon-peri_lon #mean anomaly
  E = KeplerEq(M, e)
  r = a*(1-e*np.cos(E)) #radius of the orbit
  V = np.array([-n*a**2*np.sin(E)/r,n*a**2*(1-e**2)**0.5*np.cos(E)/r, 0.0])  #Velocity in the Laplace RF
  V = np.dot(Rz(omega),V)   #Rotation around the z-axis of th eargument of pericenter
  V = np.dot(Rx(i),V)   #Inclination of the orbital plane
  V = np.dot(Rz(node_lon),V)  #Rotation to the nodal line longitude
  return V

''' Calculating JulianDate, T and LMST '''

def JulianDate(Y, M, D, H): #Julian date for Y>1582
  C = np.round((M - 14)/12)
  JD0 = D - 32075 + np.round(1461 * (Y + 4800 + C)/4)+ np.round(367 * (M - 2 - C * 12)/12) - np.round(3 * (np.round(Y + 4900 + C) / 100) /4)
  JD = JD0 + H/24 - 0.5 #It's calculated at midnight while JD0 is calculated at noon, so 0.5 is subtracted
  return JD

def getT(Y, M, D, H):
  return (JulianDate(Y, M, D, H)-2451545.0)/36525

def GMST(Y, M, D, H):  #Greenwich sideral time
  T = getT(Y, M, D, 0.0) #Calculated at midnight, so H=0
  GMST0 = 24110.54841 + 8640184.812866*T + 0.093104*T**2 - 6.2e-6*T**3  #GMST0 is in seconds, so it will be converted in days and then in radians
  GMST = GMST0/(24*3600) + 1.002737909350795 * H/24 #GMST at any hour, in days
  GMST = 2*np.pi*GMST     #GMST is converted in radians, as 24h = 1 d = 360 deg = 2 pi
  return GMST

