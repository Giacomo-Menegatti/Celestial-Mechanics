#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 23:02:15 2023
@author: menegattig
"""
import numpy as np


def RK(f,y0,t0,t,h):
    #Runge-Kutta 4th order integrator, taking as parameters the function f, the initial state y0
    # the initial time t0 end the final time t, and the step h
    
    #Struct to handle the data
    class state_vector:
        pass
    
    Y = state_vector()
    
    # Number of steps and variables in the state vector
    N = round((t-t0)/h)
    n = len(y0)
    
    # Array of times and array of state vectors
    Y.t = np.linspace(t0, t, N)
    Y.y = np.zeros((N,n))
    
    #Initialisation of the state vector array
    Y.y[0] = y0
    
    # Integrating form t0 to t
    for i in range(N-1):
        
        y = Y.y[i]
        t = Y.t[i]
        
        k1 = h*f(t,y)
        k2 = h*f(t + 0.5*h, y + 0.5*k1)
        k3 = h*f(t + 0.5*h, y + 0.5*k2)
        k4 = h*f(t + h, y + k3)
        
        Y.y[i+1] = y + (k1 + 2*k2 + 2*k3 + k4)/6.0
    
    #The state vector array is transposed to make it equal to the one given by the scipy integrators
    Y.y = np.transpose(Y.y)    
    return Y





def RK_Adapt(f,y0,t0,t,h0,hf):
    
    #Runge Kutta method with adaptive steps. Here h0 is the initial step value 
    #ad hf a function to estimate the step size from the state vector
    
    class state_vector:
        pass
    
    Y = state_vector()
    
    N = round((t-t0)/h0)
    n = len(y0)
    
    Y.t = np.linspace(t0, t, N)
    Y.y = np.zeros((N,n))
    Y.y[0] = y0
    
    #The step size is calculated as h = k/hf(y). To estimate the value of k
    # the law is reversed and calculated for the first step
    
    k = h0*hf(y0)   
    for i in range(N-1):
        
        y = Y.y[i]
        t = Y.t[i]
        
        # k/hf is the adaptive step size
        # if k/hf is smaller than the initially set step
        # the integration step is divided in smaller parts
        
        h_adapt = k/hf(y)   #Value suggested by the function
              
        n = round(h0/h_adapt)+1 #Divisions of the step
        h = h0/n                #Step used in this integration
        
        # Integration over the subdivisions
        for j in range(n):
            k1 = h*f(t,y)
            k2 = h*f(t + 0.5*h, y + 0.5*k1)
            k3 = h*f(t + 0.5*h, y + 0.5*k2)
            k4 = h*f(t + h, y + k3)
            y = y + (k1 + 2*k2 + 2*k3 + k4)/6.0
            
        Y.y[i+1] = y
        
    Y.y = np.transpose(Y.y)    
    return Y