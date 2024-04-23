from utilities import *
from orbitalElements import *
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')


def f(a,T): #functions to calculate the perturbations from the coefficients
  return a[0]+a[1]*T + a[2]*T**2 + a[3]*T**3

eps = np.array([23.439291,-0.0130042,-0.00059,0.001813]) #Earth axis inclination over the Ecliptic plane
long = 11.87167  #Latitude and longitude of Padua
lat = 45.400
R = 6378.137 #Earth radius in km
time_zone = +1
'''Conversions '''

eps = np.radians(eps) #eps is converted in radians
lat = np.radians(lat)
long = np.radians(long)

for p in [sun,mercury,venus,mars,earth]: #Conversions in kms and radians
  p.a = p.a*149597870.7 #AU in km
  p.i = np.radians(p.i)
  p.L = np.radians(p.L)
  p.Omega = np.radians(p.Omega)
  p.pi = np.radians(p.pi)
  p.RA = [] #Adding the Right Ascension and declination lists
  p.dec = []
  p.A = [] #Adding the Azimuth and elevation lists
  p.h = []

'''Cartesian coordinates in the geocentric equatorial RF'''

def equatorialCoord(planet,Y,M,D,H):
  T = getT(Y,M,D,H)
  #Calculating the coordinates in the heliocentric system
  planet.X = X(f(planet.a,T),f(planet.e,T),f(planet.i,T),f(planet.Omega,T),f(planet.pi,T),f(planet.L,T))
  earth.X = X(f(earth.a,T),f(earth.e,T),f(earth.i,T),f(earth.Omega,T),f(earth.pi,T),f(earth.L,T))

  planet.X = planet.X - earth.X # Going in the geocentric system
  planet.X = np.dot(Rx(f(eps,T)),planet.X)  #Rotation to the equatorial RF

  #print(planet.name, ' x:',planet.X[0], ' y:',planet.X[1], ' z: ',planet.X[2], '  date  ',Y,M,D,'hour',H, '\n')
  return planet.X



'''Right Ascensions and declination'''

def RAdec(planet,Y,M,D,H):

  X = equatorialCoord(planet, Y, M, D, H)
  RA = np.arctan2(X[1],X[0])  #RA is in radians and from -pi to pi
  RA = RA if RA>0 else RA+2*np.pi   #RA is mapped to 0, 2 pi
  dec = np.arctan2(X[2],(X[0]**2+X[1]**2)**0.5)
  #print(planet.name, ' RA:',np.degrees(RA), ' dec:',np.degrees(dec),'  date  ',Y,M,D,'hour',H, '\n')
  return np.degrees(RA), np.degrees(dec) #Angles are converted in degrees to make reading them easier



''' Azimuth and elevation'''
def Ah(planet,Y,M,D,H, lat, long):
  l = GMST(Y,M,D,H) + long #Mean sideral time, given by long+GMST
  X = equatorialCoord(planet, Y, M, D, H)

  X = np.dot(Rz(-np.pi/2-l),X) #Rotation around the z-axis by the longitude plus the GMST0
  X = np.dot(Rx(-np.pi/2+lat),X) #Rotation around x of the colatitude

  #In this new reference frame, y is pointing north and x is pointing east
  #The calculation is done using the not negligible size of earth

  X = X - np.array([0.0,0.0,0.0]) #(0,0,R) is the vector to the surface in this RF

  A = np.arctan2(X[0],X[1]) #The Azimuth is given by atan(east/north)
  A = A if A>0 else A+2*np.pi   #It's mapped from -pi, pi
  h =  np.arctan2(X[2],(X[0]**2+X[1]**2)**0.5)
  #print(planet.name, ' A:',np.degrees(A), ' h:',np.degrees(h),'  date  ',Y,M,D,'hour',H, '\n')
  return np.degrees(A), np.degrees(h)

'''Ephemeridis plot'''

for d in range(1,366): #for all days from 1 to 365
  for p in [sun, mercury, venus, mars]:
    RA,dec = RAdec(p, 2024, 1, d, 23-time_zone) #Padua is in the UTC + 1 time zone
    A,h = Ah(p, 2024, 1, d, 23-time_zone, lat, long)
    p.RA.append(RA)
    p.dec.append(dec)
    p.A.append(A)
    p.h.append(h)
  if(venus.h[-1]>0):
    print(d)


fig = plt.figure(label = 'Ephemeris for all planets at 11 pm')
ax1 = fig.add_subplot(1,2,1)
ax1.grid()
ax1.set(xlim = [0,360], xlabel='RA', ylim= [-90,90], ylabel='dec',title = 'Right Ascension and declination')
ax1.plot(mercury.RA, mercury.dec, '.', label=mercury.name, markersize=2)
ax1.plot(venus.RA, venus.dec, '.', label=venus.name, markersize=2)
ax1.plot(mars.RA, mars.dec, '.',  label=mars.name, markersize=2)
ax1.plot(sun.RA, sun.dec, '.', label=sun.name, markersize=2)
ax1.legend()

ax2 = fig.add_subplot(1,2,2)
ax2.grid()
ax2.set(xlim = [0,360], xlabel='A', ylim= [-90,90], ylabel='h',title = 'Azimuth and elevation')
ax2.plot(mercury.A, mercury.h, '.', label=mercury.name, markersize=2)
ax2.plot(venus.A, venus.h, '.', label=venus.name, markersize=2)
ax2.plot(mars.A, mars.h, '.',  label=mars.name, markersize=2)
ax2.plot(sun.A, sun.h, '.', label=sun.name, markersize=2)
ax2.legend()
plt.show()