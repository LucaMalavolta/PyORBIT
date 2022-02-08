#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import numpy as np # array


# radiants, degrees conversions etc.
pi = 4.*np.arctan(1.)
dpi = 2.*pi
deg2rad = pi/180.
rad2deg = 180./pi

# various
TOLERANCE = np.finfo(np.float64(1.0)).eps
d2s = 86400. # seconds in a day  =  24h  =  86400 s
d2m = 1440.  # min in a day  =  1440. min

# masses conversions
Msmer = 6.0236e6 # Msun to Mmer
Mmers = 1./Msmer # Mmer to Msun
Msven = 4.08523719e5 # Msun to Mven
Mvens = 1./Msven   #  Mven to Msun
Msear = 332946.0487 # Msun to Mear
Mears = 1./Msear  #  Mear to Msun
Msmar = 3.09870359e6 # Msun to Mmar
Mmars = 1./Msmar   #  Mmar to Msun
Msjup = 1.047348644e3 # Msun to Mjup
Mjups = 1./Msjup    #  Mjup to Msun
Mssat = 3.4979018e3 # Msun to Msat
Msats = 1./Mssat  #  Msat to Msun
Msura = 2.290298e4 # Msun to Mura
Muras = 1./Msura #  Mura to Msun
Msnep = 1.941226e4 # Msun to Mnep
Mneps = 1./Msnep #  Mnep to Msun

Mejup = Mears * Msjup # Mear to Mjup
Mjear = Mjups * Msear # Mjup to Mear

# masses of Solar System objects
Msun = 1.9884e30 # Sun mass in kg
Mmer = Msun*Mmers #  Mercury mass in kg
Mven = Msun*Mvens #  Venus mass in kg
Mear = 5.9722e24 # Earth mass in kg
Mmar = Msun*Mmars #  Mars mass in kg
Mjup = Msun*Mjups #  Jupiter mass in kg
Msat = Msun*Msats #  Saturn mass in kg
Mura = Msun*Muras #  Uranus mass in kg
Mnep = Msun*Mneps #  Neptune mass in kg

# radii of Solar System objects
Rsun = 696000. #  Sun radius in km
Rmer = 2439.7  #  Mercury radius in km
Rven = 6051.8  #  Venus radius in km
Rear = 6378.1366 # Earth radius in km
Rmar = 3396.19 #  Mars radius in km
Rjup = 71492.  #  Jupiter radius in km
Rsat = 60268.  #  Saturn radius in km
Rura = 25559.  #  Uranus radius in km
Rnep = 24764.  #  Neptune radius in km
Rplu = 1195.   #  Pluto radius in km
#
Rsjup = Rsun/Rjup # Rsun to Rjup
Rjups = Rjup/Rsun # Rjup to Rsun
Rsear = Rsun/Rear # Rsun to Rear
Rears = Rear/Rsun # Rear to Rsun
Rsnep = Rsun/Rnep # Rsun to Rnep
Rneps = Rnep/Rsun # Rnep to Rsun
#
Rejup = Rear/Rjup # Rearth to Rjupiter
Rjear = Rjup/Rear # Rjupiter to Rearth


# astronomical constants
AU = 149597870700.  # Astronomical Unit in meters
kappa = 0.01720209895  # Gaussian gravitational constant
Giau = kappa*kappa  # G [AU^3/Msun/d^2]
Gsi = 6.67428e-11  # Gravitational Constant in SI system [m^3/kg/s^2]
Gaumjd = Gsi*d2s*d2s*Mjup/(AU**3)  # G in [AU,Mjup,day]
speed = 299792458.  # speed of light (c) in [m/s]
speedaud = speed*d2s/AU  # speed of light in [AU/d]
pc2AU = 206264.806


# others
RsunAU = (Rsun*1.e3)/AU #Sun radius in AU
RjupAU = (Rjup*1.e3)/AU #Jupiter radius in AU

MJD = 2400000.5 # MJD ref time to convert to JD

rho_Sun = Msun / (4./3.*np.pi* (Rsun*1000.)**3 )

# Used in the computation of the stellar insolation
Sun_constant = 1367.
Sun_temperature = 5777.

sigma2FWHM = 2. * np.sqrt(2.* np.log(2,))