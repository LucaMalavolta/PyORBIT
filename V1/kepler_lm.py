import numpy as np

G_grav=6.67398e-11
M_sun =1.98892e30

#OLD definition now avoided
#def kepler_E(M,ec):
#  # brute-force solution of kepler equation
#  E0 = np.zeros(np.size(M),dtype=np.double)
#  E1 = np.asarray(M,dtype=np.double)
#  ecc     = np.asarray(ec,dtype=np.double)
#  eccanom = np.zeros(np.size(M),dtype=np.double)
#  
#  while np.amax(np.abs((E0-E1)/E1)) > 1e-7:
#    E0[:] = E1[:]
#    E1 = E0-(E0-ecc*np.sin(M)-M)/(1.-ecc*np.cos(E0))
#
#  return eccanom

def kepler_E(M_in,e_in):
  # calculation of the eccentric anomaly
  # using the Newton method
  # M = mean anomaly
  # e = eccentricity
  
  # Converting to numpy arrays
  M = np.asarray(M_in,dtype=np.double)
  e = np.asarray(e_in,dtype=np.double)
  
  # Fn = initial value
  EccAn = M + e*np.sin(M) + e**2/2*np.sin(2*M) # eccentric anomaly
  EccAnIt  = EccAn*2.0
  while np.amax(np.abs((EccAn-EccAnIt)/EccAnIt)) > 1e-7:
    EccAnIt[:]  = EccAn[:]
    Mn = EccAnIt-e*np.sin(EccAnIt)
    EccAn = EccAnIt+(M-Mn)/(1.-e*np.cos(EccAnIt))
  
  return EccAn

def kepler_K1(M_star1,M_star2,Period,i,e0):
  # M_star1, M_star2 in solar masses
  # P in days -> Period is converted in seconds in the routine
  # i in degrees
  # Gravitational constant is given in m^3 kg^-1 s^-2 
  # output in m/s
  K1 = (2.*np.pi*G_grav*M_sun/86400.0)**(1.0/3.0) * (np.sin(i*np.pi/180.0)/ np.sqrt(1.0-e0**2.0)) * (Period) ** (-1.0/3.0) * ( M_star2*(M_star1+M_star2)**(-2.0/3.0))
  return K1

def kepler_RV(BJD, TPeri, Period, gamma, K, e0, omega0):
  # omega = argument of pericenter
  #Mean Anomaly
  #
  MeAn = 2.*np.pi*(1. + (BJD - TPeri)/Period % 1.)

  if abs(e0)<1e-3:
     TrAn  = np.asarray(MeAn,dtype=np.double)
     e     = np.asarray(0.,dtype=np.double)
     omega = np.asarray(0.,dtype=np.double)
  else:
    if e0 < 0.:
      e     = np.asarray(-e0,dtype=np.double)
      omega = np.asarray(omega0,dtype=np.double) + np.pi
    else:
      e     = np.asarray(e0,dtype=np.double)
      omega = np.asarray(omega0,dtype=np.double)
  
    #Eccentric Anomaly
    EccAn = kepler_E(MeAn,e)
    TrAn  = 2.* np.arctan(np.sqrt((1.0+e)/(1.0-e))*np.tan(EccAn/2.0))

  rv = K*(np.cos(TrAn+omega) + e*np.cos(omega)) + gamma
  
  return rv
  
def kepler_RV_T0P(BJD0, phase, Period, K, e0, omega0):
  #BJD0 is given as BJD-T0, where T0 is arbitrarily defined by the user
  #_Tperi_ is subsituted by _phase_, which is the phase of the orbit where
  #        BJD0+T0+phase*Period = Tperi  
  # omega = argument of pericenter
  #Mean Anomaly
  #
  #MeAn = 2.*np.pi*(1. + (BJD - TPeri)/Period % 1.)
  #MeAn = 2.*np.pi*(1. + (BJD0 - phase*Period)/Period % 1.)
  MeAn = 2.*np.pi*(1. + ((BJD0/Period)-phase) % 1.)
  
  if abs(e0)<1e-3:  
     TrAn  = np.asarray(MeAn,dtype=np.double)
     e     = np.asarray(0.,dtype=np.double)
     omega = np.asarray(0.,dtype=np.double)
  else:
    if e0 < 0.:
      e     = np.asarray(-e0,dtype=np.double)
      omega = np.asarray(omega0,dtype=np.double) + np.pi
    else:
      e     = np.asarray(e0,dtype=np.double)
      omega = np.asarray(omega0,dtype=np.double)
  
    #Eccentric Anomaly
    EccAn = kepler_E(MeAn,e)
    TrAn  = 2.* np.arctan(np.sqrt((1.0+e)/(1.0-e))*np.tan(EccAn/2.0))
  
  rv = K*(np.cos(TrAn+omega) + e*np.cos(omega))

  return rv
