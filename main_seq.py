import numpy as np
from constants import Msun,Lsun,Rsun,stefan

def main_seq(M):
    """ main sequence star (empirical relation)
    M = mass / g
    return R,L,T where
      R = radius / cm
      L = luminosity / erg/s
      T = surface termperature / K
    reference: A. N. Cox, et al
      "Allen's Astrophysical Quantities" 4th edition, p382
    """
    lM = np.log(M/Msun)
    L = Lsun*10**(3.8*lM + 0.08)
    R = Rsun*10**(np.where(lM>0.12, 0.64*lM + 0.011, 0.917*lM - 0.020))
    T = (L/4/np.pi/R**2/stefan)**0.25
    return R,L,T
