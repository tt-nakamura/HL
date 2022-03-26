import numpy as np
import matplotlib.pyplot as plt
from constants import Grav,Msun,stefan,mH,kB
from polytrope import polytrope
from main_seq import main_seq

M = Msun
X,Y = 0.71,0.27
a,b = 1,3
kappa0 = 0.229/8e4**a/6e3**b
a1 = a+1
mu = 2/(1 + 3*X + Y/2)

R = np.geomspace(1e7, 1e13, 100)

for n in [1.5, 1.7, 1.9, 2.1]:
    p = polytrope(n)
    rho_c, P_c = p.central_state(M,R)
    K = P_c/rho_c**(1+1/n)
    an = a1*(n+1)
    T = (2*Grav*M*K**(a1*n)
         /(3*kappa0*R**2*(kB/mH/mu)**an))**(1/(an+b))

    L = 4*np.pi*R**2*stefan*T**4
    plt.loglog(T, L, label='n='+str(n))

M = np.geomspace(0.1*Msun, 20*Msun, 20)
_,L,T = main_seq(M)
plt.loglog(T, L, 'k--', label='main sequence') 

plt.axis([3e4, 3e2, 1e25, 1e34])
plt.xlabel('T = surface temperature / K', fontsize=14)
plt.ylabel('L = luminosity / erg/s', fontsize=14)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('fig2.eps')
plt.show()
