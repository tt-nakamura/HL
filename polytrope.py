# reference: R. Kippenhahn and A. Weigert
#   "Stellar Structure and Evolution" chapter 19

import numpy as np
from scipy.integrate import solve_bvp
from constants import Grav

class polytrope:
    def __init__(self, n, N=128, init=None, **kw):
        """ polytropic stellar model
        n = polytropic index (P = K rho^{1+1/n})
        N = number of mesh points
        init = polytrope object for initial guess
        kw = keyword arguments passed to solve_bvp
        """
        if n<0 or n>=5:
            print('bad index in polytrope')
            return

        def diff_eq(x,y,p):# Lane-Emden equation
            zn2 = p[0]**2
            d2w = -2*y[1,1:]/x[1:] - y[0,1:]**n*zn2
            return [y[1], np.r_[-zn2/3, d2w]]

        def bc(ya,yb,p):
            return [ya[0]-1, ya[1], yb[0]]

        if not isinstance(init, polytrope):
            p = [2.45]# outer radius for n=0
            x = np.linspace(0,1,N)
            z = x*p[0]
            y = [1 - z**2/6, -z/3]# exact solution for n=0
            for n in np.linspace(0, n, int(2*n)+1):
                s = solve_bvp(diff_eq, bc, x,y,p, **kw)
                x,y,p = s.x, s.y, s.p
        else:
            x = init.z/init.zn
            y = [init.w, init.wp*init.zn]
            p = [init.zn]
            s = solve_bvp(diff_eq, bc, x,y,p, **kw)

        self.n = n
        self.zn = s.p[0] # outer radius (w(zn)==0)
        self.z = s.x * s.p[0] # independent variable
        self.w = s.y[0] # dependent variable
        self.wp = s.y[1]/s.p[0] # dw/dz

    def central_state(self,M,R):# given mass and radius
        rho_c = M/4/np.pi/R**3*self.zn/(-self.wp[-1])
        P_c = 4*np.pi*Grav*(R*rho_c/self.zn)**2/(self.n + 1)
        return rho_c, P_c # density and pressure at r=0


class PolytropicStar(polytrope):
    def __init__(self, n, M, R, N=128):
        """
        n = polytropic index (P = K rho^{1+1/n})
        M,R = mass and radius of star
        N = number of mesh points
        """
        super().__init__(n,N)
        rho_c, P_c = self.central_state(M,R)
        self.r = self.z/self.zn*R # radial coordinate
        self.rho = rho_c*self.w**n # density profile
        self.P = P_c*self.w**(n+1) # pressure profile
        self.m = np.r_[0, 4*np.pi*rho_c*self.r[1:]**3
                       *(-self.wp[1:])/self.z[1:]] # mass profile
