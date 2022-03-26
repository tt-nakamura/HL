#include<cmath>
#include "HL_shoot.h"
#include "constants.h"

double HL_shoot::xf = 0.2;

static HL_shoot *hs;
static void hs_shoot(const Vec_DP& x, Vec_DP& y) { hs->shoot(x,y); }
static void hs_diff_eq(double x, const Vec_DP& y, Vec_DP& f) { hs->diff_eq(x,y,f); }
static double hs_P, hs_T, hs_L;
static double scale[3];

HL_shoot::HL_shoot(double M, // mass / gram
       double R, // radius / cm                       
       double X1,// hydrogen mass fraction
       double Y1,// helium mass fraction
       int N1     // number of data points
    ) : N(N1), mass(M), radius(R), X(X1), Y(Y1),
        r(N), data(3,N), param(3)
{
    double mm(M/Msun), rr(R/Rsun), dr;
    int j;
    Vec_DP y(3);

    param[P_shoot] = Psun*pow(mm,2)/pow(rr,4);
    param[T_shoot] = Tsun*mm/rr;
    param[L_shoot] = Lsun*mm*pow(rr, 2.5);
    scale[m_difeq] = mass;
    scale[P_difeq] = central_pressure();
    scale[T_difeq] = central_temperature();
    rf = radius*xf;
    hs = this;
    if(newt(param, hs_shoot))
        nrerror("shoot failed: bad initial guess");

    dr = radius/(N-1);
    for(j=0; j<N; j++) r[j] = j*dr;
    center(y, param);// boundary condition at center
    odeint(y, hs_diff_eq, r, data);// integrate outward
}

void HL_shoot::diff_eq(double r, const Vec_DP& y, Vec_DP& f)
// r = distance from center / cm
// y = dependent variables [m,P,T]
// output: f = dy/dr
{
    const double& m = y[m_difeq];// mass within r / g
    const double& P = y[P_difeq];// pressure at r / dyne/cm^2
    const double& T = y[T_difeq];// temperature at r / K
    double& dm_dr = f[m_difeq];
    double& dP_dr = f[P_difeq];
    double& dT_dr = f[T_difeq];
    double r2(r*r), rho, kappa;
    double delta, del, del_rad, del_ad, g;
    static double PSG(3./64./PI/Stef/Grav);
    eos(rho, del_ad, delta, kappa, P,T,X,Y);
    g = (r==0 ? 0 : Grav*m/r2);
    if(r <= rf) del = del_ad;
    else {// use mixing-length theory
        del_rad = hs_L/m*PSG*kappa*P/pow(T,4);
        mixlen(del, del_ad, del_rad, delta, P, T, rho, kappa, g);
    }
    dm_dr = 4*PI*rho*r2;
    dP_dr = -g*rho;
    dT_dr = dP_dr*T/P*del;
}

void HL_shoot::shoot(const Vec_DP& x, Vec_DP& y) {
    int i;
    Vec_DP w(3);
    center(y,x);// boundary condition at center
    odeint(y, hs_diff_eq, 0, rf);// integrate outward
    surface(w,x);// boundary condition at surface
    odeint(w, hs_diff_eq, radius, rf);// integrate inward
    for(i=0; i<3; i++) y[i] -= w[i];
    for(i=0; i<3; i++) y[i] /= scale[i];
    for(i=0; i<3; i++) std::cout << y[i] << ' ';
    std::cout << std::endl;
}

void HL_shoot::center(Vec_DP& y, const Vec_DP& x) {
    y[m_difeq] = 0;         // mass
    y[P_difeq] = x[P_shoot];// pressure
    y[T_difeq] = x[T_shoot];// temperature
}

static double hs_surface_eq(double x) { return hs->surface_eq(x); }

void HL_shoot::surface(Vec_DP& y, const Vec_DP& x) {
    double R2(radius*radius), x1, x2, y1, y2, err(1.e-12);
    hs_L = x[L_shoot];
    hs_T = pow(hs_L/(4*PI*R2*Stef), 0.25);
    hs_P = 2./3.*Grav*mass/R2;
    for(y1 = surface_eq(x1=2);
        y1 * (y2 = surface_eq(x2=x1*2)) > 0;
        x1=x2, y1=y2);
    y[P_difeq] = zbrent(hs_surface_eq, x1, x2, err);
    y[T_difeq] = hs_T;
    y[m_difeq] = mass;
}

double HL_shoot::surface_eq(double P) {
    double dum,kappa;
    eos(dum,dum,dum,kappa, P, hs_T, X,Y);
    return P - hs_P/kappa;
}