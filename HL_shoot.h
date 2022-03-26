#ifndef __HL_shoot_h__
#define __HL_shoot_h__

#include "HL.h"

struct HL_shoot {// fully convective star on Hayashi Line
    static double xf;// fitting point / stellar radius
    int N;// number of data points
    double mass;// gram
    double radius;// cm
    double X,Y;// mass fraction of H,He
    double rf;// fitting point (radius*xf)
    Vec_DP r;// radial coordinate / cm (shape(N,))
    Mat_DP data;// (m,P,T) in cgs (shape(3,N))
    Vec_DP param;// R, P_c, T_c / cgs
    enum { P_shoot, T_shoot, L_shoot };
    enum { m_difeq, P_difeq, T_difeq };
    HL_shoot(double M, // mass / gram
             double R, // radius / cm
             double X, // hydrogen mass fraction
             double Y, // helium mass fraction
             int n);   // number of data points
    void diff_eq(double, const Vec_DP&, Vec_DP&);
    void shoot(const Vec_DP&, Vec_DP&);
    void center(Vec_DP&, const Vec_DP&);
    void surface(Vec_DP&, const Vec_DP&);
    double surface_eq(double);
    void save(const char*);
    inline double& central_pressure() { return param[P_shoot]; }
    inline double& central_temperature() { return param[T_shoot]; }
    inline double& luminosity() { return param[L_shoot]; }
};

bool newt(Vec_DP&, void(const Vec_DP&, Vec_DP&));
void odeint(Vec_DP&, void(double, const Vec_DP&, Vec_DP&),
            double, double, double=1e-9);
void odeint(const Vec_DP&, void(double, const Vec_DP&, Vec_DP&),
            const Vec_DP&, Mat_DP&, double=1e-9);
double zbrent(double(double), double, double, double);

#endif // __HL_shoot_h__