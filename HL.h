#ifndef __HL_h__
#define __HL_h__

#include<cmath>
#include "nr.h"
#include "constants.h"

struct HL {// fully convective star on hayashi line
    static const char *init_file;
    static double xf;// min radius for mixing-length / stellar radius
    static double alpha;// mixing-length / scale height
    int N;// number of data points
    double mass;// gram
    double radius;// cm
    double X,Y;// mass fraction of H,He
    enum { r_idx, m_idx, P_idx, T_idx, A_idx };
    Mat_DP data;// (r,m,P,T,A) in cgs (shape(5,N))
    HL(double M, // mass / gram
       double R, // radius / cm
       double X, // hydrogen mass fraction
       double Y);// helium mass fraction
    HL(const char* =0);
    void diff_eq(const Vec_DP&, Vec_DP&);
    void center(const Vec_DP&, Vec_DP&);
    void surface(const Vec_DP&, Vec_DP&);
    void build(double M, double R, double X, double Y, double slowc=1);
    void save(const char*);
    void load(const char*);
    void draw(const char*, double, double, int);
    inline double& central_pressure() { return data[P_idx][0]; }
    inline double& surface_pressure() { return data[P_idx][N-1]; }
    inline double& central_temperature() { return data[T_idx][0]; }
    inline double& surface_temperature() { return data[T_idx][N-1]; }
    inline double luminosity() { return 4*PI*Stef*pow(radius
                                    *pow(surface_temperature(),2),2); }
};

void solvde(Mat_DP &y, double a, double b, const Vec_DP& scalv,
            void diff_eq(double, const Vec_DP&, Vec_DP&),
            void boundary1(const Vec_DP&, Vec_DP&),
            void boundary2(const Vec_DP&, Vec_DP&),
            int nb, double slowc=1,
            double eps=1.e-8, int maxiter=10000);
    
void eos(double& rho, double& del_ad, double& delta, double& kappa,
         double P, double T, double X, double Y);
void mixlen(double& del, double del_ad, double del_rad, double delta,
            double P, double T, double rho, double kappa, double g);

#endif // __HL_h__