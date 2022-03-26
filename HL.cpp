#include "HL.h"

const char* HL::init_file = "HL_M1R1X71Y27.dat";
double HL::xf = 0.2;
double HL::alpha = 1.5;

HL::HL(double M,// mass / gram
       double R,// radius / cm
       double X,// hydrogen mass fraction
       double Y // helium mass fraction
    ) {
    load(init_file);
    build(M,R,X,Y);
}

HL::HL(const char* fname) { load(fname ? fname : init_file); }

static HL *hl;
static void hl_center(const Vec_DP& y, Vec_DP& f) { hl->center(y,f); }
static void hl_surface(const Vec_DP& y, Vec_DP& f) { hl->surface(y,f); }
static void hl_diff_eq(double x, const Vec_DP& y, Vec_DP& f) { hl->diff_eq(y,f); }

void HL::build(double M, double R, double X1, double Y1, double slowc)
// slowc = control parameter of relaxation method (0<slowc<=1)
//   (slowc must be small when initial guess is bad)
{
    static double dM(1.01), dR(1.01);// gradual change in M,R
    Vec_DP scale(1.,5);
    X=X1; Y=Y1;
    hl = this;
    scale[r_idx] = radius;
    do {// adjust mass gradually
        if(mass<M) { if((mass*=dM)>M) mass=M; }
        else       { if((mass/=dM)<M) mass=M; }
        scale[m_idx] = mass;
        scale[P_idx] = central_pressure();
        scale[T_idx] = central_temperature();
        solvde(data, 0, 1, scale, hl_diff_eq, hl_center, hl_surface, 2, slowc);
    } while(mass!=M);
    scale[m_idx] = mass;
    while(radius!=R) {// adjust radius gradually
        if(radius<R) { if((radius*=dR)>R) radius=R; }
        else         { if((radius/=dR)<R) radius=R; }
        scale[r_idx] = radius;
        scale[P_idx] = central_pressure();
        scale[T_idx] = central_temperature();
        solvde(data, 0, 1, scale, hl_diff_eq, hl_center, hl_surface, 2, slowc);
    }
}

void HL::diff_eq(const Vec_DP& y, Vec_DP& f)
// y = dependent variable (r,m,P,T,A)
// ouput: f = dy/dx
{
    const double& r = y[r_idx];// distance from center / cm
    const double& m = y[m_idx];// mass within r / g
    const double& P = y[P_idx];// pressure at r / dyne/cm^2
    const double& T = y[T_idx];// temperature at r / K
    const double& A = y[A_idx];// adaptive mesh parameter
    double& dr_dx = f[r_idx];
    double& dm_dx = f[m_idx];
    double& dP_dx = f[P_idx];
    double& dT_dx = f[T_idx];
    double& dA_dx = f[A_idx];
    double r2(r*r), rho, kappa;
    double delta, del, del_rad, del_ad, g;
    double c1(1/radius);
    double c2(1/log(central_pressure()/surface_pressure()));
    static double PSG = 3./64./PI/Stef/Grav;
    eos(rho, del_ad, delta, kappa, P,T,X,Y);
    g = Grav*m/r2;
    if(r <= radius*xf) del = del_ad;
    else {// use mixing-length theory
        del_rad = luminosity()/m*PSG*kappa*P/pow(T,4);
        mixlen(del, del_ad, del_rad, delta, P, T, rho, kappa, g);
    }
    double dm_dr = 4*PI*rho*r2;
    double dP_dr = -g*rho;
    double dlnP_dr = dP_dr/P;
    double dT_dr = dlnP_dr*T*del;
    dA_dx = 0;
    dr_dx = A/(c1 - c2*dlnP_dr);// adaptive mesh
    dm_dx = dm_dr*dr_dx;
    dP_dx = dP_dr*dr_dx;
    dT_dx = dT_dr*dr_dx;
}

void HL::center(const Vec_DP& y, Vec_DP& f)
// boundary condition at r=0
{
    f[0] = y[r_idx];// r=0 at x=0
    f[1] = y[m_idx];// m=0 at x=0
}

void HL::surface(const Vec_DP& y, Vec_DP& f)
// boundary condition at r=R
{
    const double& P = y[P_idx];// pressure at R / dyne/cm^2
    const double& T = y[T_idx];// temperature at R / K
    double kappa, dum;
    eos(dum,dum,dum,kappa, P,T,X,Y);
    f[0] = y[r_idx] - radius;// r=R at x=1
    f[1] = y[m_idx] - mass;// m=M at x=1
    f[2] = y[P_idx] - 2/3.*Grav*mass/pow(radius,2)/kappa;
    // P = (2/3)(GM/R^2)/kappa at x=1
}

#include<fstream>

void HL::save(const char *fname) {
    int i,j;
    double A(data[A_idx][0]);
    std::ofstream s(fname, std::ofstream::binary);
    s.write((const char*)&X, sizeof(X));
    s.write((const char*)&Y, sizeof(Y));
    s.write((const char*)&N, sizeof(N));
    s.write((const char*)&A, sizeof(A));
    for(j=0; j<N; j++)
        for(i=0; i<4; i++)
            s.write((const char*)&data[i][j], sizeof(double));
}

void HL::load(const char *fname) {
    int i,j;
    double A;
    std::ifstream s(fname, std::ifstream::binary);
    s.read((char *)&X, sizeof(X));
    s.read((char *)&Y, sizeof(Y));
    s.read((char *)&N, sizeof(N));
    s.read((char *)&A, sizeof(A));
    data = Mat_DP(5,N);
    for(j=0; j<N; j++) {
        for(i=0; i<4; i++)
            s.read((char *)&data[i][j], sizeof(double));
        data[A_idx][j] = A;
    }
    mass = data[m_idx][N-1];
    radius = data[r_idx][N-1];
}

void HL::draw(const char *fname, double R1, double R2, int n)
// draw Hayashi-line on L-T diagram by changing R in [R1,R2]
// n = number of rows in file fname
// each row consists of (R,L,T)
{
    double R,dR(pow(R2/R1, 1./n));
    std::ofstream f(fname);
    for(int i=0; i<=n; i++) {
        R = R1*pow(dR,i);
        build(mass, R, X, Y);
        f << R << ' ';
        f << surface_temperature() << ' ';
        f << luminosity() << '\n';
    }    
}