#include<cmath>
#include<fstream>
#include "constants.h"

void main_seq(double& L, double& R, double& T, double M)
// main sequence star (empirical relation)
// input:
//   M = mass / g
// output:
//   R = radius / cm
//   L = luminosity / erg/s
//   T = surface termperature / K
// reference: A. N. Cox, et al
//  "Allen's Astrophysical Quantities" 4th edition, p382
{
    double lM = log(M/Msun);
    L = Lsun*pow(10, 3.8*lM + 0.08);
    R = Rsun*pow(10, lM>0.12 ? 0.64*lM + 0.011 : 0.917*lM - 0.020);
    T = pow(L/4/PI/R/R/Stef, 0.25);
}

main() {
    int n(30),i;
    std::ofstream f("main_seq.txt");
    double M1(1e33), M2(3e33), dM(pow(M2/M1, 1./n));
    double M,L,R,T;
    for(i=0; i<=n; i++) {
        M = M1*pow(dM,i);
        main_seq(L,R,T,M);
        f << M << ' ';
        f << L << ' ';
        f << R << ' ';
        f << T << '\n';
    }
}