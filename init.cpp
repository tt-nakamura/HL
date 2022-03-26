#include<fstream>
#include "HL_shoot.h"
#include "constants.h"

main() {
    int i,j,N(1024);
    double X(0.71), Y(0.27);
    const char *fname("HL_M1R1X71Y27.dat");
    double A(1.7);
    HL_shoot hs(Msun, Rsun, X, Y, N);
    // save initial guess
    std::ofstream s(fname, std::ofstream::binary);
    s.write((const char*)&hs.X, sizeof(hs.X));
    s.write((const char*)&hs.Y, sizeof(hs.Y));
    s.write((const char*)&hs.N, sizeof(hs.N));
    s.write((const char*)&A, sizeof(A));
    for(j=0; j<hs.N; j++) {
        s.write((const char*)&hs.r[j], sizeof(double));
        for(i=0; i<3; i++)
            s.write((const char*)&hs.data[i][j], sizeof(double));
    }
    s.close();
    HL h(fname);// load initial guess
    // slowc = 1e-4 to avoid singular matrix
    h.build(Msun, Rsun, X, Y, 1e-4);
    h.save(fname);
}