#include<cmath>
#include<fstream>
#include "HL.h"
#include "constants.h"

main() {
    double X(0.71),Y(0.27),a(0.1*pow(10,0.1));
    HL h;
    h.build(Msun, 0.1*Rsun, X, Y);
    h.draw("fig3a.txt", 0.1*Rsun, 100*Rsun, 30);
    h.build(2*Msun, 0.1*Rsun, X, Y);
    h.draw("fig3b.txt", 0.1*Rsun, 100*Rsun, 30);
    h.build(4*Msun, 0.1*Rsun, X, Y);
    h.draw("fig3c.txt", 0.1*Rsun, 100*Rsun, 30);
    h.build(Msun/2, 0.1*Rsun, X, Y);
    h.draw("fig3d.txt", 0.1*Rsun, 100*Rsun, 30);
    h.build(Msun/4, a*Rsun, X, Y);
    h.draw("fig3e.txt", a*Rsun, 100*Rsun, 29);
}