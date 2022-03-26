#include<cmath>
#include<fstream>
#include "HL.h"
#include "constants.h"

main() {
    double Y(0.27);
    HL h;
    h.build(Msun, 0.5*Rsun, 0.7, Y);
    h.draw("fig4a.txt", 0.5*Rsun, 2.5*Rsun, 20);
    h.build(Msun, 0.5*Rsun, 0.71, Y);
    h.draw("fig4b.txt", 0.5*Rsun, 2.5*Rsun, 20);
    h.build(Msun, 0.5*Rsun, 0.72, Y);
    h.draw("fig4c.txt", 0.5*Rsun, 2.5*Rsun, 20);
    h.build(Msun, 0.5*Rsun, 0.73, Y);
    h.draw("fig4d.txt", 0.5*Rsun, 2.5*Rsun, 20);
}