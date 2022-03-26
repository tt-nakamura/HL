#include<cmath>
#include<fstream>
#include "HL.h"
#include "constants.h"

main() {
    double X(0.71), Y(0.27);
    extern double mix_alpha;    
    HL h;
    mix_alpha = 1;
    h.build(Msun, Rsun, X, Y);
    h.draw("fig5a.txt", 0.5*Rsun, 5*Rsun, 20);
    mix_alpha = 1.5;
    h.build(Msun, Rsun, X, Y);
    h.draw("fig5b.txt", 0.5*Rsun, 5*Rsun, 20);
    mix_alpha = 2;
    h.build(Msun, Rsun, X, Y);
    h.draw("fig5c.txt", 0.5*Rsun, 5*Rsun, 20);
}