#include <iostream>
#include "ECAL.h"
#include <fstream>
using namespace std;

int main(){

ECAL TheEcal(5,-7.125,7.125,5,-7.125,7.125);

Double_t xx= 4;
Double_t yy= 2;
Double_t ww= 1;

TheEcal.AddHit(xx,yy,ww);
TheEcal.Draw_ECAL(ww);
    

return 0;
}