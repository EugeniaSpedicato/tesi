#ifndef ECAL_h
#define ECAL_h

#include <TH2.h>
#include <cmath>
#include <iostream>

using namespace std;

class ECAL : public TH2
{
public:

//costruttore
ECAL(Int_t nbinsx, 
    Double_t xlow, 
    Double_t xup, 
    Int_t nbinsy, 
    Double_t ylow, 
    Double_t yup);
//distruttore
~ ECAL(){}

TH2F  *EcalGrid;
double_t radlen;


void AddHit(Double_t x,Double_t y,Double_t w);
void Draw_ECAL(Double_t w);


private:
const Int_t nbinX;
const Int_t nbinY;
const Int_t Xlow;
const Int_t Xup;
const Int_t Ylow;
const Int_t Yup;



};
#endif