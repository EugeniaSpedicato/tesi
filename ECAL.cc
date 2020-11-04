#include <TH2.h>
#include <cmath>
#include <iostream>
#include "ECAL.h"

using namespace std;

ECAL::ECAL(Int_t nbinsx, 
    Double_t xlow, 
    Double_t xup, 
    Int_t nbinsy, 
    Double_t ylow, 
    Double_t yup)
    :
    nbinX(nbinsx),nbinY(nbinsy),Xlow(xlow),Xup(xup),Ylow(ylow),Yup(yup) {

TH2F  *EcalGrid = new TH2F("EcalGrid" , "EM Calorimeter",nbinsx,xlow,xup,nbinsy,ylow,yup);
double_t radlen=0.;}

void ECAL::AddHit(Double_t x, Double_t y, Double_t w){EcalGrid->Fill(x,y,w);};
void ECAL::Draw_ECAL(Double_t w){

Int_t nx = EcalGrid->GetNbinsX();
Int_t ny = EcalGrid->GetNbinsY();
for (Int_t i=1; i<nx+1; i+=w) {
for (Int_t j=1; j<ny+1; j+=w) {
    if (EcalGrid->GetBinContent(i,j)<1) EcalGrid->SetBinContent(i,j,0);}; }       
    
TCanvas * Ecal_= new TCanvas("Ecal_","Ecal_",1000,100,2500,2000); 
EcalGrid->Draw();

};



