#define atree_cxx
#include "next.h"
#include <TH2.h>
#include <TH1.h>
#include <TGraph.h>
#include <TTree.h>
#include <cmath>
#include <TMath.h>
using namespace std;

#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{
    TH1::SetDefaultSumw2();
    
Int_t n_cell; //numero di cella in cui cade l'ELETTRONE
Int_t n_cell_ph; //numero di cella in cui cade il fotone
Int_t n_tot=0;

Double_t same_cell=0;
Double_t different_cell=0;
Double_t E_CAL;
Double_t Rm = 1.959 ; //raggio di Moliere in centimetri 
    
TF2 *fxy = new TF2("fxy","(1/(2*Pi()*[0]*[1]))*(Exp(((x-[2])*(x-[2])+(y-[3])*(y-[3]))/(2*[0]*[1])",-7.125,7.125,-7.125,7.125 );
 

    fxy->SetParameter(0,Rm);
    fxy->SetParameter(1,Rm);
    fxy->SetParameter(2,1);
    fxy->SetParameter(3,3);
          fxy->Draw(); 
}