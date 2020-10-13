#define atree_cxx
#include "next.h"
#include <TH2.h>
#include <TH1.h>
#include <TGraph.h>
#include <TTree.h>
#include <cmath>
using namespace std;

#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{
    TH1::SetDefaultSumw2();
    
Int_t n_cell; //numero di cella in cui cade l'ELETTRONE
Int_t n_cell_ph; //numero di cella in cui cade il fotone
Int_t n_tot=0;

Double_t same_cell=0.;
Double_t different_cell=0.;
Double_t E_CAL;
Double_t Rm = 1.959 ; //raggio di Moliere in centimetri    

TF2 *fxy_e = new TF2("fxy_e","(1/(2*3.14159265358979323846*[0]*[1]))*(exp(((x-[2])*(x-[2])+(y-[3])*(y-[3]))/(2*[0]*[1])))",-7.125,7.125,-7.125,7.125 );
    
TF2 *fxy_ph = new TF2("fxy_ph","(1/(2*3.14159265358979323846*[0]*[1]))*(exp(((x-[2])*(x-[2])+(y-[3])*(y-[3]))/(2*[0]*[1])))",-7.125,7.125,-7.125,7.125 );
  
Double_t E1;
Double_t E2;
Double_t E3;
Double_t E4;
Double_t E5;
Double_t E6;
Double_t E7;
Double_t E8;
Double_t E9;
Double_t E10;
Double_t E11;
Double_t E12;
Double_t E13;
Double_t E14;
Double_t E15;
Double_t E16;
Double_t E17;
Double_t E18;
Double_t E19;
Double_t E20;
Double_t E21;
Double_t E22;
Double_t E23;
Double_t E24;
Double_t E25;
    
Double_t phE1;
Double_t phE2;
Double_t phE3;
Double_t phE4;
Double_t phE5;
Double_t phE6;
Double_t phE7;
Double_t phE8;
Double_t phE9;
Double_t phE10;
Double_t phE11;
Double_t phE12;
Double_t phE13;
Double_t phE14;
Double_t phE15;
Double_t phE16;
Double_t phE17;
Double_t phE18;
Double_t phE19;
Double_t phE20;
Double_t phE21;
Double_t phE22;
Double_t phE23;
Double_t phE24;
Double_t phE25;
    
    

    
    if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();


    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
        detKinBeamRot_cooXe=detKinBeamRot_cooXe*100;
        detKinBeamRot_cooYe=detKinBeamRot_cooYe*100;
        detKinBeamRot_cooXmu=detKinBeamRot_cooXmu*100;
        detKinBeamRot_cooYmu=detKinBeamRot_cooYmu*100;
        photon_coox=photon_coox*100;
        photon_cooy=photon_cooy*100;
        if (photon_coox!=-1 && photon_cooy!=-1)
    {  
     E_CAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {  E_CAL=detKinBeamRot_Ee;}
       
       
    Double_t d_e_ph=sqrt( (detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy) ); 
       
    Double_t d_e_mu=sqrt( (detKinBeamRot_cooXe-detKinBeamRot_cooXmu)*(detKinBeamRot_cooXe-detKinBeamRot_cooXmu)+(detKinBeamRot_cooYe-detKinBeamRot_cooYmu)*(detKinBeamRot_cooYe-detKinBeamRot_cooYmu) );
       
   
       
       
//if (E_CAL>1){
if (abs(detKinBeamRot_cooXe)<7.125 && abs(detKinBeamRot_cooYe)<7.125)
    
{ 
    fxy_e->SetParameter(0,Rm);
    fxy_e->SetParameter(1,Rm);
    fxy_e->SetParameter(2,detKinBeamRot_cooXe);
    fxy_e->SetParameter(3,detKinBeamRot_cooYe);
       
    
    
    if (detKinBeamRot_cooXe>-7.125 && detKinBeamRot_cooXe<-4.275 && detKinBeamRot_cooYe<7.125 && detKinBeamRot_cooYe>4.275) {
        n_cell=1; 
        Double_t E1 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,7.125,4.275);
        Double_t E2 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,7.125,4.275);
        Double_t E6 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,4.275,1.425);
        Double_t E7 = detKinBeamRot_Ee*fxy_e->Integer(4.275,1.425,4.275,1.425);
        
                }
    
    if (detKinBeamRot_cooXe>-4.275 && detKinBeamRot_cooXe<-1.425 && detKinBeamRot_cooYe<7.125 && detKinBeamRot_cooYe>4.275) {
        n_cell=2;
         E2 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,7.125,4.275);
         E1 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,7.125,4.275);
         E6 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,4.275,1.425);
         E7 = detKinBeamRot_Ee*fxy_e->Integer(4.275,1.425,4.275,1.425);
         E8 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,4.275,1.425);
         E3 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,7.125,4.275);
        
    }

    if (detKinBeamRot_cooXe>-1.425 && detKinBeamRot_cooXe<1.425 && detKinBeamRot_cooYe<7.125 && detKinBeamRot_cooYe>4.275) {
        n_cell=3;
         E3 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,7.125,4.275);
         E2 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,7.125,4.275);
         E7 = detKinBeamRot_Ee*fxy_e->Integer(4.275,1.425,4.275,1.425);
         E8 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,4.275,1.425);
         E9 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,4.275,1.425);
         E4 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,7.125,4.275);
    
    }

    if (detKinBeamRot_cooXe>1.425 && detKinBeamRot_cooXe<4.275 && detKinBeamRot_cooYe<7.125 && detKinBeamRot_cooYe>4.275) {
        n_cell=4;        
         E4 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,7.125,4.275);
         E3 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,7.125,4.275);
         E8 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,4.275,1.425);
         E9 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,4.275,1.425);
         E10 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,4.275,1.425);
         E5 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,7.125,4.275);
    }

    if (detKinBeamRot_cooXe>4.275 && detKinBeamRot_cooXe<7.125 && detKinBeamRot_cooYe<7.125 && detKinBeamRot_cooYe>4.275) {
        n_cell=5;
         E5 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,7.125,4.275);
         E4 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,7.125,4.275);
         E9 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,4.275,1.425);
         E10 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,4.275,1.425);
    }
    
    
    if (detKinBeamRot_cooXe>-7.125 && detKinBeamRot_cooXe<-4.275 && detKinBeamRot_cooYe<4.275 && detKinBeamRot_cooYe>1.425) {
        n_cell=6;
         E6 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,7.125,4.275);
         E1 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,7.125,4.275);
         E2 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,4.275,1.425);
         E7 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,4.275,1.425);
         E12 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,4.275,1.425);
         E11 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,7.125,4.275);
    }
    
    if (detKinBeamRot_cooXe>-4.275 && detKinBeamRot_cooXe<-1.425 && detKinBeamRot_cooYe<4.275 && detKinBeamRot_cooYe>1.425) {
        n_cell=7;
         E7 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,4.275,1.425);
         E1 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,7.125,4.275);
         E2 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,4.275,1.425);
         E3 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,7.125,4.275);
         E8 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,4.275,1.425);
         E13 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,7.125,4.275);
         E12 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,4.275,1.425);
         E11 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,7.125,4.275);
         E6 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,7.125,4.275);
        
    }

    if (detKinBeamRot_cooXe>-1.425 && detKinBeamRot_cooXe<1.425 && detKinBeamRot_cooYe<4.275 && detKinBeamRot_cooYe>1.425) {
        n_cell=8;
         E8 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,4.275,1.425);
         E2 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,4.275,1.425);
         E3 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,7.125,4.275);
         E4 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,7.125,4.275);
         E9 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,4.275,1.425);
         E14 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,1.425,-1.425);
         E13 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,7.125,4.275);
         E12 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,4.275,1.425);
         E7 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,4.275,1.425);
    }

    if (detKinBeamRot_cooXe>1.425 && detKinBeamRot_cooXe<4.275 && detKinBeamRot_cooYe<4.275 && detKinBeamRot_cooYe>1.425) {
        n_cell=9;
         E9 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,4.275,1.425);
         E3 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,7.125,4.275);
         E4 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,7.125,4.275);
         E5 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,7.125,4.275);
         E10 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,4.275,1.425);
         E15 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,1.425,-1.425);
         E14 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,1.425,-1.425);
         E13 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,7.125,4.275);
         E8 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,4.275,1.425);
    }

    if (detKinBeamRot_cooXe>4.275 && detKinBeamRot_cooXe<7.125 && detKinBeamRot_cooYe<4.275 && detKinBeamRot_cooYe>1.425) {
        n_cell=10;
         E10 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,4.275,1.425);
         E5 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,7.125,4.275);
         E4 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,7.125,4.275);
         E9 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,4.275,1.425);
         E14 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,1.425,-1.425);
         E15 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,1.425,-1.425);
        
    }
    

    
    if (detKinBeamRot_cooXe>-7.125 && detKinBeamRot_cooXe<-4.275 && detKinBeamRot_cooYe<1.425 && detKinBeamRot_cooYe>-1.425) {
        n_cell=11;
         E11 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,1.425,-1.425);
         E6 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,7.125,4.275);
         E7 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,4.275,1.425);
         E12 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,4.275,1.425);
         E17 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-1.425,-4.275);
         E16 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,-1.425,-4.275);
    }
    
    if (detKinBeamRot_cooXe>-4.275 && detKinBeamRot_cooXe<-1.425 && detKinBeamRot_cooYe<1.425 && detKinBeamRot_cooYe>-1.425) {
        n_cell=12;
         E12 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-1.425,1.425);
         E6 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,7.125,4.275);
         E7 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,4.275,1.425);
         E8 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,4.275,1.425);
         E13 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,7.125,4.275);
         E18 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-1.425,-4.275);
         E17 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-1.425,-4.275);
         E16 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,-1.425,-4.275);
         E11 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,1.425,-1.425);
    }

    if (detKinBeamRot_cooXe>-1.425 && detKinBeamRot_cooXe<1.425 && detKinBeamRot_cooYe<1.425 && detKinBeamRot_cooYe>-1.425) {
        n_cell=13;
         E13 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,7.125,4.275);
         E7 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,4.275,1.425);
         E8 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,4.275,1.425);
         E9 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,4.275,1.425);
         E14 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,1.425,-1.425);
         E19 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,-1.425,-4.275);
         E18 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-1.425,-4.275);
         E17 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-1.425,-4.275);
         E12 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-1.425,1.425);
    }

    if (detKinBeamRot_cooXe>1.425 && detKinBeamRot_cooXe<4.275 && detKinBeamRot_cooYe<1.425 && detKinBeamRot_cooYe>-1.425) {
        n_cell=14;
         E14 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,1.425,-1.425);
         E8 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,4.275,1.425);
         E9 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,4.275,1.425);
         E10 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,4.275,1.425);
         E15 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,1.425,-1.425);
         E20 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,-1.425,-4.275);
         E19 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,-1.425,-4.275);
         E18 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-1.425,-4.275);
         E13 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,7.125,4.275);
    }

    if (detKinBeamRot_cooXe>4.275 && detKinBeamRot_cooXe<7.125 && detKinBeamRot_cooYe<1.425 && detKinBeamRot_cooYe>-1.425) {
        n_cell=15;
         E15 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,1.425,-1.425);
         E10 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,4.275,1.425);
         E9 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,4.275,1.425);
         E14 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,1.425,-1.425);
         E19 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,-1.425,-4.275);
         E20 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,-1.425,-4.275);
    
    }
    
    
    
    
    if (detKinBeamRot_cooXe>-7.125 && detKinBeamRot_cooXe<-4.275 && detKinBeamRot_cooYe<-1.425 && detKinBeamRot_cooYe>-4.275) {
        n_cell=16;
         E16 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,-1.425,-4.275);
         E11 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,1.425,-1.425);
         E12 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-1.425,1.425);
         E17 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-1.425,-4.275);
         E22 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-4.275,-7.125);
         E21 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,-4.275,-7.125);
}
    
    if (detKinBeamRot_cooXe>-4.275 && detKinBeamRot_cooXe<-1.425 && detKinBeamRot_cooYe<-1.425 && detKinBeamRot_cooYe>-4.275) {
        n_cell=17;
         E17 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-1.425,-4.275);
         E11 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,1.425,-1.425);
         E12 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-1.425,1.425);
         E13 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,7.125,4.275);
         E18 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-1.425,-4.275);
         E23 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-4.275,-7.125);
         E22 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-4.275,-7.125);
         E21 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,-4.275,-7.125);
         E16 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,-1.425,-4.275);
    }

    if (detKinBeamRot_cooXe>-1.425 && detKinBeamRot_cooXe<1.425 && detKinBeamRot_cooYe<-1.425 && detKinBeamRot_cooYe>-4.275) {
        n_cell=18;
         E18 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-1.425,-4.275);
         E12 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-1.425,1.425);
         E13 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,7.125,4.275);
         E14 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,1.425,-1.425);
         E19 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,-1.425,-4.275);
         E24 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,-4.275,-7.125);
         E23 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-4.275,-7.125);
         E22 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-4.275,-7.125);
         E17 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-1.425,-4.275); 
    }

    if (detKinBeamRot_cooXe>1.425 && detKinBeamRot_cooXe<4.275 && detKinBeamRot_cooYe<-1.425 && detKinBeamRot_cooYe>-4.275) {
        n_cell=19;
         E19 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,-1.425,-4.275);
         E13 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,7.125,4.275);
         E14 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,1.425,-1.425);
         E15 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,1.425,-1.425);
         E20 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,-1.425,-4.275);
         E25 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,-4.275,-7.125);
         E24 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,-4.275,-7.125);
         E23 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-4.275,-7.125);
         E18 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-1.425,-4.275); 
    }

    if (detKinBeamRot_cooXe>4.275 && detKinBeamRot_cooXe<7.125 && detKinBeamRot_cooYe<-1.425 && detKinBeamRot_cooYe>-4.275) {
        n_cell=20;
         E20 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,-1.425,-4.275);
         E15 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,1.425,-1.425);
         E14 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,1.425,-1.425);
         E19 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,-1.425,-4.275);
         E24 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,-4.275,-7.125);
         E25 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,-4.275,-7.125);
        
    }
    
    
    
    if (detKinBeamRot_cooXe>-7.125 && detKinBeamRot_cooXe<-4.275 && detKinBeamRot_cooYe<-4.275 && detKinBeamRot_cooYe>-7.125) {
        n_cell=21;
         E21 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,-4.275,-7.125);
         E14 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,1.425,-1.425);
         E17 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-1.425,-4.275);
         E22 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-4.275,-7.125);
    
    }
    
    if (detKinBeamRot_cooXe>-4.275 && detKinBeamRot_cooXe<-1.425 && detKinBeamRot_cooYe<-4.275 && detKinBeamRot_cooYe>-7.125) {
        n_cell=22;
         E22 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-4.275,-7.125);
         E21 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,-4.275,-7.125);
         E16 = detKinBeamRot_Ee*fxy_e->Integer(-7.125,-4.275,-1.425,-4.275);
         E17 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-1.425,-4.275);
         E18 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-1.425,-4.275);
         E23 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-4.275,-7.125);
    
    }

    if (detKinBeamRot_cooXe>-1.425 && detKinBeamRot_cooXe<1.425 && detKinBeamRot_cooYe<-4.275 && detKinBeamRot_cooYe>-7.125) {
        n_cell=23;
         E23 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-4.275,-7.125);
         E22 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-4.275,-7.125);
         E17 = detKinBeamRot_Ee*fxy_e->Integer(-4.275,-1.425,-1.425,-4.275);
         E18 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-1.425,-4.275);
         E19 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,-1.425,-4.275);
         E24 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,-4.275,-7.125);
    }

    if (detKinBeamRot_cooXe>1.425 && detKinBeamRot_cooXe<4.275 && detKinBeamRot_cooYe<-4.275 && detKinBeamRot_cooYe>-7.125) {
        n_cell=24;
         E24 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,-4.275,-7.125);
         E23 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-4.275,-7.125);
         E18 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-1.425,-4.275);
         E19 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-1.425,-4.275);
         E20 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,-1.425,-4.275);
         E25 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,-4.275,-7.125);
    }

    if (detKinBeamRot_cooXe>4.275 && detKinBeamRot_cooXe<7.125 && detKinBeamRot_cooYe<-4.275 && detKinBeamRot_cooYe>-7.125) {
        n_cell=25;
         E25 = detKinBeamRot_Ee*fxy_e->Integer(4.275,7.125,-4.275,-7.125);
         E24 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,-4.275,-7.125);
         E19 = detKinBeamRot_Ee*fxy_e->Integer(-1.425,1.425,-1.425,-4.275);
         E20 = detKinBeamRot_Ee*fxy_e->Integer(1.425,4.275,-1.425,-4.275);
    }
    
    
  
    
//}
       
       

       
       
       
    //&& d_e_ph>2*Rm
       
       
if (abs(photon_coox)<7.125 && abs(photon_cooy)<7.125 && photon_energy>0.2 && d_e_ph>2*Rm)
    
{ n_tot++;
    
    fxy_ph->SetParameter(0,Rm);
    fxy_ph->SetParameter(1,Rm);
    fxy_ph->SetParameter(2,photon_coox);
    fxy_ph->SetParameter(3,photon_cooy);
 
    if (photon_coox>-7.125 && photon_coox<-4.275 && photon_cooy<7.125 && photon_cooy>4.275) {
        n_cell_ph=1;
         phE1 = photon_energy*fxy_ph->Integer(-7.125,-4.275,7.125,4.275);
         phE2 = photon_energy*fxy_ph->Integer(-4.275,-1.425,7.125,4.275);
         phE6 = photon_energy*fxy_ph->Integer(-7.125,-4.275,4.275,1.425);
         phE7 = photon_energy*fxy_ph->Integer(4.275,1.425,4.275,1.425);

}
    
    if (photon_coox>-4.275 && photon_coox<-1.425 && photon_cooy<7.125 && photon_cooy>4.275) {
        n_cell_ph=2;
         phE2 = photon_energy*fxy_ph->Integer(-4.275,-1.425,7.125,4.275);
         phE1 = photon_energy*fxy_ph->Integer(-7.125,-4.275,7.125,4.275);
         phE6 = photon_energy*fxy_ph->Integer(-7.125,-4.275,4.275,1.425);
         phE7 = photon_energy*fxy_ph->Integer(4.275,1.425,4.275,1.425);
         phE8 = photon_energy*fxy_ph->Integer(-1.425,1.425,4.275,1.425);
         phE3 = photon_energy*fxy_ph->Integer(-1.425,1.425,7.125,4.275);
    }

    if (photon_coox>-1.425 && photon_coox<1.425 && photon_cooy<7.125 && photon_cooy>4.275) {
        n_cell_ph=3;
         phE3 = photon_energy*fxy_ph->Integer(-1.425,1.425,7.125,4.275);
         phE2 = photon_energy*fxy_ph->Integer(-4.275,-1.425,7.125,4.275);
         phE7 = photon_energy*fxy_ph->Integer(4.275,1.425,4.275,1.425);
         phE8 = photon_energy*fxy_ph->Integer(-1.425,1.425,4.275,1.425);
         phE9 = photon_energy*fxy_ph->Integer(1.425,4.275,4.275,1.425);
         phE4 = photon_energy*fxy_ph->Integer(1.425,4.275,7.125,4.275);
    }

    if (photon_coox>1.425 && photon_coox<4.275 && photon_cooy<7.125 && photon_cooy>4.275) {
        n_cell_ph=4;
         phE4 = photon_energy*fxy_ph->Integer(1.425,4.275,7.125,4.275);
         phE3 = photon_energy*fxy_ph->Integer(-1.425,1.425,7.125,4.275);
         phE8 = photon_energy*fxy_ph->Integer(-1.425,1.425,4.275,1.425);
         phE9 = photon_energy*fxy_ph->Integer(1.425,4.275,4.275,1.425);
         phE10 = photon_energy*fxy_ph->Integer(4.275,7.125,4.275,1.425);
         phE5 = photon_energy*fxy_ph->Integer(4.275,7.125,7.125,4.275);
    }

    if (photon_coox>4.275 && photon_coox<7.125 && photon_cooy<7.125 && photon_cooy>4.275) {
        n_cell_ph=5;
         phE5 = photon_energy*fxy_ph->Integer(4.275,7.125,7.125,4.275);
         phE4 = photon_energy*fxy_ph->Integer(1.425,4.275,7.125,4.275);
         phE9 = photon_energy*fxy_ph->Integer(1.425,4.275,4.275,1.425);
         phE10 = photon_energy*fxy_ph->Integer(4.275,7.125,4.275,1.425);
    }
    
    
    if (photon_coox>-7.125 && photon_coox<-4.275 && photon_cooy<4.275 && photon_cooy>1.425) {
        n_cell_ph=6;
         phE6 = photon_energy*fxy_ph->Integer(1.425,4.275,7.125,4.275);
         phE1 = photon_energy*fxy_ph->Integer(-1.425,1.425,7.125,4.275);
         phE2 = photon_energy*fxy_ph->Integer(-1.425,1.425,4.275,1.425);
         phE7 = photon_energy*fxy_ph->Integer(1.425,4.275,4.275,1.425);
         phE12 = photon_energy*fxy_ph->Integer(4.275,7.125,4.275,1.425);
         phE11 = photon_energy*fxy_ph->Integer(4.275,7.125,7.125,4.275);
    }
    
    if (photon_coox>-4.275 && photon_coox<-1.425 && photon_cooy<4.275 && photon_cooy>1.425) {
        n_cell_ph=7;
         phE7 = photon_energy*fxy_ph->Integer(1.425,4.275,4.275,1.425);
         phE1 = photon_energy*fxy_ph->Integer(-1.425,1.425,7.125,4.275);
         phE2 = photon_energy*fxy_ph->Integer(-1.425,1.425,4.275,1.425);
         phE3 = photon_energy*fxy_ph->Integer(-1.425,1.425,7.125,4.275);
         phE8 = photon_energy*fxy_ph->Integer(-1.425,1.425,4.275,1.425);
         phE13 = photon_energy*fxy_ph->Integer(4.275,7.125,7.125,4.275);
         phE12 = photon_energy*fxy_ph->Integer(4.275,7.125,4.275,1.425);
         phE11 = photon_energy*fxy_ph->Integer(4.275,7.125,7.125,4.275);
         phE6 = photon_energy*fxy_ph->Integer(1.425,4.275,7.125,4.275);
    }

    if (photon_coox>-1.425 && photon_coox<1.425 && photon_cooy<4.275 && photon_cooy>1.425) {
        n_cell_ph=8;
         phE8 = photon_energy*fxy_ph->Integer(-1.425,1.425,4.275,1.425);
         phE2 = photon_energy*fxy_ph->Integer(-1.425,1.425,4.275,1.425);
         phE3 = photon_energy*fxy_ph->Integer(-1.425,1.425,7.125,4.275);
         phE4 = photon_energy*fxy_ph->Integer(1.425,4.275,7.125,4.275);
         phE9 = photon_energy*fxy_ph->Integer(1.425,4.275,4.275,1.425);
         phE14 = photon_energy*fxy_ph->Integer(1.425,4.275,1.425,-1.425);
         phE13 = photon_energy*fxy_ph->Integer(4.275,7.125,7.125,4.275);
         phE12 = photon_energy*fxy_ph->Integer(4.275,7.125,4.275,1.425);
         phE7 = photon_energy*fxy_ph->Integer(1.425,4.275,4.275,1.425);
    }

    if (photon_coox>1.425 && photon_coox<4.275 && photon_cooy<4.275 && photon_cooy>1.425) {
        n_cell_ph=9;
         phE9 = photon_energy*fxy_ph->Integer(1.425,4.275,4.275,1.425);
         phE3 = photon_energy*fxy_ph->Integer(-1.425,1.425,7.125,4.275);
         phE4 = photon_energy*fxy_ph->Integer(1.425,4.275,7.125,4.275);
         phE5 = photon_energy*fxy_ph->Integer(4.275,7.125,7.125,4.275);
         phE10 = photon_energy*fxy_ph->Integer(4.275,7.125,4.275,1.425);
         phE15 = photon_energy*fxy_ph->Integer(4.275,7.125,1.425,-1.425);
         phE14 = photon_energy*fxy_ph->Integer(1.425,4.275,1.425,-1.425);
         phE13 = photon_energy*fxy_ph->Integer(4.275,7.125,7.125,4.275);
         phE8 = photon_energy*fxy_ph->Integer(-1.425,1.425,4.275,1.425);
    }

    if (photon_coox>4.275 && photon_coox<7.125 && photon_cooy<4.275 && photon_cooy>1.425) {
        n_cell_ph=10;
         phE10 = photon_energy*fxy_ph->Integer(4.275,7.125,4.275,1.425);
         phE5 = photon_energy*fxy_ph->Integer(4.275,7.125,7.125,4.275);
         phE4 = photon_energy*fxy_ph->Integer(1.425,4.275,7.125,4.275);
         phE9 = photon_energy*fxy_ph->Integer(1.425,4.275,4.275,1.425);
         phE14 = photon_energy*fxy_ph->Integer(1.425,4.275,1.425,-1.425);
         phE15 = photon_energy*fxy_ph->Integer(4.275,7.125,1.425,-1.425);
    }
    

    
    if (photon_coox>-7.125 && photon_coox<-4.275 && photon_cooy<1.425 && photon_cooy>-1.425) {
        n_cell_ph=11;
         phE11 = photon_energy*fxy_ph->Integer(-7.125,-4.275,1.425,-1.425);
         phE6 = photon_energy*fxy_ph->Integer(1.425,4.275,7.125,4.275);
         phE7 = photon_energy*fxy_ph->Integer(1.425,4.275,4.275,1.425);
         phE12 = photon_energy*fxy_ph->Integer(4.275,7.125,4.275,1.425);
         phE17 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-1.425,-4.275);
         phE16 = photon_energy*fxy_ph->Integer(-7.125,-4.275,-1.425,-4.275);
    }
    
    if (photon_coox>-4.275 && photon_coox<-1.425 && photon_cooy<1.425 && photon_cooy>-1.425) {
        n_cell_ph=12;
         phE12 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-1.425,1.425);
         phE6 = photon_energy*fxy_ph->Integer(1.425,4.275,7.125,4.275);
         phE7 = photon_energy*fxy_ph->Integer(1.425,4.275,4.275,1.425);
         phE8 = photon_energy*fxy_ph->Integer(-1.425,1.425,4.275,1.425);
         phE13 = photon_energy*fxy_ph->Integer(4.275,7.125,7.125,4.275);
         phE18 = photon_energy*fxy_ph->Integer(-1.425,1.425,-1.425,-4.275);
         phE17 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-1.425,-4.275);
         phE16 = photon_energy*fxy_ph->Integer(-7.125,-4.275,-1.425,-4.275);
         phE11 = photon_energy*fxy_ph->Integer(-7.125,-4.275,1.425,-1.425);
    }

    if (photon_coox>-1.425 && photon_coox<1.425 && photon_cooy<1.425 && photon_cooy>-1.425) {
        n_cell_ph=13;
         phE13 = photon_energy*fxy_ph->Integer(4.275,7.125,7.125,4.275);
         phE7 = photon_energy*fxy_ph->Integer(1.425,4.275,4.275,1.425);
         phE8 = photon_energy*fxy_ph->Integer(-1.425,1.425,4.275,1.425);
         phE9 = photon_energy*fxy_ph->Integer(1.425,4.275,4.275,1.425);
         phE14 = photon_energy*fxy_ph->Integer(1.425,4.275,1.425,-1.425);
         phE19 = photon_energy*fxy_ph->Integer(1.425,4.275,-1.425,-4.275);
         phE18 = photon_energy*fxy_ph->Integer(-1.425,1.425,-1.425,-4.275);
         phE17 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-1.425,-4.275);
         phE12 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-1.425,1.425);
    }

    if (photon_coox>1.425 && photon_coox<4.275 && photon_cooy<1.425 && photon_cooy>-1.425) {
        n_cell_ph=14;
         phE14 = photon_energy*fxy_ph->Integer(1.425,4.275,1.425,-1.425);
         phE8 = photon_energy*fxy_ph->Integer(-1.425,1.425,4.275,1.425);
         phE9 = photon_energy*fxy_ph->Integer(1.425,4.275,4.275,1.425);
         phE10 = photon_energy*fxy_ph->Integer(4.275,7.125,4.275,1.425);
         phE15 = photon_energy*fxy_ph->Integer(4.275,7.125,1.425,-1.425);
         phE20 = photon_energy*fxy_ph->Integer(4.275,7.125,-1.425,-4.275);
         phE19 = photon_energy*fxy_ph->Integer(1.425,4.275,-1.425,-4.275);
         phE18 = photon_energy*fxy_ph->Integer(-1.425,1.425,-1.425,-4.275);
         phE13 = photon_energy*fxy_ph->Integer(4.275,7.125,7.125,4.275);
    }

    if (photon_coox>4.275 && photon_coox<7.125 && photon_cooy<1.425 && photon_cooy>-1.425) {
        n_cell_ph=15;
         phE15 = photon_energy*fxy_ph->Integer(4.275,7.125,1.425,-1.425);
         phE10 = photon_energy*fxy_ph->Integer(4.275,7.125,4.275,1.425);
         phE9 = photon_energy*fxy_ph->Integer(1.425,4.275,4.275,1.425);
         phE14 = photon_energy*fxy_ph->Integer(1.425,4.275,1.425,-1.425);
         phE19 = photon_energy*fxy_ph->Integer(1.425,4.275,-1.425,-4.275);
         phE20 = photon_energy*fxy_ph->Integer(4.275,7.125,-1.425,-4.275);
    }
    
    
    
    
    if (photon_coox>-7.125 && photon_coox<-4.275 && photon_cooy<-1.425 && photon_cooy>-4.275) {
        n_cell_ph=16;
         phE16 = photon_energy*fxy_ph->Integer(-7.125,-4.275,-1.425,-4.275);
         phE11 = photon_energy*fxy_ph->Integer(-7.125,-4.275,1.425,-1.425);
         phE12 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-1.425,1.425);
         phE17 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-1.425,-4.275);
         phE22 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-4.275,-7.125);
         phE21 = photon_energy*fxy_ph->Integer(-7.125,-4.275,-4.275,-7.125);
    }
    
    if (photon_coox>-4.275 && photon_coox<-1.425 && photon_cooy<-1.425 && photon_cooy>-4.275) {
        n_cell_ph=17;
         phE17 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-1.425,-4.275);
         phE11 = photon_energy*fxy_ph->Integer(-7.125,-4.275,1.425,-1.425);
         phE12 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-1.425,1.425);
         phE13 = photon_energy*fxy_ph->Integer(4.275,7.125,7.125,4.275);
         phE18 = photon_energy*fxy_ph->Integer(-1.425,1.425,-1.425,-4.275);
         phE23 = photon_energy*fxy_ph->Integer(-1.425,1.425,-4.275,-7.125);
         phE22 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-4.275,-7.125);
         phE21 = photon_energy*fxy_ph->Integer(-7.125,-4.275,-4.275,-7.125);
         phE16 = photon_energy*fxy_ph->Integer(-7.125,-4.275,-1.425,-4.275);
    }

    if (photon_coox>-1.425 && photon_coox<1.425 && photon_cooy<-1.425 && photon_cooy>-4.275) {
        n_cell_ph=18;
         phE18 = photon_energy*fxy_ph->Integer(-1.425,1.425,-1.425,-4.275);
         phE12 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-1.425,1.425);
         phE13 = photon_energy*fxy_ph->Integer(4.275,7.125,7.125,4.275);
         phE14 = photon_energy*fxy_ph->Integer(1.425,4.275,1.425,-1.425);
         phE19 = photon_energy*fxy_ph->Integer(1.425,4.275,-1.425,-4.275);
         phE24 = photon_energy*fxy_ph->Integer(1.425,4.275,-4.275,-7.125);
         phE23 = photon_energy*fxy_ph->Integer(-1.425,1.425,-4.275,-7.125);
         phE22 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-4.275,-7.125);
         phE17 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-1.425,-4.275);
    }

    if (photon_coox>1.425 && photon_coox<4.275 && photon_cooy<-1.425 && photon_cooy>-4.275) {
        n_cell_ph=19;
         phE19 = photon_energy*fxy_ph->Integer(1.425,4.275,-1.425,-4.275);
         phE13 = photon_energy*fxy_ph->Integer(4.275,7.125,7.125,4.275);
         phE14 = photon_energy*fxy_ph->Integer(1.425,4.275,1.425,-1.425);
         phE15 = photon_energy*fxy_ph->Integer(4.275,7.125,1.425,-1.425);
         phE20 = photon_energy*fxy_ph->Integer(4.275,7.125,-1.425,-4.275);
         phE25 = photon_energy*fxy_ph->Integer(4.275,7.125,-4.275,-7.125);
         phE24 = photon_energy*fxy_ph->Integer(1.425,4.275,-4.275,-7.125);
         phE23 = photon_energy*fxy_ph->Integer(-1.425,1.425,-4.275,-7.125);
         phE18 = photon_energy*fxy_ph->Integer(-1.425,1.425,-1.425,-4.275); 
    }

    if (photon_coox>4.275 && photon_coox<7.125 && photon_cooy<-1.425 && photon_cooy>-4.275) {
        n_cell_ph=20;
         phE20 = photon_energy*fxy_ph->Integer(4.275,7.125,-1.425,-4.275);
         phE15 = photon_energy*fxy_ph->Integer(4.275,7.125,1.425,-1.425);
         phE14 = photon_energy*fxy_ph->Integer(1.425,4.275,1.425,-1.425);
         phE19 = photon_energy*fxy_ph->Integer(1.425,4.275,-1.425,-4.275);
         phE24 = photon_energy*fxy_ph->Integer(1.425,4.275,-4.275,-7.125);
         phE25 = photon_energy*fxy_ph->Integer(4.275,7.125,-4.275,-7.125);
    }
    
    
    
    if (photon_coox>-7.125 && photon_coox<-4.275 && photon_cooy<-4.275 && photon_cooy>-7.125) {
        n_cell_ph=21;
         phE21 = photon_energy*fxy_ph->Integer(-7.125,-4.275,-4.275,-7.125);
         phE14 = photon_energy*fxy_ph->Integer(1.425,4.275,1.425,-1.425);
         phE17 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-1.425,-4.275);
         phE22 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-4.275,-7.125);
    }
    
    if (photon_coox>-4.275 && photon_coox<-1.425 && photon_cooy<-4.275 && photon_cooy>-7.125) {
        n_cell_ph=22;
         phE22 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-4.275,-7.125);
         phE21 = photon_energy*fxy_ph->Integer(-7.125,-4.275,-4.275,-7.125);
         phE16 = photon_energy*fxy_ph->Integer(-7.125,-4.275,-1.425,-4.275);
         phE17 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-1.425,-4.275);
         phE18 = photon_energy*fxy_ph->Integer(-1.425,1.425,-1.425,-4.275);
         phE23 = photon_energy*fxy_ph->Integer(-1.425,1.425,-4.275,-7.125);
    }

    if (photon_coox>-1.425 && photon_coox<1.425 && photon_cooy<-4.275 && photon_cooy>-7.125) {
        n_cell_ph=23;
         phE23 = photon_energy*fxy_ph->Integer(-1.425,1.425,-4.275,-7.125);
         phE22 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-4.275,-7.125);
         phE17 = photon_energy*fxy_ph->Integer(-4.275,-1.425,-1.425,-4.275);
         phE18 = photon_energy*fxy_ph->Integer(-1.425,1.425,-1.425,-4.275);
         phE19 = photon_energy*fxy_ph->Integer(1.425,4.275,-1.425,-4.275);
         phE24 = photon_energy*fxy_ph->Integer(1.425,4.275,-4.275,-7.125);
    }

    if (photon_coox>1.425 && photon_coox<4.275 && photon_cooy<-4.275 && photon_cooy>-7.125) {
        n_cell_ph=24;
         phE24 = photon_energy*fxy_ph->Integer(1.425,4.275,-4.275,-7.125);
         phE23 = photon_energy*fxy_ph->Integer(-1.425,1.425,-4.275,-7.125);
         phE18 = photon_energy*fxy_ph->Integer(-1.425,1.425,-1.425,-4.275);
         phE19 = photon_energy*fxy_ph->Integer(-1.425,1.425,-1.425,-4.275);
         phE20 = photon_energy*fxy_ph->Integer(1.425,4.275,-1.425,-4.275);
         phE25 = photon_energy*fxy_ph->Integer(4.275,7.125,-4.275,-7.125);
    }

    if (photon_coox>4.275 && photon_coox<7.125 && photon_cooy<-4.275 && photon_cooy>-7.125) {
        n_cell_ph=25;
         phE25 = photon_energy*fxy_ph->Integer(4.275,7.125,-4.275,-7.125);
         phE24 = photon_energy*fxy_ph->Integer(1.425,4.275,-4.275,-7.125);
         phE19 = photon_energy*fxy_ph->Integer(-1.425,1.425,-1.425,-4.275);
         phE20 = photon_energy*fxy_ph->Integer(1.425,4.275,-1.425,-4.275);
    }
 
if (n_cell_ph==n_cell)
{same_cell+=wgt_full;}
else {different_cell+=wgt_full;}  // contatori double 
} }
}

}
