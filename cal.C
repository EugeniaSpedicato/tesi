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
Double_t n_tot=0.;

Double_t same_cell=0.;
Double_t different_cell=0.;
Double_t E_CAL;
Double_t Rm = 1.959 ; //raggio di Moliere in centimetri    
 

    
    if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();


    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
        detKinBeamRot_cooXe=detKinBeamRot_cooXe*100;
        detKinBeamRot_cooYe=detKinBeamRot_cooYe*100;
           photon_coox=photon_coox*100;
        photon_cooy=photon_cooy*100;
        if (photon_coox!=-1 && photon_cooy!=-1)
    {  
     E_CAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {  E_CAL=detKinBeamRot_Ee;}
       
       
    Double_t d_e_ph=sqrt( (detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy) ); 
       
    Double_t d_e_mu=sqrt( (detKinBeamRot_cooXe-detKinBeamRot_cooXmu)*(detKinBeamRot_cooXe-detKinBeamRot_cooXmu)+(detKinBeamRot_cooYe-detKinBeamRot_cooYmu)*(detKinBeamRot_cooYe-detKinBeamRot_cooYmu) ); 
       
       
if (d_e_ph>1*Rm)
{
    if (photon_n_cell_ph==detKinBeamRot_n_cell_e && photon_n_cell_ph!=0)
{same_cell+=wgt_full;}
else {different_cell+=wgt_full;}  // contatori double 
}

 
}

    
 cout << "Elettroni e fotoni nella stessa cella: " << same_cell << endl;
 cout << "Elettroni e fotoni in una diversa cella: " << different_cell << endl;   
cout << "n tot con fotoni= " << n_tot << endl;
    
}