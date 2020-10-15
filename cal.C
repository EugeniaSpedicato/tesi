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
Double_t n_one=0.;
Double_t n_two=0.;

        

//Double_t same_cell=0.;
//Double_t different_cell=0.;
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
       
  /*if (photon_n_cell_ph!=0 && detKinBeamRot_n_cell_e!=0){ 
      
if (d_e_ph>1*Rm)
{   n_tot+=wgt_full;
 if (photon_n_cell_ph==detKinBeamRot_n_cell_e)
    {same_cell+=wgt_full;}
    
    else {different_cell+=wgt_full;}  // contatori double 
}
  }*/
       
       if (detKinBeamRot_n_cell_e!=0)
    {n_tot+=wgt_full;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (photon_n_cell_ph!=0)    
           {
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=N*RM DALL'ELETTRONE             
                if (d_e_ph>1*Rm)
                {
                    //if (photon_energy>0.2) {
                    n_two+=wgt_full;
                }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { 
                //if (photon_energy>0.2) {
                n_one+=wgt_full;
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one+=wgt_full;}
    } 

       
       
 
}

    
 /*cout << "Elettroni e fotoni nella stessa cella: " << same_cell << endl;
 cout << "Elettroni e fotoni in una diversa cella: " << different_cell << endl;   
cout << "n tot con fotoni= " << n_tot << endl;*/
    
}