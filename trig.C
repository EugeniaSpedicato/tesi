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
    
Double_t Rm= 0.022; //raggio di Moliere in metri
Double_t d; //distanza elettronce fotone.
Double_t z  =2.025; //distanza origine-calorimetro in metri.
Double_t nEl=0; // numero di elettroni con un fotone
Double_t nElNO=0; // numero di elettroni senzs un fotone
    
 
    
    
     if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
    
   //TGraph *energyThEl= new TGraph(nentries); 

    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
   
       //energyThEl->SetPoint(jentry,detKinBeamRot_the,detKinBeamRot_Ee);
       if (detKinBeamRot_cooXe < 0.07 && detKinBeamRot_cooYe < 0.07 && detKinBeamRot_cooXe > -0.07 && detKinBeamRot_cooYe > -0.07)
        {
       nElNO++;
       
double_t posEl = sqrt(detKinBeamRot_cooXe*detKinBeamRot_cooXe+detKinBeamRot_cooYe*detKinBeamRot_cooYe);
double_t posPh = sqrt(photon_coox*photon_coox+photon_cooy*photon_cooy);
d=posEl-posPh;   
  
       if (photon_coox < 0.07 && photon_cooy < 0.07 && photon_coox > -0.07 && photon_cooy > -0.07 && abs(d)>2*Rm)
       {
          nEl++; 
       }
   
   
           }
   }
    
    
    

    /*energyThEl->SetTitle("Energy_e(theta_e)");
    energyThEl->SetMarkerColor(50);
    energyThEl->SetMarkerStyle(8);
    energyThEl->SetLineColor(9);
    energyThEl->SetMarkerStyle(kFullDotSmall);
    energyThEl->GetXaxis()->SetTitle("ThetaEl(mrad)");
    energyThEl->GetYaxis()->SetTitle("Energy(GeV)");   
    
    TCanvas * theE= new TCanvas("theE","theE",1000,100,2500,2000);
    energyThEl->Draw("AP");
    energyThEl->SaveAs("thetaEn.png");*/
    cout << "numero di elettroni senza fotoni: " << nElNO << endl;
    cout << "numero di elettroni con fotoni: " << nEl << endl;
    
    
    
}