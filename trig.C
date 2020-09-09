#define atree_cxx
#include "next.h"
#include <TH2.h>
#include <TH1.h>
#include <TGraph.h>
#include <TTree.h>

#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{
    TH1::SetDefaultSumw2();
    

    
    
    
     if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
    Double_t Ee[nentries];
    Double_t The[nentries]; 
   
    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
 
       
       Ee[jentry]=detKinBeamRot_Ee;
       The[jentry]=detKinBeamRot_the;
  
       
   }
    
 TGraph *energyThEl= new TGraph(4,The,Ee); 
    energyThEl->SetTitle("Energy_e(theta_e)");
    energyThEl->SetMarkerColor(50);
    energyThEl->SetMarkerStyle(8);
    energyThEl->SetLineColor(9);
    energyThEl->GetXaxis()->SetTitle("ThetaEl(mrad)");
    energyThEl->GetYaxis()->SetTitle("Energy(GeV)");   
    
    TCanvas * theE= new TCanvas("theE","theE",1000,100,2500,2000);
    energyThEl->Draw();
    energyThEl->SaveAs("thetaEn.png")
    
    
}