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
    
    TH1F* E_el=new TH1F("h1aN", "Energy electron tar 0", 200,0,160);
    TH1F* E_elLO=new TH1F("h1aN", "Energy electron tar 0", 200,0,160);
    
    
    
  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
    
   
    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
 
       if (detKinBeamRot_tar==0)
       {
        E_el->Fill(detKinBeamRot_Ee,wgt_full);
        E_elLO->Fill(detKinBeamRot_Ee,wgt_LO);
       }
       
       
}
 TCanvas * e= new TCanvas("e","e",400,10,1500,1000);
    e->Divide(3,2);
    e->cd(1);
    E_el->SetLineColor(46);
    E_el->GetXaxis()->SetTitle("E [GeV] log scale");
    E_el->Draw("HIST");
    E_elLO->SetLineColor(32);
    E_elLO->Draw("HIST same");
    gPad->SetLogx();
    
    
    e->SaveAs("energy_el0.png");
    