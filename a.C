#define atree_cxx
#include "next.h"
#include <TH2.h>
#include <TH1.h>

#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{
    TH1::SetDefaultSumw2();
TH1F* px_mu=new TH1F("h1", "pX_in muon without divergence", 190,-0.2,0.2);
TH1F* py_mu=new TH1F("h2", "pY_in muon without divergence", 190,-0.2,0.2);
TH1F* pz_mu=new TH1F("h3", "pZ_in muon without divergence", 190,120,180);    
    
if (fChain == 0) return;

Long64_t nentries = fChain->GetEntriesFast();
TGraph* E3x3 = new TGraph(nentries);
TGraph* E3x3noph = new TGraph(nentries);
    
    


    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
       
px_mu->Fill(detKinBeamRot_pXmu,wgt_full);
py_mu->Fill(detKinBeamRot_pYmu,wgt_full);
pz_mu->Fill(detKinBeamRot_pZmu,wgt_full);
}
  TCanvas * Pin= new TCanvas("Pin","Pin",1500,1500,1000,500);
    Pin->Divide(3,1);
    Pin->cd(1);
    px_mu->SetLineWidth(3);
    px_mu->SetLineColor(kRed);
    px_mu->Draw("HIST");
    px_mu->GetXaxis()->SetTitle("Px [GeV]");
    Pin->cd(2);
    py_mu->SetLineWidth(3);
    py_mu->SetLineColor(kRed);
    py_mu->Draw("HIST");
    py_mu->GetXaxis()->SetTitle("Py [GeV]");
    Pin->cd(3);
    pz_mu->SetLineWidth(3);
    pz_mu->SetLineColor(kRed);
    pz_mu->Draw("HIST");
    pz_mu->GetXaxis()->SetTitle("Pz [GeV]");
    
   Pin->SaveAs("/home/LHCB-T3/espedicato/tesi/p_inNODIV.png");    
    
    
    
    
    
}
