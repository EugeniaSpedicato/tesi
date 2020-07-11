#define atree_cxx
#include "atree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{
TH1F* E_mu=new TH1F("h1", "Energy muon", 150,0,200);
TH1F* E_e=new TH1F("h2", "Energy electron", 150,0,200);

TH1F* px_mu=new TH1F("h1", "pX_in muon", 190,-0.3,0.3);
TH1F* py_mu=new TH1F("h2", "pY_in muon", 190,-0.3,0.3);
TH1F* pz_mu=new TH1F("h3", "pZ_in muon", 190,50,180);
    
TH1F* px_mu_out=new TH1F("h1a", "pX_out muon", 150,-0.3,0.3);
TH1F* py_mu_out=new TH1F("h2a", "pY_out muon", 150,-0.3,0.3);
TH1F* pz_mu_out=new TH1F("h3a", "pZ_out muon", 150,0,180);
    
TH1F* px_e_out=new TH1F("h1b", "pX_out electron", 150,-0.3,0.3);
TH1F* py_e_out=new TH1F("h2b", "pY_out electron", 150,-0.3,0.3);
TH1F* pz_e_out=new TH1F("h3b", "pZ_out electron", 150,0,5);
 

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

       E_mu->Fill(detKinBeamRot_Emu);
       E_e->Fill(detKinBeamRot_Ee);
    
       
       px_mu->Fill(detKinBeamRot_pXmu);
       py_mu->Fill(detKinBeamRot_pYmu);
       pz_mu->Fill(detKinBeamRot_pZmu);
       
       px_mu_out->Fill(detKinBeamRot_pXmu_out);
       py_mu_out->Fill(detKinBeamRot_pYmu_out);
       pz_mu_out->Fill(detKinBeamRot_pZmu_out);
       
       px_e_out->Fill(detKinBeamRot_pXe_out);
       py_e_out->Fill(detKinBeamRot_pYe_out);
       pz_e_out->Fill(detKinBeamRot_pZe_out);
       
       
       
   }

    px_mu->Fit("gaus");
    py_mu->Fit("gaus");
    pz_mu->Fit("gaus");
    
    
    TCanvas * E= new TCanvas("E","E",400,10,1500,1000);
    E->Divide(2,1);
    E->cd(1);
    E_mu->Draw();
    E->cd(2);
    E_e->SetLineColor(kRed);
    E_e->Draw();
    
    TCanvas * Pin= new TCanvas("Pin","Pin",400,10,1500,1000);
    Pin->Divide(1,3);
    Pin->cd(1);
    px_mu->Draw();
    Pin->cd(2);
    py_mu->Draw();
    Pin->cd(3);
    pz_mu->Draw();
        
    TCanvas * diff= new TCanvas("diff","diff",400,10,1500,1000);
    diff->Divide(2,2);
    diff->cd(1);
    px_e_out->SetLineColor(kRed);
    px_e_out->Draw();
    px_mu_out->Draw("same");
    
    
    diff->cd(2);
    py_e_out->SetLineColor(kRed);
    py_e_out->Draw();
    py_mu_out->Draw("same");
    
    
    diff->cd(3);
    pz_mu_out->Draw();
    
    diff->cd(4);
    pz_e_out->Draw();
    
}
