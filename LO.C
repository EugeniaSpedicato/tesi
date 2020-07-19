#define atree_cxx
#include "LO.h"
#include <TH2.h>
#include <TH1.h>

#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{
    TH1::SetDefaultSumw2();
    
     
           
TH1F* px_mu_out=new TH1F("h1a", "pX_out muon", 150,-0.3,0.3);
TH1F* py_mu_out=new TH1F("h2a", "pY_out muon", 150,-0.3,0.3);
TH1F* pz_mu_out=new TH1F("h3a", "pZ_out muon", 150,0,180);
TH1F* px_mu_outNO=new TH1F("h1aN", "pX_out muon NO DIV", 150,-0.3,0.3);
TH1F* py_mu_outNO=new TH1F("h2aN", "pY_out muon NO DIV", 150,-0.3,0.3);
TH1F* pz_mu_outNO=new TH1F("h3aN", "pZ_out muon NO DIV", 150,0,180);
    
TH1F* px_e_out=new TH1F("h1b", "pX_out electron", 150,-0.3,0.3);
TH1F* py_e_out=new TH1F("h2b", "pY_out electron", 150,-0.3,0.3);
TH1F* pz_e_out=new TH1F("h3b", "pZ_out electron", 150,0,5);
TH1F* px_e_outNO=new TH1F("h1bN", "pX_out electron NO DIV", 150,-0.3,0.3);
TH1F* py_e_outNO=new TH1F("h2bN", "pY_out electron NO DIV", 150,-0.3,0.3);
TH1F* pz_e_outNO=new TH1F("h3bN", "pZ_out electron NO DIV", 150,0,5);
      
    if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
       
       
       px_mu_out->Fill(detKinBeamRot_pXmu_out,wgt_LO);
       py_mu_out->Fill(detKinBeamRot_pYmu_out,wgt_LO);
       pz_mu_out->Fill(detKinBeamRot_pZmu_out,wgt_LO);
       px_mu_outNO->Fill(detKinBeamRot_pXmu_out,wgt_full);
       py_mu_outNO->Fill(detKinBeamRot_pYmu_out,wgt_full);
       pz_mu_outNO->Fill(detKinBeamRot_pZmu_out,wgt_full);
       
       px_e_out->Fill(detKinBeamRot_pXe_out,wgt_LO);
       py_e_out->Fill(detKinBeamRot_pYe_out,wgt_LO);
       pz_e_out->Fill(detKinBeamRot_pZe_out,wgt_LO);
       px_e_outNO->Fill(detKinBeamRot_pXe_out,wgt_full);
       py_e_outNO->Fill(detKinBeamRot_pYe_out,wgt_full);
       pz_e_outNO->Fill(detKinBeamRot_pZe_out,wgt_full);
       
       
   }
    
    TCanvas * diffP= new TCanvas("diff","diff",400,10,1500,1000);
    diffP->Divide(3,2);
    diffP->cd(1);
    px_e_outNO->SetLineColor(kRed);
    px_e_outNO->Draw("HIST");
    px_e_out->SetLineColor(kBlack);
    px_e_out->Draw("HIST same");

    
    
    
    diffP->cd(2);
    py_e_outNO->SetLineColor(kRed);
    py_e_outNO->Draw("HIST");
    py_e_out->SetLineColor(kBlack);
    py_e_out->Draw("HIST same");

    
    
    
    
    diffP->cd(3);
    pz_e_outNO->SetLineColor(kRed);
    pz_e_outNO->Draw("HIST");
    pz_e_out->SetLineColor(kBlack);
    pz_e_out->Draw("HIST same");

    
    
    
    diffP->cd(4);
    px_mu_outNO->SetLineColor(40);
    px_mu_outNO->Draw("HIST");
    px_mu_out->SetLineColor(46);
    px_mu_out->Draw("HIST same");

    
    
    diffP->cd(5);
    py_mu_outNO->SetLineColor(40);
    py_mu_outNO->Draw("HIST");
    py_mu_out->SetLineColor(46);
    py_mu_out->Draw("HIST same");

    
    
    diffP->cd(6);
    pz_mu_outNO->SetLineColor(40);
    pz_mu_outNO->Draw("HIST");
    pz_mu_out->SetLineColor(46);
    pz_mu_out->Draw("HIST same");

    
    
    diffP->SaveAs("LO-NLO.png");
    
        
    
}