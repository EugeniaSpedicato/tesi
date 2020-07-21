#define atree_cxx
#include "atree.h"
#include <TH2.h>
#include <TH1.h>

#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{
    TH1::SetDefaultSumw2();
    
    
    
TH1F* px_e_out=new TH1F("h1bq", "pX_out NO DIV", 150,-0.3,0.3);
TH1F* py_e_out=new TH1F("h2bq", "pY_out NO DIV", 150,-0.3,0.3);
TH1F* pz_e_out=new TH1F("h3bq", "pZ_out NO DIV", 150,0,5);
TH1F* px_mu_out=new TH1F("h1aq", "pX_out NO DIV", 150,-0.3,0.3);
TH1F* py_mu_out=new TH1F("h2aq", "pY_out NO DIV", 150,-0.3,0.3);
TH1F* pz_mu_out=new TH1F("h3aq", "pZ_out NO DIV", 150,0,180);
    
TH1F* GKpx_e_out=new TH1F("h1", "pX_out NO DIV", 150,-0.3,0.3);
TH1F* GKpy_e_out=new TH1F("h2", "pY_out NO DIV", 150,-0.3,0.3);
TH1F* GKpz_e_out=new TH1F("h3", "pZ_out NO DIV", 150,0,5);
TH1F* GKpx_mu_out=new TH1F("h1a", "pX_out NO DIV", 150,-0.3,0.3);
TH1F* GKpy_mu_out=new TH1F("h2a", "pY_out NO DIV", 150,-0.3,0.3);
TH1F* GKpz_mu_out=new TH1F("h3a", "pZ_out NO DIV", 150,0,180);
    
TH1F* thXZmu=new TH1F("a", "theta XZ mu", 150,-0.005,0.005);
TH1F* thYZmu=new TH1F("c", "theta YZ mu", 150,-0.005,0.005);
    
TH1F* thXZe=new TH1F("v", "theta XZ e", 150,-0.05,0.05);
TH1F* thYZe=new TH1F("b", "theta YZ e", 150,-0.05,0.05);
    
     if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
       
       
       px_mu_out->Fill(detKin_pXmu_out,wgt_full);
       py_mu_out->Fill(detKin_pYmu_out,wgt_full);
       pz_mu_out->Fill(detKin_pZmu_out,wgt_full);
       px_e_out->Fill(detKin_pXe_out,wgt_full);
       py_e_out->Fill(detKin_pYe_out,wgt_full);
       pz_e_out->Fill(detKin_pZe_out,wgt_full);
       
       GKpx_mu_out->Fill(genKin_pXmu_out,wgt_full);
       GKpy_mu_out->Fill(genKin_pYmu_out,wgt_full);
       GKpz_mu_out->Fill(genKin_pZmu_out,wgt_full);
       GKpx_e_out->Fill(genKin_pXe_out,wgt_full);
       GKpy_e_out->Fill(genKin_pYe_out,wgt_full);
       GKpz_e_out->Fill(genKin_pZe_out,wgt_full);
       
       
    Double_t anglex_mu = atan2(detKin_pXmu_out, detKin_pZmu_out);
    Double_t angley_mu = atan2(detKin_pYmu_out, detKin_pZmu_out);       
    Double_t anglex_e = atan2(detKin_pXe_out, detKin_pZe_out);
    Double_t angley_e = atan2(detKin_pYe_out, detKin_pZe_out);
       
    thXZmu->Fill(anglex_mu,wgt_full);
    thYZmu->Fill(angley_mu,wgt_full);
    thXZe->Fill(anglex_e,wgt_full);
    thYZe->Fill(angley_e,wgt_full);
       
  
   }
   
    TCanvas * c= new TCanvas("c","c",400,10,1500,1000);  
    c->Divide(2,3);
    c->cd(1);
    px_mu_out->Draw("HIST");
    GKpx_mu_out->SetLineColor(30);
    GKpx_mu_out->Draw("HIST same");
    
    c->cd(2);
    py_mu_out->Draw("HIST");
    GKpy_mu_out->SetLineColor(30);
    GKpy_mu_out->Draw("HIST same");
    
    c->cd(3);
    pz_mu_out->Draw("HIST");
    GKpz_mu_out->SetLineColor(30);
    GKpz_mu_out->Draw("HIST same");
    
    c->cd(4);
    px_e_out->Draw("HIST");
    GKpx_e_out->SetLineColor(30);
    GKpx_e_out->Draw("HIST same");
    c->cd(5);
    py_e_out->Draw("HIST");
    GKpy_e_out->SetLineColor(30);
    GKpy_e_out->Draw("HIST same");
    c->cd(6);
    pz_e_out->Draw("HIST");
    GKpz_e_out->SetLineColor(30);
    GKpz_e_out->Draw("HIST same");
    
    c->SaveAs("Poriginal.png");
    
    TCanvas * theC= new TCanvas("tar","tar",400,10,1500,1000);
    theC->Divide(2,2);
    theC->cd(1);
    thXZmu->SetLineColor(30);
    thXZmu->Draw("HIST");
    theC->cd(2);
    thXZe->SetLineColor(38);
    thXZe->Draw("HIST");
    
    theC->cd(3);
    thYZmu->SetLineColor(30);
    thYZmu->Draw("HIST");
    theC->cd(4);
    thYZe->SetLineColor(38);
    thYZe->Draw("HIST");
    
    
  theC->SaveAs("ThOriginal.png");
    
    
    
}
