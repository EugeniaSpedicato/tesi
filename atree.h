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
    
    
TH1F* coox_mu=new TH1F("h1", "Coo X mu", 140,-0.1,0.1);
TH1F* cooy_mu=new TH1F("h2", "Coo Y mu", 140,-0.1,0.1);
 TH1F* coox_e=new TH1F("h1e", "Coo X e", 140,-0.1,0.1);
TH1F* cooy_e=new TH1F("h2e", "Coo Y e", 140,-0.1,0.1);
 
TH1F* diffX_mue=new TH1F("h", "DiffCoo X mu and e-", 140,-0.3,0.3);
TH1F* diffY_mue=new TH1F("h", "DiffCoo Y mu and e-", 140,-0.3,0.3);
 
TH2F  *X_Y_mu  = new TH2F("h2d" , " X  Vs. y of the muon",140,-0.3,-0.3,100,0,40);
TH2F  *X_Y_e  = new TH2F("h2da" , " X  Vs. y of the electron",140,-0.3,-0.3,100,0,40);
 
 

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
       
        coox_mu->Fill(detKinBeamRot_cooXmu);
       cooy_mu->Fill(detKinBeamRot_cooYmu);
       
       coox_e->Fill(detKinBeamRot_cooXe);
       cooy_e->Fill(detKinBeamRot_cooYe);
       
       Double_t diffX=detKinBeamRot_cooXe-detKinBeamRot_cooXmu;
       diffX_mue->Fill(diffX);
       Double_t diffY=detKinBeamRot_cooYe-detKinBeamRot_cooYmu;
       diffY_mue->Fill(diffY);
       
       
     X_Y_mu ->Fill(detKinBeamRot_cooXmu, detKinBeamRot_cooYmu);
     X_Y_e ->Fill(detKinBeamRot_cooXe, detKinBeamRot_cooYe);
       
       
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
    
        TCanvas * cooX= new TCanvas("cooX","cooX",400,10,600,400);
    cooX->Divide(2,2);
    cooX->cd(1);
    coox_mu->Draw();
    coox_e->SetLineColor(kRed);
    coox_e->Draw("same");

 
    cooX->cd(2);
    cooy_mu->Draw();
    cooy_e->SetLineColor(kRed);
    cooy_e->Draw("same");
 
    cooX->cd(3);
    diffX_mue->Draw();
      
    cooX->cd(4);
    diffY_mue->SetLineColor(kRed);
    diffY_mue->Draw();
  cooX->SaveAs("coo.png");
 
    TCanvas * dued= new TCanvas("dued","dued",400,10,600,400);
  X_Y_mu->SetMarkerColor(kBlack);
    X_Y_mu->Draw();
  X_Y_e->SetMarkerColor(kRed);
    X_Y_e->Draw("same");
  dued->SaveAs("duedcoo.png");
 
    
}
