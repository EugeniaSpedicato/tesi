#define atree_cxx
#include "atreecoo.h"
#include <TH2.h>
#include <TH3.h>

#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{TH1F* coox_mu=new TH1F("h1", "Coo X mu", 140,-0.1,0.1);
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
 
 
    TCanvas * dued= new TCanvas("dued","dued",400,10,600,400);
  X_Y_mu->SetMarkerColor(kBlack);
    X_Y_mu->Draw();
  X_Y_e->SetMarkerColor(kRed);
    X_Y_e->Draw("same");
 
 /*THStack *hs = new THStack("hs","Stacked histograms");
     hs->Add(X_Y_mu);
     X_Y_e->SetFillColor(kRed);
     hs->Add(X_Y_e);
     hs->Draw();*/
 
}