#define atree_cxx
#include "next.h"
#include <TH2.h>
#include <TH1.h>

#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{
    TH1::SetDefaultSumw2();
    
TH1F* px_mu_out=new TH1F("h1a", "pX_out muon", 150,-10,10);
TH1F* py_mu_out=new TH1F("h2a", "pY_out muon", 150,-10,10);
TH1F* pz_mu_out=new TH1F("h3a", "pZ_out muon", 150,0,180);
TH1F* Emuin=new TH1F("h1aN", "Energy in", 150,0,160);
TH1F* Emuout=new TH1F("h1aN", "Energy out", 150,0,160);
    
TH1F* px_e_out=new TH1F("h1b", "pX_out electron", 150,-0.3,0.3);
TH1F* py_e_out=new TH1F("h2b", "pY_out electron", 150,-0.3,0.3);
TH1F* pz_e_out=new TH1F("h3b", "pZ_out electron", 150,0,5);
TH1F* Eein=new TH1F("h2aN", "Energ iny", 150,0,10);
TH1F* Eeout=new TH1F("h2aN", "Energy out", 150,0,10);
    
TH1F* thmu=new TH1F("h3bNj", "theta", 150,0,0.3);  
TH1F* the=new TH1F("h3bNj", "theta", 150,0,0.5);

TH1F* thXZmu=new TH1F("a", "theta XZ mu", 150,-0.2,0.2);
TH1F* thYZmu=new TH1F("c", "theta YZ mu", 150,-0.2,0.2);
    
TH1F* thXZe=new TH1F("v", "theta XZ e", 150,-0.5,0.5);
TH1F* thYZe=new TH1F("b", "theta YZ e", 150,-0.5,0.5);
    
TH1F* tarONEXmu=new TH1F("h1a", "Coo X mu tar1  ", 140,-0.4,0.4);
TH1F* tarONEYmu=new TH1F("h2a", "Coo Y mu tar1 ", 140,-0.4,0.4);
TH1F* tarONEXe=new TH1F("h1ea", "Coo X e tar1 ", 140,-0.8,0.8);
TH1F* tarONEYe=new TH1F("h2ea", "Coo Y e tar1 ", 140,-0.8,0.8);
    
TH1F* tarTWOXmu=new TH1F("h1a", "Coo X mu tar2  ", 140,-0.4,0.4);
TH1F* tarTWOYmu=new TH1F("h2a", "Coo Y mu tar2 ", 140,-0.4,0.4);
TH1F* tarTWOXe=new TH1F("h1ea", "Coo X e tar2 ", 140,-0.8,0.8);
TH1F* tarTWOYe=new TH1F("h2ea", "Coo Y e tar2 ", 140,-0.8,0.8);
    
TH2F  *X_Y_mu  = new TH2F("h2d" , " X  Vs. y of the muon",140,-0.5,-0.5,140,-0.5,0.5);
TH2F  *X_Y_e  = new TH2F("h2da" , " X  Vs. y of the electron",140,-0.5,-0.5,140,-0.5,0.5);
    
    
     if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
       if (detKinBeamRot_tar==0){
       px_mu_out->Fill(detKinBeamRot_pXmu_out,wgt_full);
       py_mu_out->Fill(detKinBeamRot_pYmu_out,wgt_full);
       pz_mu_out->Fill(detKinBeamRot_pZmu_out,wgt_full);
       px_e_out->Fill(detKinBeamRot_pXe_out,wgt_full);
       py_e_out->Fill(detKinBeamRot_pYe_out,wgt_full);
       pz_e_out->Fill(detKinBeamRot_pZe_out,wgt_full);       
       
Double_t mmu= 105.6583745 *0.001;
Double_t me= 0.5109989461 *0.001;
       
    Double_t Pmu=sqrt(detKinBeamRot_pXmu_out*detKinBeamRot_pXmu_out+detKinBeamRot_pYmu_out*detKinBeamRot_pYmu_out+detKinBeamRot_pZmu_out*detKinBeamRot_pZmu_out);
    Double_t Pe=sqrt(detKinBeamRot_pXe_out*detKinBeamRot_pXe_out+detKinBeamRot_pYe_out*detKinBeamRot_pYe_out+detKinBeamRot_pZe_out*detKinBeamRot_pZe_out);
       
    Double_t Emu_out=sqrt(Pmu*Pmu+mmu*mmu);
    Double_t Ee_out=sqrt(Pe*Pe+me*me);
       
    Emuout->Fill(Emu_out,wgt_full);   
    Eeout->Fill(Ee_out,wgt_full);
       
    Emuin->Fill(detKinBeamRot_Emu,wgt_full);   
    Eein->Fill(detKinBeamRot_Ee,wgt_full);
       
    thmu->Fill(detKinBeamRot_thmu*0.001,wgt_full);
    the->Fill(detKinBeamRot_the*0.001,wgt_full);  
        
    Double_t anglex_mu = atan2(detKinBeamRot_pXmu_out, detKinBeamRot_pZmu_out);
    Double_t angley_mu = atan2(detKinBeamRot_pYmu_out, detKinBeamRot_pZmu_out);       
    Double_t anglex_e = atan2(detKinBeamRot_pXe_out, detKinBeamRot_pZe_out);
    Double_t angley_e = atan2(detKinBeamRot_pYe_out, detKinBeamRot_pZe_out);
       
    thXZmu->Fill(anglex_mu,wgt_full);
    thYZmu->Fill(angley_mu,wgt_full);
    thXZe->Fill(anglex_e,wgt_full);
    thYZe->Fill(angley_e,wgt_full); 
       
      /*     if (detKinBeamRot_tar==0)
       {*/
         tarONEXmu->Fill(detKinBeamRot_cooXmu,wgt_full);
         tarONEYmu->Fill(detKinBeamRot_cooYmu,wgt_full);
         tarONEXe->Fill(detKinBeamRot_cooXe,wgt_full);
         tarONEYe->Fill(detKinBeamRot_cooYe,wgt_full);
      // }
       
      /* if (detKinBeamRot_tar==1)
       {*/
         tarTWOXmu->Fill(detKinBeamRot_cooXmu,wgt_full);
         tarTWOYmu->Fill(detKinBeamRot_cooYmu,wgt_full);
         tarTWOXe->Fill(detKinBeamRot_cooXe,wgt_full);
         tarTWOYe->Fill(detKinBeamRot_cooYe,wgt_full);
           
      // }
       
    X_Y_mu ->Fill(detKinBeamRot_cooXmu, detKinBeamRot_cooYmu,wgt_full);
     X_Y_e ->Fill(detKinBeamRot_cooXe, detKinBeamRot_cooYe,wgt_full);
      
      }}
      
    TCanvas * p= new TCanvas("p","p",400,10,1500,1000);
    p->Divide(2,3);
    p->cd(1);
    px_mu_out->SetLineColor(46);
    px_mu_out->Draw("HIST");
    p->cd(2);    
    py_mu_out->SetLineColor(46);
    py_mu_out->Draw("HIST");
    p->cd(3);
    pz_mu_out->SetLineColor(46);
    pz_mu_out->Draw("HIST");  
    p->cd(4);
    px_e_out->Draw("HIST");  
    p->cd(5);
    py_e_out->Draw("HIST");   
    p->cd(6);
    pz_e_out->Draw("HIST");
    p->SaveAs("Pemu.png");
    
    TCanvas * e= new TCanvas("e","e",400,10,1500,1000);
    e->Divide(2,2);
    e->cd(1);
    Emuin->SetLineColor(46);
    Emuin->Draw("HIST");
    e->cd(2);
    Emuout->SetLineColor(46);
    Emuout->Draw("HIST");
    e->cd(3);
    Eein->Draw("HIST");
    e->cd(4);
    Eeout->Draw("HIST");
    e->SaveAs("energy.png");

        TCanvas * t= new TCanvas("t","t",400,10,1500,1000);
    t->Divide(1,2);
    t->cd(1);
    thmu->SetLineColor(46);
    thmu->Draw("HIST");

    t->cd(2);
    the->Draw("HIST");

    t->SaveAs("THpolar.png");
    
    TCanvas * theC= new TCanvas("tar","tar",400,10,1500,1000);
    theC->Divide(2,2);
    theC->cd(1);
    thXZmu->SetLineColor(46);
    thXZmu->Draw("HIST");
    theC->cd(2);
    thXZe->Draw("HIST");
    theC->cd(3);
    thYZmu->SetLineColor(46);
    thYZmu->Draw("HIST");
    theC->cd(4);
    thYZe->Draw("HIST");
    
    
  theC->SaveAs("ThOriginal.png");
    
    
    
        TCanvas * tar= new TCanvas("tar","tar",400,10,1500,1000);
    tar->Divide(2,4);
    tar->cd(1);
    tarONEXmu->SetLineColor(46);
    tarONEXmu->Draw("HIST");
    tar->cd(2);
    tarONEXe->Draw("HIST");
    tar->cd(3);
    tarONEYmu->SetLineColor(46);
    tarONEYmu->Draw("HIST");
    tar->cd(4);
    tarONEYe->Draw("HIST");
    tar->cd(5);
    tarTWOXmu->SetLineColor(46);
    tarTWOXmu->Draw("HIST");
    tar->cd(6);
    tarTWOXe->Draw("HIST");
    tar->cd(7);
    tarTWOYmu->SetLineColor(46);
    tarTWOYmu->Draw("HIST");
    tar->cd(8);
    tarTWOYe->Draw("HIST");
    
  tar->SaveAs("tar.png");
    
    TCanvas * dued= new TCanvas("dued","dued",1000,100,2500,2000);
    dued->Divide(1,2);
    dued->cd(1);
    X_Y_mu->SetMarkerColor(46);
    X_Y_mu->Draw("HIST");
    dued->cd(2);
    X_Y_e->Draw("HIST");
  dued->SaveAs("dued.png");
    

    
      }