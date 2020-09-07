#define atree_cxx
#include "next.h"
#include <TH2.h>
#include <TH1.h>

#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{
    TH1::SetDefaultSumw2();
    
TH1F* px_mu=new TH1F("h1", "pX_in muon", 190,-0.3,0.3);
TH1F* py_mu=new TH1F("h2", "pY_in muon", 190,-0.3,0.3);
TH1F* pz_mu=new TH1F("h3", "pZ_in muon", 190,50,180);
    
TH1F* px_mu_out=new TH1F("h1a", "pX_out muon", 150,-0.3,0.3);
TH1F* py_mu_out=new TH1F("h2a", "pY_out muon", 150,-0.3,0.3);
TH1F* pz_mu_out=new TH1F("h3a", "pZ_out muon", 150,0,180);
TH1F* px_mu_outLO=new TH1F("h1a", "pX_out muon LO", 150,-0.3,0.3);
TH1F* py_mu_outLO=new TH1F("h2a", "pY_out muon LO", 150,-0.3,0.3);
TH1F* pz_mu_outLO=new TH1F("h3a", "pZ_out muon LO", 150,0,180);
TH1F* Emuin=new TH1F("h1aN", "Energy in mu", 150,0,160);
TH1F* Emuout=new TH1F("h1aN", "Energy out mu", 150,0,160);
TH1F* Ep=new TH1F("h1aN", "Energy out p", 150,0,1);
    
    
TH1F* px_e_out=new TH1F("h1b", "pX_out electron", 150,-0.3,0.3);
TH1F* py_e_out=new TH1F("h2b", "pY_out electron", 150,-0.3,0.3);
TH1F* pz_e_out=new TH1F("h3b", "pZ_out electron", 150,0,5);
TH1F* px_e_outLO=new TH1F("h1b", "pX_out electron LO", 150,-0.3,0.3);
TH1F* py_e_outLO=new TH1F("h2b", "pY_out electron LO", 150,-0.3,0.3);
TH1F* pz_e_outLO=new TH1F("h3b", "pZ_out electron LO", 150,0,5);
TH1F* Eein=new TH1F("h2aN", "Energ in e", 150,0,10);
TH1F* Eeout=new TH1F("h2aN", "Energy out e", 150,0,10);
    
TH1F* thmu=new TH1F("h3bNj", "theta", 180,0,0.002);  
TH1F* the=new TH1F("h3bNj", "theta", 180,0,0.1);

TH1F* thXZmu=new TH1F("a", "theta XZ mu", 150,-0.002,0.002);
TH1F* thYZmu=new TH1F("c", "theta YZ mu", 150,-0.002,0.002);
    
TH1F* thXZe=new TH1F("v", "theta XZ e", 150,-0.1,0.1);
TH1F* thYZe=new TH1F("b", "theta YZ e", 150,-0.1,0.1);
    
TH1F* tarONEXmu=new TH1F("h1a", "Coo X mu tar1  ", 140,-0.15,0.15);
TH1F* tarONEYmu=new TH1F("h2a", "Coo Y mu tar1 ", 140,-0.15,0.15);
TH1F* tarONEXe=new TH1F("h1ea", "Coo X e tar1 ", 140,-0.15,0.15);
TH1F* tarONEYe=new TH1F("h2ea", "Coo Y e tar1 ", 140,-0.15,0.15);
TH1F* tarONEXp=new TH1F("h1ea", "Coo X ph tar2 ", 100,-0.5,0.5);
TH1F* tarONEYp=new TH1F("h2ea", "Coo Y ph tar2 ", 100,-0.5,0.5);
    
TH1F* tarTWOXmu=new TH1F("h1a", "Coo X mu tar2  ", 140,-0.15,0.15);
TH1F* tarTWOYmu=new TH1F("h2a", "Coo Y mu tar2 ", 140,-0.15,0.15);
TH1F* tarTWOXe=new TH1F("h1ea", "Coo X e tar2 ", 140,-0.15,0.15);
TH1F* tarTWOYe=new TH1F("h2ea", "Coo Y e tar2 ", 140,-0.15,0.15);
TH1F* tarTWOXp=new TH1F("h1ea", "Coo X ph tar2 ", 100,-0.5,0.5);
TH1F* tarTWOYp=new TH1F("h2ea", "Coo Y ph tar2 ", 100,-0.5,0.5);
    
TH1F* dx=new TH1F("h2ea", "diff coo X e and mu", 140,-0.1,0.1);
TH1F* dy=new TH1F("h2ea", "diff coo Y e and mu ", 140,-0.1,0.1);
TH1F* dxmp=new TH1F("h2ea", "diff coo X photon and mu", 140,-0.1,0.1);
TH1F* dymp=new TH1F("h2ea", "diff coo Y photon and mu ", 140,-0.1,0.1);
TH1F* dxep=new TH1F("h2ea", "diff coo X e and photon", 140,-0.1,0.1);
TH1F* dyep=new TH1F("h2ea", "diff coo Y e and photon ", 140,-0.1,0.1);
    
TH2F  *X_Y_mu  = new TH2F("h2d" , " X  Vs. y of the muon",140,-0.5,-0.5,140,-0.5,0.5);
TH2F  *X_Y_e  = new TH2F("h2da" , " X  Vs. y of the electron",140,-0.5,-0.5,140,-0.5,0.5);
TH2F  *X_Y_p  = new TH2F("h2da" , " X  Vs. y of the photon",140,-0.5,-0.5,140,-0.5,0.5);
    
    
     if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
        px_mu->Fill(detKinBeamRot_pXmu,wgt_full);
       py_mu->Fill(detKinBeamRot_pYmu,wgt_full);
       pz_mu->Fill(detKinBeamRot_pZmu,wgt_full);
       px_mu_out->Fill(detKinBeamRot_pXmu_out,wgt_full);
       py_mu_out->Fill(detKinBeamRot_pYmu_out,wgt_full);
       pz_mu_out->Fill(detKinBeamRot_pZmu_out,wgt_full);
       px_e_out->Fill(detKinBeamRot_pXe_out,wgt_LO);
       py_e_out->Fill(detKinBeamRot_pYe_out,wgt_LO);
       pz_e_out->Fill(detKinBeamRot_pZe_out,wgt_LO);  
       px_mu_outLO->Fill(detKinBeamRot_pXmu_out,wgt_full);
       py_mu_outLO->Fill(detKinBeamRot_pYmu_out,wgt_full);
       pz_mu_outLO->Fill(detKinBeamRot_pZmu_out,wgt_full);
       px_e_outLO->Fill(detKinBeamRot_pXe_out,wgt_LO);
       py_e_outLO->Fill(detKinBeamRot_pYe_out,wgt_LO);
       pz_e_outLO->Fill(detKinBeamRot_pZe_out,wgt_LO);  
       
Double_t mmu= 105.6583745 *0.001;
Double_t me= 0.5109989461 *0.001;
       
    Double_t Pmu=sqrt(detKinBeamRot_pXmu_out*detKinBeamRot_pXmu_out+detKinBeamRot_pYmu_out*detKinBeamRot_pYmu_out+detKinBeamRot_pZmu_out*detKinBeamRot_pZmu_out);
    Double_t Pe=sqrt(detKinBeamRot_pXe_out*detKinBeamRot_pXe_out+detKinBeamRot_pYe_out*detKinBeamRot_pYe_out+detKinBeamRot_pZe_out*detKinBeamRot_pZe_out);
       
    Double_t Emu_out=sqrt(Pmu*Pmu+mmu*mmu);
    Double_t Ee_out=sqrt(Pe*Pe+me*me);
       
    Emuout->Fill(Emu_out,wgt_full);   
    Eeout->Fill(Ee_out,wgt_full);
    Ep->Fill(photon_energy,wgt_full);

       
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
       
           if (detKinBeamRot_tar==0)
       
           { 
        tarONEXmu->Fill(detKinBeamRot_cooXmu,wgt_full);
         tarONEYmu->Fill(detKinBeamRot_cooYmu,wgt_full);
         tarONEXe->Fill(detKinBeamRot_cooXe,wgt_full);
         tarONEYe->Fill(detKinBeamRot_cooYe,wgt_full);
             if (photon_coox!=-1 && photon_cooy!=-1)
        {tarONEXp->Fill(photon_coox,wgt_full);
         tarONEYp->Fill(photon_cooy,wgt_full);}
      }
       
       if (detKinBeamRot_tar==1)
       {  
         tarTWOXmu->Fill(detKinBeamRot_cooXmu,wgt_full);
         tarTWOYmu->Fill(detKinBeamRot_cooYmu,wgt_full);
         tarTWOXe->Fill(detKinBeamRot_cooXe,wgt_full);
         tarTWOYe->Fill(detKinBeamRot_cooYe,wgt_full);
        if (photon_coox!=-1)
        {tarTWOXp->Fill(photon_coox,wgt_full);
         tarTWOYp->Fill(photon_cooy,wgt_full);}
       }

    X_Y_mu ->Fill(detKinBeamRot_cooXmu, detKinBeamRot_cooYmu,wgt_full);
     X_Y_e ->Fill(detKinBeamRot_cooXe, detKinBeamRot_cooYe,wgt_full);
     if (photon_coox!=-1 && photon_cooy!=-1)
        {X_Y_p ->Fill(photon_coox, photon_cooy,wgt_full);}
       
       
    Double_t Dx = detKinBeamRot_cooXe-detKinBeamRot_cooXmu;
    Double_t Dy = detKinBeamRot_cooYe-detKinBeamRot_cooXmu;
    Double_t Dxep = detKinBeamRot_cooXe-photon_coox;
    Double_t Dyep = detKinBeamRot_cooYe-photon_cooy;    
    Double_t Dxmp = photon_coox-detKinBeamRot_cooXmu;
    Double_t Dymp = photon_cooy-detKinBeamRot_cooXmu;

    dx->Fill(Dx,wgt_full);
    dy->Fill(Dy,wgt_full);
    dxep->Fill(Dxep,wgt_full);
    dyep->Fill(Dyep,wgt_full);    
    dxmp->Fill(Dxmp,wgt_full);
    dymp->Fill(Dymp,wgt_full);
    
      
      }
      
    TCanvas * p= new TCanvas("p","p",400,10,1500,1000);
    p->Divide(2,3);
    p->cd(1);
    px_mu_out->SetLineColor(46);
    px_mu_out->Draw("HIST");
    px_mu_outLO->SetLineColor(kBlack);
    px_mu_outLO->Draw("HIST same");
    p->cd(2);    
    py_mu_out->SetLineColor(46);
    py_mu_out->Draw("HIST");
    py_mu_outLO->SetLineColor(kBlack);
    py_mu_outLO->Draw("HIST same");
    p->cd(3);
    pz_mu_out->SetLineColor(46);
    pz_mu_out->Draw("HIST");  
    pz_mu_outLO->SetLineColor(kBlack);
    pz_mu_outLO->Draw("HIST same");  
    p->cd(4);
    px_e_out->Draw("HIST");
    px_e_outLO->SetLineColor(kBlack);
    px_e_outLO->Draw("HIST same");  
    p->cd(5);
    py_e_out->Draw("HIST"); 
    py_e_outLO->SetLineColor(kBlack);
    py_e_outLO->Draw("HIST same");
    p->cd(6);
    pz_e_out->Draw("HIST");
    pz_e_outLO->SetLineColor(kBlack);
    pz_e_outLO->Draw("HIST same");
    p->SaveAs("Pemu.png");
    
    TCanvas * e= new TCanvas("e","e",400,10,1500,1000);
    e->Divide(3,2);
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
    e->cd(5);
    Ep->SetLineColor(30);
    Ep->Draw("HIST");
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
    tar->Divide(3,4);
    tar->cd(1);
    tarONEXmu->SetLineColor(46);
    tarONEXmu->Draw("HIST");
    tar->cd(2);    
    tarONEYmu->SetLineColor(46);
    tarONEYmu->Draw("HIST");
    tar->cd(3);
    tarONEXe->Draw("HIST");
    tar->cd(4);
    tarONEYe->Draw("HIST");
    tar->cd(5);
    tarTWOXmu->SetLineColor(46);
    tarTWOXmu->Draw("HIST");
    tar->cd(6);
    tarTWOYmu->SetLineColor(46);
    tarTWOYmu->Draw("HIST");
    tar->cd(7);
    tarTWOXe->Draw("HIST");
    tar->cd(8);
    tarTWOYe->Draw("HIST");
    tar->cd(9);
    tarONEXp->SetLineColor(30);
    tarONEXp->Draw("HIST");
    tar->cd(10);
    tarONEYp->SetLineColor(30);
    tarONEYp->Draw("HIST");
    tar->cd(11);
    tarTWOXp->SetLineColor(30);
    tarTWOXp->Draw("HIST");
    tar->cd(12);
    tarTWOYp->SetLineColor(30);
    tarTWOYp->Draw("HIST");
    
    
  tar->SaveAs("tar.png");
    
    TCanvas * dued= new TCanvas("dued","dued",1000,100,2500,2000);
    dued->Divide(1,3);
    dued->cd(1);
    X_Y_mu->SetMarkerColor(46);
    X_Y_mu->Draw("HIST");
    dued->cd(2);
    X_Y_e->SetMarkerColor(30);
    X_Y_e->Draw("HIST");
    dued->cd(3);
    X_Y_p->Draw("HIST");
  dued->SaveAs("dued.png");
    
        TCanvas * Pin= new TCanvas("Pin","Pin",400,10,1500,1000);
    Pin->Divide(1,3);
    Pin->cd(1);
    px_mu->Draw("HIST");
    Pin->cd(2);
    py_mu->Draw("HIST");
    Pin->cd(3);
    pz_mu->Draw("HIST");
    
   Pin->SaveAs("p_in.png");
    
    TCanvas * d= new TCanvas("d","d",1000,100,2500,2000);
    d->Divide(3,2);
    d->cd(1);
    dx->Draw("HIST");
    d->cd(2);
    dy->Draw("HIST");
    d->cd(3);
    dxep->Draw("HIST");
    d->cd(4);
    dyep->Draw("HIST");
    d->cd(5);
    dxmp->Draw("HIST");
    d->cd(6);
    dymp->Draw("HIST");
  d->SaveAs("diffCoo.png");
    

    
      }