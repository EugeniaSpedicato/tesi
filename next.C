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
    
Double_t mmu= 105.6583745 *0.001;
Double_t me= 0.5109989461 *0.001;
Double_t Rm = 0.01959; //raggio di Moliere in metri
Double_t E_ECAL;

    
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
TH1F* Ep=new TH1F("h1aN", "Energy out p", 200,0,160);

    
    
TH1F* px_e_out=new TH1F("h1b", "pX_out electron", 150,-0.3,0.3);
TH1F* py_e_out=new TH1F("h2b", "pY_out electron", 150,-0.3,0.3);
TH1F* pz_e_out=new TH1F("h3b", "pZ_out electron", 150,0,5);
TH1F* px_e_outLO=new TH1F("h1b", "pX_out electron LO", 150,-0.3,0.3);
TH1F* py_e_outLO=new TH1F("h2b", "pY_out electron LO", 150,-0.3,0.3);
TH1F* pz_e_outLO=new TH1F("h3b", "pZ_out electron LO", 150,0,5);
TH1F* Eein=new TH1F("h2aN", "Energ in e", 200,0,160);
TH1F* Eeout=new TH1F("h2aN", "Energy out e", 200,0,160);
    
TH1F* thmu=new TH1F("h3bNj", "theta", 180,0,5);  
TH1F* the=new TH1F("h3bNj", "theta", 180,0,100);

TH1F* thXZmu=new TH1F("a", "theta XZ mu", 150,-0.002,0.002);
TH1F* thYZmu=new TH1F("c", "theta YZ mu", 150,-0.002,0.002);
    
TH1F* thXZe=new TH1F("v", "theta XZ e", 150,-0.1,0.1);
TH1F* thYZe=new TH1F("b", "theta YZ e", 150,-0.1,0.1);
    
TH1F* tarONEXmu=new TH1F("h1a", "Coo X mu tar1  ", 070,-0.15,0.15);
TH1F* tarONEYmu=new TH1F("h2a", "Coo Y mu tar1 ", 070,-0.15,0.15);
TH1F* tarONEXe=new TH1F("h1ea", "Coo X e tar1 ", 070,-0.15,0.15);
TH1F* tarONEYe=new TH1F("h2ea", "Coo Y e tar1 ", 070,-0.15,0.15);
TH1F* tarONEXp=new TH1F("h1ea", "Coo X ph tar2 ", 100,-0.5,0.5);
TH1F* tarONEYp=new TH1F("h2ea", "Coo Y ph tar2 ", 100,-0.5,0.5);
    
TH1F* tarTWOXmu=new TH1F("h1a", "Coo X mu tar2  ", 070,-0.15,0.15);
TH1F* tarTWOYmu=new TH1F("h2a", "Coo Y mu tar2 ", 070,-0.15,0.15);
TH1F* tarTWOXe=new TH1F("h1ea", "Coo X e tar2 ", 070,-0.15,0.15);
TH1F* tarTWOYe=new TH1F("h2ea", "Coo Y e tar2 ", 070,-0.15,0.15);
TH1F* tarTWOXp=new TH1F("h1ea", "Coo X ph tar2 ", 100,-0.5,0.5);
TH1F* tarTWOYp=new TH1F("h2ea", "Coo Y ph tar2 ", 100,-0.5,0.5);
    
TH1F* dx=new TH1F("h2ea", "diff coo X e and mu", 070,-0.1,0.1);
TH1F* dy=new TH1F("h2ea", "diff coo Y e and mu ", 070,-0.1,0.1);
TH1F* dxmp=new TH1F("h2ea", "diff coo X photon and mu", 070,-0.1,0.1);
TH1F* dymp=new TH1F("h2ea", "diff coo Y photon and mu ", 070,-0.1,0.1);
TH1F* dxep=new TH1F("h2ea", "diff coo X e and photon", 070,-0.1,0.1);
TH1F* dyep=new TH1F("h2ea", "diff coo Y e and photon ", 070,-0.1,0.1);
    
TH1F* DR=new TH1F("DR", "Distanza elettrone-fotone", 70,0,0.14);
TH1F* DR_cut=new TH1F("DR", "Distanza elettrone-fotone cut Ee>1Gev", 70,0,0.14);
TH1F* DRmu=new TH1F("DR", "Distanza elettrone-muone", 70,0,0.14);
TH1F* DRmu_cut=new TH1F("DR", "Distanza elettrone-muone cut Ee>1Gev", 70,0,0.14);


    
TH2F  *X_Y_mu  = new TH2F("h2d" , " X  Vs. y of the muon",70,-0.1,0.1,70,-0.1,0.1);
TH2F  *X_Y_e  = new TH2F("h2da" , " X  Vs. y of the electron",70,-0.1,0.1,70,-0.1,0.1);
TH2F  *X_Y_p  = new TH2F("h2da" , " X  Vs. y of the photon",70,-0.1,0.1,70,-0.1,0.1);
TH2F  *E_r  = new TH2F("h2da" , " X  Vs. y of the photon",70,-0.1,0.1,70,-0.1,0.1);

 
    
     if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
    
   
    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
 
       if (detKinBeamRot_tar==1){
       px_mu->Fill(detKinBeamRot_pXmu,wgt_full);
       py_mu->Fill(detKinBeamRot_pYmu,wgt_full);
       pz_mu->Fill(detKinBeamRot_pZmu,wgt_full);
       
           if (photon_coox!=-1)
    {  
     E_ECAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {  E_ECAL=detKinBeamRot_Ee;}
       
        if (abs(detKinBeamRot_cooXmu) <0.07125 && abs(detKinBeamRot_cooYmu) <0.07125)
        {
       px_mu_out->Fill(detKinBeamRot_pXmu_out,wgt_full);
       py_mu_out->Fill(detKinBeamRot_pYmu_out,wgt_full);
       pz_mu_out->Fill(detKinBeamRot_pZmu_out,wgt_full); 
       px_mu_outLO->Fill(detKinBeamRot_pXmu_out,wgt_LO);
       py_mu_outLO->Fill(detKinBeamRot_pYmu_out,wgt_LO);
       pz_mu_outLO->Fill(detKinBeamRot_pZmu_out,wgt_LO);
            
    Double_t Pmu=sqrt(detKinBeamRot_pXmu_out*detKinBeamRot_pXmu_out+detKinBeamRot_pYmu_out*detKinBeamRot_pYmu_out+detKinBeamRot_pZmu_out*detKinBeamRot_pZmu_out);
            
    Double_t Emu_out=sqrt(Pmu*Pmu+mmu*mmu);
    Emuout->Fill(Emu_out,wgt_full);  
    Emuin->Fill(detKinBeamRot_Emu,wgt_LO);   
            
    thmu->Fill(detKinBeamRot_thmu,wgt_full);
    Double_t anglex_mu = atan2(detKinBeamRot_pXmu_out, detKinBeamRot_pZmu_out);
    Double_t angley_mu = atan2(detKinBeamRot_pYmu_out, detKinBeamRot_pZmu_out);
            
           
    thXZmu->Fill(anglex_mu,wgt_full);
    thYZmu->Fill(angley_mu,wgt_full);
            
              if (detKinBeamRot_tar==0)
       
              {
         tarONEXmu->Fill(detKinBeamRot_cooXmu,wgt_full);
         tarONEYmu->Fill(detKinBeamRot_cooYmu,wgt_full);  
               }
            
               if (detKinBeamRot_tar==1)
       
              {
         tarTWOXmu->Fill(detKinBeamRot_cooXmu,wgt_full);
         tarTWOYmu->Fill(detKinBeamRot_cooYmu,wgt_full);    
              }
if(E_ECAL>1)
{X_Y_mu ->Fill(detKinBeamRot_cooXmu, detKinBeamRot_cooYmu,wgt_full);}
        
        }
        
       
       
       if (abs(detKinBeamRot_cooXe) <0.07125 && abs(detKinBeamRot_cooYe) <0.07125)
        {
       px_e_out->Fill(detKinBeamRot_pXe_out,wgt_full);
       py_e_out->Fill(detKinBeamRot_pYe_out,wgt_full);
       pz_e_out->Fill(detKinBeamRot_pZe_out,wgt_full);  
       px_e_outLO->Fill(detKinBeamRot_pXe_out,wgt_LO);
       py_e_outLO->Fill(detKinBeamRot_pYe_out,wgt_LO);
       pz_e_outLO->Fill(detKinBeamRot_pZe_out,wgt_LO);
      
       Double_t Pe=sqrt(detKinBeamRot_pXe_out*detKinBeamRot_pXe_out+detKinBeamRot_pYe_out*detKinBeamRot_pYe_out+detKinBeamRot_pZe_out*detKinBeamRot_pZe_out);
   
       Double_t Ee_out=sqrt(Pe*Pe+me*me);
       Eeout->Fill(Ee_out,wgt_full);
       Eein->Fill(detKinBeamRot_Ee,wgt_LO);
           
       the->Fill(detKinBeamRot_the,wgt_full);  
        Double_t anglex_e = atan2(detKinBeamRot_pXe_out, detKinBeamRot_pZe_out);
        Double_t angley_e = atan2(detKinBeamRot_pYe_out, detKinBeamRot_pZe_out);
        thXZe->Fill(anglex_e,wgt_full);
        thYZe->Fill(angley_e,wgt_full);            
               if (detKinBeamRot_tar==0)
       
               {
             tarONEXe->Fill(detKinBeamRot_cooXe,wgt_full);
             tarONEYe->Fill(detKinBeamRot_cooYe,wgt_full);               
               }
            
               if (detKinBeamRot_tar==1)
       
               {
             tarTWOXe->Fill(detKinBeamRot_cooXe,wgt_full);
             tarTWOYe->Fill(detKinBeamRot_cooYe,wgt_full);  
               } 
           
           if(E_ECAL>1)
           {X_Y_e ->Fill(detKinBeamRot_cooXe, detKinBeamRot_cooYe,wgt_full);}
           
        }
        
       
       
       if(abs(photon_coox) < 0.07125 && abs(photon_cooy) < 0.07125 && photon_energy>0.2)
       {
        Ep->Fill(photon_energy,wgt_full);            
               if (detKinBeamRot_tar==0)
       
               {
            tarONEXp->Fill(photon_coox,wgt_full);
         tarONEYp->Fill(photon_cooy,wgt_full);
               }
            
               if (detKinBeamRot_tar==1)
       
               {
         tarTWOXp->Fill(photon_coox,wgt_full);
         tarTWOYp->Fill(photon_cooy,wgt_full);
               }                  
        
                
        if(E_ECAL>1)   
        {X_Y_p ->Fill(photon_coox, photon_cooy,wgt_full);}
       }
       
       



       
       if (abs(detKinBeamRot_cooXmu) <0.07125 && abs(detKinBeamRot_cooYmu) <0.07125 && abs(detKinBeamRot_cooXe) <0.07125 && abs(detKinBeamRot_cooYe) <0.07125 && abs(photon_coox) < 0.07125 && abs(photon_cooy) < 0.07125)
       {Double_t Dx = detKinBeamRot_cooXe-detKinBeamRot_cooXmu;
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
    dymp->Fill(Dymp,wgt_full);}
       
       
       if (abs(detKinBeamRot_cooXe) <0.07125 && abs(detKinBeamRot_cooYe) <0.07125 && abs(photon_coox) < 0.07125 && abs(photon_cooy) < 0.07125 && E_ECAL>1 && photon_energy>0.2)
       {
           Double_t dist=sqrt((detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy));
           
           DR->Fill(dist,wgt_full);
           
           if (detKinBeamRot_Ee>10)
           { DR_cut->Fill(dist,wgt_full);}
       }
              if (abs(detKinBeamRot_cooXe) <0.07125 && abs(detKinBeamRot_cooYe) <0.07125 && abs(detKinBeamRot_cooXmu) < 0.07125 && abs(detKinBeamRot_cooYmu) < 0.07125 && E_ECAL>1)
       {
           Double_t dist=sqrt((detKinBeamRot_cooXe-detKinBeamRot_cooXmu)*(detKinBeamRot_cooXe-detKinBeamRot_cooXmu)+(detKinBeamRot_cooYe-detKinBeamRot_cooYmu)*(detKinBeamRot_cooYe-detKinBeamRot_cooYmu));
           
           DRmu->Fill(dist,wgt_full);

        if (detKinBeamRot_Ee>10)
           { DRmu_cut->Fill(dist,wgt_full);}
       }

       }
       
}

    

    /* TCanvas * p= new TCanvas("p","p",400,10,1500,1000);
    p->Divide(2,3);
    p->cd(1);
    px_mu_out->SetLineColor(46);
    px_mu_out->Draw("HIST");
    px_mu_outLO->SetLineColor(kBlack);
    px_mu_outLO->Draw("HIST same");
    px_mu_out->GetXaxis()->SetTitle("Px [GeV]");
    
    p->cd(2);    
    py_mu_out->SetLineColor(46);
    py_mu_out->Draw("HIST");
    py_mu_outLO->SetLineColor(kBlack);
    py_mu_outLO->Draw("HIST same");
    py_mu_out->GetXaxis()->SetTitle("Py [GeV]");
    p->cd(3);
    pz_mu_out->SetLineColor(46);
    pz_mu_out->Draw("HIST");  
    pz_mu_outLO->SetLineColor(kBlack);
    pz_mu_outLO->Draw("HIST same");
    pz_mu_out->GetXaxis()->SetTitle("Pz [GeV]");
    p->cd(4);
    px_e_out->Draw("HIST");
    px_e_outLO->SetLineColor(kBlack);
    px_e_outLO->Draw("HIST same"); 
    px_e_out->GetXaxis()->SetTitle("Px [GeV]");
    p->cd(5);
    py_e_out->Draw("HIST"); 
    py_e_outLO->SetLineColor(kBlack);
    py_e_outLO->Draw("HIST same");
    py_e_out->GetXaxis()->SetTitle("Py [GeV]");
    p->cd(6);
    pz_e_out->Draw("HIST");
    pz_e_outLO->SetLineColor(kBlack);
    pz_e_outLO->Draw("HIST same");
    pz_e_out->GetXaxis()->SetTitle("Pz [GeV]");
    p->SaveAs("Pemu.png");
    
    
    TCanvas * e= new TCanvas("e","e",400,10,1500,1000);
    e->Divide(3,2);
    e->cd(1);
    Emuin->SetLineColor(46);
    Emuin->GetXaxis()->SetTitle("E [GeV]");
    Emuin->Draw("HIST");
    Emuout->SetLineColor(kBlack);
    Emuout->Draw("HIST same");
    
    
    e->cd(2);
    Emuout->GetXaxis()->SetTitle("E [GeV]");
    Emuout->Draw("HIST");
    
    e->cd(3);
    Eein->Draw("HIST");
    Eein->GetXaxis()->SetTitle("E [GeV] (log scale)");
    gPad->SetLogx();
    Eeout->SetLineColor(kBlack);
    Eeout->Draw("HIST same");
    
    e->cd(4);
    Eeout->Draw("HIST");
    Eeout->GetXaxis()->SetTitle("E [GeV] (log scale)");
    gPad->SetLogx();
    
    e->cd(5);
    Ep->SetLineColor(30);
    Ep->Draw("HIST");
    Ep->GetXaxis()->SetTitle("E [GeV] (log scale)");
    gPad->SetLogx();
    
    e->SaveAs("energy.png");
  

    TCanvas * t= new TCanvas("t","t",400,10,1500,1000);
    t->Divide(1,2);
    t->cd(1);
    thmu->SetLineColor(46);
    thmu->Draw("HIST");
    thmu->GetXaxis()->SetTitle("Polar Theta [mrad]");

    t->cd(2);
    the->Draw("HIST");
    the->GetXaxis()->SetTitle("Polar Theta [mrad]");

    t->SaveAs("THpolar.png");
    
    TCanvas * theC= new TCanvas("tar","tar",400,10,1500,1000);
    theC->Divide(2,2);
    theC->cd(1);
    thXZmu->SetLineColor(46);
    thXZmu->Draw("HIST");
    thXZmu->GetXaxis()->SetTitle("Theta XZ [rad]");
    theC->cd(2);
    thXZe->Draw("HIST");
    thXZe->GetXaxis()->SetTitle("Theta XZ [rad]");
    theC->cd(3);
    thYZmu->SetLineColor(46);
    thYZmu->Draw("HIST");
    thYZmu->GetXaxis()->SetTitle("Theta YZ [rad]");
    theC->cd(4);
    thYZe->Draw("HIST");
    thYZe->GetXaxis()->SetTitle("Theta YZ [rad]");
    
  theC->SaveAs("ThOriginal.png");
    
    
    
        TCanvas * tar= new TCanvas("tar","tar",400,10,1500,1000);
    tar->Divide(3,4);
    tar->cd(1);
    tarONEXmu->SetLineColor(46);
    tarONEXmu->Draw("HIST");
    tarONEXmu->GetXaxis()->SetTitle("x [m]");
    tar->cd(2);    
    tarONEYmu->SetLineColor(46);
    tarONEYmu->Draw("HIST");
    tarONEYmu->GetXaxis()->SetTitle("y [m]");
    tar->cd(3);
    tarONEXe->Draw("HIST");
    tarONEXe->GetXaxis()->SetTitle("x [m]");
    tar->cd(4);
    tarONEYe->Draw("HIST");
    tarONEYe->GetXaxis()->SetTitle("y [m]");
    tar->cd(5);
    tarTWOXmu->SetLineColor(46);
    tarTWOXmu->Draw("HIST");
    tarTWOXmu->GetXaxis()->SetTitle("x [m]");
    tar->cd(6);
    tarTWOYmu->SetLineColor(46);
    tarTWOYmu->Draw("HIST");
    tarTWOYmu->GetXaxis()->SetTitle("y [m]");
    tar->cd(7);
    tarTWOXe->Draw("HIST");
    tarTWOXe->GetXaxis()->SetTitle("x [m]");
    tar->cd(8);
    tarTWOYe->Draw("HIST");
    tarTWOYe->GetXaxis()->SetTitle("y [m]");
    tar->cd(9);
    tarONEXp->SetLineColor(30);
    tarONEXp->Draw("HIST");
    tarONEXp->GetXaxis()->SetTitle("x [m]");
    tar->cd(10);
    tarONEYp->SetLineColor(30);
    tarONEYp->Draw("HIST");
    tarONEYp->GetXaxis()->SetTitle("y [m]");
    tar->cd(11);
    tarTWOXp->SetLineColor(30);
    tarTWOXp->Draw("HIST");
    tarTWOXp->GetXaxis()->SetTitle("x [m]");
    tar->cd(12);
    tarTWOYp->SetLineColor(30);
    tarTWOYp->Draw("HIST");
    tarTWOYp->GetXaxis()->SetTitle("y [m]");
    
    
  tar->SaveAs("tar.png");*/
        
Int_t nx1 = X_Y_mu->GetNbinsX();
Int_t ny1 = X_Y_mu->GetNbinsY();
for (Int_t i=1; i<nx1+1; i++) {
for (Int_t j=1; j<ny1+1; j++) {
    if (X_Y_mu->GetBinContent(i,j)<1) X_Y_mu->SetBinContent(i,j,0);}} 
        
Int_t nx2 = X_Y_e->GetNbinsX();
Int_t ny2 = X_Y_e->GetNbinsY();
for (Int_t i=1; i<nx2+1; i++) {
for (Int_t j=1; j<ny2+1; j++) {
    if (X_Y_e->GetBinContent(i,j)<1) X_Y_e->SetBinContent(i,j,0);}} 
        
Int_t nx3 = X_Y_p->GetNbinsX();
Int_t ny3 = X_Y_p->GetNbinsY();
for (Int_t i=1; i<nx3+1; i++) {
for (Int_t j=1; j<ny3+1; j++) {
    if (X_Y_p->GetBinContent(i,j)<1) X_Y_p->SetBinContent(i,j,0);}} 
    
    
    
    TCanvas * duedmu= new TCanvas("duedmu","duedmu",1000,100,2500,2000);
    gStyle->SetPalette(kCherry);
    TColor::InvertPalette();
    X_Y_mu->SetMarkerColor(46);
    X_Y_mu->Draw("COLZ");
    X_Y_mu->GetXaxis()->SetTitle("x [m]");
    X_Y_mu->GetYaxis()->SetTitle("y [m]");
  duedmu->SaveAs("duedmu.png");
    
    TCanvas * duede= new TCanvas("duede","duede",1000,100,2500,2000);
    X_Y_e->SetMarkerColor(30);
    X_Y_e->Draw("COLZ");
    X_Y_e->GetXaxis()->SetTitle("x [m]");
    X_Y_e->GetYaxis()->SetTitle("y [m]");
  duede->SaveAs("duede.png");
    

    TCanvas * duedp= new TCanvas("duedp","duedp",1000,100,2500,2000);
    X_Y_p->Draw("COLZ");
    X_Y_p->GetXaxis()->SetTitle("x [m]");
    X_Y_p->GetYaxis()->SetTitle("y [m]");
  duedp->SaveAs("duedph.png");
  /*  
        TCanvas * Pin= new TCanvas("Pin","Pin",400,10,1500,1000);
    Pin->Divide(1,3);
    Pin->cd(1);
    px_mu->Draw("HIST");
    px_mu->GetXaxis()->SetTitle("Px [GeV]");
    Pin->cd(2);
    py_mu->Draw("HIST");
    py_mu->GetXaxis()->SetTitle("Py [GeV]");
    Pin->cd(3);
    pz_mu->Draw("HIST");
    pz_mu->GetXaxis()->SetTitle("Pz [GeV]");
    
   Pin->SaveAs("p_in.png");
    
    TCanvas * d= new TCanvas("d","d",1000,100,2500,2000);
    d->Divide(3,2);
    d->cd(1);
    dx->Draw("HIST");
    dx->GetXaxis()->SetTitle("Delta_x [m]");
    d->cd(2);
    dy->Draw("HIST");
    dy->GetXaxis()->SetTitle("Delta_y [m]");
    d->cd(3);
    dxep->Draw("HIST");
    dxep->GetXaxis()->SetTitle("Delta_x [m]");
    d->cd(4);
    dyep->Draw("HIST");
    dyep->GetXaxis()->SetTitle("Delta_y [m]");
    d->cd(5);
    dxmp->Draw("HIST");
    dxmp->GetXaxis()->SetTitle("Delta_x [m]");
    d->cd(6);
    dymp->Draw("HIST");
    dymp->GetXaxis()->SetTitle("Delta_y [m]");
  d->SaveAs("diffCoo.png");
    
TCanvas * DRR= new TCanvas("d","d",1000,100,2500,2000);
DRR->Divide(2,1);
DRR->cd(1);
DRmu->Draw("HIST");
DRmu_cut->SetLineColor(kRed);
DRmu_cut->Draw("HIST same");
DRR->cd(2);
DR->Draw("HIST");
DR_cut->SetLineColor(kRed);
DR_cut->Draw("HIST same");
    
DRR->SaveAs("DRphe.png");*/
      }




void EMShower::compute() {
  double t = 0.;
  double dt = 0.;
  if (!stepsCalculated)
    prepareSteps();

  // Prepare the grids in EcalHitMaker
  // theGrid->setInnerAndOuterDepth(innerDepth,outerDepth);
  float pstot = 0.;
  float ps2tot = 0.;
  float ps1tot = 0.;
  bool status = false;
  //  double E1 = 0.;  // Energy layer 1
  //  double E2 = 0.;  // Energy layer 2
  //  double n1 = 0.;  // #mips layer 1
  //  double n2 = 0.;  // #mips layer 2
  //  double E9 = 0.;  // Energy ECAL

  // Loop over all segments for the longitudinal development
  double totECalc = 0;

  for (unsigned iStep = 0; iStep < nSteps; ++iStep) {
    // The length of the shower in this segment
    dt = steps[iStep].second;

    //    std::cout << " Detector " << steps[iStep].first << " t " << t << " " << dt << std::endl;

    // The elapsed length
    t += dt;

    // In what detector are we ?
    unsigned detector = steps[iStep].first;

    bool presh1 = detector == 0;
    bool presh2 = detector == 1;
    bool ecal = detector == 2;
    bool hcal = detector == 3;
    bool vfcal = detector == 4;
    bool gap = detector == 5;

    // Temporary. Will be removed
    if (theHCAL == nullptr)
      hcal = false;

    // Keep only ECAL for now
    if (vfcal)
      continue;

    // Nothing to do in the gap
    if (gap)
      continue;

    //    cout << " t = " << t << endl;
    // Build the grid of crystals at this ECAL depth
    // Actually, it might be useful to check if this grid is empty or not.
    // If it is empty (because no crystal at this depth), it is of no use
    // (and time consuming) to generate the spots

    // middle of the step
    double tt = t - 0.5 * dt;

    double realTotalEnergy = 0.;
    for (unsigned int i = 0; i < nPart; ++i) {
      realTotalEnergy += depositedEnergy[iStep][i] * E[i];
    }

    //    std::cout << " Step " << tt << std::endl;
    //    std::cout << "ecal " << ecal << " hcal "  << hcal <<std::endl;

    // If the amount of energy is greater than 1 MeV, make a new grid
    // otherwise put in the previous one.
    bool usePreviousGrid = (realTotalEnergy < 0.001);

    // If the amount of energy is greater than 1 MeV, make a new grid
    // otherwise put in the previous one.

    // If less than 1 kEV. Just skip
    if (iStep > 2 && realTotalEnergy < 0.000001)
      continue;

    if (ecal && !usePreviousGrid) {
      status = theGrid->getPads(meanDepth[iStep]);
    }
    if (hcal) {
      status = theHcalHitMaker->setDepth(tt);
    }
    if ((ecal || hcal) && !status)
      continue;

    bool detailedShowerTail = false;
    // check if a detailed treatment of the rear leakage should be applied
    if (ecal && !usePreviousGrid) {
      detailedShowerTail = (t - dt > theGrid->getX0back());
    }

    // The particles of the shower are processed in parallel
    for (unsigned int i = 0; i < nPart; ++i) {
      //      double Edepo=deposit(t,a[i],b[i],dt);

      //  integration of the shower profile between t-dt and t
      double dE = (!hcal) ? depositedEnergy[iStep][i] : 1. - deposit(a[i], b[i], t - dt);

      // no need to do the full machinery if there is ~nothing to distribute)
      if (dE * E[i] < 0.000001)
        continue;

      if (ecal && !theECAL->isHom()) {
        double mean = dE * E[i];
        double sigma = theECAL->resE() * sqrt(mean);

        /*
	  double meanLn = log(mean);
	  double kLn = sigma/mean+1;
	  double sigmaLn = log(kLn);
	*/

        double dE0 = dE;

        //	  std::cout << "dE before shoot = " << dE << std::endl;
        dE = random->gaussShoot(mean, sigma) / E[i];

        //	  myGammaGenerator->setParameters(aSam,bSam,0);
        //	  dE = myGammaGenerator->shoot()/E[i];
        //	  std::cout << "dE shooted = " << dE << " E[i] = " << E[i] << std::endl;
        if (dE * E[i] < 0.000001)
          continue;
        photos[i] = photos[i] * dE / dE0;
      }

      /*
      if (ecal && !theParam->ecalProperties()->isHom()){
	double cSquare = TMath::Power(theParam->ecalProperties()->resE(),2);
	double aSam = dE/cSquare;
	double bSam = 1./cSquare;
	//	dE = dE*gam(bSam*dE, aSam)/tgamma(aSam);
      }
      */

      totECalc += dE;

      // The number of energy spots (or mips)
      double nS = 0;

      // ECAL case : Account for photostatistics and long'al non-uniformity
      if (ecal) {
        //	double aSam = E[i]*dE*one_over_resoSquare;
        //	double bSam = one_over_resoSquare;

        dE = random->poissonShoot(dE * photos[i]) / photos[i];
        double z0 = random->gaussShoot(0., 1.);
        dE *= 1. + z0 * theECAL->lightCollectionUniformity();

        // Expected spot number
        nS = (theNumberOfSpots[i] * gam(bSpot[i] * tt, aSpot[i]) * bSpot[i] * dt / tgamma(aSpot[i]));

        // Preshower : Expected number of mips + fluctuation
      } else if (hcal) {
        nS = (theNumberOfSpots[i] * gam(bSpot[i] * tt, aSpot[i]) * bSpot[i] * dt / tgamma(aSpot[i])) *
             theHCAL->spotFraction();
        double nSo = nS;
        nS = random->poissonShoot(nS);
        // 'Quick and dirty' fix (but this line should be better removed):
        if (nSo > 0. && nS / nSo < 10.)
          dE *= nS / nSo;

        //	if(true)
        //	  {
        //	    std::cout << " theHCAL->spotFraction = " <<theHCAL->spotFraction() <<std::endl;
        //	    std::cout << " nSpot Ecal : " << nSo/theHCAL->spotFraction() << " Final " << nS << std::endl;
        //	  }
      } else if (presh1) {
        nS = random->poissonShoot(dE * E[i] * theLayer1->mipsPerGeV());
        //	std::cout << " dE *E[i] (1)" << dE*E[i] << " " << dE*E[i]*theLayer1->mipsPerGeV() << "  "<< nS << std::endl;
        pstot += dE * E[i];
        ps1tot += dE * E[i];
        dE = nS / (E[i] * theLayer1->mipsPerGeV());

        //        E1 += dE*E[i];
        //	n1 += nS;
        //	if (presh2) { E2 += SpotEnergy; ++n2; }

      } else if (presh2) {
        nS = random->poissonShoot(dE * E[i] * theLayer2->mipsPerGeV());
        //	std::cout << " dE *E[i] (2) " << dE*E[i] << " " << dE*E[i]*theLayer2->mipsPerGeV() << "  "<< nS << std::endl;
        pstot += dE * E[i];
        ps2tot += dE * E[i];
        dE = nS / (E[i] * theLayer2->mipsPerGeV());

        //        E2 += dE*E[i];
        //	n2 += nS;
      }

      if (detailedShowerTail)
        myGammaGenerator->setParameters(floor(a[i] + 0.5), b[i], t - dt);

      //    myHistos->fill("h100",t,dE);

      // The lateral development parameters

      // Energy of the spots
      double eSpot = (nS > 0.) ? dE / nS : 0.;
      double SpotEnergy = eSpot * E[i];

      if (hasPreshower && (presh1 || presh2))
        thePreshower->setSpotEnergy(0.00009);
      if (hcal) {
        SpotEnergy *= theHCAL->hOverPi();
        theHcalHitMaker->setSpotEnergy(SpotEnergy);
      }
      // Poissonian fluctuations for the number of spots
      //    int nSpot = random->poissonShoot(nS);
      int nSpot = (int)(nS + 0.5);

      // Fig. 11 (right) *** Does not match.
      //    myHistos->fill("h101",t,(double)nSpot/theNumberOfSpots);

      //double taui = t/T;
      double taui = tt / Ti[i];
      double proba = theParam->p(taui, E[i]);
      double theRC = theParam->rC(taui, E[i]);
      double theRT = theParam->rT(taui, E[i]);

      // Fig. 10
      //    myHistos->fill("h300",taui,theRC);
      //    myHistos->fill("h301",taui,theRT);
      //    myHistos->fill("h302",taui,proba);

      double dSpotsCore = random->gaussShoot(proba * nSpot, std::sqrt(proba * (1. - proba) * nSpot));

      if (dSpotsCore < 0)
        dSpotsCore = 0;

      unsigned nSpots_core = (unsigned)(dSpotsCore + 0.5);
      unsigned nSpots_tail = ((unsigned)nSpot > nSpots_core) ? nSpot - nSpots_core : 0;

      for (unsigned icomp = 0; icomp < 2; ++icomp) {
        double theR = (icomp == 0) ? theRC : theRT;
        unsigned ncompspots = (icomp == 0) ? nSpots_core : nSpots_tail;

        RadialInterval radInterval(theR, ncompspots, SpotEnergy, random);
        if (ecal) {
          if (icomp == 0) {
            setIntervals(icomp, radInterval);
          } else {
            setIntervals(icomp, radInterval);
          }
        } else {
          radInterval.addInterval(100., 1.);  // 100% of the spots
        }

        radInterval.compute();
        // irad = 0 : central circle; irad=1 : outside

        unsigned nrad = radInterval.nIntervals();

        for (unsigned irad = 0; irad < nrad; ++irad) {
          double spote = radInterval.getSpotEnergy(irad);
          if (ecal)
            theGrid->setSpotEnergy(spote);
          if (hcal)
            theHcalHitMaker->setSpotEnergy(spote);
          unsigned nradspots = radInterval.getNumberOfSpots(irad);
          double umin = radInterval.getUmin(irad);
          double umax = radInterval.getUmax(irad);
          // Go for the lateral development
          //	       std::cout << "Couche " << iStep << " irad = " << irad << " Ene = " << E[i] << " eSpot = " << eSpot << " spote = " << spote << " nSpot = " << nS << std::endl;

          for (unsigned ispot = 0; ispot < nradspots; ++ispot) {
            double z3 = random->flatShoot(umin, umax);
            double ri = theR * std::sqrt(z3 / (1. - z3));

            // Generate phi
            double phi = 2. * M_PI * random->flatShoot();

            // Add the hit in the crystal
            //	if( ecal ) theGrid->addHit(ri*theECAL->moliereRadius(),phi);
            // Now the *moliereRadius is done in EcalHitMaker
            if (ecal) {
              if (detailedShowerTail) {
                //			   std::cout << "About to call addHitDepth " << std::endl;
                double depth;
                do {
                  depth = myGammaGenerator->shoot(random);
                } while (depth > t);
                theGrid->addHitDepth(ri, phi, depth);
                //			   std::cout << " Done " << std::endl;
              } else
                theGrid->addHit(ri, phi);
            } else if (hasPreshower && presh1)
              thePreshower->addHit(ri, phi, 1);
            else if (hasPreshower && presh2)
              thePreshower->addHit(ri, phi, 2);
            else if (hcal) {
              //		       std::cout << " About to add a spot in the HCAL" << status << std::endl;
              theHcalHitMaker->addHit(ri, phi);
              //		       std::cout << " Added a spot in the HCAL" << status << std::endl;
            }
            //	if (ecal) E9 += SpotEnergy;
            //	if (presh1) { E1 += SpotEnergy; ++n1; }
            //	if (presh2) { E2 += SpotEnergy; ++n2; }

            Etot[i] += spote;
          }
        }
      }
      //      std::cout << " Done with the step " << std::endl;
      // The shower!
      //myHistos->fill("h500",theSpot.z(),theSpot.perp());
    }
    //    std::cout << " nPart " << nPart << std::endl;
  }
  //  std::cout << " Finshed the step loop " << std::endl;
  //  myHistos->fill("h500",E1+0.7*E2,E9);
  //  myHistos->fill("h501",n1+0.7*n2,E9);
  //  myHistos->fill("h400",n1);
  //  myHistos->fill("h401",n2);
  //  myHistos->fill("h402",E9+E1+0.7*E2);
  //  if(!standalone)theGrid->printGrid();
  double Etotal = 0.;
  for (unsigned i = 0; i < nPart; ++i) {
    //      myHistos->fill("h10",Etot[i]);
    Etotal += Etot[i];
  }

  //  std::cout << "Etotal = " << Etotal << " nPart = "<< nPart << std::endl;
  //  std::cout << "totECalc = " << totECalc << std::endl;

  //  myHistos->fill("h20",Etotal);
  //  if(thePreshower)
  //    std::cout << " PS " << thePreshower->layer1Calibrated() << " " << thePreshower->layer2Calibrated() << " " << thePreshower->totalCalibrated() << " " << ps1tot << " " <<ps2tot << " " << pstot << std::endl;
}

