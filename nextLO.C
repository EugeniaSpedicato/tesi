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
Double_t E_ECAL;

TH1F* Emuout=new TH1F("EnergyMU", "Energy Mu out", 75,0.2,150); 
TH1F* Emuout_E =new TH1F("EnergyMU_E", "Energy Mu out", 75,0.2,150); 
TH1F* Emuout1=new TH1F("EnergyMU1", "Energy Mu out Tar 1", 75,0.2,150); 
TH1F* Emuout_E1 =new TH1F("EnergyMU_E1", "Energy Mu out Tar 1", 75,0.2,150); 
TH1F* Emuout2=new TH1F("EnergyMU2", "Energy Mu out Tar 2", 75,0.2,150); 
TH1F* Emuout_E2 =new TH1F("EnergyMU_E2", "Energy Mu out Tar 2", 75,0.2,150); 
    
TH1F* Eelout=new TH1F("EnergyEL", "Energy El out", 70,0.2,140); 
TH1F* Eelout_E =new TH1F("EnergyEL_E", "Energy El out", 70,0.2,140);
TH1F* Eelout1=new TH1F("EnergyEL1", "Energy El out Tar 1", 70,0.2,140); 
TH1F* Eelout_E1 =new TH1F("EnergyEL_E1", "Energy El out Tar 1", 70,0.2,140);
TH1F* Eelout2=new TH1F("EnergyEL1", "Energy El out Tar 1", 70,0.2,140); 
TH1F* Eelout_E2 =new TH1F("EnergyEL_E1", "Energy El out Tar 1", 70,0.2,140);

TH1F* Ephout=new TH1F("EnergyPH", "Energy Ph out", 75,0.2,150); 

TH1F* Ephout1=new TH1F("EnergyPH1", "Energy Ph out Tar 1", 75,0.2,150); 
 
TH1F* Ephout2=new TH1F("EnergyPH2", "Energy Ph out Tar 2", 75,0.2,150); 


    
TH1F* thmu=new TH1F("thetaMU", "Muon Polar Angle", 180,0,0.005);  
TH1F* the=new TH1F("thetaEL", "Electron Polar Angle", 180,0,0.1);

TH1F* thXZmu=new TH1F("thetaXZ", "theta XZ plane mu", 150,-0.002,0.002);
TH1F* thYZmu=new TH1F("thetaYZ", "theta YZ plane mu", 150,-0.002,0.002);
    
TH1F* thXZe=new TH1F("thetaXZ", "theta XZ plane e", 150,-0.1,0.1);
TH1F* thYZe=new TH1F("thetaYZ", "theta YZ plane e", 150,-0.1,0.1);
    
TH1F* thXZmu1=new TH1F("thetaXZ1", "theta XZ plane mu1", 150,-0.002,0.002);
TH1F* thYZmu1=new TH1F("thetaYZ1", "theta YZ plane mu1", 150,-0.002,0.002);
    
TH1F* thXZe1=new TH1F("thetaXZ1", "theta XZ plane e1", 150,-0.1,0.1);
TH1F* thYZe1=new TH1F("thetaYZ1", "theta YZ plane e1", 150,-0.1,0.1);
    
TH1F* thXZmu2=new TH1F("thetaXZ2", "theta XZ plane mu2", 150,-0.002,0.002);
TH1F* thYZmu2=new TH1F("thetaYZ2", "theta YZ plane mu2", 150,-0.002,0.002);
    
TH1F* thXZe2=new TH1F("thetaXZ2", "theta XZ plane e2", 150,-0.1,0.1);
TH1F* thYZe2=new TH1F("thetaYZ2", "theta YZ plane e2", 150,-0.1,0.1);
    
TH1F* tarONEthe=new TH1F("thetaEL1", "Electron Polar Angle Tar 1", 180,0,0.1);
TH1F* tarTWOthe=new TH1F("thetaEL2", "Electron Polar Angle Tar 2", 180,0,0.1);
TH1F* tarONEthmu=new TH1F("thetaMU1", "Muon Polar Angle Tar 1", 180,0,0.005);
TH1F* tarTWOthmu=new TH1F("thetaMU2", "Muon Polar Angle Tar 2", 180,0,0.005);
    
TH2F  *X_Y_mu  = new TH2F("CooMU" , " X  Vs. Y of the muon",100,-0.1,0.1,100,-0.1,0.1);
TH2F  *X_Y_e  = new TH2F("CooEL" , " X  Vs. Y of the electron",100,-0.1,0.1,100,-0.1,0.1);
TH2F  *X_Y_p  = new TH2F("CooPH" , " X  Vs. Y of the photon",100,-0.1,0.1,100,-0.1,0.1);
    
TH2F  *X_Y_mu1  = new TH2F("CooMU1" , " X  Vs. Y of the muon TAR 1",100,-0.1,0.1,100,-0.1,0.1);
TH2F  *X_Y_e1  = new TH2F("CooEL1" , " X  Vs. Y of the electron TAR 1",100,-0.1,0.1,100,-0.1,0.1);
TH2F  *X_Y_p1  = new TH2F("CooPH1" , " X  Vs. Y of the photon TAR 1",100,-0.1,0.1,100,-0.1,0.1);
    
TH2F  *X_Y_mu2  = new TH2F("CooMU2" , " X  Vs. Y of the muon TAR 2",100,-0.1,0.1,100,-0.1,0.1);
TH2F  *X_Y_e2  = new TH2F("CooEL2" , " X  Vs. Y of the electron TAR 2",100,-0.1,0.1,100,-0.1,0.1);
TH2F  *X_Y_p2  = new TH2F("CooPH2" , " X  Vs. Y of the photon TAR 2",100,-0.1,0.1,100,-0.1,0.1);
 
TH2F  *Th_E_el  = new TH2F("ThEel" , " Theta el Vs. E_ECAL",180,0,0.1,70,0.2,140);
TH2F  *Th_E_el1  = new TH2F("ThEel1" , " Theta el Vs. E_ECAL TAR 1",180,0,0.1,70,0.2,140);
TH2F  *Th_E_el2  = new TH2F("ThEel2" , " Theta el Vs. E_ECAL TAR 2",180,0,0.1,70,0.2,140);
    
     if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
    
   
    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
       
    if (photon_coox!=-1)
    {  E_ECAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {  E_ECAL=detKinBeamRot_Ee;}

      if (E_ECAL>1){
        if (abs(detKinBeamRot_cooXmu) < 0.07125 && abs(detKinBeamRot_cooYmu) < 0.07125)
        {
            
    Double_t Pmu=sqrt(detKinBeamRot_pXmu_out*detKinBeamRot_pXmu_out+detKinBeamRot_pYmu_out*detKinBeamRot_pYmu_out+detKinBeamRot_pZmu_out*detKinBeamRot_pZmu_out);
            
    Double_t Emu_out=sqrt(Pmu*Pmu+mmu*mmu);
    Emuout->Fill(Emu_out,wgt_full);  
    Emuout_E->Fill(detKinBeamRot_Emu,wgt_full);   
            
    thmu->Fill(detKinBeamRot_thmu*0.001,wgt_full);
    Double_t anglex_mu = atan2(detKinBeamRot_pXmu_out, detKinBeamRot_pZmu_out);
    Double_t angley_mu = atan2(detKinBeamRot_pYmu_out, detKinBeamRot_pZmu_out);
            
           
    thXZmu->Fill(anglex_mu,wgt_full);
    thYZmu->Fill(angley_mu,wgt_full);
            
               if (detKinBeamRot_tar==0)
       
              {
        Emuout1->Fill(Emu_out,wgt_full); 
        Emuout_E1->Fill(detKinBeamRot_Emu,wgt_full); 
        tarONEthmu->Fill(detKinBeamRot_thmu*0.001,wgt_full);
        X_Y_mu1->Fill(detKinBeamRot_cooXmu, detKinBeamRot_cooYmu,wgt_full);
        thXZmu1->Fill(anglex_mu,wgt_full);
        thYZmu1->Fill(angley_mu,wgt_full);
              }
            
                if (detKinBeamRot_tar==1)
       
               {
        Emuout2->Fill(Emu_out,wgt_full); 
        Emuout_E2->Fill(detKinBeamRot_Emu,wgt_full); 
        tarTWOthmu->Fill(detKinBeamRot_thmu*0.001,wgt_full);
        X_Y_mu2->Fill(detKinBeamRot_cooXmu, detKinBeamRot_cooYmu,wgt_full); 
        thXZmu2->Fill(anglex_mu,wgt_full);
        thYZmu2->Fill(angley_mu,wgt_full);
               }
    X_Y_mu ->Fill(detKinBeamRot_cooXmu, detKinBeamRot_cooYmu,wgt_full);
        
        }
        
       
       
       if (abs(detKinBeamRot_cooXe) < 0.07125 && abs(detKinBeamRot_cooYe) < 0.07125)
        {
       Double_t Pe=sqrt(detKinBeamRot_pXe_out*detKinBeamRot_pXe_out+detKinBeamRot_pYe_out*detKinBeamRot_pYe_out+detKinBeamRot_pZe_out*detKinBeamRot_pZe_out);
   
       Double_t Ee_out=sqrt(Pe*Pe+me*me);
       Eelout->Fill(Ee_out,wgt_full);
        Eelout_E->Fill(detKinBeamRot_Ee,wgt_full);
        
       the->Fill(detKinBeamRot_the*0.001,wgt_full);  
        Double_t anglex_e = atan2(detKinBeamRot_pXe_out, detKinBeamRot_pZe_out);
        Double_t angley_e = atan2(detKinBeamRot_pYe_out, detKinBeamRot_pZe_out);
        thXZe->Fill(anglex_e,wgt_full);
        thYZe->Fill(angley_e,wgt_full);            
         
        Th_E_el->Fill(detKinBeamRot_the*0.001,E_ECAL,wgt_full);
           
           if (detKinBeamRot_tar==0)
              {
        Eelout1->Fill(Ee_out,wgt_full); 
        Eelout_E1->Fill(detKinBeamRot_Ee,wgt_full); 
        tarONEthe->Fill(detKinBeamRot_the*0.001,wgt_full);
        X_Y_e1->Fill(detKinBeamRot_cooXe, detKinBeamRot_cooYe,wgt_full);
        thXZe1->Fill(anglex_e,wgt_full);
        thYZe1->Fill(angley_e,wgt_full);
        Th_E_el1->Fill(detKinBeamRot_the*0.001,E_ECAL,wgt_full);
              }
            
            if (detKinBeamRot_tar==1)
            {
        Eelout2->Fill(Ee_out,wgt_full); 
        Eelout_E2->Fill(detKinBeamRot_Ee,wgt_full); 
        tarTWOthe->Fill(detKinBeamRot_the*0.001,wgt_full);
        X_Y_e2->Fill(detKinBeamRot_cooXe, detKinBeamRot_cooYe,wgt_full);  
        thXZe2->Fill(anglex_e,wgt_full);
        thYZe2->Fill(angley_e,wgt_full);
        Th_E_el2->Fill(detKinBeamRot_the*0.001,E_ECAL,wgt_full);    
             } 
           
            X_Y_e ->Fill(detKinBeamRot_cooXe, detKinBeamRot_cooYe,wgt_full);
           
        }
        
       
       
       if(abs(photon_coox) < 0.07125 && abs(photon_cooy) < 0.07125 && photon_energy>0.2)
       {
        Ephout->Fill(photon_energy,wgt_full);            
            
        if (detKinBeamRot_tar==0)
        {X_Y_p1->Fill(photon_coox, photon_cooy,wgt_full);
        Ephout1->Fill(photon_energy,wgt_full); }

            
           if (detKinBeamRot_tar==1)
           {X_Y_p2->Fill(photon_coox, photon_cooy,wgt_full);
           Ephout2->Fill(photon_energy,wgt_full);}                  
           
        X_Y_p->Fill(photon_coox, photon_cooy,wgt_full);
      }       
   }
   }
    
    
    /*TCanvas * e= new TCanvas("e","e",1500,1000,3500,2000);
    e->Divide(2,3);
    e->cd(1);
    Emuout->GetXaxis()->SetTitle("E [GeV]");
    Emuout->SetLineColor(kRed);
    Emuout->SetLineWidth(2);
    Emuout->Draw("HIST");
    gPad->SetLogy();
    e->cd(2);
    Emuout_E->GetXaxis()->SetTitle("E [GeV]");
    Emuout_E->SetLineColor(kRed);
    Emuout_E->SetLineWidth(2);
    Emuout_E->Draw("HIST");
    gPad->SetLogy();
    e->cd(3);
    Emuout1->GetXaxis()->SetTitle("E [GeV]");
    Emuout1->SetLineColor(8);
    Emuout1->SetLineWidth(2);
    Emuout1->Draw("HIST");
    gPad->SetLogy();
    e->cd(4);
    Emuout_E1->GetXaxis()->SetTitle("E [GeV]");
    Emuout_E1->SetLineColor(8);
    Emuout_E1->SetLineWidth(2);
    Emuout_E1->Draw("HIST");
    gPad->SetLogy();
    e->cd(5);
    Emuout2->GetXaxis()->SetTitle("E [GeV]");
    Emuout2->SetLineColor(kBlack);
    Emuout2->SetLineWidth(2);
    Emuout2->Draw("HIST");
    gPad->SetLogy();
    e->cd(6);
    Emuout_E2->GetXaxis()->SetTitle("E [GeV]");
    Emuout_E2->SetLineColor(kBlack);
    Emuout_E2->SetLineWidth(2);
    Emuout_E2->Draw("HIST");
    gPad->SetLogy();
    
    e->SaveAs("energyMU.png");
    
    TCanvas * ee= new TCanvas("e","e",1500,1000,3500,2000);
    ee->Divide(2,3);
    ee->cd(1);
    Eelout->GetXaxis()->SetTitle("E [GeV]");
    Eelout->SetLineWidth(2);
    Eelout->Draw("HIST");
    gPad->SetLogy();
    ee->cd(2);
    Eelout_E->GetXaxis()->SetTitle("E [GeV]");
    Eelout_E->SetLineWidth(2);
    Eelout_E->Draw("HIST");
    gPad->SetLogy();
    ee->cd(3);
    Eelout1->GetXaxis()->SetTitle("E [GeV]");
    Eelout1->SetLineWidth(2);
    Eelout1->SetLineColor(8);
    Eelout1->Draw("HIST");
    gPad->SetLogy();
    ee->cd(4);
    Eelout_E1->GetXaxis()->SetTitle("E [GeV]");
    Eelout_E1->SetLineWidth(2);
    Eelout_E1->SetLineColor(8);
    Eelout_E1->Draw("HIST");
    gPad->SetLogy();
    ee->cd(5);
    Eelout2->GetXaxis()->SetTitle("E [GeV]");
    Eelout2->SetLineWidth(2);
    Eelout2->SetLineColor(kBlack);
    Eelout2->Draw("HIST");
    gPad->SetLogy();
    ee->cd(6);
    Eelout_E2->GetXaxis()->SetTitle("E [GeV]");
    Eelout_E2->SetLineWidth(2);
    Eelout_E2->SetLineColor(kBlack);
    Eelout_E2->Draw("HIST");
    gPad->SetLogy();
    
    ee->SaveAs("energyEL.png");

    TCanvas * eew= new TCanvas("ew","ew",1500,1000,3500,2000);
    eew->Divide(2,3);
    eew->cd(1);
    Ephout->GetXaxis()->SetTitle("E [GeV]");
    Ephout->SetLineColor(9);
    Ephout->SetLineWidth(2);
    Ephout->Draw("HIST");
    gPad->SetLogy();
    eew->cd(3);
    Ephout1->GetXaxis()->SetTitle("E [GeV]");
    Ephout1->SetLineColor(9);
    Ephout1->SetLineWidth(2);
    Ephout1->SetLineColor(8);
    Ephout1->Draw("HIST");
    gPad->SetLogy();
    eew->cd(5);
    Ephout2->GetXaxis()->SetTitle("E [GeV]");
    Ephout2->SetLineColor(9);
    Ephout2->SetLineWidth(2);
    Ephout2->SetLineColor(kBlack);
    Ephout2->Draw("HIST");
    gPad->SetLogy();    
    
    eew->SaveAs("energyPh.png");*/


    
    TCanvas * t= new TCanvas("t","t",1500,1000,3500,2000);
    t->Divide(1,2);
    t->cd(1);
    thmu->SetLineColor(46);
    thmu->SetLineWidth(2);
    tarONEthmu->SetLineColor(8);
    tarONEthmu->SetLineWidth(2);
    tarTWOthmu->SetLineColor(kBlack);
    tarTWOthmu->SetLineWidth(2);
    thmu->Draw("HIST");
    tarONEthmu->Draw("HIST SAME");
    tarTWOthmu->Draw("HIST SAME");
    gPad->SetLogy();    
    thmu->GetXaxis()->SetTitle("Polar Theta [rad]");

    t->cd(2);
    tarONEthe->SetLineColor(8);
    tarONEthe->SetLineWidth(2);
    tarTWOthe->SetLineColor(kBlack);
    tarTWOthe->SetLineWidth(2);
    the->SetLineWidth(2);
    the->Draw("HIST");
    tarONEthe->Draw("HIST SAME");
    tarTWOthe->Draw("HIST SAME");
    the->GetXaxis()->SetTitle("Polar Theta [rad]");

    t->SaveAs("THpolar.png");    


TCanvas * theC= new TCanvas("tar","tar",1500,1000,3500,2000);
    theC->Divide(2,2);
    theC->cd(1);
    thXZmu->SetLineColor(46);
    thXZmu->SetLineWidth(2);
    thXZmu->Draw("HIST");
    thXZmu1->SetLineColor(8);
    thXZmu1->SetLineWidth(2);
    thXZmu1->Draw("HIST SAME");
    thXZmu2->SetLineColor(kBlack);
    thXZmu2->SetLineWidth(2);
    thXZmu2->Draw("HIST SAME");
    thXZmu->GetXaxis()->SetTitle("Theta XZ [rad]");
    theC->cd(2);
    thXZe->SetLineWidth(2);
    thXZe->Draw("HIST");
    thXZe1->SetLineColor(8);
    thXZe1->SetLineWidth(2);
    thXZe1->Draw("HIST SAME");
    thXZe2->SetLineColor(kBlack);
    thXZe2->SetLineWidth(2);
    thXZe2->Draw("HIST SAME");
    thXZe->GetXaxis()->SetTitle("Theta XZ [rad]");
    theC->cd(3);
    thYZmu->SetLineColor(46);
    thYZmu->SetLineWidth(2);
    thYZmu->Draw("HIST");
    thYZmu1->SetLineColor(8);
    thYZmu1->SetLineWidth(2);
    thYZmu1->Draw("HIST SAME");
    thYZmu2->SetLineColor(kBlack);
    thYZmu2->SetLineWidth(2);
    thYZmu2->Draw("HIST SAME");
    thYZmu->GetXaxis()->SetTitle("Theta YZ [rad]");
    theC->cd(4);
    thYZe->SetLineWidth(2);
    thYZe->Draw("HIST");
    thYZe1->SetLineColor(8);
    thYZe1->SetLineWidth(2);
    thYZe1->Draw("HIST SAME");
    thYZe2->SetLineColor(kBlack);
    thYZe2->SetLineWidth(2);
    thYZe2->Draw("HIST SAME");
    thYZe->GetXaxis()->SetTitle("Theta YZ [rad]");
    

  theC->SaveAs("ThXZYZ.png");


/*Int_t nx13_cut = X_Y_mu->GetNbinsX();
Int_t ny13_cut = X_Y_mu->GetNbinsY();
for (Int_t i=1; i<nx13_cut+1; i++) {
for (Int_t j=1; j<ny13_cut+1; j++) {
    if (X_Y_mu->GetBinContent(i,j)<1) X_Y_mu->SetBinContent(i,j,0);}} 
Int_t nx14_cut = X_Y_e->GetNbinsX();
Int_t ny14_cut = X_Y_e->GetNbinsY();
for (Int_t i=1; i<nx14_cut+1; i++) {
for (Int_t j=1; j<ny14_cut+1; j++) {
    if (X_Y_e->GetBinContent(i,j)<1) X_Y_e->SetBinContent(i,j,0);}} 
Int_t nx15_cut = X_Y_p->GetNbinsX();
Int_t ny15_cut = X_Y_p->GetNbinsY();
for (Int_t i=1; i<nx15_cut+1; i++) {
for (Int_t j=1; j<ny15_cut+1; j++) {
if (X_Y_p->GetBinContent(i,j)<1) X_Y_p->SetBinContent(i,j,0);}} 

Int_t nx13 = X_Y_mu1->GetNbinsX();
Int_t ny13 = X_Y_mu1->GetNbinsY();
for (Int_t i=1; i<nx13+1; i++) {
for (Int_t j=1; j<ny13+1; j++) {
    if (X_Y_mu1->GetBinContent(i,j)<1) X_Y_mu1->SetBinContent(i,j,0);}} 
Int_t nx14 = X_Y_e1->GetNbinsX();
Int_t ny14 = X_Y_e1->GetNbinsY();
for (Int_t i=1; i<nx14+1; i++) {
for (Int_t j=1; j<ny14+1; j++) {
    if (X_Y_e1->GetBinContent(i,j)<1) X_Y_e1->SetBinContent(i,j,0);}} 
Int_t nx15 = X_Y_p1->GetNbinsX();
Int_t ny15 = X_Y_p1->GetNbinsY();
for (Int_t i=1; i<nx15+1; i++) {
for (Int_t j=1; j<ny15+1; j++) {
if (X_Y_p1->GetBinContent(i,j)<1) X_Y_p1->SetBinContent(i,j,0);}} 


Int_t nx1 = X_Y_mu2->GetNbinsX();
Int_t ny1 = X_Y_mu2->GetNbinsY();
for (Int_t i=1; i<nx1+1; i++) {
for (Int_t j=1; j<ny1+1; j++) {
    if (X_Y_mu2->GetBinContent(i,j)<1) X_Y_mu2->SetBinContent(i,j,0);}} 
Int_t nx4 = X_Y_e2->GetNbinsX();
Int_t ny4 = X_Y_e2->GetNbinsY();
for (Int_t i=1; i<nx4+1; i++) {
for (Int_t j=1; j<ny4+1; j++) {
    if (X_Y_e2->GetBinContent(i,j)<1) X_Y_e2->SetBinContent(i,j,0);}} 
Int_t nx5 = X_Y_p2->GetNbinsX();
Int_t ny5 = X_Y_p2->GetNbinsY();
for (Int_t i=1; i<nx5+1; i++) {
for (Int_t j=1; j<ny5+1; j++) {
if (X_Y_p2->GetBinContent(i,j)<1) X_Y_p2->SetBinContent(i,j,0);}} 
    
    
    TCanvas * duedmu= new TCanvas("duedmu","duedmu",1000,100,2500,2000);
    duedmu->Divide(3,3);
    duedmu->cd(1);
    gStyle->SetPalette(kLake);
    TColor::InvertPalette(); 
    X_Y_mu->Draw("COLZ");
    X_Y_mu->GetXaxis()->SetTitle("x [m]");
    X_Y_mu->GetYaxis()->SetTitle("y [m]");
    duedmu->cd(2);
    X_Y_mu1->Draw("COLZ");
    X_Y_mu1->GetXaxis()->SetTitle("x [m]");
    X_Y_mu1->GetYaxis()->SetTitle("y [m]");
    duedmu->cd(3);
    X_Y_mu2->Draw("COLZ");
    X_Y_mu2->GetXaxis()->SetTitle("x [m]");
    X_Y_mu2->GetYaxis()->SetTitle("y [m]");
    
    duedmu->cd(4);
    X_Y_e->Draw("COLZ");
    X_Y_e->GetXaxis()->SetTitle("x [m]");
    X_Y_e->GetYaxis()->SetTitle("y [m]");
    duedmu->cd(5);
    X_Y_e1->Draw("COLZ");
    X_Y_e1->GetXaxis()->SetTitle("x [m]");
    X_Y_e1->GetYaxis()->SetTitle("y [m]");
    duedmu->cd(6);
    X_Y_e2->Draw("COLZ");
    X_Y_e2->GetXaxis()->SetTitle("x [m]");
    X_Y_e2->GetYaxis()->SetTitle("y [m]");
    
    duedmu->cd(7);
    X_Y_p->Draw("COLZ");
    X_Y_p->GetXaxis()->SetTitle("x [m]");
    X_Y_p->GetYaxis()->SetTitle("y [m]");
    duedmu->cd(8);
    X_Y_p1->Draw("COLZ");
    X_Y_p1->GetXaxis()->SetTitle("x [m]");
    X_Y_p1->GetYaxis()->SetTitle("y [m]");
    duedmu->cd(9);
    X_Y_p2->Draw("COLZ");
    X_Y_p2->GetXaxis()->SetTitle("x [m]");
    X_Y_p2->GetYaxis()->SetTitle("y [m]");
    
  duedmu->SaveAs("dued.png");  
    
    
Int_t nx7 = Th_E_el->GetNbinsX();
Int_t ny7 = Th_E_el->GetNbinsY();
for (Int_t i=1; i<nx7+1; i++) {
for (Int_t j=1; j<ny7+1; j++) {
    if (Th_E_el->GetBinContent(i,j)<1) Th_E_el->SetBinContent(i,j,0);}}     
    
TCanvas * th_en= new TCanvas("th_en","th_en",1000,100,2500,2000);      
Th_E_el->Draw("COLZ");  
th_en->SaveAs("theta-energy-electron.png");
    
Int_t nx9 = Th_E_el1->GetNbinsX();
Int_t ny9 = Th_E_el1->GetNbinsY();
for (Int_t i=1; i<nx9+1; i++) {
for (Int_t j=1; j<ny9+1; j++) {
    if (Th_E_el1->GetBinContent(i,j)<1) Th_E_el1->SetBinContent(i,j,0);}}  

TCanvas * th_en1= new TCanvas("th_en1","th_en1",1000,100,2500,2000);    
Th_E_el1->Draw("COLZ");  
th_en1->SaveAs("theta-energy-electron1.png"); 
    
Int_t nx11 = Th_E_el2->GetNbinsX();
Int_t ny11 = Th_E_el2->GetNbinsY();
for (Int_t i=1; i<nx11+1; i++) {
for (Int_t j=1; j<ny11+1; j++) {
    if (Th_E_el2->GetBinContent(i,j)<1) Th_E_el2->SetBinContent(i,j,0);}} 
    
TCanvas * th_en2= new TCanvas("th_en0","th_en0",1000,100,2500,2000);   
Th_E_el2->Draw("COLZ");  
th_en2->SaveAs("theta-energy-electron2.png");  */

    
      }