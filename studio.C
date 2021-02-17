#define atree_cxx
#include "atree.h"
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>

#include <TGraph.h>
#include <TTree.h>
#include <cmath>
using namespace std;

#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{
    TH1::SetDefaultSumw2();


Double_t E_CAL=0.;
Double_t E9=0.;
 
TH1F* TheCUT=new TH1F("th", "th El core CUT", 120,0,30); 
TH1F* The1CUT=new TH1F("th", "th El TAR 1 core CUT", 120,0,30); 
TH1F* The2CUT=new TH1F("th", "th El TAR 2 core CUT", 120,0,30);


TH2F  *E3x31CUT  = new TH2F("ThEel1" , " Th_el Vs. E_E3x3 core NLO CUT",90,0,30,280,0,140);
TH2F  *E3x32CUT  = new TH2F("ThEel2" , " Th_el Vs. E_E3x3 core LO CUT",90,0,30,280,0,140);


TH1F* hist_E9_e=new TH1F("E9e", "E9 e- tot", 500,0,1);
TH1F* hist_E9_eLO=new TH1F("E9eLO", "E9 e- tot LO", 500,0,1);

    
if (fChain == 0) return;

Long64_t nentries = fChain->GetEntriesFast();



    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

        detKinBeamRot_cooXe=detKinBeamRot_cooXe*100; // cm
        detKinBeamRot_cooYe=detKinBeamRot_cooYe*100; // cm
        detKinBeamRot_x_in=detKinBeamRot_x_in*100; // cm coo muone entrante
        detKinBeamRot_y_in=detKinBeamRot_y_in*100; // cm coo muone entrante
        photon_coox=photon_coox*100; // cm
        photon_cooy=photon_cooy*100; // cm
               
       double r_mu=sqrt((detKinBeamRot_x_in*detKinBeamRot_x_in)+(detKinBeamRot_y_in*detKinBeamRot_y_in));
              
        double r=sqrt((detKinBeamRot_cooXe*detKinBeamRot_cooXe)+(detKinBeamRot_cooYe*detKinBeamRot_cooYe));
              
       
        E9=detKinBeamRot_E_1/detKinBeamRot_E_clus3x3;
       
if (detKinBeamRot_n_cell_e!=0 && detKinBeamRot_E_clus3x3>0.2)
{  
    if(abs(detKinBeamRot_cooXe)<4.275 && abs(detKinBeamRot_cooYe)<4.275) 
    {
        if(r_mu<1.7 && detKinBeamRot_def_angle_mu>0.2 && detKinBeamRot_E_clus3x3>1)
        {TheCUT->Fill(detKinBeamRot_def_angle_e,wgt_full);
         hist_E9_e->Fill(E9,wgt_full);
         hist_E9_eLO->Fill(E9,wgt_LO);}
    
            if (detKinBeamRot_E_clus3x3!=0) 
            {if(r_mu<1.7 && detKinBeamRot_def_angle_mu>0.2 && detKinBeamRot_E_clus3x3>1) E3x31CUT->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);
            } 
        
            if (detKinBeamRot_E_clus3x3!=0) 
            {if(r_mu<1.7 && detKinBeamRot_def_angle_mu>0.2 && detKinBeamRot_E_clus3x3>1) E3x32CUT->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_LO);
            } 
        
        if (detKinBeamRot_tar==0)
        {if(r_mu<1.7 && detKinBeamRot_def_angle_mu>0.2 && detKinBeamRot_E_clus3x3>1) The1CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);

          /*  if (detKinBeamRot_E_clus3x3!=0) 
            {if(r_mu<1.7 && detKinBeamRot_def_angle_mu>0.2 && detKinBeamRot_E_clus3x3>1) E3x31CUT->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);
            }*/
        }   
    
        if (detKinBeamRot_tar==1)
        {if(r_mu<1.7 && detKinBeamRot_def_angle_mu>0.2 && detKinBeamRot_E_clus3x3>1) The2CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);

    /*        if (detKinBeamRot_E_clus3x3!=0) 
            {if(r_mu<1.7 && detKinBeamRot_def_angle_mu>0.2 && detKinBeamRot_E_clus3x3>1) E3x32CUT->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);
            } */
        } 
    }
}
}
      
    
Int_t nx1CUT = E3x31CUT->GetNbinsX();
Int_t ny1CUT = E3x31CUT->GetNbinsY();
for (Int_t i=1; i<nx1CUT+1; i++) {
for (Int_t j=1; j<ny1CUT+1; j++) {
if (E3x31CUT->GetBinContent(i,j)<10) E3x31CUT->SetBinContent(i,j,0);}}
    
Int_t nx2CUT = E3x32CUT->GetNbinsX();
Int_t ny2CUT = E3x32CUT->GetNbinsY();
for (Int_t i=1; i<nx2CUT+1; i++) {
for (Int_t j=1; j<ny2CUT+1; j++) {
if (E3x32CUT->GetBinContent(i,j)<10) E3x32CUT->SetBinContent(i,j,0);}}
    
TCanvas * c4a= new TCanvas("c4a","c4a",100,100,2500,2000);
c4a->Divide(1,2);
gStyle->SetPalette(kLake);
TColor::InvertPalette(); 
c4a->cd(1);   
E3x31CUT->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x31CUT->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x31CUT->Draw("COLZ");
c4a->cd(2);   
E3x32CUT->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x32CUT->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x32CUT->Draw("COLZ");


c4a->SaveAs("/home/LHCB-T3/espedicato/tesi/studio/thE_cut.png");    
    
TCanvas * c9= new TCanvas("c9","c9",1000,100,2500,2000);
hist_E9_e->GetXaxis()->SetTitle("Ecentral/E3x3");
hist_E9_e->SetLineWidth(3);
hist_E9_e->Draw("HIST"); 

hist_E9_eLO->GetXaxis()->SetTitle("Ecentral/E3x3");
hist_E9_eLO->SetLineWidth(3);
hist_E9_eLO->SetLineColor(kRed);
hist_E9_eLO->Draw("HIST same"); 
gPad->BuildLegend(0.25,0.15,0.25,0.15);

c9->SaveAs("/home/LHCB-T3/espedicato/studio/tesi/E9.png");




}