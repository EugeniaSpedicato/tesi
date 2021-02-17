#define atree_cxx
#include "atree.h"
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>
#include <iostream>
#include <sstream>
#include <TGraph.h>
#include <TTree.h>
#include <cmath>

using namespace std;

#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{
    TH1::SetDefaultSumw2();

typedef map<int, double>  energy_cell; 
energy_cell number;
energy_cell Rev_number;
energy_cell en_c; 
Double_t E_CAL=0.;
Double_t E9=0.;
 
TH1F* TheCUT=new TH1F("th", "th El core CUT", 120,0,30); 
TH1F* The1CUT=new TH1F("th", "th El TAR 1 core CUT", 120,0,30); 
TH1F* The2CUT=new TH1F("th", "th El TAR 2 core CUT", 120,0,30);


TH2F  *E3x31CUT  = new TH2F("ThEel1" , " Th_el Vs. E_E3x3 core NLO CUT",90,0,30,280,0,140);
TH2F  *E3x32CUT  = new TH2F("ThEel2" , " Th_el Vs. E_E3x3 core LO CUT",90,0,30,280,0,140);


TH1F* hist_E9_e=new TH1F("E9e", "E9 e- tot", 500,0,1);
TH1F* hist_E9_eLO=new TH1F("E9eLO", "E9 e- tot LO", 500,0,1);

    

number[36]=1; number[37]=2; number[38]=3; number[39]=4; number[40]=5;
number[29]=6; number[30]=7; number[31]=8; number[32]=9; number[33]=10;
number[22]=11; number[23]=12; number[24]=13; number[25]=14; number[26]=15;
number[15]=16; number[16]=17; number[17]=18; number[18]=19; number[19]=20;
number[8]=21; number[9]=22; number[10]=23; number[11]=24; number[12]=25;    
    
Rev_number[1]=36; Rev_number[2]=37; Rev_number[3]=38; Rev_number[4]=39; Rev_number[5]=40;
Rev_number[6]=29; Rev_number[7]=30; Rev_number[8]=31; Rev_number[9]=32; Rev_number[10]=33;
Rev_number[11]=22; Rev_number[12]=23; Rev_number[13]=24; Rev_number[14]=25; Rev_number[15]=26;
Rev_number[16]=15; Rev_number[17]=16; Rev_number[18]=17; Rev_number[19]=18; Rev_number[20]=19;
Rev_number[21]=8; Rev_number[22]=9; Rev_number[23]=10; Rev_number[24]=11; Rev_number[25]=12;
    
int *Array9=0;
    
if (fChain == 0) return;

Long64_t nentries = fChain->GetEntriesFast();



    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
Double_t E_1=0.;
Double_t E_clus3x3=0.;
Double_t Etotcal=0.;
              
        detKinBeamRot_cooXe=detKinBeamRot_cooXe*100; // cm
        detKinBeamRot_cooYe=detKinBeamRot_cooYe*100; // cm
        detKinBeamRot_x_in=detKinBeamRot_x_in*100; // cm coo muone entrante
        detKinBeamRot_y_in=detKinBeamRot_y_in*100; // cm coo muone entrante
        photon_coox=photon_coox*100; // cm
        photon_cooy=photon_cooy*100; // cm
  
    
       double r_mu=sqrt((detKinBeamRot_x_in*detKinBeamRot_x_in)+(detKinBeamRot_y_in*detKinBeamRot_y_in));
       
       if (r_mu<1.7 && detKinBeamRot_def_angle_mu>0.2 && abs(detKinBeamRot_cooXe)<4.275 && abs(detKinBeamRot_cooYe)<4.275){
    TH2F* myGrid= new TH2F("myGrid" , "EM Calorimeter with E in GeV",5,-7.125,7.125,5,-7.125,7.125);
       en_c[1]=detKinBeamRot_Ecell1; en_c[2]=detKinBeamRot_Ecell2; en_c[3]=detKinBeamRot_Ecell3; en_c[4]=detKinBeamRot_Ecell4; en_c[5]=detKinBeamRot_Ecell5;
        en_c[6]=detKinBeamRot_Ecell6; en_c[7]=detKinBeamRot_Ecell7; en_c[8]=detKinBeamRot_Ecell8; en_c[9]=detKinBeamRot_Ecell9; en_c[10]=detKinBeamRot_Ecell10;
        en_c[11]=detKinBeamRot_Ecell11; en_c[12]=detKinBeamRot_Ecell12; en_c[13]=detKinBeamRot_Ecell13; en_c[14]=detKinBeamRot_Ecell14; en_c[15]=detKinBeamRot_Ecell15;
        en_c[16]=detKinBeamRot_Ecell16; en_c[17]=detKinBeamRot_Ecell17; en_c[18]=detKinBeamRot_Ecell18; en_c[19]=detKinBeamRot_Ecell19; en_c[20]=detKinBeamRot_Ecell20;
        en_c[21]=detKinBeamRot_Ecell21; en_c[22]=detKinBeamRot_Ecell22; en_c[23]=detKinBeamRot_Ecell23; en_c[24]=detKinBeamRot_Ecell24; en_c[25]=detKinBeamRot_Ecell25;
       
    for(int i=1;i<26;++i){myGrid->SetBinContent(Rev_number[i],en_c[i]);}
    
       
    int binx_e = myGrid->GetXaxis()->FindBin(detKinBeamRot_cooXe);
    int biny_e = myGrid->GetYaxis()->FindBin(detKinBeamRot_cooYe);
    int nbin_e = myGrid->GetBin(binx_e,biny_e);
    int n_cell_e=number[nbin_e];
    int binx_ph = myGrid->GetXaxis()->FindBin(photon_coox);
    int biny_ph = myGrid->GetYaxis()->FindBin(photon_cooy);
    int nbin_ph = myGrid->GetBin(binx_ph,biny_ph);
    int n_cell_ph=number[nbin_ph];
       
       
int binMax=myGrid->GetMaximumBin();  
int CentralCell=number[binMax];

E_1=myGrid->GetBinContent(binMax);
    if (CentralCell==1) {Array9= new int[9]{1,6,7,2,0,0,0,0,0};}
    if (CentralCell==2) {Array9= new int[9]{1,2,6,7,8,3,0,0,0};}
    if (CentralCell==3) {Array9= new int[9]{2,3,7,8,9,4,0,0,0};}
    if (CentralCell==4) {Array9= new int[9]{3,4,8,9,10,0,0,0};}
    if (CentralCell==5) {Array9= new int[9]{4,5,9,10,0,0,0};}
    if (CentralCell==6) {Array9= new int[9]{1,2,6,7,12,11,0,0,0};}
    if (CentralCell==7) {Array9= new int[9]{1,2,3,6,7,8,11,12,13};}
    if (CentralCell==8) {Array9= new int[9]{2,3,4,7,8,9,12,13,14};}
    if (CentralCell==9) {Array9= new int[9]{3,4,5,8,9,10,13,14,15};}
    if (CentralCell==10) {Array9= new int[9]{4,5,9,10,14,15,0,0,0};}
    if (CentralCell==11) {Array9= new int[9]{6,7,11,12,16,17,0,0,0};}
    if (CentralCell==12) {Array9= new int[9]{6,7,8,11,12,13,16,17,18};}
    if (CentralCell==13) {Array9= new int[9]{7,8,9,12,13,14,17,18,19};}
    if (CentralCell==14) {Array9= new int[9]{8,9,10,13,14,15,18,19,20};}
    if (CentralCell==15) {Array9= new int[9]{9,10,14,15,19,20,0,0,0};}
    if (CentralCell==16) {Array9= new int[9]{11,12,16,17,21,22,0,0,0};}
    if (CentralCell==17) {Array9= new int[9]{11,12,13,16,17,18,21,22,23};}
    if (CentralCell==18) {Array9= new int[9]{12,13,14,17,18,19,22,23,24};}
    if (CentralCell==19) {Array9= new int[9]{13,14,15,18,19,20,23,24,25};}
    if (CentralCell==20) {Array9= new int[9]{14,15,19,20,24,25,0,0,0};}
    if (CentralCell==21) {Array9= new int[9]{16,17,21,22,0,0,0,0,0};}
    if (CentralCell==22) {Array9= new int[9]{16,17,18,21,22,23,0,0,0};}
    if (CentralCell==23) {Array9= new int[9]{17,18,19,22,23,24,0,0,0};}
    if (CentralCell==24) {Array9= new int[9]{18,19,20,23,24,25,0,0,0};}
    if (CentralCell==25) {Array9= new int[9]{19,20,24,25,0,0,0,0,0};}  

for(int i=0; i<9; ++i)
{
    if (Array9[i]>0 && Array9[i]<26) E_clus3x3+=myGrid->GetBinContent(Rev_number[Array9[i]]);

}  

    double r=sqrt((detKinBeamRot_cooXe*detKinBeamRot_cooXe)+(detKinBeamRot_cooYe*detKinBeamRot_cooYe));
    
    E9=E_1/E_clus3x3;
    
    for(int i=1;i<26;++i)
    {Etotcal+=en_c[i];}
           
    double Eout=(Etotcal-E_clus3x3)/E_clus3x3;

        if(E_clus3x3>1)
        {TheCUT->Fill(detKinBeamRot_def_angle_e,wgt_full);
         if(n_cell_ph!=0)hist_E9_e->Fill(Eout,wgt_full);
         if(n_cell_ph==0)hist_E9_eLO->Fill(Eout,wgt_full);
    
        if (E_clus3x3!=0){E3x31CUT->Fill(detKinBeamRot_def_angle_e,E_clus3x3,wgt_full);} 
        
        if (E_clus3x3!=0){E3x32CUT->Fill(detKinBeamRot_def_angle_e,E_clus3x3,wgt_LO);} 
        
        if (detKinBeamRot_tar==0)
        {The1CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);}   
    
        if (detKinBeamRot_tar==1)
        {The2CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);} 
    }
           
           
/*TCanvas * Ecal_= new TCanvas("Ecal_","Ecal_",1500,100,3500,2000);
Ecal_->Divide(2,1);
Ecal_->cd(1);
gStyle->SetPalette(kAquamarine);
//TColor::InvertPalette();
myGrid->SetXTitle("x (cm)");
myGrid->SetYTitle("y (cm)");
myGrid->Draw("COL");
myGrid->Draw("TEXT SAME");
Ecal_->cd(2);
myGrid->Draw("LEGO");
std::ostringstream name1;
name1 <<"/home/LHCB-T3/espedicato/tesi/studio/Ecal"<< jentry << ".png";
TString name =name1.str();
Ecal_->SaveAs(name);   */        
 
delete myGrid; 
}}
      
    
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

c9->SaveAs("/home/LHCB-T3/espedicato/tesi/studio/E9.png");



}