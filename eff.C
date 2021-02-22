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
int n_cell_e=0;
int n_cell_ph=0;
    
    
typedef map<int, double>  energy_cell; 
energy_cell number;
energy_cell Rev_number;
energy_cell Rev_numberX;
energy_cell Rev_numberY;
energy_cell en_c; 
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
Rev_numberX[1]=1; Rev_numberX[2]=2; Rev_numberX[3]=3; Rev_numberX[4]=4; Rev_numberX[5]=5;
Rev_numberX[6]=1; Rev_numberX[7]=2; Rev_numberX[8]=3; Rev_numberX[9]=4; Rev_numberX[10]=5;
Rev_numberX[11]=1; Rev_numberX[12]=2; Rev_numberX[13]=3; Rev_numberX[14]=4; Rev_numberX[15]=5;
Rev_numberX[16]=1; Rev_numberX[17]=2; Rev_numberX[18]=3; Rev_numberX[19]=4; Rev_numberX[20]=5;
Rev_numberX[21]=1; Rev_numberX[22]=2; Rev_numberX[23]=3; Rev_numberX[24]=4; Rev_numberX[25]=5;
    
    
Rev_numberY[1]=5; Rev_numberY[2]=5; Rev_numberY[3]=5; Rev_numberY[4]=5; Rev_numberY[5]=5;
Rev_numberY[6]=4; Rev_numberY[7]=4; Rev_numberY[8]=4; Rev_numberY[9]=4; Rev_numberY[10]=4;
Rev_numberY[11]=3; Rev_numberY[12]=3; Rev_numberY[13]=3; Rev_numberY[14]=3; Rev_numberY[15]=3;
Rev_numberY[16]=2; Rev_numberY[17]=2; Rev_numberY[18]=2; Rev_numberY[19]=2; Rev_numberY[20]=2;
Rev_numberY[21]=1; Rev_numberY[22]=1; Rev_numberY[23]=1; Rev_numberY[24]=1; Rev_numberY[25]=1;
    
int *Array9=0;
Double_t E_1=0.;
Double_t E_clus3x3=0.;

    

TH1F* The_trueCUT=new TH1F("th", "th El true CUT", 120,0,30);
TH1F* The_trueCUT1=new TH1F("th", "th El true TAR 1 CUT", 120,0,30);
TH1F* The_trueCUT2=new TH1F("th", "th El true TAR 2 CUT", 120,0,30);

TH1F* The_trueCUTmu=new TH1F("th", "th El  true CUT th_mu>0.2mrad", 120,0,30);
TH1F* The_trueCUT1mu=new TH1F("th", "th El true TAR 1 CUT th_mu>0.2mrad", 120,0,30);
TH1F* The_trueCUT2mu=new TH1F("th", "th El true TAR 2 CUT th_mu>0.2mrad", 120,0,30);
    
TH1F* The_trueCUTEe=new TH1F("th", "th El true CUT on E_e", 120,0,30);
TH1F* The_trueCUT1Ee=new TH1F("th", "th El true TAR 1 CUT on E_e", 120,0,30);
TH1F* The_trueCUT2Ee=new TH1F("th", "th El true TAR 2 CUT on E_e", 120,0,30);
    
TH1F* The_trueCUTtot=new TH1F("th", "th El true CUT th_mu+Ee", 120,0,30);
TH1F* The_trueCUT1tot=new TH1F("th", "th El true TAR 1 CUT th_mu+Ee", 120,0,30);
TH1F* The_trueCUT2tot=new TH1F("th", "th El true TAR 2 CUT th_mu+Ee", 120,0,30);

TH1F* The_true=new TH1F("th", "th El true", 120,0,30);    
TH1F* The=new TH1F("th", "th El core", 120,0,30); 
/*TH1F* TheBIG=new TH1F("th", "th El 5X5", 120,0,30); 
TH1F* The2P=new TH1F("th", "th El crown", 120,0,30); */

TH1F* The_true1=new TH1F("th", "th El true TAR 1", 120,0,30);    
TH1F* The1=new TH1F("th", "th El core TAR 1", 120,0,30); 
/*TH1F* TheBIG1=new TH1F("th", "th El 5X5 TAR 1", 120,0,30); 
TH1F* The2P1=new TH1F("th", "th El crown TAR 1", 120,0,30); */

    
TH1F* The_true2=new TH1F("th", "th El true TAR 2", 120,0,30);     
TH1F* The2=new TH1F("th", "th El core TAR 2", 120,0,30); 
/*TH1F* TheBIG2=new TH1F("th", "th El 5X5 TAR 2", 120,0,30); 
TH1F* The2P2=new TH1F("th", "th El crown TAR 2", 120,0,30);*/ 
    
    
TH1F* TheCUT=new TH1F("th", "th El core CUT", 120,0,30); 
TH1F* The1CUT=new TH1F("th", "th El TAR 1 core CUT", 120,0,30); 
TH1F* The2CUT=new TH1F("th", "th El TAR 2 core CUT", 120,0,30); 

/*TH1F* TheBIGCUT=new TH1F("th", "th El 5X5 CUT", 120,0,30); 
TH1F* The2PCUT=new TH1F("th", "th El crown CUT", 120,0,30);*/ 
    
    
TH1F* TheCUTmu=new TH1F("th", "th El core CUT th_mu>0.2mrad", 120,0,30); 
TH1F* The1CUTmu=new TH1F("th", "th El TAR 1 core CUT th_mu>0.2mrad", 120,0,30); 
TH1F* The2CUTmu=new TH1F("th", "th El TAR 2 core CUT th_mu>0.2mrad", 120,0,30); 

/*TH1F* TheBIG1CUT=new TH1F("th", "th El 5X5 TAR 1 CUT", 120,0,30); 
TH1F* The2P1CUT=new TH1F("th", "th El crown TAR 1 CUT", 120,0,30); */
    
    
TH1F* TheCUTEe=new TH1F("th", "th El core CUT on E_e", 120,0,30); 
TH1F* The1CUTEe=new TH1F("th", "th El core TAR 1 CUT on E_e", 120,0,30); 
TH1F* The2CUTEe=new TH1F("th", "th El core TAR 2 CUT on E_e", 120,0,30); 
    
/*TH1F* TheBIG2CUT=new TH1F("th", "th El 5X5 TAR 2 CUT", 120,0,30);
TH1F* The2P2CUT=new TH1F("th", "th El crown TAR 2 CUT", 120,0,30); */

TH1F* TheCUTtot=new TH1F("th", "th El core CUT th_mu+Ee", 120,0,30); 
TH1F* The1CUTtot=new TH1F("th", "th El core TAR 1 CUT th_mu+Ee", 120,0,30); 
TH1F* The2CUTtot=new TH1F("th", "th El core TAR 2 CUT th_mu+Ee", 120,0,30); 
    
TH1F* rmu=new TH1F("rmu", "impact point", 100,0,10); 


TH2F  *E3x31CUT  = new TH2F("ThEel1" , " Th_el Vs. E_E3x3 core TAR1 CUT",90,0,30,280,0,140);
TH2F  *E3x32CUT  = new TH2F("ThEel2" , " Th_el Vs. E_E3x3 core TAR2 CUT",90,0,30,280,0,140);
    
TH2F  *E3x31CUTmu  = new TH2F("ThEel1" , " Th_el Vs. E_E3x3 core TAR1 CUT th_mu>0.2mrad",90,0,30,280,0,140);
TH2F  *E3x32CUTmu  = new TH2F("ThEel2" , " Th_el Vs. E_E3x3 core TAR2 CUT th_mu>0.2mrad",90,0,30,280,0,140);
    
TH2F  *E3x31CUTEe  = new TH2F("ThEel1" , " Th_el Vs. E_E3x3 core TAR1 CUT on E_e",90,0,30,280,0,140);
TH2F  *E3x32CUTEe  = new TH2F("ThEel2" , " Th_el Vs. E_E3x3 core TAR2 CUT on E_e",90,0,30,280,0,140);
    
TH2F  *E3x31CUTtot  = new TH2F("ThEel1" , " Th_el Vs. E_E3x3 core TAR1 CUT th_mu+Ee",90,0,30,280,0,140);
TH2F  *E3x32CUTtot  = new TH2F("ThEel2" , " Th_el Vs. E_E3x3 core TAR2 CUT th_mu+Ee",90,0,30,280,0,140);

TH1F* hist_E9_e=new TH1F("E9e", "E9 e- tot", 500,0,1);
TH1F* hist_E9_eLO=new TH1F("E9eLO", "E9 e- tot LO", 500,0,1);
    
    
if (fChain == 0) return;

Long64_t nentries = fChain->GetEntriesFast();



    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
   
       
TH2F* myGrid= new TH2F("myGrid" , "EM Calorimeter with E in GeV",5,-7.125,7.125,5,-7.125,7.125);
       en_c[1]=detKinBeamRot_Ecell1; en_c[2]=detKinBeamRot_Ecell2; en_c[3]=detKinBeamRot_Ecell3; en_c[4]=detKinBeamRot_Ecell4; en_c[5]=detKinBeamRot_Ecell5;
        en_c[6]=detKinBeamRot_Ecell6; en_c[7]=detKinBeamRot_Ecell7; en_c[8]=detKinBeamRot_Ecell8; en_c[9]=detKinBeamRot_Ecell9; en_c[10]=detKinBeamRot_Ecell10;
        en_c[11]=detKinBeamRot_Ecell11; en_c[12]=detKinBeamRot_Ecell12; en_c[13]=detKinBeamRot_Ecell13; en_c[14]=detKinBeamRot_Ecell14; en_c[15]=detKinBeamRot_Ecell15;
        en_c[16]=detKinBeamRot_Ecell16; en_c[17]=detKinBeamRot_Ecell17; en_c[18]=detKinBeamRot_Ecell18; en_c[19]=detKinBeamRot_Ecell19; en_c[20]=detKinBeamRot_Ecell20;
        en_c[21]=detKinBeamRot_Ecell21; en_c[22]=detKinBeamRot_Ecell22; en_c[23]=detKinBeamRot_Ecell23; en_c[24]=detKinBeamRot_Ecell24; en_c[25]=detKinBeamRot_Ecell25;
       
    for(int i=1;i<26;++i){myGrid->SetBinContent(Rev_number[i],en_c[i]);}    
       

int binMax=myGrid->GetMaximumBin();  
int CentralCell=number[binMax];
       

E_1=myGrid->GetBinContent(binMax);

// con if E_1!=0 è gia imposto r<1.7 perchè gia in fastsim         
if(E_1!=0){           
    if (detKinBeamRot_Ee>0.2){int binx_e = myGrid->GetXaxis()->FindBin(detKinBeamRot_cooXe);
    int biny_e = myGrid->GetYaxis()->FindBin(detKinBeamRot_cooYe);
    int nbin_e = myGrid->GetBin(binx_e,biny_e);
     n_cell_e=number[nbin_e];} else n_cell_e=0;
    if(photon_coox!=-100 && photon_energy>0.2)
    {int binx_ph = myGrid->GetXaxis()->FindBin(photon_coox);
    int biny_ph = myGrid->GetYaxis()->FindBin(photon_cooy);
    int nbin_ph = myGrid->GetBin(binx_ph,biny_ph);
     n_cell_ph=number[nbin_ph];}
    else n_cell_ph=0;}  

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
           
       
       
       
       
              
    if (photon_coox!=-1 && photon_cooy!=-1)
    {  
     E_CAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {E_CAL=detKinBeamRot_Ee;}

        detKinBeamRot_cooXe=detKinBeamRot_cooXe*100; // cm
        detKinBeamRot_cooYe=detKinBeamRot_cooYe*100; // cm
        detKinBeamRot_x_in=detKinBeamRot_x_in*100; // cm coo muone entrante
        detKinBeamRot_y_in=detKinBeamRot_y_in*100; // cm coo muone entrante
        photon_coox=photon_coox*100; // cm
        photon_cooy=photon_cooy*100; // cm
               
       double r_mu=sqrt((detKinBeamRot_x_in*detKinBeamRot_x_in)+(detKinBeamRot_y_in*detKinBeamRot_y_in));
              
        E9=E_1/E_clus3x3;

       
The_true->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<5) The_trueCUT->Fill(detKinBeamRot_def_angle_e,wgt_full);
//if(r_mu<2 && detKinBeamRot_def_angle_mu>0.2) The_trueCUTmu->Fill(detKinBeamRot_def_angle_e,wgt_full);
//if(r_mu<1.7 && E_clus3x3>2.5) The_trueCUTEe->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<1.7 && E_clus3x3>2.5) The_trueCUTtot->Fill(detKinBeamRot_def_angle_e,wgt_full);
    

    if (detKinBeamRot_tar==0)
    {The_true1->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<5) The_trueCUT1->Fill(detKinBeamRot_def_angle_e,wgt_full);
//if(r_mu<2 && detKinBeamRot_def_angle_mu>0.2) The_trueCUT1mu->Fill(detKinBeamRot_def_angle_e,wgt_full);
//if(r_mu<1.7 && E_clus3x3>2.5) The_trueCUT1Ee->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<1.7 && E_clus3x3>2.5) The_trueCUT1tot->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    if (detKinBeamRot_tar==1)
{The_true2->Fill(detKinBeamRot_def_angle_e,wgt_full);
     
if(r_mu<5) The_trueCUT2->Fill(detKinBeamRot_def_angle_e,wgt_full);
//if(r_mu<2 && detKinBeamRot_def_angle_mu>0.2) The_trueCUT2mu->Fill(detKinBeamRot_def_angle_e,wgt_full);
//if(r_mu<1.7 && E_clus3x3>2.5) The_trueCUT2Ee->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<1.7 && E_clus3x3>2.5) The_trueCUT2tot->Fill(detKinBeamRot_def_angle_e,wgt_full);}
       
    if (detKinBeamRot_tar==0) rmu->Fill(r_mu,wgt_full);
       
if (n_cell_e!=0)
{  

    /*TheBIG->Fill(detKinBeamRot_def_angle_e,wgt_full);
    TheBIGCUT->Fill(detKinBeamRot_def_angle_e,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {TheBIG1->Fill(detKinBeamRot_def_angle_e,wgt_full);
    TheBIG1CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    if (detKinBeamRot_tar==1)
    {TheBIG2->Fill(detKinBeamRot_def_angle_e,wgt_full);
    TheBIG2CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);}*/

    
if(abs(detKinBeamRot_cooXe)<4.275 && abs(detKinBeamRot_cooYe)<4.275) {

The->Fill(detKinBeamRot_def_angle_e,wgt_full);
    
if(r_mu<5) TheCUT->Fill(detKinBeamRot_def_angle_e,wgt_full);
//if(r_mu<2 && detKinBeamRot_def_angle_mu>0.2) TheCUTmu->Fill(detKinBeamRot_def_angle_e,wgt_full);
//if(r_mu<1.7 && E_clus3x3>2.5) TheCUTEe->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<1.7 && E_clus3x3>2.5)
{TheCUTtot->Fill(detKinBeamRot_def_angle_e,wgt_full);
     hist_E9_e->Fill(E9,wgt_full);
     hist_E9_eLO->Fill(E9,wgt_LO);}
    
    if (detKinBeamRot_tar==0)
    {The1->Fill(detKinBeamRot_def_angle_e,wgt_full);
     
if(r_mu<5) The1CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);
//if(r_mu<2 && detKinBeamRot_def_angle_mu>0.2) The1CUTmu->Fill(detKinBeamRot_def_angle_e,wgt_full);
//if(r_mu<1.7 && E_clus3x3>2.5) The1CUTEe->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<1.7 && E_clus3x3>2.5)The1CUTtot->Fill(detKinBeamRot_def_angle_e,wgt_full);

     if (E_clus3x3!=0) 
     {
if(r_mu<5) E3x31CUT->Fill(detKinBeamRot_def_angle_e,E_clus3x3,wgt_full);
//if(r_mu<2 && detKinBeamRot_def_angle_mu>0.2) E3x31CUTmu->Fill(detKinBeamRot_def_angle_e,E_clus3x3,wgt_full);
//if(r_mu<1.7 && E_clus3x3>2.5) E3x31CUTEe->Fill(detKinBeamRot_def_angle_e,E_clus3x3,wgt_full);
if(r_mu<1.7 && E_clus3x3>2.5) E3x31CUTtot->Fill(detKinBeamRot_def_angle_e,E_clus3x3,wgt_full);
} }
    if (detKinBeamRot_tar==1)
    {The2->Fill(detKinBeamRot_def_angle_e,wgt_full);
     
if(r_mu<5) The2CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);
//if(r_mu<2 && detKinBeamRot_def_angle_mu>0.2) The2CUTmu->Fill(detKinBeamRot_def_angle_e,wgt_full);
//if(r_mu<1.7 && E_clus3x3>2.5) The2CUTEe->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<1.7 && E_clus3x3>2.5)The2CUTtot->Fill(detKinBeamRot_def_angle_e,wgt_full);
     if (E_clus3x3!=0) 
     {
if(r_mu<5) E3x32CUT->Fill(detKinBeamRot_def_angle_e,E_clus3x3,wgt_full);
//if(r_mu<2 && detKinBeamRot_def_angle_mu>0.2) E3x32CUTmu->Fill(detKinBeamRot_def_angle_e,E_clus3x3,wgt_full);
//if(r_mu<1.7 && E_clus3x3>2.5) E3x32CUTEe->Fill(detKinBeamRot_def_angle_e,E_clus3x3,wgt_full);
if(r_mu<1.7 && E_clus3x3>2.5) E3x32CUTtot->Fill(detKinBeamRot_def_angle_e,E_clus3x3,wgt_full);
} 
    
    } /*else { 

    The2P->Fill(detKinBeamRot_def_angle_e,wgt_full);
    The2PCUT->Fill(detKinBeamRot_def_angle_e,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {The2P1->Fill(detKinBeamRot_def_angle_e,wgt_full);
    The2P1CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    if (detKinBeamRot_tar==1)
    {The2P2->Fill(detKinBeamRot_def_angle_e,wgt_full);
    The2P2CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);}
}*/
     
}}delete myGrid;}
    
TH1F *Eff1CUT = new TH1F("ef1cut", "Eff Tar1 cut", 120,0,30);
TH1F *Eff2CUT = new TH1F("ef2cut", "Eff Tar2 cut", 120,0,30);
Eff1CUT->Divide(The1CUT,The_trueCUT1,1,1,"B");
Eff2CUT->Divide(The2CUT,The_trueCUT2,1,1,"B");

TH1F *Eff1CUTmu = new TH1F("ef1cutmu", "Eff Tar1 cut th_mu>0.2mrad", 120,0,30);
TH1F *Eff2CUTmu = new TH1F("ef2cutmu", "Eff Tar2 cut th_mu>0.2mrad", 120,0,30);
Eff1CUTmu->Divide(The1CUTmu,The_trueCUT1mu,1,1,"B");
Eff2CUTmu->Divide(The2CUTmu,The_trueCUT2mu,1,1,"B");
    
TH1F *Eff1CUTEe = new TH1F("ef1cutEe", "Eff Tar1 cut E_e", 120,0,30);
TH1F *Eff2CUTEe = new TH1F("ef2cutEe", "Eff Tar2 cut E_e", 120,0,30);
Eff1CUTEe->Divide(The1CUTEe,The_trueCUT1Ee,1,1,"B");
Eff2CUTEe->Divide(The2CUTEe,The_trueCUT2Ee,1,1,"B");

TH1F *Eff1CUTtot = new TH1F("ef1cuttot", "Eff Tar1 cut th_mu+Ee", 120,0,30);
TH1F *Eff2CUTtot = new TH1F("ef2cuttot", "Eff Tar2 cut th_mu+Ee", 120,0,30);
Eff1CUTtot->Divide(The1CUTtot,The_trueCUT1tot,1,1,"B");
Eff2CUTtot->Divide(The2CUTtot,The_trueCUT2tot,1,1,"B");

       

TCanvas * ef= new TCanvas("ef","ef",1000,100,2500,2000);
ef->Divide(1,2);
ef->cd(1); 
Eff1CUT->GetXaxis()->SetTitle("Theta el[mrad]");
Eff1CUT->GetYaxis()->SetTitle("Efficency");
Eff1CUT->SetLineWidth(3);
Eff1CUT->SetLineColor(kRed);
Eff1CUT->SetMinimum(0);
Eff1CUT->Draw(); 

ef->cd(2);
Eff2CUT->GetXaxis()->SetTitle("Theta el[mrad]");
Eff2CUT->GetYaxis()->SetTitle("Efficency");
Eff2CUT->SetLineWidth(3);
Eff2CUT->SetLineColor(kRed);
Eff2CUT->SetMinimum(0.7);
Eff2CUT->Draw();   
    
gPad->BuildLegend(0.25,0.15,0.25,0.15);

ef->SaveAs("/home/LHCB-T3/espedicato/tesi/eff/Effrmu5.png");   
    
    
    TCanvas * ef1= new TCanvas("ef","ef",1000,100,2500,2000);
ef1->Divide(1,2);
ef1->cd(1);
Eff1CUT->Draw(); 
Eff1CUTtot->GetXaxis()->SetTitle("Theta el[mrad]");
Eff1CUTtot->GetYaxis()->SetTitle("Efficency");
Eff1CUTtot->SetLineWidth(3);
Eff1CUTtot->SetLineColor(kBlue);
Eff1CUTtot->SetMinimum(0);
Eff1CUTtot->Draw(); 
gPad->BuildLegend(0.25,0.15,0.25,0.15);
ef1->cd(2);  
Eff2CUTtot->GetXaxis()->SetTitle("Theta el[mrad]");
Eff2CUTtot->GetYaxis()->SetTitle("Efficency");
Eff2CUTtot->SetLineWidth(3);
Eff2CUTtot->SetLineColor(kBlue);
Eff2CUTtot->SetMinimum(0.7);
Eff2CUTtot->Draw(); 
gPad->BuildLegend(0.25,0.15,0.25,0.15);

ef1->SaveAs("/home/LHCB-T3/espedicato/tesi/eff/Effcut.png");


    
/*TCanvas * efCUT= new TCanvas("efcut","efcut",1000,100,2500,2000);
efCUT->Divide(1,2);
efCUT->cd(1);
Eff1CUT->GetXaxis()->SetTitle("Theta el[mrad]");
Eff1CUT->GetYaxis()->SetTitle("Efficency");
Eff1CUT->SetLineWidth(3);
Eff1CUT->SetLineColor(kBlack);
Eff1CUT->SetMaximum(1);
Eff1CUT->SetMinimum(0);
Eff1CUT->Draw();  
gPad->BuildLegend(0.25,0.15,0.25,0.15);
efCUT->cd(2);
Eff2CUT->GetXaxis()->SetTitle("Theta el[mrad]");
Eff2CUT->GetYaxis()->SetTitle("Efficency");
Eff2CUT->SetLineWidth(3);
Eff2CUT->SetLineColor(kRed);
Eff2CUT->SetMaximum(1);
Eff2CUT->SetMinimum(0);
Eff2CUT->Draw();   
gPad->BuildLegend(0.25,0.15,0.25,0.15);
efCUT->SaveAs("/home/LHCB-T3/espedicato/tesi/eff/EffCUT.png");*/
    
    
TCanvas * c5= new TCanvas("c5","c5",1000,1000,2500,2000);
c5->Divide(2,1);
/*c5->cd(1);
The_true->SetLineColor(kBlack);
The_true->SetLineWidth(3);
The_true->Draw("HIST"); 
The->SetLineWidth(3);
The->Draw("HIST same");
gPad->BuildLegend(0.3,0.21,0.3,0.21);*/
c5->cd(1);
The_trueCUT1->SetLineColor(kBlack);
The_trueCUT1->SetLineWidth(3);
The_trueCUT1->Draw("HIST");
The1CUT->SetLineWidth(3);
The1CUT->Draw("HIST same");  

gPad->BuildLegend(0.3,0.21,0.3,0.21);
c5->cd(2);
The_trueCUT2->SetLineColor(kBlack);
The_trueCUT2->SetLineWidth(3);
The_trueCUT2->Draw("HIST");   
The2CUT->SetLineWidth(3);
The2CUT->Draw("HIST same"); 

gPad->BuildLegend(0.3,0.21,0.3,0.21);
    
c5->SaveAs("/home/LHCB-T3/espedicato/tesi/eff/th_el.png"); 
    
TCanvas * c5MCS= new TCanvas("c5MCS","c5MCS",1000,100,2500,2000);
c5MCS->Divide(2,4);
c5MCS->cd(1);
The_trueCUT1->SetLineColor(kRed);
The_trueCUT1->SetLineWidth(3);
The_trueCUT1->Draw("HIST");
The1CUT->SetLineWidth(3);
The1CUT->Draw("HIST same");  

gPad->BuildLegend(0.3,0.21,0.3,0.21);
c5MCS->cd(2);
The_trueCUT2->SetLineColor(kRed);
The_trueCUT2->SetLineWidth(3);
The_trueCUT2->Draw("HIST");     

The2CUT->SetLineWidth(3);
The2CUT->Draw("HIST same"); 

gPad->BuildLegend(0.3,0.21,0.3,0.21);
    
    
c5MCS->cd(3);
The_trueCUT1mu->SetLineColor(kOrange);
The_trueCUT1mu->SetLineWidth(3);
The_trueCUT1mu->Draw("HIST");
The1CUTmu->SetLineWidth(3);
The1CUTmu->Draw("HIST same");  

gPad->BuildLegend(0.3,0.21,0.3,0.21);
c5MCS->cd(4);
The_trueCUT2mu->SetLineColor(kOrange);
The_trueCUT2mu->SetLineWidth(3);
The_trueCUT2mu->Draw("HIST");     

The2CUTmu->SetLineWidth(3);
The2CUTmu->Draw("HIST same"); 

gPad->BuildLegend(0.3,0.21,0.3,0.21);
    
    
c5MCS->cd(5);
The_trueCUT1Ee->SetLineColor(kGreen);
The_trueCUT1Ee->SetLineWidth(3);
The_trueCUT1Ee->Draw("HIST");
The1CUTEe->SetLineWidth(3);
The1CUTEe->Draw("HIST same");  

gPad->BuildLegend(0.3,0.21,0.3,0.21);
c5MCS->cd(6);
The_trueCUT2Ee->SetLineColor(kGreen);
The_trueCUT2Ee->SetLineWidth(3);
The_trueCUT2Ee->Draw("HIST");     

The2CUTEe->SetLineWidth(3);
The2CUTEe->Draw("HIST same"); 

gPad->BuildLegend(0.3,0.21,0.3,0.21);
    
c5MCS->cd(7);
The_trueCUT1tot->SetLineColor(kBlack);
The_trueCUT1tot->SetLineWidth(3);
The_trueCUT1tot->Draw("HIST");
The1CUTtot->SetLineWidth(3);
The1CUTtot->Draw("HIST same");  

gPad->BuildLegend(0.3,0.21,0.3,0.21);
c5MCS->cd(8);
The_trueCUT2tot->SetLineColor(kBlack);
The_trueCUT2tot->SetLineWidth(3);
The_trueCUT2tot->Draw("HIST");     

The2CUTtot->SetLineWidth(3);
The2CUTtot->Draw("HIST same"); 

gPad->BuildLegend(0.3,0.21,0.3,0.21);
        
    
c5MCS->SaveAs("/home/LHCB-T3/espedicato/tesi/eff/th_elCUT.png"); 
  
TCanvas * c= new TCanvas("c5MCS","c5MCS",1000,100,2500,2000);
rmu->Draw();
c->SaveAs("/home/LHCB-T3/espedicato/tesi/eff/rmu.png"); 
    
    
Int_t nx1CUT = E3x31CUT->GetNbinsX();
Int_t ny1CUT = E3x31CUT->GetNbinsY();
for (Int_t i=1; i<nx1CUT+1; i++) {
for (Int_t j=1; j<ny1CUT+1; j++) {
if (E3x31CUT->GetBinContent(i,j)<5) E3x31CUT->SetBinContent(i,j,0);}}
    
Int_t nx2CUT = E3x32CUT->GetNbinsX();
Int_t ny2CUT = E3x32CUT->GetNbinsY();
for (Int_t i=1; i<nx2CUT+1; i++) {
for (Int_t j=1; j<ny2CUT+1; j++) {
if (E3x32CUT->GetBinContent(i,j)<5) E3x32CUT->SetBinContent(i,j,0);}}
    
    
    
Int_t nx1CUTmu = E3x31CUTmu->GetNbinsX();
Int_t ny1CUTmu = E3x31CUTmu->GetNbinsY();
for (Int_t i=1; i<nx1CUTmu+1; i++) {
for (Int_t j=1; j<ny1CUTmu+1; j++) {
if (E3x31CUTmu->GetBinContent(i,j)<5) E3x31CUTmu->SetBinContent(i,j,0);}}
    
Int_t nx2CUTmu = E3x32CUTmu->GetNbinsX();
Int_t ny2CUTmu = E3x32CUTmu->GetNbinsY();
for (Int_t i=1; i<nx2CUTmu+1; i++) {
for (Int_t j=1; j<ny2CUTmu+1; j++) {
if (E3x32CUTmu->GetBinContent(i,j)<5) E3x32CUTmu->SetBinContent(i,j,0);}} 
  

Int_t nx1CUTEe = E3x31CUTEe->GetNbinsX();
Int_t ny1CUTEe = E3x31CUTEe->GetNbinsY();
for (Int_t i=1; i<nx1CUTEe+1; i++) {
for (Int_t j=1; j<ny1CUTEe+1; j++) {
if (E3x31CUTEe->GetBinContent(i,j)<5) E3x31CUTEe->SetBinContent(i,j,0);}}
    
Int_t nx2CUTEe = E3x32CUTEe->GetNbinsX();
Int_t ny2CUTEe = E3x32CUTEe->GetNbinsY();
for (Int_t i=1; i<nx2CUTEe+1; i++) {
for (Int_t j=1; j<ny2CUTEe+1; j++) {
if (E3x32CUTEe->GetBinContent(i,j)<5) E3x32CUTEe->SetBinContent(i,j,0);}}
    
    
Int_t nx1CUTtot = E3x31CUTtot->GetNbinsX();
Int_t ny1CUTtot = E3x31CUTtot->GetNbinsY();
for (Int_t i=1; i<nx1CUT+1; i++) {
for (Int_t j=1; j<ny1CUT+1; j++) {
if (E3x31CUTtot->GetBinContent(i,j)<5) E3x31CUTtot->SetBinContent(i,j,0);}}
    
Int_t nx2CUTtot = E3x32CUTtot->GetNbinsX();
Int_t ny2CUTtot = E3x32CUTtot->GetNbinsY();
for (Int_t i=1; i<nx2CUTtot+1; i++) {
for (Int_t j=1; j<ny2CUTtot+1; j++) {
if (E3x32CUTtot->GetBinContent(i,j)<5) E3x32CUTtot->SetBinContent(i,j,0);}}
    
TCanvas * c4a= new TCanvas("c4a","c4a",100,100,2500,2000);
c4a->Divide(2,4);
c4a->cd(1);
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
c4a->cd(3);   
E3x31CUTmu->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x31CUTmu->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x31CUTmu->Draw("COLZ");
c4a->cd(4);   
E3x32CUTmu->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x32CUTmu->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x32CUTmu->Draw("COLZ");
c4a->cd(5);   
E3x31CUTEe->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x31CUTEe->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x31CUTEe->Draw("COLZ");
c4a->cd(6);   
E3x32CUTEe->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x32CUTEe->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x32CUTEe->Draw("COLZ");
c4a->cd(7);   
E3x31CUTtot->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x31CUTtot->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x31CUTtot->Draw("COLZ");
c4a->cd(8);   
E3x32CUTtot->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x32CUTtot->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x32CUTtot->Draw("COLZ");


c4a->SaveAs("/home/LHCB-T3/espedicato/tesi/eff/thE_cut.png");    
    
TCanvas * c9= new TCanvas("c9","c9",1000,100,2500,2000);
hist_E9_e->GetXaxis()->SetTitle("Ecentral/E3x3");
hist_E9_e->SetLineWidth(3);
hist_E9_e->Draw("HIST"); 

hist_E9_eLO->GetXaxis()->SetTitle("Ecentral/E3x3");
hist_E9_eLO->SetLineWidth(3);
hist_E9_eLO->SetLineColor(kRed);
hist_E9_eLO->Draw("HIST same"); 
gPad->BuildLegend(0.25,0.15,0.25,0.15);

c9->SaveAs("/home/LHCB-T3/espedicato/tesi/eff/E9.png");

}