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
energy_cell Rev_numberX;
energy_cell Rev_numberY;
    
int n=0;
energy_cell en_c; 
int n_cell_e;
int n_cell_ph;
Double_t E_CAL=0.;
Double_t E9=0.;
Double_t Eout=0.;
Double_t Emean_out=0.;
Double_t Eres_in=0.;
   
const Double_t m_e = 0.5109989461 *0.001; //GeV
const Double_t m_mu =105.6583745 *0.001; //GeV
    
    
TH1F* hist_dist=new TH1F("dist1", "Dist e-centroide", 200,0,2);
TH1F* hist_distCUT=new TH1F("dist2", "Dist e-centroide CUT", 200,0,2);


TH2F  *E3x31CUT  = new TH2F("Eel1" , " Th_el Vs. E_3x3 ",180,0,0.030,380,0,140);
TH2F  *E3x32CUT  = new TH2F("Eel2" , " Th_el Vs. E_3x3 (Fiducial cut) ",180,0,0.030,380,0,140);
    
    
TH2F  *Th1  = new TH2F("ThEel1" , " Th_el Vs. Th_mu core Tar 2 ",180,0,0.030,250,0,0.005);
TH2F  *Th2  = new TH2F("ThEel2" , " Th_el Vs. Th_mu core Tar 2  (Fiducial cut) ",180,0,0.030,250,0,0.005);    

 

TH1F* hist_E9_e=new TH1F("E9e", "E9", 100,0.,1);
TH1F* hist_E9_eCUT=new TH1F("E9eCUT", "E9 cut", 100,0.,1);
    
TH1F* hist_Eout_e=new TH1F("en", "Eout", 100,0.,0.1);
TH1F* hist_Eout_eCUT=new TH1F("encut", "Eout cut", 100,0.,0.1);
    
TH1F* hist_E92_e=new TH1F("Emean", "Mean E out", 100,0.,0.15);
TH1F* hist_E92_eCUT=new TH1F("Emeancut", "Mean E out cut", 100,0.,0.15);
    
TH1F* hist_E3x3_e=new TH1F("E3x3", "Energy Reco 3x3", 70,0.,140);
TH1F* hist_E3x3_eCUT=new TH1F("E3x3cut", "Energy Reco 3x3 cut", 70,0.,140); 


TH1F* hist_DE=new TH1F("E3x3", "DE", 100,-0.7,0.7);
TH1F* hist_DE5=new TH1F("E3x3", "DE", 100,-0.4,0.4);
    
    
    
//caratteristiche fotoni
    TH1F* Ephout=new TH1F("EnergyPH", "Energy Ph tot", 75,0.2,150); 
    TH1F* Thphout=new TH1F("thetaPH", "Theta gen Ph tot", 110,0.,100); 
    TH1F* diff_th_phe=new TH1F("thetaPH", "Diff Th_e-Th_ph tot", 50,-30,30); 
    TH1F* diff_r_phe=new TH1F("thetaPH", "Diff r_e-r_ph tot", 75,0,10); 
    
    
    TH1F* EphoutCUT=new TH1F("EnergyPH1", "Energy Ph cut off", 75,0.2,150); 
    TH1F* ThphoutCUT=new TH1F("thetaPH1", "Theta gen Ph cut off", 110,0.,100); 
    TH1F* diff_th_pheCUT=new TH1F("thetaPH1", "Diff Th_e-Th_ph cut off", 50,-30,30); 
    TH1F* diff_r_pheCUT =new TH1F("thetaPH", "Diff r_e-r_ph cut off", 75,0,10); 
    

    TH1F* EphoutCUTafter=new TH1F("EnergyPH1", "Energy Ph after cut off", 75,0.2,150); 
    TH1F* ThphoutCUTafter=new TH1F("thetaPH1", "Theta gen Ph after cut off", 110,0.,100); 
    TH1F* diff_th_pheCUTafter=new TH1F("thetaPH1", "Diff Th_e-Th_ph after cut off", 50,-30,30); 
    TH1F* diff_r_pheCUTafter =new TH1F("thetaPH", "Diff r_e-r_ph after cut off", 75,0,10); 

    
    TH1F* Eeout=new TH1F("EnergyPH", "Energy el- tot", 75,0.2,150); 
    TH1F* Theout=new TH1F("thetaPH", "Theta gen el- tot", 75,0.,40); 
TH1F* EeoutCUT=new TH1F("EnergyPH1", "Energy el- cut off", 75,0.2,150); 
    TH1F* TheoutCUT=new TH1F("thetaPH1", "Theta gen el- cut off", 75,0.,40); 
 TH1F* EeoutCUTafter=new TH1F("EnergyPH1", "Energy el- after cut off", 75,0.2,150); 
    TH1F* TheoutCUTafter=new TH1F("thetaPH1", "Theta gen el- after cut off", 75,0.,40); 
    
    
    TH1F* Emuout=new TH1F("EnergyPH", "Energy mu tot", 75,0.2,160); 
    TH1F* Thmuout=new TH1F("thetaPH", "Theta gen mu tot", 100,0.,5); 
TH1F* EmuoutCUT=new TH1F("EnergyPH1", "Energy mu cut off", 75,0.2,160); 
    TH1F* ThmuoutCUT=new TH1F("thetaPH1", "Theta gen mu cut off", 100,0.,5); 
 TH1F* EmuoutCUTafter=new TH1F("EnergyPH1", "Energy mu after cut off", 75,0.2,160); 
    TH1F* ThmuoutCUTafter=new TH1F("thetaPH1", "Theta gen mu after cut off", 100,0.,5); 
    

    
    
    TH1F* DeltaR=new TH1F("res", "r_cal-r_trak", 100,0,2);
    TH1F* DeltaRCUT=new TH1F("resLO", "r_cal-r_trak LO", 100,0,2);
    
    TH1F* residuoX=new TH1F("res1", "Residual X_cal-X_trak", 30,-0.8,0.8);
    TH1F* residuoY=new TH1F("res2", "Residual Y_cal-Y_trak", 30,-0.8,0.8);

    
    TH1F* diff_r_mue=new TH1F("thetaPH", "Diff r_e-r_mu tot", 75,0,10); 
    TH1F* diff_r_mueCUT=new TH1F("thetaPH", "Diff r_e-r_mu tot", 75,0,10); 
    


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
    
if (fChain == 0) return;

Long64_t nentries = fChain->GetEntriesFast();



    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
Double_t E_1=0.;
Double_t E2=0.;
Double_t E3=0.;
Double_t E4=0.;
       
Double_t E_clus3x3=0.;
Double_t Etotcal=0.;
              
        detKinBeamRot_cooXe=detKinBeamRot_cooXe*100; // cm
        detKinBeamRot_cooYe=detKinBeamRot_cooYe*100; // cm
        detKinBeamRot_x_in=detKinBeamRot_x_in*100; // cm coo muone entrante
        detKinBeamRot_y_in=detKinBeamRot_y_in*100; // cm coo muone entrante
        photon_coox=photon_coox*100; // cm
        photon_cooy=photon_cooy*100; // cm
         
       double r_mu=sqrt((detKinBeamRot_x_in*detKinBeamRot_x_in)+(detKinBeamRot_y_in*detKinBeamRot_y_in));
       
       double d_e_ph=sqrt((detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy)); 
       
       double diffTh=detKinBeamRot_def_angle_e-photon_def_angle_ph;

       if (r_mu<5){
           
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

// con if E_1!=0 è gia imposto r<1 perchè gia in fastsim         
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
    else n_cell_ph=0;           

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
    
for(int i=1;i<26;++i)
    {Etotcal+=en_c[i];}
           
if(CentralCell==7 || CentralCell==8 || CentralCell==9 || CentralCell==12 || CentralCell==13 || CentralCell==14 || CentralCell==17 || CentralCell==18 || CentralCell==19)
{ 

//seconda    
double en_Maxcell=0.;
int Maxcell=0.;

for(int i=1;i<26;++i)
{
    if(en_c[i]>en_Maxcell && i!=CentralCell){en_Maxcell=en_c[i]; Maxcell=i;}
    else continue;
}
int SecondCentralCell=Maxcell;
E2=en_Maxcell ;
    
int SecondCentralCell_in9=0;
 for(int i=0; i<9; ++i)
{          
if (SecondCentralCell==Array9[i] && SecondCentralCell!=0)
    {//++sec_9;
    SecondCentralCell_in9=SecondCentralCell;
    SecondCentralCell=0;
    break;} 
else continue;
}
//terza
double en_Maxcell1=0.;
int Maxcell1=0.;

for(int i=0; i<9; ++i)
{
    if(en_c[Array9[i]]>en_Maxcell1 && Array9[i]!=CentralCell && Array9[i]!=SecondCentralCell_in9){en_Maxcell1=en_c[Array9[i]]; Maxcell1=Array9[i];}
    else continue;
}
int ThirdCentralCell=Maxcell1;
E3=en_Maxcell1 ;

//quarta
double en_Maxcell2=0.;
int Maxcell2=0.;

for(int i=0; i<9; ++i)
{
    if(en_c[Array9[i]]>en_Maxcell2 && Array9[i]!=CentralCell && Array9[i]!=SecondCentralCell_in9 && Array9[i]!=ThirdCentralCell){en_Maxcell2=en_c[Array9[i]]; Maxcell2=Array9[i];}
    else continue;
}
int FourthCentralCell=Maxcell2;
E4=en_Maxcell2 ;

    
E9=E_1/E_clus3x3;           
Eout=(Etotcal-E_clus3x3)/E_clus3x3;
Emean_out=(Etotcal-E_clus3x3)/16;    
    
double Ex=0.;
double Ey=0.;
double wtot=0.;

for(int i=0; i<9; ++i)
{
double xx = myGrid->GetXaxis()->GetBinCenter(Rev_numberX[Array9[i]]);
double yy = myGrid->GetYaxis()->GetBinCenter(Rev_numberY[Array9[i]]);

double w=4.0+log(en_c[Array9[i]]/E_clus3x3);
double wi=(w>0)?w:0; 
    
Ex+=wi*xx;
Ey+=wi*yy;
wtot+=wi;
}    
    
double centroidX=(Ex)/wtot;
double centroidY=(Ey)/wtot;         

double r_mue=sqrt((detKinBeamRot_cooXe-detKinBeamRot_cooXmu)*(detKinBeamRot_cooXe-detKinBeamRot_cooXmu)+(detKinBeamRot_cooYe-detKinBeamRot_cooYmu)*(detKinBeamRot_cooYe-detKinBeamRot_cooYmu)); 
 
double ddd=sqrt((centroidX-detKinBeamRot_cooXe)*(centroidX-detKinBeamRot_cooXe)+(centroidY-detKinBeamRot_cooYe)*(centroidY-detKinBeamRot_cooYe)); 
    
double x=detKinBeamRot_def_angle_e*0.001;    
double beta=(sqrt(150*150-(m_mu*m_mu))/(150+m_e));
double Ethe=m_e*((1+(beta*beta*cos(x)*cos(x)))/(1-(beta*beta*cos(x)*cos(x))));
//cout << " E teorico " << Ethe  << " E reco " << E_clus3x3 << endl;
//double der_Ee=-4*m_e*beta*beta*( (cos(x)*sin(x))/((1-(beta*beta*cos(x)*cos(x)))*(1-(beta*beta*cos(x)*cos(x)))) );
double DE= (E_clus3x3/Ethe-0.955);
//double DE5= (E_clus3x3-(95.4/100)*Ethe)/(sqrt(der_Ee*der_Ee+1));
// /(sqrt(der_Ee*der_Ee+1));
// cout << " DE " << DE << " E_clus3x3-Ethe " << E_clus3x3-Ethe << endl;
    
/*Elastic->SetPoint(n,x,Ethe);
++n;*/

    
if(r_mu<5 && detKinBeamRot_tar==1) 
{   E3x31CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
Th1->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);}

if(r_mu<1.7 && detKinBeamRot_tar==1 ){ //E_clus3x3>1

E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);

   /* if (DE>-0.08 && DE<0.08 && detKinBeamRot_def_angle_e<=5) 
    {E3x31CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
    Th1->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);}   
    else if (DE>-0.2 && DE<0.2 && detKinBeamRot_def_angle_e>5 && detKinBeamRot_def_angle_e<=10) 
    {E3x31CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th1->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);} else if (DE>-0.3 && DE<0.3 && detKinBeamRot_def_angle_e>10 && detKinBeamRot_def_angle_e<=15) {E3x31CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th1->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);} else if ( DE>-0.4 && DE<0.4 && detKinBeamRot_def_angle_e>15 && detKinBeamRot_def_angle_e<=20) 
    {E3x31CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th1->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);}     
    else if (detKinBeamRot_def_angle_e>20) {E3x31CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th1->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);} */

    
    
    
/*//ZONA 1
    
if(E_clus3x3<=10){

    if(ddd>0.1*9)
    {
 
        DeltaRCUT->Fill(ddd,wgt_full);
        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full); diff_r_pheCUT->Fill(d_e_ph,wgt_full);}
                
            EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
            TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
            EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
            ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); 
            diff_r_mueCUT->Fill(r_mue,wgt_full);  
            
            
    } else {    
    
           if (DE>-0.08 && DE<0.08 && detKinBeamRot_def_angle_e<=5) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
    Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }   
    else if (DE>-0.2 && DE<0.2 && detKinBeamRot_def_angle_e>5 && detKinBeamRot_def_angle_e<=10) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if (DE>-0.3 && DE<0.3 && detKinBeamRot_def_angle_e>10 && detKinBeamRot_def_angle_e<=15) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if ( DE>-0.4 && DE<0.4 && detKinBeamRot_def_angle_e>15 && detKinBeamRot_def_angle_e<=20) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }     
    else if (detKinBeamRot_def_angle_e>20) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }
        else{        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full); diff_r_pheCUT->Fill(d_e_ph,wgt_full);}
                
            EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
            TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
            EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
            ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); 
            diff_r_mueCUT->Fill(r_mue,wgt_full);  }
 }    
                   
              }

 

       
//ZONA 2
    
if(E_clus3x3>10 && E_clus3x3<=20){
if(ddd>0.1*7.3)
    {
 
        DeltaRCUT->Fill(ddd,wgt_full);
        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full);
            diff_r_pheCUT->Fill(d_e_ph,wgt_full);}
                EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
               TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
                EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
               ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); 
    diff_r_mueCUT->Fill(r_mue,wgt_full); 
            

    } else { 
    
              if (DE>-0.08 && DE<0.08 && detKinBeamRot_def_angle_e<=5) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
    Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }   
    else if (DE>-0.2 && DE<0.2 && detKinBeamRot_def_angle_e>5 && detKinBeamRot_def_angle_e<=10) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if (DE>-0.3 && DE<0.3 && detKinBeamRot_def_angle_e>10 && detKinBeamRot_def_angle_e<=15) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if ( DE>-0.4 && DE<0.4 && detKinBeamRot_def_angle_e>15 && detKinBeamRot_def_angle_e<=20) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }     
    else if (detKinBeamRot_def_angle_e>20) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }
            else{        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full); diff_r_pheCUT->Fill(d_e_ph,wgt_full);}
                
            EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
            TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
            EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
            ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); 
            diff_r_mueCUT->Fill(r_mue,wgt_full);  }
                   
}
}
//ZONA 3
    
if(E_clus3x3>20 && E_clus3x3<=30){
if(ddd>0.1*6.4)
    {
 
        DeltaRCUT->Fill(ddd,wgt_full);
        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full);
            diff_r_pheCUT->Fill(d_e_ph,wgt_full);} 
                EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
               TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
                EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
               ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); diff_r_mueCUT->Fill(r_mue,wgt_full); 

    }else {    
    
                 if (DE>-0.08 && DE<0.08 && detKinBeamRot_def_angle_e<=5) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
    Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }   
    else if (DE>-0.2 && DE<0.2 && detKinBeamRot_def_angle_e>5 && detKinBeamRot_def_angle_e<=10) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if (DE>-0.3 && DE<0.3 && detKinBeamRot_def_angle_e>10 && detKinBeamRot_def_angle_e<=15) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if ( DE>-0.4 && DE<0.4 && detKinBeamRot_def_angle_e>15 && detKinBeamRot_def_angle_e<=20) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }     
    else if (detKinBeamRot_def_angle_e>20) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } 
            else{        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full); diff_r_pheCUT->Fill(d_e_ph,wgt_full);}
                
            EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
            TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
            EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
            ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); 
            diff_r_mueCUT->Fill(r_mue,wgt_full);  }
    

}
}
    
//ZONA 4
    
if(E_clus3x3>30 && E_clus3x3<=40){
if(ddd>0.1*6)
    {
 
        DeltaRCUT->Fill(ddd,wgt_full);
        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full);
            diff_r_pheCUT->Fill(d_e_ph,wgt_full);} 
                EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
               TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
                EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
               ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); diff_r_mueCUT->Fill(r_mue,wgt_full); 

    } else {   
    
             if (DE>-0.08 && DE<0.08 && detKinBeamRot_def_angle_e<=5) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
    Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }   
    else if (DE>-0.2 && DE<0.2 && detKinBeamRot_def_angle_e>5 && detKinBeamRot_def_angle_e<=10) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if (DE>-0.3 && DE<0.3 && detKinBeamRot_def_angle_e>10 && detKinBeamRot_def_angle_e<=15) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if ( DE>-0.4 && DE<0.4 && detKinBeamRot_def_angle_e>15 && detKinBeamRot_def_angle_e<=20) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }     
    else if (detKinBeamRot_def_angle_e>20) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } 
            else{        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full); diff_r_pheCUT->Fill(d_e_ph,wgt_full);}
                
            EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
            TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
            EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
            ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); 
            diff_r_mueCUT->Fill(r_mue,wgt_full);  }

}
}
    
//ZONA 5
    
if(E_clus3x3>40 && E_clus3x3<=50){
 if(ddd>0.1*5.8)
    {
 
        DeltaRCUT->Fill(ddd,wgt_full);
        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full);
            diff_r_pheCUT->Fill(d_e_ph,wgt_full);} 
                    EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
               TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
                EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
               ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); diff_r_mueCUT->Fill(r_mue,wgt_full); 
            

    } else {    
    
              if (DE>-0.08 && DE<0.08 && detKinBeamRot_def_angle_e<=5) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
    Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }   
    else if (DE>-0.2 && DE<0.2 && detKinBeamRot_def_angle_e>5 && detKinBeamRot_def_angle_e<=10) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if (DE>-0.3 && DE<0.3 && detKinBeamRot_def_angle_e>10 && detKinBeamRot_def_angle_e<=15) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if ( DE>-0.4 && DE<0.4 && detKinBeamRot_def_angle_e>15 && detKinBeamRot_def_angle_e<=20) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }     
    else if (detKinBeamRot_def_angle_e>20) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } 
        else{        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full); diff_r_pheCUT->Fill(d_e_ph,wgt_full);}
                
            EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
            TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
            EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
            ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); 
            diff_r_mueCUT->Fill(r_mue,wgt_full);  }
     
 }
} 

//ZONA 6
    
if(E_clus3x3>50 && E_clus3x3<=60){
 if(ddd>0.1*5.5)
    {
 
        DeltaRCUT->Fill(ddd,wgt_full);
        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full);
            diff_r_pheCUT->Fill(d_e_ph,wgt_full);} 
                EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
               TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
                EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
               ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); diff_r_mueCUT->Fill(r_mue,wgt_full); 
            

    } else{
     
    if (DE>-0.08 && DE<0.08 && detKinBeamRot_def_angle_e<=5) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
    Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }   
    else if (DE>-0.2 && DE<0.2 && detKinBeamRot_def_angle_e>5 && detKinBeamRot_def_angle_e<=10) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if (DE>-0.3 && DE<0.3 && detKinBeamRot_def_angle_e>10 && detKinBeamRot_def_angle_e<=15) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if ( DE>-0.4 && DE<0.4 && detKinBeamRot_def_angle_e>15 && detKinBeamRot_def_angle_e<=20) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }     
    else if (detKinBeamRot_def_angle_e>20) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } 
             else{        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full); diff_r_pheCUT->Fill(d_e_ph,wgt_full);}
                
            EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
            TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
            EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
            ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); 
            diff_r_mueCUT->Fill(r_mue,wgt_full);  }
     
    
}
}
    
//ZONA 7
    
if(E_clus3x3>60 && E_clus3x3<=75){
  if(ddd>0.1*5.12)
    {
 
        DeltaRCUT->Fill(ddd,wgt_full);
        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full);
            diff_r_pheCUT->Fill(d_e_ph,wgt_full);} 
                    EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
               TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
                EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
               ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); diff_r_mueCUT->Fill(r_mue,wgt_full); 
            

    } else {   
    
             if (DE>-0.08 && DE<0.08 && detKinBeamRot_def_angle_e<=5) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
    Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }   
    else if (DE>-0.2 && DE<0.2 && detKinBeamRot_def_angle_e>5 && detKinBeamRot_def_angle_e<=10) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if (DE>-0.3 && DE<0.3 && detKinBeamRot_def_angle_e>10 && detKinBeamRot_def_angle_e<=15) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if ( DE>-0.4 && DE<0.4 && detKinBeamRot_def_angle_e>15 && detKinBeamRot_def_angle_e<=20) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }     
    else if (detKinBeamRot_def_angle_e>20) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } 
              else{        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full); diff_r_pheCUT->Fill(d_e_ph,wgt_full);}
                
            EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
            TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
            EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
            ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); 
            diff_r_mueCUT->Fill(r_mue,wgt_full);  }
  
} 
}
//ZONA 8
    
if(E_clus3x3>75){
  if(ddd>0.1*4.8)
    {
 
        DeltaRCUT->Fill(ddd,wgt_full);
        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full);
            diff_r_pheCUT->Fill(d_e_ph,wgt_full);} 
                EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
               TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
                EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
               ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); diff_r_mueCUT->Fill(r_mue,wgt_full); 

    } else {  
    
         if (DE>-0.08 && DE<0.08 && detKinBeamRot_def_angle_e<=5) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
    Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }   
    else if (DE>-0.2 && DE<0.2 && detKinBeamRot_def_angle_e>5 && detKinBeamRot_def_angle_e<=10) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if (DE>-0.3 && DE<0.3 && detKinBeamRot_def_angle_e>10 && detKinBeamRot_def_angle_e<=15) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } else if ( DE>-0.4 && DE<0.4 && detKinBeamRot_def_angle_e>15 && detKinBeamRot_def_angle_e<=20) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     }     
    else if (detKinBeamRot_def_angle_e>20) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);  if(photon_energy!=-1 && n_cell_ph!=0){             EphoutCUTafter->Fill(photon_energy,wgt_full);             ThphoutCUTafter->Fill(photon_def_angle_ph,wgt_full);             diff_th_pheCUTafter->Fill(diffTh,wgt_full);             diff_r_pheCUTafter->Fill(d_e_ph,wgt_full);}                          EeoutCUTafter->Fill(detKinBeamRot_Ee,wgt_full);              TheoutCUTafter->Fill(detKinBeamRot_def_angle_e,wgt_full); EmuoutCUTafter->Fill(detKinBeamRot_Emu,wgt_full);                 ThmuoutCUTafter->Fill(detKinBeamRot_def_angle_mu,wgt_full);     } 
              else{        if(photon_energy!=-1 && n_cell_ph!=0){
            EphoutCUT->Fill(photon_energy,wgt_full);
            ThphoutCUT->Fill(photon_def_angle_ph,wgt_full);
            diff_th_pheCUT->Fill(diffTh,wgt_full); diff_r_pheCUT->Fill(d_e_ph,wgt_full);}
                
            EeoutCUT->Fill(detKinBeamRot_Ee,wgt_full); 
            TheoutCUT->Fill(detKinBeamRot_def_angle_e,wgt_full); 
            EmuoutCUT->Fill(detKinBeamRot_Emu,wgt_full); 
            ThmuoutCUT->Fill(detKinBeamRot_def_angle_mu,wgt_full); 
            diff_r_mueCUT->Fill(r_mue,wgt_full);  }
       }  
   

}  */
    
    
 
 /*if(photon_energy==-1 && n_cell_ph==0)   
 {if (DE>abs(0.08) && detKinBeamRot_def_angle_e<=5) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
    Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);}   
    else if ( DE>abs(0.2) && detKinBeamRot_def_angle_e>5 && detKinBeamRot_def_angle_e<=10) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);} else if ( DE>abs(0.3) && detKinBeamRot_def_angle_e>10 && detKinBeamRot_def_angle_e<=15) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);} else if (  DE>abs(0.4) && detKinBeamRot_def_angle_e>15 && detKinBeamRot_def_angle_e<=20) 
    {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
     Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);}     } */
  
 }      
    
    
/*if (photon_energy>1 || photon_energy<0) {E3x32CUT->Fill(detKinBeamRot_def_angle_e*0.001,E_clus3x3,wgt_full);
Th2->Fill(detKinBeamRot_def_angle_e*0.001,detKinBeamRot_def_angle_mu*0.001,wgt_full);}
Ephout->Fill(photon_energy,wgt_full);
Thphout->Fill(photon_def_angle_ph,wgt_full);
diff_th_phe->Fill(diffTh,wgt_full);
diff_r_phe->Fill(d_e_ph,wgt_full);
 hist_distCUT->Fill(ddd,wgt_full);
 double eq=3*sqrt((5.78/sqrt(E_clus3x3))*(5.78/sqrt(E_clus3x3))+1.095*1.095)*0.1;*/

   
}

}       

delete myGrid; 
}}


TCanvas * el= new TCanvas("el","el",1000,100,2500,2000);
TF1 *Elastic = new TF1("Elastic","0.5109989461*0.001*((1+(sqrt(150*150-(105.6583745 *0.001*105.6583745 *0.001))/(150+0.5109989461*0.001))*(sqrt(150*150-(105.6583745 *0.001*105.6583745 *0.001))/(150+0.5109989461*0.001))*cos(x)*cos(x))/(1-(sqrt(150*150-(105.6583745 *0.001*105.6583745 *0.001))/(150+0.5109989461*0.001))*(sqrt(150*150-(105.6583745 *0.001*105.6583745 *0.001))/(150+0.5109989461*0.001))*cos(x)*cos(x)))",0,0.030); 
Elastic->SetLineColor(kViolet);

TF1 *Elastic2 = new TF1("Elastic2","asin( (sin(x)*sqrt(Elastic(x)*Elastic(x)-0.5109989461*0.001*0.5109989461*0.001))/sqrt( (150+0.5109989461*0.001-Elastic(x))*(150+0.5109989461*0.001-Elastic(x))-105.6583745 *0.001*105.6583745 *0.001 ) )",0,0.030); 
Elastic2->SetLineColor(kViolet);
    
Int_t nx1CUT = E3x31CUT->GetNbinsX();
Int_t ny1CUT = E3x31CUT->GetNbinsY();
for (Int_t i=1; i<nx1CUT+1; i++) {
for (Int_t j=1; j<ny1CUT+1; j++) {
if (E3x31CUT->GetBinContent(i,j)<1) E3x31CUT->SetBinContent(i,j,0);}}
    
Int_t nx2CUT = E3x32CUT->GetNbinsX();
Int_t ny2CUT = E3x32CUT->GetNbinsY();
for (Int_t i=1; i<nx2CUT+1; i++) {
for (Int_t j=1; j<ny2CUT+1; j++) {
if (E3x32CUT->GetBinContent(i,j)<1) E3x32CUT->SetBinContent(i,j,0);}}
    
Int_t nxTh = Th1->GetNbinsX();
Int_t nyTh = Th1->GetNbinsY();
for (Int_t i=1; i<nxTh+1; i++) {
for (Int_t j=1; j<nyTh+1; j++) {
if (Th1->GetBinContent(i,j)<1) Th1->SetBinContent(i,j,0);}}
    
Int_t nx2thcut = Th2->GetNbinsX();
Int_t ny2thcut = Th2->GetNbinsY();
for (Int_t i=1; i<nx2thcut+1; i++) {
for (Int_t j=1; j<ny2thcut+1; j++) {
if (Th2->GetBinContent(i,j)<1) Th2->SetBinContent(i,j,0);}}
           
    
    
TCanvas * thu= new TCanvas("thu","thu",1000,1000,2500,2000);
thu->Divide(2,2);
gStyle->SetPalette(kRainBow);
 
thu->cd(1);   
E3x31CUT->GetXaxis()->SetTitle("Theta_el[rad]");
E3x31CUT->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x31CUT->Draw("COLZ");
gStyle->SetOptStat("ne"); 
Elastic->Draw("same");
thu->cd(2);   
E3x32CUT->GetXaxis()->SetTitle("Theta_el[rad]");
E3x32CUT->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x32CUT->Draw("COLZ");
gStyle->SetOptStat("ne"); 
Elastic->Draw("same");
/*c4a->SaveAs("/home/LHCB-T3/espedicato/tesi/studio4/thE_cut.png");
c4a->SaveAs("/home/LHCB-T3/espedicato/tesi/studio4/thE_cut.C");*/

    
/*TCanvas * thu= new TCanvas("c4a","c4a",100,100,2500,2000);
thu->Divide(2,1);
gStyle->SetPalette(kRainBow);*/
 
thu->cd(3);   
Th1->GetXaxis()->SetTitle("Theta_el[rad]");
Th1->GetYaxis()->SetTitle("Theta_mu[rad]");
Th1->LabelsDeflate("Y");
Th1->Draw("COLZ");
Elastic2->Draw("same");
gStyle->SetOptStat("ne"); 
thu->cd(4);   
Th2->GetXaxis()->SetTitle("Theta_el[rad]");
Th2->GetYaxis()->SetTitle("Theta_mu[rad]");
Th2->LabelsDeflate("Y");
Th2->Draw("COLZ");
Elastic2->Draw("same");
gStyle->SetOptStat("ne"); 
thu->SaveAs("/home/LHCB-T3/espedicato/tesi/studio4/thu.png");
thu->SaveAs("/home/LHCB-T3/espedicato/tesi/studio4/thu.C");

/*
    
TCanvas * cres= new TCanvas("cres","cres",1000,100,2500,2000);  
cres->Divide(1,2);
cres->cd(1);    
diff_r_mue->GetXaxis()->SetTitle("r [cm]");
diff_r_mue->SetLineWidth(3);
diff_r_mue->Draw("HIST");
diff_r_mueCUT->SetLineWidth(3);
diff_r_mueCUT->SetLineColor(kRed-4);
diff_r_mueCUT->Draw("HIST same"); 
cres->cd(2); 
DeltaR->GetXaxis()->SetTitle("r [cm]");
DeltaR->SetLineWidth(3);
DeltaR->Draw("HIST");
DeltaRCUT->SetLineWidth(3);
DeltaRCUT->SetLineColor(kRed-4);
DeltaRCUT->Draw("HIST same"); 
cres->SaveAs("/home/LHCB-T3/espedicato/tesi/studio4/res.png");
   
TCanvas * c= new TCanvas("c","c",1000,100,2500,2000);

hist_E92_e->GetXaxis()->SetTitle("<Eres> [GeV]");
hist_E92_e->SetLineWidth(3);
hist_E92_e->Draw("HIST"); 
hist_E92_e->SetMinimum(1);
gPad->SetLogy();
    
hist_E92_eCUT->GetXaxis()->SetTitle("<Eres> [GeV]");
hist_E92_eCUT->SetLineWidth(3);
hist_E92_eCUT->SetLineColor(kRed);
hist_E92_eCUT->Draw("HIST same");   
gPad->BuildLegend(0.25,0.15,0.25,0.15);
c->SaveAs("/home/LHCB-T3/espedicato/tesi/studio4/meanout.png");*/

} 