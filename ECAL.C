#define atree_cxx
#include "next.h"
#include <TH2.h>
#include <TH1.h>
#include <TGraph.h>
#include <TTree.h>
#include <cmath>
using namespace std;

#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{
    TH1::SetDefaultSumw2();
    
Int_t n_cell; //numero di cella in cui cade l'ELETTRONE
Int_t n_cell_ph; //numero di cella in cui cade il fotone
Int_t n_tot=0;

Double_t same_cell=0.;
Double_t different_cell=0.;
Double_t E_CAL;
Double_t Rm = 1.959 ; //raggio di Moliere in centimetri    

    
TF2 *fxy_e = new TF2("fxy_e","(1/(2*pi*[0]*[1]))*(exp(-((x-[2])*(x-[2])+(y-[3])*(y-[3]))/(2*[0]*[1])))",-7.125,7.125,-7.125,7.125 );
    
TF2 *fxy_ph = new TF2("fxy_ph","(1/(2*pi*[0]*[1]))*(exp(-((x-[2])*(x-[2])+(y-[3])*(y-[3]))/(2*[0]*[1])))",-7.125,7.125,-7.125,7.125 );

    
TMatrixD E(3,25);

    
    if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();


    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
        detKinBeamRot_cooXe=detKinBeamRot_cooXe*100;
        detKinBeamRot_cooYe=detKinBeamRot_cooYe*100;
        detKinBeamRot_cooXmu=detKinBeamRot_cooXmu*100;
        detKinBeamRot_cooYmu=detKinBeamRot_cooYmu*100;
        photon_coox=photon_coox*100;
        photon_cooy=photon_cooy*100;
       
        if (photon_coox!=-1 && photon_cooy!=-1)
    {  
     E_CAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {  E_CAL=detKinBeamRot_Ee;}
       
       
    Double_t d_e_ph=sqrt( (detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy) ); 
       
    Double_t d_e_mu=sqrt( (detKinBeamRot_cooXe-detKinBeamRot_cooXmu)*(detKinBeamRot_cooXe-detKinBeamRot_cooXmu)+(detKinBeamRot_cooYe-detKinBeamRot_cooYmu)*(detKinBeamRot_cooYe-detKinBeamRot_cooYmu) );
       
   
       
       
//if (E_CAL>1){
if (abs(detKinBeamRot_cooXe)<7.125 && abs(detKinBeamRot_cooYe)<7.125)
    
{ 
    fxy_e->SetParameter(0,Rm);
    fxy_e->SetParameter(1,Rm);
    fxy_e->SetParameter(2,detKinBeamRot_cooXe);
    fxy_e->SetParameter(3,detKinBeamRot_cooYe);
    
    if (detKinBeamRot_cooXe>-7.125 && detKinBeamRot_cooXe<-4.275 && detKinBeamRot_cooYe<7.125 && detKinBeamRot_cooYe>4.275) {n_cell=1;}
    
    if (detKinBeamRot_cooXe>-4.275 && detKinBeamRot_cooXe<-1.425 && detKinBeamRot_cooYe<7.125 && detKinBeamRot_cooYe>4.275) {n_cell=2;}

    if (detKinBeamRot_cooXe>-1.425 && detKinBeamRot_cooXe<1.425 && detKinBeamRot_cooYe<7.125 && detKinBeamRot_cooYe>4.275) {n_cell=3;}

    if (detKinBeamRot_cooXe>1.425 && detKinBeamRot_cooXe<4.275 && detKinBeamRot_cooYe<7.125 && detKinBeamRot_cooYe>4.275) {n_cell=4;}

    if (detKinBeamRot_cooXe>4.275 && detKinBeamRot_cooXe<7.125 && detKinBeamRot_cooYe<7.125 && detKinBeamRot_cooYe>4.275) {n_cell=5;}
    
    
    if (detKinBeamRot_cooXe>-7.125 && detKinBeamRot_cooXe<-4.275 && detKinBeamRot_cooYe<4.275 && detKinBeamRot_cooYe>1.425) {n_cell=6;}
    
    if (detKinBeamRot_cooXe>-4.275 && detKinBeamRot_cooXe<-1.425 && detKinBeamRot_cooYe<4.275 && detKinBeamRot_cooYe>1.425) {n_cell=7;}

    if (detKinBeamRot_cooXe>-1.425 && detKinBeamRot_cooXe<1.425 && detKinBeamRot_cooYe<4.275 && detKinBeamRot_cooYe>1.425) {n_cell=8;}

    if (detKinBeamRot_cooXe>1.425 && detKinBeamRot_cooXe<4.275 && detKinBeamRot_cooYe<4.275 && detKinBeamRot_cooYe>1.425) {n_cell=9;}

    if (detKinBeamRot_cooXe>4.275 && detKinBeamRot_cooXe<7.125 && detKinBeamRot_cooYe<4.275 && detKinBeamRot_cooYe>1.425) {n_cell=10;}
    

    
    if (detKinBeamRot_cooXe>-7.125 && detKinBeamRot_cooXe<-4.275 && detKinBeamRot_cooYe<1.425 && detKinBeamRot_cooYe>-1.425) {n_cell=11;}
    
    if (detKinBeamRot_cooXe>-4.275 && detKinBeamRot_cooXe<-1.425 && detKinBeamRot_cooYe<1.425 && detKinBeamRot_cooYe>-1.425) {n_cell=12;}

    if (detKinBeamRot_cooXe>-1.425 && detKinBeamRot_cooXe<1.425 && detKinBeamRot_cooYe<1.425 && detKinBeamRot_cooYe>-1.425) {n_cell=13;}

    if (detKinBeamRot_cooXe>1.425 && detKinBeamRot_cooXe<4.275 && detKinBeamRot_cooYe<1.425 && detKinBeamRot_cooYe>-1.425) {n_cell=14;}

    if (detKinBeamRot_cooXe>4.275 && detKinBeamRot_cooXe<7.125 && detKinBeamRot_cooYe<1.425 && detKinBeamRot_cooYe>-1.425) {n_cell=15;}
    
    
    
    
    if (detKinBeamRot_cooXe>-7.125 && detKinBeamRot_cooXe<-4.275 && detKinBeamRot_cooYe<-1.425 && detKinBeamRot_cooYe>-4.275) {n_cell=16;}
    
    if (detKinBeamRot_cooXe>-4.275 && detKinBeamRot_cooXe<-1.425 && detKinBeamRot_cooYe<-1.425 && detKinBeamRot_cooYe>-4.275) {n_cell=17;}

    if (detKinBeamRot_cooXe>-1.425 && detKinBeamRot_cooXe<1.425 && detKinBeamRot_cooYe<-1.425 && detKinBeamRot_cooYe>-4.275) {n_cell=18;}

    if (detKinBeamRot_cooXe>1.425 && detKinBeamRot_cooXe<4.275 && detKinBeamRot_cooYe<-1.425 && detKinBeamRot_cooYe>-4.275) {n_cell=19;}

    if (detKinBeamRot_cooXe>4.275 && detKinBeamRot_cooXe<7.125 && detKinBeamRot_cooYe<-1.425 && detKinBeamRot_cooYe>-4.275) {n_cell=20;}
    
    
    
    if (detKinBeamRot_cooXe>-7.125 && detKinBeamRot_cooXe<-4.275 && detKinBeamRot_cooYe<-4.275 && detKinBeamRot_cooYe>-7.125) {n_cell=21;}
    
    if (detKinBeamRot_cooXe>-4.275 && detKinBeamRot_cooXe<-1.425 && detKinBeamRot_cooYe<-4.275 && detKinBeamRot_cooYe>-7.125) {n_cell=22;}

    if (detKinBeamRot_cooXe>-1.425 && detKinBeamRot_cooXe<1.425 && detKinBeamRot_cooYe<-4.275 && detKinBeamRot_cooYe>-7.125) {n_cell=23;}

    if (detKinBeamRot_cooXe>1.425 && detKinBeamRot_cooXe<4.275 && detKinBeamRot_cooYe<-4.275 && detKinBeamRot_cooYe>-7.125) {n_cell=24;}

    if (detKinBeamRot_cooXe>4.275 && detKinBeamRot_cooXe<7.125 && detKinBeamRot_cooYe<-4.275 && detKinBeamRot_cooYe>-7.125) {n_cell=25;}
    
    cout << "cella elettrone:" << n_cell << endl; 
       
E[0][0]= detKinBeamRot_Ee*fxy_e->Integral(-7.125,-4.275,4.275,7.125); //cella 1
E[0][1]= detKinBeamRot_Ee*fxy_e->Integral(-4.275,-1.425,4.275,7.125); //cella 2
E[0][2]= detKinBeamRot_Ee*fxy_e->Integral(-1.425,1.425,4.275,7.125); //cella 3
E[0][3]= detKinBeamRot_Ee*fxy_e->Integral(1.425,4.275,4.275,7.125); //cella 4
E[0][4]= detKinBeamRot_Ee*fxy_e->Integral(4.275,7.125,4.275,7.125); //cella 5
E[0][5]= detKinBeamRot_Ee*fxy_e->Integral(-7.125,-4.275,1.425,4.275); //cella 6
E[0][6]= detKinBeamRot_Ee*fxy_e->Integral(-4.275,-1.425,1.425,4.275); //cella 7
E[0][7]= detKinBeamRot_Ee*fxy_e->Integral(-1.425,1.425,1.425,4.275); //cella 8
E[0][8]= detKinBeamRot_Ee*fxy_e->Integral(1.425,4.275,1.425,4.275); //cella 9
E[0][9]= detKinBeamRot_Ee*fxy_e->Integral(4.275,7.125,1.425,4.275); //cella 10
E[0][10]= detKinBeamRot_Ee*fxy_e->Integral(-7.125,-4.275,-1.425,1.425); //cella 11
E[0][11]= detKinBeamRot_Ee*fxy_e->Integral(-4.275,-1.425,-1.425,1.425); //cella 12
E[0][12]= detKinBeamRot_Ee*fxy_e->Integral(-1.425,1.425,-1.425,1.425); //cella 13
E[0][13]= detKinBeamRot_Ee*fxy_e->Integral(1.425,4.275,-1.425,1.425); //cella 14
E[0][14]= detKinBeamRot_Ee*fxy_e->Integral(4.275,7.125,-1.425,1.425); //cella 15
E[0][15]= detKinBeamRot_Ee*fxy_e->Integral(-7.125,-4.275,-4.275,-1.425); //cella 16
E[0][16]= detKinBeamRot_Ee*fxy_e->Integral(-4.275,-1.425,-4.275,-1.425); //cella 17
E[0][17]= detKinBeamRot_Ee*fxy_e->Integral(-1.425,1.425,-4.275,-1.425); //cella 18
E[0][18]= detKinBeamRot_Ee*fxy_e->Integral(1.425,4.275,-4.275,-1.425); //cella 19
E[0][19]= detKinBeamRot_Ee*fxy_e->Integral(4.275,7.125,-4.275,-1.425); //cella 20
E[0][20]= detKinBeamRot_Ee*fxy_e->Integral(-7.125,-4.275,-7.125,-4.275); //cella 21
E[0][21]= detKinBeamRot_Ee*fxy_e->Integral(-4.275,-1.425,-7.125,-4.275); //cella 22
E[0][22]= detKinBeamRot_Ee*fxy_e->Integral(-1.425,1.425,-7.125,-4.275); //cella 23
E[0][23]= detKinBeamRot_Ee*fxy_e->Integral(1.425,4.275,-7.125,-4.275); //cella 24
E[0][24]= detKinBeamRot_Ee*fxy_e->Integral(4.275,7.125,-7.125,-4.275); //cella 25
       
if (abs(photon_coox)<7.125 && abs(photon_cooy)<7.125 && photon_energy>0.2)
    
{   fxy_ph->SetParameter(0,Rm);
    fxy_ph->SetParameter(1,Rm);
    fxy_ph->SetParameter(2,photon_coox);
    fxy_ph->SetParameter(3,photon_cooy);
    
    if (photon_coox>-7.125 && photon_coox<-4.275 && photon_cooy<7.125 && photon_cooy>4.275) {n_cell_ph=1;}
    
    if (photon_coox>-4.275 && photon_coox<-1.425 && photon_cooy<7.125 && photon_cooy>4.275) {n_cell_ph=2;}

    if (photon_coox>-1.425 && photon_coox<1.425 && photon_cooy<7.125 && photon_cooy>4.275) {n_cell_ph=3;}

    if (photon_coox>1.425 && photon_coox<4.275 && photon_cooy<7.125 && photon_cooy>4.275) {n_cell_ph=4;}

    if (photon_coox>4.275 && photon_coox<7.125 && photon_cooy<7.125 && photon_cooy>4.275) {n_cell_ph=5;}
    
    
    if (photon_coox>-7.125 && photon_coox<-4.275 && photon_cooy<4.275 && photon_cooy>1.425) {n_cell_ph=6;}
    
    if (photon_coox>-4.275 && photon_coox<-1.425 && photon_cooy<4.275 && photon_cooy>1.425) {n_cell_ph=7;}

    if (photon_coox>-1.425 && photon_coox<1.425 && photon_cooy<4.275 && photon_cooy>1.425) {n_cell_ph=8;}

    if (photon_coox>1.425 && photon_coox<4.275 && photon_cooy<4.275 && photon_cooy>1.425) {n_cell_ph=9;}

    if (photon_coox>4.275 && photon_coox<7.125 && photon_cooy<4.275 && photon_cooy>1.425) {n_cell_ph=10;}
    

    
    if (photon_coox>-7.125 && photon_coox<-4.275 && photon_cooy<1.425 && photon_cooy>-1.425) {n_cell_ph=11;}
    
    if (photon_coox>-4.275 && photon_coox<-1.425 && photon_cooy<1.425 && photon_cooy>-1.425) {n_cell_ph=12;}

    if (photon_coox>-1.425 && photon_coox<1.425 && photon_cooy<1.425 && photon_cooy>-1.425) {n_cell_ph=13;}

    if (photon_coox>1.425 && photon_coox<4.275 && photon_cooy<1.425 && photon_cooy>-1.425) {n_cell_ph=14;}

    if (photon_coox>4.275 && photon_coox<7.125 && photon_cooy<1.425 && photon_cooy>-1.425) {n_cell_ph=15;}
    
    
    
    
    if (photon_coox>-7.125 && photon_coox<-4.275 && photon_cooy<-1.425 && photon_cooy>-4.275) {n_cell_ph=16;}
    
    if (photon_coox>-4.275 && photon_coox<-1.425 && photon_cooy<-1.425 && photon_cooy>-4.275) {n_cell_ph=17;}

    if (photon_coox>-1.425 && photon_coox<1.425 && photon_cooy<-1.425 && photon_cooy>-4.275) {n_cell_ph=18;}

    if (photon_coox>1.425 && photon_coox<4.275 && photon_cooy<-1.425 && photon_cooy>-4.275) {n_cell_ph=19;}

    if (photon_coox>4.275 && photon_coox<7.125 && photon_cooy<-1.425 && photon_cooy>-4.275) {n_cell_ph=20;}
    
    
    
    if (photon_coox>-7.125 && photon_coox<-4.275 && photon_cooy<-4.275 && photon_cooy>-7.125) {n_cell_ph=21;}
    
    if (photon_coox>-4.275 && photon_coox<-1.425 && photon_cooy<-4.275 && photon_cooy>-7.125) {n_cell_ph=22;}

    if (photon_coox>-1.425 && photon_coox<1.425 && photon_cooy<-4.275 && photon_cooy>-7.125) {n_cell_ph=23;}

    if (photon_coox>1.425 && photon_coox<4.275 && photon_cooy<-4.275 && photon_cooy>-7.125) {n_cell_ph=24;}

    if (photon_coox>4.275 && photon_coox<7.125 && photon_cooy<-4.275 && photon_cooy>-7.125) {n_cell_ph=25;} 
 
 cout << "cella fotone:" << n_cell_ph << endl; 
    
E[1][0]= photon_energy*fxy_ph->Integral(-7.125,-4.275,4.275,7.125); //cella 1
E[1][1]= photon_energy*fxy_ph->Integral(-4.275,-1.425,4.275,7.125); //cella 2
E[1][2]= photon_energy*fxy_ph->Integral(-1.425,1.425,4.275,7.125); //cella 3
E[1][3]= photon_energy*fxy_ph->Integral(1.425,4.275,4.275,7.125); //cella 4
E[1][4]= photon_energy*fxy_ph->Integral(4.275,7.125,4.275,7.125); //cella 5
E[1][5]= photon_energy*fxy_ph->Integral(-7.125,-4.275,1.425,4.275); //cella 6
E[1][6]= photon_energy*fxy_ph->Integral(-4.275,-1.425,1.425,4.275); //cella 7
E[1][7]= photon_energy*fxy_ph->Integral(-1.425,1.425,1.425,4.275); //cella 8
E[1][8]= photon_energy*fxy_ph->Integral(1.425,4.275,1.425,4.275); //cella 9
E[1][9]= photon_energy*fxy_ph->Integral(4.275,7.125,1.425,4.275); //cella 10
E[1][10]= photon_energy*fxy_ph->Integral(-7.125,-4.275,-1.425,1.425); //cella 11
E[1][11]= photon_energy*fxy_ph->Integral(-4.275,-1.425,-1.425,1.425); //cella 12
E[1][12]= photon_energy*fxy_ph->Integral(-1.425,1.425,-1.425,1.425); //cella 13
E[1][13]= photon_energy*fxy_ph->Integral(1.425,4.275,-1.425,1.425); //cella 14
E[1][14]= photon_energy*fxy_ph->Integral(4.275,7.125,-1.425,1.425); //cella 15
E[1][15]= photon_energy*fxy_ph->Integral(-7.125,-4.275,-4.275,-1.425); //cella 16
E[1][16]= photon_energy*fxy_ph->Integral(-4.275,-1.425,-4.275,-1.425); //cella 17
E[1][17]= photon_energy*fxy_ph->Integral(-1.425,1.425,-4.275,-1.425); //cella 18
E[1][18]= photon_energy*fxy_ph->Integral(1.425,4.275,-4.275,-1.425); //cella 19
E[1][19]= photon_energy*fxy_ph->Integral(4.275,7.125,-4.275,-1.425); //cella 20
E[1][20]= photon_energy*fxy_ph->Integral(-7.125,-4.275,-7.125,-4.275); //cella 21
E[1][21]= photon_energy*fxy_ph->Integral(-4.275,-1.425,-7.125,-4.275); //cella 22
E[1][22]= photon_energy*fxy_ph->Integral(-1.425,1.425,-7.125,-4.275); //cella 23
E[1][23]= photon_energy*fxy_ph->Integral(1.425,4.275,-7.125,-4.275); //cella 24
E[1][24]= photon_energy*fxy_ph->Integral(4.275,7.125,-7.125,-4.275); //cella 25 
  
// ENERGIA TOTALE NELLE CELLE QUANDO CI SONO ELETTRONE E FOTONE
E[2][0]= E[0][0]+E[1][0]; //cella 1
E[2][1]= E[0][1]+E[1][1]; //cella 2
E[2][2]= E[0][2]+E[1][2]; //cella 3
E[2][3]= E[0][3]+E[1][3]; //cella 4
E[2][4]= E[0][4]+E[1][4]; //cella 5
E[2][5]= E[0][5]+E[1][5]; //cella 6
E[2][6]= E[0][6]+E[1][6]; //cella 7
E[2][7]= E[0][7]+E[1][7]; //cella 8
E[2][8]= E[0][8]+E[1][8]; //cella 9
E[2][9]= E[0][9]+E[1][9]; //cella 10
E[2][10]= E[0][10]+E[1][10]; //cella 11
E[2][11]= E[0][11]+E[1][11]; //cella 12
E[2][12]= E[0][12]+E[1][12]; //cella 13
E[2][13]= E[0][13]+E[1][13]; //cella 14
E[2][14]= E[0][14]+E[1][14]; //cella 15
E[2][15]= E[0][15]+E[1][15]; //cella 16
E[2][16]= E[0][16]+E[1][16]; //cella 17
E[2][17]= E[0][17]+E[1][17]; //cella 18
E[2][18]= E[0][18]+E[1][18]; //cella 19
E[2][19]= E[0][19]+E[1][19]; //cella 20
E[2][20]= E[0][20]+E[1][20]; //cella 21
E[2][21]= E[0][21]+E[1][21]; //cella 22
E[2][22]= E[0][22]+E[1][22]; //cella 23
E[2][23]= E[0][23]+E[1][23]; //cella 24
E[2][24]= E[0][24]+E[1][24]; //cella 25 
            

            } 
    else { 
// ENERGIA TOTALE NELLE CELLE QUANDO C'E' SOLO L'ELETTRONE
E[2][0]= E[0][0]; //cella 1
E[2][1]= E[0][1]; //cella 2
E[2][2]= E[0][2]; //cella 3
E[2][3]= E[0][3]; //cella 4
E[2][4]= E[0][4]; //cella 5
E[2][5]= E[0][5]; //cella 6
E[2][6]= E[0][6]; //cella 7
E[2][7]= E[0][7]; //cella 8
E[2][8]= E[0][8]; //cella 9
E[2][9]= E[0][9]; //cella 10
E[2][10]= E[0][10]; //cella 11
E[2][11]= E[0][11]; //cella 12
E[2][12]= E[0][12]; //cella 13
E[2][13]= E[0][13]; //cella 14
E[2][14]= E[0][14]; //cella 15
E[2][15]= E[0][15]; //cella 16
E[2][16]= E[0][16]; //cella 17
E[2][17]= E[0][17]; //cella 18
E[2][18]= E[0][18]; //cella 19
E[2][19]= E[0][19]; //cella 20
E[2][20]= E[0][20]; //cella 21
E[2][21]= E[0][21]; //cella 22
E[2][22]= E[0][22]; //cella 23
E[2][23]= E[0][23]; //cella 24
E[2][24]= E[0][24]; //cella 25 
        }
    if (jentry==132941) {for (Int_t i=0,i<25,i++){cout << "Energia nella cella " << i <<" : " << E[3][i] << endl;}}
    
   }
       
  
       
       
}
    
}