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
    
Double_t Rm= 0.022; //raggio di Moliere in metri
Double_t d; //distanza elettronce fotone.
Double_t z  =2.025; //distanza origine-calorimetro in metri.
Double_t nEl=0; // numero di elettroni con un fotone
Double_t nElPh=0; // numero elecctorni con fotone in generale
Double_t nElNO=0; // numero di elettroni senzs un fotone
Double_t nElNORm=0; // numero di elettroni con un fotone ma senza richesta Rm
    
Double_t nEl0=0; // numero di elettroni con un fotone dal target 0
Double_t nElPh0=0; // numero elecctorni con fotone in generale dal target 0
Double_t nElNO0=0; // numero di elettroni senzs un fotone dal target 0
Double_t nElNORm0=0; // numero di elettroni con un fotone ma senza richesta Rm dal target 0
    
Double_t nEl1=0; // numero di elettroni con un fotone dal target 1
Double_t nElPh1=0; // numero elecctorni con fotone in generale dal target 1
Double_t nElNO1=0; // numero di elettroni senzs un fotone dal target 1
Double_t nElNORm1=0; // numero di elettroni con un fotone ma senza richesta Rm dal target 1

TH2F  *X_Y_e  = new TH2F("h2da" , " X  Vs. y of the electron",140,-0.5,-0.5,140,-0.5,0.5);
TH2F  *X_Y_p  = new TH2F("h2daa" , " X  Vs. y of the photon",140,-0.5,-0.5,140,-0.5,0.5); 
 
TH2F  *X_Y_e0  = new TH2F("h2da0" , " X  Vs. y of the electron target 0",140,-0.5,-0.5,140,-0.5,0.5);
TH2F  *X_Y_p0  = new TH2F("h2da00" , " X  Vs. y of the photon target 0",140,-0.5,-0.5,140,-0.5,0.5);
    
TH2F  *X_Y_e1  = new TH2F("h2da1" , " X  Vs. y of the electron target 1",140,-0.5,-0.5,140,-0.5,0.5);
TH2F  *X_Y_p1  = new TH2F("h2da11" , " X  Vs. y of the photon  target 1",140,-0.5,-0.5,140,-0.5,0.5);
    
    
     if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
    
   //TGraph *energyThEl= new TGraph(nentries); 

    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
   
       //energyThEl->SetPoint(jentry,detKinBeamRot_the,detKinBeamRot_Ee);
       if (detKinBeamRot_cooXe < 0.07 && detKinBeamRot_cooYe < 0.07 && detKinBeamRot_cooXe > -0.07 && detKinBeamRot_cooYe > -0.07)
        {
       nElNO++;
       if (photon_coox != -1)
       {
           nElPh++;
double_t posEl = sqrt(detKinBeamRot_cooXe*detKinBeamRot_cooXe+detKinBeamRot_cooYe*detKinBeamRot_cooYe);
double_t posPh = sqrt(photon_coox*photon_coox+photon_cooy*photon_cooy);
d=posEl-posPh;   
  
  
      if (photon_coox < 0.07 && photon_cooy < 0.07 && photon_coox > -0.07 && photon_cooy > -0.07)
       {
          nElNORm++; 
           if (abs(d)>2*Rm)
           {
               nEl++;
           }
       }
   }
   
           }
       
              if (detKinBeamRot_cooXe < 0.07 && detKinBeamRot_cooYe < 0.07 && detKinBeamRot_cooXe > -0.07 && detKinBeamRot_cooYe > -0.07 && detKinBeamRot_tar==0 && detKinBeamRot_Ee>1)
        {
       nElNO0++;
              if (photon_coox != -1)
       {
           nElPh0++;
double_t posEl = sqrt(detKinBeamRot_cooXe*detKinBeamRot_cooXe+detKinBeamRot_cooYe*detKinBeamRot_cooYe);
double_t posPh = sqrt(photon_coox*photon_coox+photon_cooy*photon_cooy);
d=posEl-posPh;   
  
       if (photon_coox < 0.07 && photon_cooy < 0.07 && photon_coox > -0.07 && photon_cooy > -0.07)
       {
          nElNORm0++; 
           if (abs(d)>2*Rm)
           {
               nEl0++;
           }
       }
   
              }
           }
       
              if (detKinBeamRot_cooXe < 0.07 && detKinBeamRot_cooYe < 0.07 && detKinBeamRot_cooXe > -0.07 && detKinBeamRot_cooYe > -0.07 && detKinBeamRot_tar==1)
        {
       nElNO1++;
              if (photon_coox != -1)
              {
           nElPh1++;
double_t posEl = sqrt(detKinBeamRot_cooXe*detKinBeamRot_cooXe+detKinBeamRot_cooYe*detKinBeamRot_cooYe);
double_t posPh = sqrt(photon_coox*photon_coox+photon_cooy*photon_cooy);
d=posEl-posPh;   
  
       if (photon_coox < 0.07 && photon_cooy < 0.07 && photon_coox > -0.07 && photon_cooy > -0.07)
       {
          nElNORm1++; 
           if (abs(d)>2*Rm)
           {
               nEl1++;
           }
       }
   
              }
           }
   }
    
    
    

    /*energyThEl->SetTitle("Energy_e(theta_e)");
    energyThEl->SetMarkerColor(50);
    energyThEl->SetMarkerStyle(8);
    energyThEl->SetLineColor(9);
    energyThEl->SetMarkerStyle(kFullDotSmall);
    energyThEl->GetXaxis()->SetTitle("ThetaEl(mrad)");
    energyThEl->GetYaxis()->SetTitle("Energy(GeV)");   
    
    TCanvas * theE= new TCanvas("theE","theE",1000,100,2500,2000);
    energyThEl->Draw("AP");
    energyThEl->SaveAs("thetaEn.png");*/
    cout << "numero di elettroni senza fotoni: " << nElNO << endl;
    cout << "numero di elettroni con fotoni: " << nElPh << endl;
    cout << "numero di elettroni con fotoni nel calorimetro: " << nEl << endl;
    cout << "numero di elettroni con fotoni senza richiesta raggio Moliere: " << nElNORm << endl;
    cout<<endl;
    cout << "numero di elettroni senza fotoni dal tar 0: " << nElNO0 << endl;
        cout << "numero di elettroni con fotoni dal tar 0: " << nElPh0 << endl;
    cout << "numero di elettroni con fotoni nel calorimetro dal tar 0: " << nEl0 << endl;
    cout << "numero di elettroni con fotoni senza richiesta raggio Moliere dal tar 0: " << nElNORm0 << endl;
    cout<<endl;
    cout << "numero di elettroni senza fotoni dal tar 1: " << nElNO1 << endl;
    cout << "numero di elettroni con fotoni dal tar 1: " << nElPh1 << endl;
    cout << "numero di elettroni con fotoni nel calorimetro dal tar 1: " << nEl1 << endl;
    cout << "numero di elettroni con fotoni senza richiesta raggio Moliere dal tar 1: " << nElNORm1 << endl;
    
    
    
    
}