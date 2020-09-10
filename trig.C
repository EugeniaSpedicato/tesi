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

TH1F* EnCalNORm=new TH1F("h2aN", "Energy e not 2 Rm distant from photons", 200,0,160);
TH1F* EnCalNORm0=new TH1F("h2aN", "Energy e not 2 Rm distant from photons TAR 0", 200,0,160);
TH1F* EnCalNORm1=new TH1F("h2aN", "Energy e not 2 Rm distant from photons TAR 1", 200,0,160);
    
    
    
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
          else EnCalNORm->Fill(detKinBeamRot_Ee,wgt_full);
       }
      else nElPh++;     
   }
   
           }
       
              if (detKinBeamRot_cooXe < 0.07 && detKinBeamRot_cooYe < 0.07 && detKinBeamRot_cooXe > -0.07 && detKinBeamRot_cooYe > -0.07 && detKinBeamRot_tar==0 && detKinBeamRot_Ee>1)
        {
       nElNO0++;
              if (photon_coox != -1)
       {
           
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
           else EnCalNORm0->Fill(detKinBeamRot_Ee,wgt_full);
       }
   else nElPh0++;
              }
           }
       
              if (detKinBeamRot_cooXe < 0.07 && detKinBeamRot_cooYe < 0.07 && detKinBeamRot_cooXe > -0.07 && detKinBeamRot_cooYe > -0.07 && detKinBeamRot_tar==1)
        {
       nElNO1++;
              if (photon_coox != -1)
              {
          
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
           else EnCalNORm1->Fill(detKinBeamRot_Ee,wgt_full);
       }
   else  nElPh1++;
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
    cout << "numero di elettroni con fotoni fuori dal calorimetro: " << nElPh << endl;
    cout << "numero di elettroni con fotoni nel calorimetro: " << nEl << endl;
    cout << "numero di elettroni con fotoni senza richiesta raggio Moliere: " << nElNORm << endl;
    cout<<endl;
    cout << "numero di elettroni senza fotoni dal tar 0: " << nElNO0 << endl;
        cout << "numero di elettroni con fotoni fuori dal calorimetro dal tar 0: " << nElPh0 << endl;
    cout << "numero di elettroni con fotoni nel calorimetro dal tar 0: " << nEl0 << endl;
    cout << "numero di elettroni con fotoni senza richiesta raggio Moliere dal tar 0: " << nElNORm0 << endl;
    cout<<endl;
    cout << "numero di elettroni senza fotoni dal tar 1: " << nElNO1 << endl;
    cout << "numero di elettroni con fotoni fuori dal calorimetro dal tar 1: " << nElPh1 << endl;
    cout << "numero di elettroni con fotoni nel calorimetro dal tar 1: " << nEl1 << endl;
    cout << "numero di elettroni con fotoni senza richiesta raggio Moliere dal tar 1: " << nElNORm1 << endl;
    
    
    TCanvas * e= new TCanvas("e","e",400,10,1500,1000);
    e->Divide(1,3);
    e->cd(1);
    EnCalNORm->Draw("HIST");
    gPad->SetLogx();
    
    e->cd(2);
    EnCalNORm0->Draw("HIST");
    gPad->SetLogx();
    
    e->cd(3);
    EnCalNORm1->Draw("HIST");
    gPad->SetLogx();
    
    e->SaveAs("EnergyElnoRm.png");
    
}