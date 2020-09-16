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
    
Double_t Rm = 0.022; //raggio di Moliere in metri
Double_t d; //distanza elettrone fotone.
Double_t z  = 2.025; //distanza origine-calorimetro in metri.
Double_t nEl=0; // numero di elettroni con un fotone
Double_t nElPh=0; // numero elettroni con fotone in generale
Double_t nElNO=0; // numero di elettroni nel cal
Double_t nElNORm=0; // numero di elettroni con un fotone ma senza richesta Rm
Double_t nElNoPh=0; // numero di elettroni senza fotone
    
    
Double_t nEl0=0; // numero di elettroni con un fotone dal target 0
Double_t nElPh0=0; // numero elettroni con fotone in generale dal target 0
Double_t nElNO0=0; // numero di elettroni nel cal dal target 0
Double_t nElNORm0=0; // numero di elettroni con un fotone ma senza richesta Rm dal target 0
Double_t nElNoPh0=0; // numero di elettroni senza fotone dal target 0
    
    
Double_t nEl1=0; // numero di elettroni con un fotone dal target 1
Double_t nElPh1=0; // numero elettroni con fotone in generale dal target 1
Double_t nElNO1=0; // numero di elettroni nel cal dal target 1
Double_t nElNORm1=0; // numero di elettroni con un fotone ma senza richesta Rm dal target 1
Double_t nElNoPh1=0; // numero di elettroni senza fotone dal target 1


TH1F* EnCalNORm=new TH1F("h2aN", "Energy e not 2 Rm distant from photons", 200,0,160);
TH1F* EnCalNORm0=new TH1F("h2aN", "Energy e not 2 Rm distant from photons TAR 0", 200,0,160);
TH1F* EnCalNORm1=new TH1F("h2aN", "Energy e not 2 Rm distant from photons TAR 1", 200,0,160);
TH1F* EnCalNoPh=new TH1F("h2aN", "Energy e wh/out ph", 200,0,160);
TH1F* EnCalNoPh0=new TH1F("h2aN", "Energy e wh/out ph TAR 0", 200,0,160);
TH1F* EnCalNoPh1=new TH1F("h2aN", "Energy e wh/out ph TAR 1", 200,0,160);   

TH1F* ThCalNORm=new TH1F("h2aN", "Theta e not 2 Rm distant from photons",180,0,0.1);
TH1F* ThCalNORm0=new TH1F("h2aN", "Theta e not 2 Rm distant from photons TAR 0", 180,0,0.1);
TH1F* ThCalNORm1=new TH1F("h2aN", "Theta e not 2 Rm distant from photons TAR 1", 180,0,0.1);
TH1F* ThCalNoPh=new TH1F("h2aN", "Theta e wh/out ph", 180,0,0.1);
TH1F* ThCalNoPh0=new TH1F("h2aN", "Theta e wh/out ph TAR 0", 180,0,0.1);
TH1F* ThCalNoPh1=new TH1F("h2aN", "Theta e wh/out ph TAR 1", 180,0,0.1);   
    
TH2F  *Th_E_noph  = new TH2F("h2da" , " Th  Vs. E of the electrons",140,0,0.1,140,0,160);
TH2F  *Th_E_noRm = new TH2F("h2da" , " Th  Vs. E of the electrons",140,0,0.1,140,0,160);
    
TH2F  *Th_E_noph0  = new TH2F("h2da0" , " Th  Vs. E of the electrons",140,0,0.1,140,0,160);
TH2F  *Th_E_noRm0 = new TH2F("h2da0" , " Th  Vs. E of the electrons",140,0,0.1,140,0,160);
    
TH2F  *Th_E_noph1  = new TH2F("h2da1" , " Th  Vs. E of the electrons",140,0,0.1,140,0,160);
TH2F  *Th_E_noRm1 = new TH2F("h2da1" , " Th  Vs. E of the electrons",140,0,0.1,140,0,160);


    
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
double_t posEl = sqrt(detKinBeamRot_cooXe*detKinBeamRot_cooXe+detKinBeamRot_cooYe*detKinBeamRot_cooYe);
double_t posPh = sqrt(photon_coox*photon_coox+photon_cooy*photon_cooy);
d=posEl-posPh;   
  
  
      if (photon_coox < 0.07 && photon_cooy < 0.07 && photon_coox > -0.07 && photon_cooy > -0.07)
       {
           
           if (abs(d)>2*Rm)
           {
               nEl++;
           }
          else {nElNORm++; EnCalNORm->Fill(detKinBeamRot_Ee,wgt_full); ThCalNORm->Fill(detKinBeamRot_the,wgt_full); Th_E_noRm->Fill(detKinBeamRot_the,detKinBeamRot_Ee,wgt_full);}
       }
      else nElPh++;     
   }
else {nElNOph++;
      EnCalNoPh->Fill(detKinBeamRot_Ee,wgt_full); ThCalNoPh->Fill(detKinBeamRot_the,wgt_full); Th_E_noph->Fill(detKinBeamRot_the, detKinBeamRot_Ee,wgt_full);}
 
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
     
           if (abs(d)>2*Rm)
           {
               nEl0++;
           }
           else {nElNORm0++; EnCalNORm0->Fill(detKinBeamRot_Ee,wgt_full); ThCalNORm0->Fill(detKinBeamRot_the,wgt_full); Th_E_noRm0->Fill(detKinBeamRot_the, detKinBeamRot_Ee,wgt_full);}
       }
   else nElPh0++;
              }
           else {nElNOph0++;
                 EnCalNoPh0->Fill(detKinBeamRot_Ee,wgt_full); ThCalNoPh0->Fill(detKinBeamRot_the,wgt_full); Th_E_noph0->Fill(detKinBeamRot_the, detKinBeamRot_Ee,wgt_full);}
    
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
           if (abs(d)>2*Rm)
           {
               nEl1++;
           }
           else {nElNORm1++; EnCalNORm1->Fill(detKinBeamRot_Ee,wgt_full); ThCalNORm1->Fill(detKinBeamRot_the,wgt_full); Th_E_noRm1->Fill(detKinBeamRot_the, detKinBeamRot_Ee,wgt_full);}
       }
   else  nElPh1++;
              }
    else {nElNOph1++; 
          EnCalNoPh1->Fill(detKinBeamRot_Ee,wgt_full); ThCalNoPh1->Fill(detKinBeamRot_the,wgt_full); Th_E_noph1->Fill(detKinBeamRot_the, detKinBeamRot_Ee,wgt_full);}

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
    cout << "numero di elettroni nel calorimetro: " << nElNO << endl;
    cout << "numero di elettroni con fotoni fuori dal calorimetro: " << nElPh << endl;
    cout << "numero di elettroni con fotoni nel calorimetro: " << nEl << endl;
    cout << "numero di elettroni con fotoni senza richiesta raggio Moliere: " << nElNORm << endl;
    cout<<endl;
    cout << "numero di elettroni senza fotoni: " << nElNoPh << endl;
    
    cout << "numero di elettroni nel calorimetro dal tar 0: " << nElNO0 << endl;
    cout << "numero di elettroni con fotoni fuori dal calorimetro dal tar 0: " << nElPh0 << endl;
    cout << "numero di elettroni con fotoni nel calorimetro dal tar 0: " << nEl0 << endl;
    cout << "numero di elettroni con fotoni senza richiesta raggio Moliere dal tar 0: " << nElNORm0 << endl;
    cout << "numero di elettroni senza fotoni dal tar 0: " << nElNoPh0 << endl;
    
    cout<<endl;
    cout << "numero di elettroni nel calorimetro dal tar 1: " << nElNO1 << endl;
    cout << "numero di elettroni con fotoni fuori dal calorimetro dal tar 1: " << nElPh1 << endl;
    cout << "numero di elettroni con fotoni nel calorimetro dal tar 1: " << nEl1 << endl;
    cout << "numero di elettroni con fotoni senza richiesta raggio Moliere dal tar 1: " << nElNORm1 << endl;
    cout << "numero di elettroni senza fotoni: " << nElNoPh1 << endl;
    
    
    
    TCanvas * e= new TCanvas("e","e",200,10,1000,1000);
    e->Divide(1,3);
    e->cd(1);
    EnCalNORm->Draw("HIST");
    EnCalNoPh->SetLineColor(kRed);
    EnCalNoPh->Draw("HIST same");
    gPad->SetLogx();
    
    e->cd(2);
    EnCalNORm0->Draw("HIST");
    EnCalNoPh0->SetLineColor(kRed);
    EnCalNoPh0->Draw("HIST same");
    gPad->SetLogx();
    
    e->cd(3);
    EnCalNORm1->Draw("HIST");
    EnCalNoPh1->SetLineColor(kRed);
    EnCalNoPh1->Draw("HIST same");
    gPad->SetLogx();
    
    e->SaveAs("EnergyElnoRm.png");
    
    TCanvas * te= new TCanvas("te","te",200,10,1000,1000);
    te->Divide(1,3);
    te->cd(1);
    ThCalNORm->Draw("HIST");
    ThCalNoPh->SetLineColor(kRed);
    ThCalNoPh->Draw("HIST same");

    
    te->cd(2);
    ThCalNORm0->Draw("HIST");
    ThCalNoPh0->SetLineColor(kRed);
    ThCalNoPh0->Draw("HIST same");
   
    
    te->cd(3);
    ThCalNORm1->Draw("HIST");
    ThCalNoPh1->SetLineColor(kRed);
    ThCalNoPh1->Draw("HIST same");
   
    
    te->SaveAs("ThElnoRm.png");
    
    TCanvas * dued1= new TCanvas("dued1","dued1",1000,100,2500,2000);
    Th_E_noph1->Draw("HIST");
    Th_E_noph1->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noph1->GetYaxis()->SetTitle("E [GeV]");
    Th_E_noRm1->SetMarkerColor(kRed);
    Th_E_noRm1->Draw("HIST same");
    Th_E_noRm1->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noRm1->GetYaxis()->SetTitle("E [GeV]");
  dued1->SaveAs("Eth1.png");
    
    TCanvas * dued= new TCanvas("dued","dued",1000,100,2500,2000);
    Th_E_noph->Draw("HIST");
    Th_E_noph->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noph->GetYaxis()->SetTitle("E [GeV]");
    Th_E_noRm->SetMarkerColor(kRed);
    Th_E_noRm->Draw("HIST same");
    Th_E_noRm->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noRm->GetYaxis()->SetTitle("E [GeV]");
  dued1->SaveAs("Eth.png");
    
    TCanvas * dued0= new TCanvas("dued0","dued0",1000,100,2500,2000);
    Th_E_noph0->Draw("HIST");
    Th_E_noph0->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noph0->GetYaxis()->SetTitle("E [GeV]");
    Th_E_noRm0->SetMarkerColor(kRed);
    Th_E_noRm0->Draw("HIST same");
    Th_E_noRm0->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noRm0->GetYaxis()->SetTitle("E [GeV]");
  dued1->SaveAs("Eth0.png");
    
}