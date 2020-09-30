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
    
    Double_t Rm = 0.01959 ; //raggio di Moliere in metri
    Double_t E_CAL; //energy in the calorimeter
    Double_t d_e_ph; //distanza elettrone-fotone
    //distanza elettrone-fotone
    Double_t n_tot=0; //numero di elettroni nel calorimetro
    Double_t n_two=0; //numero di elettroni che hanno una distanza maggiore di 2RM dal fotone, quindi 2 clusters
    Double_t n_one=0; //casi rimanenti che formano 1 cluster

    TH2F  *Th_emu = new TH2F("h2da1" , " Th e Vs. Th mu one cluster",500,0,100,500,0,5);
    TH1F* DR=new TH1F("DR", "Distanza elettrone-fotone", 70,0,0.14);
    
    
    
    
    if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();


    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
   
    if (photon_coox!=-1 && photon_cooy!=-1 && photon_energy>0.2)
    {  
     E_CAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {  E_CAL=detKinBeamRot_Ee;}
       
    d_e_ph=sqrt( (detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy) ); 
    
       //if(E_CAL>1){
    if (abs(detKinBeamRot_cooXe) < 0.07 && abs(detKinBeamRot_cooYe) < 0.07)
    {
        n_tot++;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (abs(photon_coox)<0.07 && abs(photon_cooy)<0.07)
               
           {    DR->Fill(d_e_ph,wgt_full);
               
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=2RM DALL'ELETTRONE             
                if (d_e_ph>2*Rm )
                {
                    n_two++;
                }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else {n_one++;
                  Th_emu->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);}
           }
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one++;
            Th_emu->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);}
    
    }
  //  }
       
       
       
   }
    
          
cout << "Elettroni totali nel calorimetro: " << n_tot << endl;
cout << "Elettroni ad una distanza 2RM dal fotone: " << n_two << endl;
cout << "Eventi in cui vedo solo un cluster: " << n_one << endl;
    
TCanvas * tmue= new TCanvas("tmue","tmue",1000,100,2500,2000); 
tmue->Divide(2,2);
tmue->cd(1);
Th_emu->Draw("LEGO");
Th_emu->SetMarkerSize(5);
tmue->cd(2);
Th_emu->ProjectionY()->DrawClone();
tmue->cd(3);
Th_emu->ProjectionX()->DrawClone();
tmue ->SaveAs("Th_emu.png");   
    
TCanvas * Drr= new TCanvas("Drr","Drr",1000,100,2500,2000);
DR->Draw();   
Drr ->SaveAs("DReph.png");   
    
    
}