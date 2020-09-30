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
    
    Double_t n_tot_cut=0; //numero di elettroni nel calorimetro
    Double_t n_two_cut=0; //numero di elettroni che hanno una distanza maggiore di 2RM dal fotone, quindi 2 clusters
    Double_t n_one_cut=0; //casi rimanenti che formano 1 cluster

    TH2F  *Th_emu = new TH2F("h2da1" , " Th e Vs. Th mu one cluster",500,0,100,500,0,5);
    TH2F  *Th_emu_cut = new TH2F("h2da1" , " Th e Vs. Th mu one cluster with cut",500,0,100,500,0,5);
    TH2F  *E_R = new TH2F("h2da1" , " R Vs. E_CAL one cluster",70,0,0.14,500,1000,160);
    
    TH1F* DR=new TH1F("DR", "Distanza elettrone-fotone", 70,0,0.14);
    TH1F* DR_cut=new TH1F("DR", "Distanza elettrone-fotone", 70,0,0.14);
    
    
    
    
    
    if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();


    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
   
    if (photon_coox!=-1 && photon_cooy!=-1)
    {  
     E_CAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {  E_CAL=detKinBeamRot_Ee;}
       
    d_e_ph=sqrt( (detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy) ); 
    
       if(E_CAL>1){
       
    if (abs(detKinBeamRot_cooXe) < 0.07 && abs(detKinBeamRot_cooYe) < 0.07)
    {
        n_tot_cut++;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (abs(photon_coox)<0.07 && abs(photon_cooy)<0.07)
               
           {    //if (photon_energy>0.2) DR->Fill(d_e_ph,wgt_full);
               
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=2RM DALL'ELETTRONE             
                if (d_e_ph>2*Rm )
                {
                    if (photon_energy>0.2) {n_two_cut++; 
                                            DR_cut->Fill(d_e_ph,wgt_full);
                                           }
                    }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { if (photon_energy>0.2) { n_one_cut++;
                                            Th_emu_cut->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
                                            DR_cut->Fill(d_e_ph,wgt_full);
                                          }
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one_cut++;
            Th_emu_cut->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);}
    
    } 
    }
       
       
       // senza TAGLIO
       
if (abs(detKinBeamRot_cooXe) < 0.07 && abs(detKinBeamRot_cooYe) < 0.07)
    {
        n_tot++;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (abs(photon_coox)<0.07 && abs(photon_cooy)<0.07)
               
           {    //if (photon_energy>0.2) DR->Fill(d_e_ph,wgt_full);
               
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=2RM DALL'ELETTRONE             
                if (d_e_ph>2*Rm )
                {
                    if (photon_energy>0.2) {n_two++; 
                                            DR->Fill(d_e_ph,wgt_full);
                                            E_R->Fill(d_e_ph,E_CAL,wgt_full);
                                           }
                    }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { if (photon_energy>0.2) { n_one++;
                                            Th_emu->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
                                            DR->Fill(d_e_ph,wgt_full);
                                            E_R->Fill(d_e_ph,E_CAL,wgt_full);
                                           
                                          }
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one++;
            Th_emu->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
             }
    
    } 
    
    
       
   }
    
          
cout << "Elettroni totali nel calorimetro: " << n_tot << endl;
cout << "Elettroni ad una distanza 2RM dal fotone: " << n_two << endl;
cout << "Eventi in cui vedo solo un cluster: " << n_one << endl;
cout << "Elettroni totali nel calorimetro CON TAGLIO: " << n_tot_cut << endl;
cout << "Elettroni ad una distanza 2RM dal fotone CON TAGLIO: " << n_two_cut << endl;
cout << "Eventi in cui vedo solo un cluster CON TAGLIO: " << n_one_cut << endl;
    
TCanvas * tmue= new TCanvas("tmue","tmue",1000,100,2500,2000); 
tmue->Divide(2,2);
tmue->cd(1);
Th_emu_cut->SetMarkerColor(kRed);
Th_emu_cut->Draw("COLZ");
tmue->cd(2);
Th_emu->Draw("COLZ");
tmue->cd(3);
Th_emu->ProjectionY()->DrawClone();
Th_emu_cut->ProjectionY()->DrawClone("same");
tmue->cd(4);
Th_emu->ProjectionX()->DrawClone();
Th_emu_cut->ProjectionX()->DrawClone("same");
tmue ->SaveAs("Th_emu.png"); 
    
TCanvas * Drr= new TCanvas("Drr","Drr",1000,100,2500,2000);
Drr->Divide(2,1);
Drr->cd(1);
DR->Draw();
DR_cut->SetLineColor(kRed);
DR_cut->Draw("same");  
Drr->cd(2);
E_R->Draw("COLZ");
Drr ->SaveAs("DReph.png");   

    
    
}