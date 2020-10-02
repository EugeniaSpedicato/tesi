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
    Double_t Re; //posizione elettrone
    Double_t n_tot=0; //numero di elettroni nel calorimetro
    Double_t n_two=0; //numero di elettroni che hanno una distanza maggiore di 2RM dal fotone, quindi 2 clusters
    Double_t n_one=0; //casi rimanenti che formano 1 cluster
    
    Double_t n_tot_cut=0; //numero di elettroni nel calorimetro
    Double_t n_two_cut=0; //numero di elettroni che hanno una distanza maggiore di 2RM dal fotone, quindi 2 clusters
    Double_t n_one_cut=0; //casi rimanenti che formano 1 cluster
    
    
    Double_t n_tot0=0; //numero di elettroni nel calorimetro
    Double_t n_two0=0; //numero di elettroni che hanno una distanza maggiore di 2RM dal fotone, quindi 2 clusters
    Double_t n_one0=0; //casi rimanenti che formano 1 cluster
    
    Double_t n_tot_cut0=0; //numero di elettroni nel calorimetro
    Double_t n_two_cut0=0; //numero di elettroni che hanno una distanza maggiore di 2RM dal fotone, quindi 2 clusters
    Double_t n_one_cut0=0; //casi rimanenti che formano 1 cluster
    
    
    Double_t n_tot1=0; //numero di elettroni nel calorimetro
    Double_t n_two1=0; //numero di elettroni che hanno una distanza maggiore di 2RM dal fotone, quindi 2 clusters
    Double_t n_one1=0; //casi rimanenti che formano 1 cluster
    
    Double_t n_tot_cut1=0; //numero di elettroni nel calorimetro
    Double_t n_two_cut1=0; //numero di elettroni che hanno una distanza maggiore di 2RM dal fotone, quindi 2 clusters
    Double_t n_one_cut1=0; //casi rimanenti che formano 1 cluster

    TH2F  *Th_emu = new TH2F("h2da1" , " Th e Vs. Th mu one cluster",500,0,100,500,0,5);
    TH2F  *Th_emu_cut = new TH2F("h2da1" , " Th e Vs. Th mu one cluster with cut",500,0,100,500,0,5);
    TH2F  *E_R = new TH2F("h2da1" , " R Vs. E_CAL one cluster",70,0,0.14,500,1000,160);
    TH2F  *Th_E_el  = new TH2F("h2da" , " Th e Vs. E_CAL one cluster",500,0,100,500,0,160);
    TH2F  *Th_E_el_cut  = new TH2F("h2da" , " Th e Vs. E_CAL one cluster with cut",500,0,100,500,0,160);
    
    TH1F* DR=new TH1F("DR", "Distanza elettrone-fotone", 70,0,0.14);
    TH1F* DR_cut=new TH1F("DR", "Distanza elettrone-fotone", 70,0,0.14);
    
    TH2F  *Th_emu0 = new TH2F("h2da1" , " Th e Vs. Th mu one cluster TAR 0",500,0,100,500,0,5);
    TH2F  *Th_emu_cut0 = new TH2F("h2da1" , " Th e Vs. Th mu one cluster with cut TAR 0",500,0,100,500,0,5);
    TH2F  *E_R0 = new TH2F("h2da1" , " R Vs. E_CAL one cluster TAR 0",70,0,0.14,500,1000,160);
    TH2F  *Th_E_el0  = new TH2F("h2da" , " Th e Vs. E_CAL one cluster TAR 0",500,0,100,500,0,160);
    TH2F  *Th_E_el_cut0  = new TH2F("h2da" , " Th e Vs. E_CAL one cluster with cut TAR 0",500,0,100,500,0,160);
    
    TH1F* DR0=new TH1F("DR", "Distanza elettrone-fotone TAR 0", 70,0,0.14);
    TH1F* DR_cut0=new TH1F("DR", "Distanza elettrone-fotone TAR 0", 70,0,0.14);
    
    TH2F  *Th_emu1 = new TH2F("h2da1" , " Th e Vs. Th mu one cluster TAR 1",500,0,100,500,0,5);
    TH2F  *Th_emu_cut1 = new TH2F("h2da1" , " Th e Vs. Th mu one cluster with cut TAR 1",500,0,100,500,0,5);
    TH2F  *E_R1 = new TH2F("h2da1" , " R Vs. E_CAL one cluster TAR 1",70,0,0.14,500,1000,160);
    TH2F  *Th_E_el1  = new TH2F("h2da" , " Th e Vs. E_CAL one cluster TAR 1",500,0,100,500,0,160);
    TH2F  *Th_E_el_cut1  = new TH2F("h2da" , " Th e Vs. E_CAL one cluster with cut TAR 1",500,0,100,500,0,160);
    
    TH1F* DR1=new TH1F("DR", "Distanza elettrone-fotone TAR 1", 70,0,0.14);
    TH1F* DR_cut1=new TH1F("DR", "Distanza elettrone-fotone TAR 1", 70,0,0.14); 
    
    
    
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
    
    Re=sqrt( (detKinBeamRot_cooXe*detKinBeamRot_cooXe)+(detKinBeamRot_cooYe*detKinBeamRot_cooYe));
    
      
       
    //TUTTI I TARGET   
       if(E_CAL>1){
       
    if (abs(detKinBeamRot_cooXe) < 0.07 && abs(detKinBeamRot_cooYe) < 0.07)
    {
        n_tot_cut++;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (abs(photon_coox)<0.07 && abs(photon_cooy)<0.07)
               
           {    //if (photon_energy>0.2) DR->Fill(d_e_ph,wgt_full);
               
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=2RM DALL'ELETTRONE             
                if (d_e_ph>2*Rm)
                {
                    if (photon_energy>0.2) {n_two_cut++; 
                                            DR_cut->Fill(d_e_ph,wgt_full);
                                           }
                    }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { if (photon_energy>0.2) { n_one_cut++;
                                            Th_emu_cut->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
                                            DR_cut->Fill(d_e_ph,wgt_full);
                                            Th_E_el_cut->Fill(detKinBeamRot_the,E_CAL,wgt_full);
                                          }
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one_cut++;
            Th_emu_cut->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
            Th_E_el_cut->Fill(detKinBeamRot_the,E_CAL,wgt_full);
}
    
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
                                            E_R->Fill(Re,E_CAL,wgt_full);
                                           }
                    }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { if (photon_energy>0.2) { n_one++;
                                            Th_emu->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
                                            DR->Fill(d_e_ph,wgt_full);
                                            E_R->Fill(Re,E_CAL,wgt_full);
                                            Th_E_el->Fill(detKinBeamRot_the,E_CAL,wgt_full);
                                          }
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one++;
            Th_emu->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
            Th_E_el->Fill(detKinBeamRot_the,E_CAL,wgt_full);
             }
    
    } 
       
       
       
       
       
       
// DAL TARGET 0 
if (detKinBeamRot_tar==0){
    if(E_CAL>5){
       
    if (abs(detKinBeamRot_cooXe) < 0.07 && abs(detKinBeamRot_cooYe) < 0.07)
    {
        n_tot_cut0++;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (abs(photon_coox)<0.07 && abs(photon_cooy)<0.07)
               
           {    //if (photon_energy>0.2) DR->Fill(d_e_ph,wgt_full);
               
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=2RM DALL'ELETTRONE             
                if (d_e_ph>2*Rm)
                {
                    if (photon_energy>0.2) {n_two_cut0++; 
                                            DR_cut0->Fill(d_e_ph,wgt_full);
                                           }
                    }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { if (photon_energy>0.2) { n_one_cut0++;
                                            Th_emu_cut0->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
                                            DR_cut0->Fill(d_e_ph,wgt_full);
                                            Th_E_el_cut0->Fill(detKinBeamRot_the,E_CAL,wgt_full);
                                          }
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one_cut0++;
            Th_emu_cut0->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
            Th_E_el_cut0->Fill(detKinBeamRot_the,E_CAL,wgt_full);
}
    
    } 
    }
       
       
       // senza TAGLIO
       
if (abs(detKinBeamRot_cooXe) < 0.07 && abs(detKinBeamRot_cooYe) < 0.07)
    {
        n_tot0++;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (abs(photon_coox)<0.07 && abs(photon_cooy)<0.07)
               
           {    //if (photon_energy>0.2) DR->Fill(d_e_ph,wgt_full);
               
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=2RM DALL'ELETTRONE             
                if (d_e_ph>2*Rm )
                {
                    if (photon_energy>0.2) {n_two0++; 
                                            DR0->Fill(d_e_ph,wgt_full);
                                            E_R0->Fill(Re,E_CAL,wgt_full);
                                           }
                    }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { if (photon_energy>0.2) { n_one0++;
                                            Th_emu0->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
                                            DR0->Fill(d_e_ph,wgt_full);
                                            E_R0->Fill(Re,E_CAL,wgt_full);
                                            Th_E_el0->Fill(detKinBeamRot_the,E_CAL,wgt_full);
                                          }
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one0++;
            Th_emu0->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
            Th_E_el0->Fill(detKinBeamRot_the,E_CAL,wgt_full);
             }
    
    } }
       
       
       
// DAL TARGET 1
if (detKinBeamRot_tar==1){
    if(E_CAL>5){
       
    if (abs(detKinBeamRot_cooXe) < 0.07 && abs(detKinBeamRot_cooYe) < 0.07)
    {
        n_tot_cut1++;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (abs(photon_coox)<0.07 && abs(photon_cooy)<0.07)
               
           {    //if (photon_energy>0.2) DR->Fill(d_e_ph,wgt_full);
               
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=2RM DALL'ELETTRONE             
                if (d_e_ph>2*Rm)
                {
                    if (photon_energy>0.2) {n_two_cut1++; 
                                            DR_cut1->Fill(d_e_ph,wgt_full);
                                           }
                    }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { if (photon_energy>0.2) { n_one_cut1++;
                                            Th_emu_cut1->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
                                            DR_cut1->Fill(d_e_ph,wgt_full);
                                            Th_E_el_cut1->Fill(detKinBeamRot_the,E_CAL,wgt_full);
                                          }
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one_cut1++;
            Th_emu_cut1->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
            Th_E_el_cut1->Fill(detKinBeamRot_the,E_CAL,wgt_full);
}
    
    } 
    }
       
       
       // senza TAGLIO
       
if (abs(detKinBeamRot_cooXe) < 0.07 && abs(detKinBeamRot_cooYe) < 0.07)
    {
        n_tot1++;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (abs(photon_coox)<0.07 && abs(photon_cooy)<0.07)
               
           {    //if (photon_energy>0.2) DR->Fill(d_e_ph,wgt_full);
               
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=2RM DALL'ELETTRONE             
                if (d_e_ph>2*Rm )
                {
                    if (photon_energy>0.2) {n_two1++; 
                                            DR1->Fill(d_e_ph,wgt_full);
                                            E_R1->Fill(Re,E_CAL,wgt_full);
                                           }
                    }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { if (photon_energy>0.2) { n_one1++;
                                            Th_emu1->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
                                            DR1->Fill(d_e_ph,wgt_full);
                                            E_R1->Fill(Re,E_CAL,wgt_full);
                                            Th_E_el1->Fill(detKinBeamRot_the,E_CAL,wgt_full);
                                          }
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one1++;
            Th_emu1->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
            Th_E_el1->Fill(detKinBeamRot_the,E_CAL,wgt_full);
             }
    
    } }
       
   }
    
          
cout << "Elettroni totali nel calorimetro: " << n_tot << endl;
cout << "Elettroni ad una distanza 2RM dal fotone: " << n_two << endl;
cout << "Eventi in cui vedo solo un cluster: " << n_one << endl;
cout << "Elettroni totali nel calorimetro CON TAGLIO: " << n_tot_cut << endl;
cout << "Elettroni ad una distanza 2RM dal fotone CON TAGLIO: " << n_two_cut << endl;
cout << "Eventi in cui vedo solo un cluster CON TAGLIO: " << n_one_cut << endl;
    
<< endl;
    
cout << "Elettroni totali nel calorimetro TAR 0: " << n_tot0 << endl;
cout << "Elettroni ad una distanza 2RM dal fotone TAR 0: " << n_two0 << endl;
cout << "Eventi in cui vedo solo un cluster TAR 0: " << n_one0 << endl;
cout << "Elettroni totali nel calorimetro CON TAGLIO TAR 0: " << n_tot_cut0 << endl;
cout << "Elettroni ad una distanza 2RM dal fotone CON TAGLIO TAR 0: " << n_two_cut0 << endl;
cout << "Eventi in cui vedo solo un cluster CON TAGLIO TAR 0: " << n_one_cut0 << endl;
    
<< endl;
    
cout << "Elettroni totali nel calorimetro TAR 1: " << n_tot1 << endl;
cout << "Elettroni ad una distanza 2RM dal fotone TAR 1: " << n_two1 << endl;
cout << "Eventi in cui vedo solo un cluster TAR 1: " << n_one1 << endl;
cout << "Elettroni totali nel calorimetro CON TAGLIO TAR 1: " << n_tot_cut1 << endl;
cout << "Elettroni ad una distanza 2RM dal fotone CON TAGLIO TAR 1: " << n_two_cut1 << endl;
cout << "Eventi in cui vedo solo un cluster CON TAGLIO TAR 1: " << n_one_cut1 << endl;
    
TCanvas * tmue= new TCanvas("tmue","tmue",1000,100,2500,2000); 
tmue->Divide(2,2);
tmue->cd(1);
Th_emu_cut->SetMarkerColor(kOrange);
Th_emu_cut->Draw("COLZ");
tmue->cd(2);
Th_emu->Draw("COLZ");
gStyle->SetPalette(kCherry);
TColor::InvertPalette();
tmue->cd(3);
Th_emu->ProjectionY()->DrawClone("HIST");
Th_emu_cut->ProjectionY()->DrawClone("HIST same");
tmue->cd(4);
Th_emu->ProjectionX()->DrawClone("HIST");
Th_emu_cut->ProjectionX()->DrawClone("HIST same");
tmue->SaveAs("Th_emu.png"); 
    
TCanvas * tmue1= new TCanvas("tmue1","tmue1",1000,100,2500,2000); 
tmue1->Divide(2,2);
tmue1->cd(1);
Th_emu_cut1->SetMarkerColor(kOrange);
Th_emu_cut1->Draw("COLZ");
tmue1->cd(2);
Th_emu1->Draw("COLZ");
tmue1->cd(3);
Th_emu1->ProjectionY()->DrawClone("HIST");
Th_emu_cut1->ProjectionY()->DrawClone("HIST same");
tmue1->cd(4);
Th_emu1->ProjectionX()->DrawClone("HIST");
Th_emu_cut1->ProjectionX()->DrawClone("HIST same");
tmue1->SaveAs("Th_emu1.png");    
    
TCanvas * tmue0= new TCanvas("tmue0","tmue0",1000,100,2500,2000); 
tmue0->Divide(2,2);
tmue0->cd(1);
Th_emu_cut0->SetMarkerColor(kOrange);
Th_emu_cut0->Draw("COLZ");
tmue0->cd(2);
Th_emu0->Draw("COLZ");
tmue0->cd(3);
Th_emu0->ProjectionY()->DrawClone("HIST");
Th_emu_cut0->ProjectionY()->DrawClone("HIST same");
tmue0->cd(4);
Th_emu0->ProjectionX()->DrawClone("HIST");
Th_emu_cut0->ProjectionX()->DrawClone("HIST same");
tmue0->SaveAs("Th_emu0.png"); 

    
TCanvas * Drr= new TCanvas("Drr","Drr",1000,100,2500,2000);
Drr->Divide(2,1);
Drr->cd(1);
DR->Draw();
DR_cut->SetLineColor(kRed);
DR_cut->Draw("same");  
Drr->cd(2);
E_R->Draw("HIST");
Drr->SaveAs("DReph.png");  
    
TCanvas * Drr0= new TCanvas("Drr0","Drr0",1000,100,2500,2000);
Drr0->Divide(2,1);
Drr0->cd(1);
DR0->Draw();
DR_cut0->SetLineColor(kRed);
DR_cut0->Draw("same");  
Drr0->cd(2);
E_R0->Draw("HIST");
Drr0->SaveAs("DReph0.png"); 
    
TCanvas * Drr1= new TCanvas("Drr1","Drr1",1000,100,2500,2000);
Drr1->Divide(2,1);
Drr1->cd(1);
DR1->Draw();
DR_cut1->SetLineColor(kRed);
DR_cut1->Draw("same");  
Drr1->cd(2);
E_R1->Draw("HIST");
Drr1->SaveAs("DReph1.png");  

TCanvas * th_en= new TCanvas("th_en","th_en",1000,100,2500,2000); 
th_en->Divide(2,1);
th_en->cd(1);
Th_E_el->SetMarkerColor(kOrange);    
Th_E_el_cut->SetMarkerColor(kOrange);    
Th_E_el->Draw("COLZ");  
th_en->cd(2);
Th_E_el_cut->Draw("COLZ");  
th_en->SaveAs("theta-energy-electron.png");   

TCanvas * th_en1= new TCanvas("th_en1","th_en1",1000,100,2500,2000); 
th_en1->Divide(2,1);
th_en1->cd(1);
Th_E_el1->SetMarkerColor(kOrange);    
Th_E_el_cut1->SetMarkerColor(kOrange);    
Th_E_el1->Draw("COLZ");  
th_en1->cd(2);
Th_E_el_cut1->Draw("COLZ");  
th_en1->SaveAs("theta-energy-electron1.png");  
    
TCanvas * th_en0= new TCanvas("th_en0","th_en0",1000,100,2500,2000); 
th_en0->Divide(2,1);
th_en0->cd(1);
Th_E_el0->SetMarkerColor(kOrange);    
Th_E_el_cut0->SetMarkerColor(kOrange);    
Th_E_el0->Draw("COLZ");  
th_en0->cd(2);
Th_E_el_cut0->Draw("COLZ");  
th_en0->SaveAs("theta-energy-electron0.png");   
}