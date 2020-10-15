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
Double_t n_tot_eph=0.;
    
    Double_t n_tot=0.;
    Double_t n_one=0.;
    Double_t n_two=0.;
    Double_t ratio=0.; // ratio of two clusters over total, events to drop
    
    Double_t n_tot_cut=0.; //numero di elettroni nel calorimetro
    Double_t n_two_cut=0.; //numero di elettroni che hanno una distanza maggiore di 2RM dal fotone, quindi 2 clusters
    Double_t n_one_cut=0.; //casi rimanenti che formano 1 cluster
    Double_t ratio_cut=0.; // ratio of two clusters over total, events to drop
    
    Double_t n_tot0=0.; //numero di elettroni nel calorimetro
    Double_t n_two0=0.; //numero di elettroni che hanno una distanza maggiore di 2RM dal fotone, quindi 2 clusters
    Double_t n_one0=0.; //casi rimanenti che formano 1 cluster
    Double_t ratio0=0.; // ratio of two clusters over total, events to drop
    
    Double_t n_tot_cut0=0.; //numero di elettroni nel calorimetro
    Double_t n_two_cut0=0.; //numero di elettroni che hanno una distanza maggiore di 2RM dal fotone, quindi 2 clusters
    Double_t n_one_cut0=0.; //casi rimanenti che formano 1 cluster
    Double_t ratio_cut0=0.;
    
    Double_t n_tot1=0.; //numero di elettroni nel calorimetro
    Double_t n_two1=0.; //numero di elettroni che hanno una distanza maggiore di 2RM dal fotone, quindi 2 clusters
    Double_t n_one1=0.; //casi rimanenti che formano 1 cluster
    Double_t ratio1=0.; // ratio of two clusters over total, events to drop
    
    Double_t n_tot_cut1=0.; //numero di elettroni nel calorimetro
    Double_t n_two_cut1=0.; //numero di elettroni che hanno una distanza maggiore di 2RM dal fotone, quindi 2 clusters
    Double_t n_one_cut1=0.; //casi rimanenti che formano 1 cluster
    Double_t ratio_cut1=0.;

        

//Double_t same_cell=0.;
//Double_t different_cell=0.;
Double_t E_CAL;
Double_t Rm = 1.959 ; //raggio di Moliere in centimetri    
 

    
    if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();


    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
        detKinBeamRot_cooXe=detKinBeamRot_cooXe*100;
        detKinBeamRot_cooYe=detKinBeamRot_cooYe*100;
           photon_coox=photon_coox*100;
        photon_cooy=photon_cooy*100;
    if (photon_coox!=-1 && photon_cooy!=-1)
    {  
     E_CAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {  E_CAL=detKinBeamRot_Ee;}
       
       
    Double_t d_e_ph=sqrt( (detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy) ); 
       
    Double_t d_e_mu=sqrt( (detKinBeamRot_cooXe-detKinBeamRot_cooXmu)*(detKinBeamRot_cooXe-detKinBeamRot_cooXmu)+(detKinBeamRot_cooYe-detKinBeamRot_cooYmu)*(detKinBeamRot_cooYe-detKinBeamRot_cooYmu) ); 
       
  if (photon_n_cell_ph!=0 && detKinBeamRot_n_cell_e!=0){ 
      
if (d_e_ph>1*Rm)
{   n_tot_eph+=wgt_full;
 if (photon_n_cell_ph==detKinBeamRot_n_cell_e)
    {same_cell+=wgt_full;}
    
    else {different_cell+=wgt_full;}  // contatori double 
}
  }
       
       
if (detKinBeamRot_n_cell_e!=0)
    {n_tot+=wgt_full;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (photon_n_cell_ph!=0)    
           {
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=N*RM DALL'ELETTRONE             
                if (d_e_ph>1*Rm)
                {
                    //if (photon_energy>0.2) {
                    n_two+=wgt_full;
                }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { 
                //if (photon_energy>0.2) {
                n_one+=wgt_full;
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one+=wgt_full;}
    } 
       
if (E_CAL>1)
{
    if (detKinBeamRot_n_cell_e!=0)
    {n_tot_cut+=wgt_full;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (photon_n_cell_ph!=0)    
           {
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=N*RM DALL'ELETTRONE             
                if (d_e_ph>1*Rm)
                {
                    n_two_cut+=wgt_full;
                }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { 
                n_one_cut+=wgt_full;
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one_cut+=wgt_full;}
    } 
}
       
       
//-------->TARGET ZERO
if (detKinBeamRot_tar==0)
{
    if (detKinBeamRot_n_cell_e!=0)
    {n_tot0+=wgt_full;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (photon_n_cell_ph!=0)    
           {
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=N*RM DALL'ELETTRONE             
                if (d_e_ph>1*Rm)
                {
                    //if (photon_energy>0.2) {
                    n_two0+=wgt_full;
                }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { 
                //if (photon_energy>0.2) {
                n_one0+=wgt_full;
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one0+=wgt_full;}
    } 
       
if (E_CAL>1)
{
    if (detKinBeamRot_n_cell_e!=0)
    {n_tot_cut0+=wgt_full;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (photon_n_cell_ph!=0)    
           {
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=N*RM DALL'ELETTRONE             
                if (d_e_ph>1*Rm)
                {
                    n_two_cut0+=wgt_full;
                }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { 
                n_one_cut0+=wgt_full;
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one_cut0+=wgt_full;}
    } 
}
}

//-------->TARGET UNO
if(detKinBeamRot_tar==1)
{
    if (detKinBeamRot_n_cell_e!=0)
    {n_tot1+=wgt_full;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (photon_n_cell_ph!=0)    
           {
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=N*RM DALL'ELETTRONE             
                if (d_e_ph>1*Rm)
                {
                    //if (photon_energy>0.2) {
                    n_two1+=wgt_full;
                }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { 
                //if (photon_energy>0.2) {
                n_one1+=wgt_full;
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one1+=wgt_full;}
    } 
       
if (E_CAL>1)
{
    if (detKinBeamRot_n_cell_e!=0)
    {n_tot_cut1+=wgt_full;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (photon_n_cell_ph!=0)    
           {
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=N*RM DALL'ELETTRONE             
                if (d_e_ph>1*Rm)
                {
                    n_two_cut1+=wgt_full;
                }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { 
                n_one_cut1+=wgt_full;
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one_cut1+=wgt_full;}
    } 
}
}
 ratio=n_two/n_tot;
ratio_cut=n_two_cut/n_tot_cut;   
ratio0=n_two0/n_tot0;
ratio_cut0=n_two_cut0/n_tot_cut0; 
ratio1=n_two1/n_tot1;
ratio_cut1=n_two_cut1/n_tot_cut1;        
 
}

    
 cout << "Elettroni e fotoni nella stessa cella: " << same_cell << endl;
 cout << "Elettroni e fotoni in una diversa cella: " << different_cell << endl;   
cout << "n tot con fotoni= " << n_tot_eph << endl;
cout << "-------------------------------------------"<<endl;

cout << "Elettroni totali nel calorimetro: " << n_tot << endl;
cout << "Elettroni ad una distanza 2RM dal fotone: " << n_two << endl;
cout << "Eventi in cui vedo solo un cluster: " << n_one << endl;
cout << "Frazione di eventi scartabili: " << ratio <<endl;
cout << "-------------------------------------------"<<endl;
cout << "Elettroni totali nel calorimetro CON TAGLIO: " << n_tot_cut << endl;
cout << "Elettroni ad una distanza 2RM dal fotone CON TAGLIO: " << n_two_cut << endl;
cout << "Eventi in cui vedo solo un cluster CON TAGLIO: " << n_one_cut << endl;
cout << "Frazione di eventi scartabili CON TAGLIO: " << ratio_cut <<endl;
    
cout << endl;
    
cout << "Elettroni totali nel calorimetro TAR 0: " << n_tot0 << endl;
cout << "Elettroni ad una distanza 2RM dal fotone TAR 0: " << n_two0 << endl;
cout << "Eventi in cui vedo solo un cluster TAR 0: " << n_one0 << endl;
cout << "Frazione di eventi scartabili TAR 0: " << ratio0 <<endl;
cout << "-------------------------------------------"<<endl;
cout << "Elettroni totali nel calorimetro CON TAGLIO TAR 0: " << n_tot_cut0 << endl;
cout << "Elettroni ad una distanza 2RM dal fotone CON TAGLIO TAR 0: " << n_two_cut0 << endl;
cout << "Eventi in cui vedo solo un cluster CON TAGLIO TAR 0: " << n_one_cut0 << endl;
cout << "Frazione di eventi scartabili CON TAGLIO TAR0: " << ratio_cut0 <<endl;
    
cout << endl;
    
cout << "Elettroni totali nel calorimetro TAR 1: " << n_tot1 << endl;
cout << "Elettroni ad una distanza 2RM dal fotone TAR 1: " << n_two1 << endl;
cout << "Eventi in cui vedo solo un cluster TAR 1: " << n_one1 << endl;
cout << "Frazione di eventi scartabili TAR 1: " << ratio1 <<endl;
cout << "-------------------------------------------"<<endl;
cout << "Elettroni totali nel calorimetro CON TAGLIO TAR 1: " << n_tot_cut1 << endl;
cout << "Elettroni ad una distanza 2RM dal fotone CON TAGLIO TAR 1: " << n_two_cut1 << endl;
cout << "Eventi in cui vedo solo un cluster CON TAGLIO TAR 1: " << n_one_cut1 << endl;
cout << "Frazione di eventi scartabili CON TAGLIO TAR 1: " << ratio_cut1 <<endl;
    
}