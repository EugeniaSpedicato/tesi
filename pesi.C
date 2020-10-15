#define atree_cxx
#include "next.h"
#include <TTree.h>
#include <cmath>
using namespace std;


void atree::Loop()
{

    
    Double_t Rm = 0.01959 ; //raggio di Moliere in metri
    Double_t E_CAL; //energy in the calorimeter
    Double_t d_e_ph; //distanza elettrone-fotone

    
    //CONSIDERANDO I PESI
    Double_t n_tot=0.; //numero di elettroni nel calorimetro
    Double_t n_two=0.; //numero di elettroni che hanno una distanza maggiore di 2RM dal fotone, quindi 2 clusters
    Double_t n_one=0.; //casi rimanenti che formano 1 cluster
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
    
    
    if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();


    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry+=wgt_full) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
   
    // Energia nel calorimetro    
       if (photon_coox!=-1 && photon_cooy!=-1)
    {  
     E_CAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {  E_CAL=detKinBeamRot_Ee;}
     
       
    //Distanza punti di impatto elettrone-fotone
    d_e_ph=sqrt( (detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy) );
       
       
//----->DA TUTTI I TARGET

//SENZA TAGLIO
if (abs(detKinBeamRot_cooXe) < 0.07125 && abs(detKinBeamRot_cooYe) < 0.07125)
    {
        n_tot+=wgt_full;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (abs(photon_coox)<0.07125 && abs(photon_cooy)<0.07125)
               
           {
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=N*RM DALL'ELETTRONE             
                if (d_e_ph>1*Rm)
                {
                    if (photon_energy>0.2) {n_two+=wgt_full;}
                }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { if (photon_energy>0.2) {n_one+=wgt_full;}
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one+=wgt_full;}
    } 


// CON TAGLIO
if (E_CAL>1) {if (abs(detKinBeamRot_cooXe) < 0.07125 && abs(detKinBeamRot_cooYe) < 0.07125)
    {
        n_tot_cut+=wgt_full;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (abs(photon_coox)<0.07125 && abs(photon_cooy)<0.07125)
               
           {
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=N*RM DALL'ELETTRONE             
                if (d_e_ph>1*Rm)
                {
                    if (photon_energy>0.2) {n_two_cut+=wgt_full;}
                }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { if (photon_energy>0.2) {n_one_cut+=wgt_full;}
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one_cut+=wgt_full;}
    } }

//---->TARGET 0

if (detKinBeamRot_tar==0)
{
    if (abs(detKinBeamRot_cooXe) < 0.07125 && abs(detKinBeamRot_cooYe) < 0.07125)
    {
        n_tot0+=wgt_full;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (abs(photon_coox)<0.07125 && abs(photon_cooy)<0.07125)
               
           {
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=N*RM DALL'ELETTRONE             
                if (d_e_ph>1*Rm)
                {
                    if (photon_energy>0.2) {n_two0+=wgt_full;}
                }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { if (photon_energy>0.2) {n_one0+=wgt_full;}
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one0+=wgt_full;}
    } 


// CON TAGLIO
if (E_CAL>1) {if (abs(detKinBeamRot_cooXe) < 0.07125 && abs(detKinBeamRot_cooYe) < 0.07125)
    {
        n_tot_cut0+=wgt_full;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (abs(photon_coox)<0.07125 && abs(photon_cooy)<0.07125)
               
           {
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=N*RM DALL'ELETTRONE             
                if (d_e_ph>1*Rm)
                {
                    if (photon_energy>0.2) {n_two_cut0+=wgt_full;}
                }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { if (photon_energy>0.2) {n_one_cut0+=wgt_full;}
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one_cut0+=wgt_full;}
    } }
}

//---->TARGET 1

if (detKinBeamRot_tar==1)
{
    if (abs(detKinBeamRot_cooXe) < 0.07125 && abs(detKinBeamRot_cooYe) < 0.07125)
    {
        n_tot1+=wgt_full;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (abs(photon_coox)<0.07125 && abs(photon_cooy)<0.07125)
               
           {
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=N*RM DALL'ELETTRONE             
                if (d_e_ph>1*Rm)
                {
                    if (photon_energy>0.2) {n_two1+=wgt_full;}
                }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { if (photon_energy>0.2) {n_one1+=wgt_full;}
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one1+=wgt_full;}
    } 


//CON TAGLIO
if (E_CAL>1) {if (abs(detKinBeamRot_cooXe) < 0.07125 && abs(detKinBeamRot_cooYe) < 0.07125)
    {
        n_tot_cut1+=wgt_full;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (abs(photon_coox)<0.07125 && abs(photon_cooy)<0.07125)
               
           {
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=N*RM DALL'ELETTRONE             
                if (d_e_ph>1*Rm)
                {
                    if (photon_energy>0.2) {n_two_cut1+=wgt_full;}
                }
// SE E' NEL CALORIMETRO MA AD UNA d<2RM
            else { if (photon_energy>0.2) {n_one_cut1+=wgt_full;}
                 }
           } 
// SE E' NON E' PRODOTTO O NON E' NEL CALORIMETRO        
        else {n_one_cut1+=wgt_full;}
    } }
}

ratio=n_two/n_tot;
ratio_cut=n_two_cut/n_tot_cut;   
ratio0=n_two0/n_tot0;
ratio_cut0=n_two_cut0/n_tot_cut0; 
ratio1=n_two1/n_tot1;
ratio_cut1=n_two_cut1/n_tot_cut1;   
       
       
   }
    
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
       