#define atree_cxx
#include "atree.h"
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
    
   // typedef map<int, double>  energy_cell; 
   // energy_cell en_c;    

Double_t n=0.;
Double_t n_tot=0.;
Double_t n_tot_e=0.;
Double_t n_tot_eBIG=0.;
    
Double_t n_tot_eph=0.;
Double_t n_tot_NOph=0.;
int i=0;
int j=0;    
   /* Double_t n_tot=0.;
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
    Double_t ratio_cut1=0.;*/

    
Double_t same_cell=0.;
Double_t different_cell=0.;
Double_t E_CAL=0.;
Double_t Rm = 2.190 ; //raggio di Moliere in centimetri    
Double_t E9=0.;

TH1F* hist_E9_e=new TH1F("E9e", "E9 e tot", 500,0,1);
TH1F* hist_thxz_e=new TH1F("thetaXZ", "theta XZ plane e", 150,-0.1,0.1);
TH1F* hist_thyz_e=new TH1F("thetaYZ", "theta YZ plane e", 150,-0.1,0.1);
TH1F* hist_thxz_e1=new TH1F("thetaXZ", "theta XZ plane e", 150,-0.1,0.1);
TH1F* hist_thyz_e1=new TH1F("thetaYZ", "theta YZ plane e", 150,-0.1,0.1);
TH1F* hist_thxz_e2=new TH1F("thetaXZ", "theta XZ plane e", 150,-0.1,0.1);
TH1F* hist_thyz_e2=new TH1F("thetaYZ", "theta YZ plane e", 150,-0.1,0.1);
    
TH1F* hist_E9_eBIG=new TH1F("E9e5X5", "E9 e tot 5X5", 500,0,1);
TH1F* hist_thxz_eBIG=new TH1F("thetaXZ5X5", "theta XZ plane e 5X5", 150,-0.1,0.1);
TH1F* hist_thyz_eBIG=new TH1F("thetaYZ5X5", "theta YZ plane e 5X5", 150,-0.1,0.1);
TH1F* hist_thxz_e1BIG=new TH1F("thetaXZ5X5", "theta XZ plane e 5X5", 150,-0.1,0.1);
TH1F* hist_thyz_e1BIG=new TH1F("thetaYZ5X5", "theta YZ plane e 5X5", 150,-0.1,0.1);
TH1F* hist_thxz_e2BIG=new TH1F("thetaXZ5X5", "theta XZ plane e 5X5", 150,-0.1,0.1);
TH1F* hist_thyz_e2BIG=new TH1F("thetaYZ5X5", "theta YZ plane e 5X5", 150,-0.1,0.1);

/*TH1F* hist_E9_NOph=new TH1F("E9noph", "E9 NO photons", 500,0,1);
TH1F* hist_E9_eph_same=new TH1F("E9eph", "E9 e+ph same cell", 500,0,1);
TH1F* hist_E9_eph_diff=new TH1F("E9eph", "E9 e+ph different cell", 500,0,1);
    
    
TH1F* hist_Eout_9_eph_same=new TH1F("E9outeph", "E_3x3/Etotcal e+ph same cell", 300,0.89,1);
TH1F* hist_Eout_9_eph_diff=new TH1F("E9outeph", "E_3x3/Etotcal e+ph diff cell",300,0.89,1);
TH1F* hist_Eout_9_e=new TH1F("E9oute", "E_3x3/Etotcal e",300,0.89,1);
TH1F* hist_Eout_9_NOph=new TH1F("E9outnoph", "E_3x3/Etotcal NO photons",300,0.89,1);

TH1F* hist_dist=new TH1F("dist", "Dist e-gamma", 400,0,4);
TH1F* hist_dist_same=new TH1F("dist", "Dist e-gamma same cell", 400,0,4);
TH1F* hist_dist_diff=new TH1F("dist", "Dist e-gamma diff cel", 400,0,4);
    
TH1F* hist_ang=new TH1F("dist", "DTheta (Thel-Thph) e-gamma", 200,-25,25);
TH1F* hist_ang_same=new TH1F("dist", "DTheta (Thel-Thph) e-gamma same cell", 200,-25,25);
TH1F* hist_ang_diff=new TH1F("dist", "DTheta (Thel-Thph) e-gamma diff cel", 200,-25,25);
    
    
TH1F* Ephout=new TH1F("EnergyPH", "Energy Ph out", 75,0.2,150); 
TH1F* Thph=new TH1F("th", "th Ph out", 75,50,50); */
TH1F* The=new TH1F("th", "th El out", 75,0,90); 
TH1F* TheBIG=new TH1F("th", "th El out BIG", 75,0,90); 
    
TH1F* The1=new TH1F("th", "th El out TAR 1", 75,0,90); 
TH1F* TheBIG1=new TH1F("th", "th El out BIG TAR 1", 75,0,90); 
    
TH1F* The2=new TH1F("th", "th El out TAR 2", 75,0,90); 
TH1F* TheBIG2=new TH1F("th", "th El out BIG TAR 2", 75,0,90); 
    
    
TH2F  *E3x3  = new TH2F("ThEel" , " Theta el Vs. E_ECAL",1400,0,70,1400,0.2,140);
//TH2F  *E3x3noph  = new TH2F("ThEel" , " Theta el Vs. E_ECAL no ph",100,0,70,70,0.2,140);

    
if (fChain == 0) return;

Long64_t nentries = fChain->GetEntriesFast();

    
    


    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
       n_tot+=wgt_full;

        /*en_c[1]=detKinBeamRot_Ecell1; en_c[2]=detKinBeamRot_Ecell2; en_c[3]=detKinBeamRot_Ecell3; en_c[4]=detKinBeamRot_Ecell4; en_c[5]=detKinBeamRot_Ecell5;
        en_c[6]=detKinBeamRot_Ecell6; en_c[7]=detKinBeamRot_Ecell7; en_c[8]=detKinBeamRot_Ecell8; en_c[9]=detKinBeamRot_Ecell9; en_c[10]=detKinBeamRot_Ecell10;
        en_c[11]=detKinBeamRot_Ecell11; en_c[12]=detKinBeamRot_Ecell12; en_c[13]=detKinBeamRot_Ecell13; en_c[14]=detKinBeamRot_Ecell14; en_c[15]=detKinBeamRot_Ecell15;
        en_c[16]=detKinBeamRot_Ecell16; en_c[17]=detKinBeamRot_Ecell17; en_c[18]=detKinBeamRot_Ecell18; en_c[19]=detKinBeamRot_Ecell19; en_c[20]=detKinBeamRot_Ecell20;
        en_c[21]=detKinBeamRot_Ecell21; en_c[22]=detKinBeamRot_Ecell22; en_c[23]=detKinBeamRot_Ecell23; en_c[24]=detKinBeamRot_Ecell24; en_c[25]=detKinBeamRot_Ecell25;*/
       
    if (photon_coox!=-1 && photon_cooy!=-1)
    {  
     E_CAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {  E_CAL=detKinBeamRot_Ee;}

        detKinBeamRot_cooXe=detKinBeamRot_cooXe*100; // cm
        detKinBeamRot_cooYe=detKinBeamRot_cooYe*100; // cm
        photon_coox=photon_coox*100; // cm
        photon_cooy=photon_cooy*100; // cm
       
       //double Eout=(Etotcal-detKinBeamRot_E_clus3x3)/detKinBeamRot_E_clus3x3;
       //double Eout_9=detKinBeamRot_E_clus3x3/Etotcal;
       E9=detKinBeamRot_E_1/detKinBeamRot_E_clus3x3;
       /*cout << detKinBeamRot_n_max_Cell << " cella impatto elettrone " << detKinBeamRot_n_cell_e << "con energia " <<detKinBeamRot_Ee << " cella impatto fotone " << photon_n_cell_ph<< "con energia " <<photon_energy <<endl;*/
        Double_t anglex_e = atan2(detKinBeamRot_pXe_out, detKinBeamRot_pZe_out);
        Double_t angley_e = atan2(detKinBeamRot_pYe_out, detKinBeamRot_pZe_out);

       
    /*Double_t d_e_mu=sqrt( (detKinBeamRot_cooXe-detKinBeamRot_cooXmu)*(detKinBeamRot_cooXe-detKinBeamRot_cooXmu)+(detKinBeamRot_cooYe-detKinBeamRot_cooYmu)*(detKinBeamRot_cooYe-detKinBeamRot_cooYmu) ); detKinBeamRot_n_cell_e!=1 && detKinBeamRot_n_cell_e!=2 && detKinBeamRot_n_cell_e!=3 && detKinBeamRot_n_cell_e!=4 && detKinBeamRot_n_cell_e!=5 && detKinBeamRot_n_cell_e!=10 && detKinBeamRot_n_cell_e!=15 && detKinBeamRot_n_cell_e!=20 && detKinBeamRot_n_cell_e!=25 && detKinBeamRot_n_cell_e!=24 && detKinBeamRot_n_cell_e!=23 && detKinBeamRot_n_cell_e!=22 && detKinBeamRot_n_cell_e!=21 && detKinBeamRot_n_cell_e!=16 && detKinBeamRot_n_cell_e!=11 && detKinBeamRot_n_cell_e!=6*/

if (detKinBeamRot_n_cell_e!=0 && abs(detKinBeamRot_cooXe)<4.275 && abs(detKinBeamRot_cooYe)<4.275)  {     
    
    n_tot_e+=wgt_full;
    hist_E9_e->Fill(E9,wgt_full);
    hist_thxz_e->Fill(anglex_e,wgt_full);
    hist_thyz_e->Fill(angley_e,wgt_full);
    The->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {hist_thxz_e1->Fill(anglex_e,wgt_full);
    hist_thyz_e1->Fill(angley_e,wgt_full);
    The1->Fill(detKinBeamRot_ThEl_interaction,wgt_full);}
    if (detKinBeamRot_tar==1)
    {hist_thxz_e2->Fill(anglex_e,wgt_full);
    hist_thyz_e2->Fill(angley_e,wgt_full);
    The2->Fill(detKinBeamRot_ThEl_interaction,wgt_full);}
    
    //hist_Eout_9_e->Fill(Eout_9,wgt_full);
    
    if (detKinBeamRot_E_clus3x3!=0) {E3x3->Fill(detKinBeamRot_ThEl_interaction,detKinBeamRot_E_clus3x3,wgt_full);}
    
  /*if (photon_n_cell_ph!=0 && photon_n_cell_ph!=1 && photon_n_cell_ph!=2 && photon_n_cell_ph!=3 && photon_n_cell_ph!=4 && photon_n_cell_ph!=5 && photon_n_cell_ph!=10 && photon_n_cell_ph!=15 && photon_n_cell_ph!=20 && photon_n_cell_ph!=25 && photon_n_cell_ph!=24 && photon_n_cell_ph!=23 && photon_n_cell_ph!=22 && photon_n_cell_ph!=21 && photon_n_cell_ph!=16 && photon_n_cell_ph!=11 && photon_n_cell_ph!=6)
  { Ephout->Fill(photon_energy,wgt_full);
    Thph->Fill(photon_theta,wgt_full);
    The->Fill(detKinBeamRot_the,wgt_full);
    
   Double_t d_e_ph=sqrt( (detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy) )/Rm; 
    double Dtheta=detKinBeamRot_the-photon_theta;

      n_tot_eph+=wgt_full; // e+gamma sul calorimetro
//      hist_E9_eph->Fill(E9,wgt_full);
      hist_dist->Fill(d_e_ph,wgt_full);
    hist_ang->Fill(Dtheta,wgt_full);
        if (photon_n_cell_ph==detKinBeamRot_n_cell_e)
        {same_cell+=wgt_full;
         hist_E9_eph_same->Fill(E9,wgt_full);
         hist_dist_same->Fill(d_e_ph,wgt_full);
        hist_ang_same->Fill(Dtheta,wgt_full);
         
      hist_Eout_9_eph_same->Fill(Eout_9,wgt_full);
        }//stessa cella
        else {different_cell+=wgt_full;
              hist_E9_eph_diff->Fill(E9,wgt_full);
              hist_dist_diff->Fill(d_e_ph,wgt_full);
                hist_Eout_9_eph_diff->Fill(Eout_9,wgt_full);
              hist_ang_diff->Fill(Dtheta,wgt_full);
             }  // cella diversa
  } else {
      n_tot_NOph+=wgt_full;   
      hist_E9_NOph->Fill(E9,wgt_full);
      hist_Eout_9_NOph->Fill(Eout_9,wgt_full);
      //if (detKinBeamRot_E_clus3x3!=0) {E3x3noph->Fill(detKinBeamRot_ThEl_interaction,detKinBeamRot_E_clus3x3);}
        }*/

  }
if (detKinBeamRot_n_cell_e!=0)  {     
    
    n_tot_eBIG+=wgt_full;
    hist_E9_eBIG->Fill(E9,wgt_full);
    hist_thxz_eBIG->Fill(anglex_e,wgt_full);
    hist_thyz_eBIG->Fill(angley_e,wgt_full);
    TheBIG->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {hist_thxz_e1BIG->Fill(anglex_e,wgt_full);
    hist_thyz_e1BIG->Fill(angley_e,wgt_full);
    TheBIG1->Fill(detKinBeamRot_ThEl_interaction,wgt_full);}
    if (detKinBeamRot_tar==1)
    {hist_thxz_e2BIG->Fill(anglex_e,wgt_full);
    hist_thyz_e2BIG->Fill(angley_e,wgt_full);
    TheBIG2->Fill(detKinBeamRot_ThEl_interaction,wgt_full);}
    
    //hist_Eout_9_e->Fill(Eout_9,wgt_full);
    
//if (detKinBeamRot_E_clus3x3!=0) {E3x3BIG->Fill(detKinBeamRot_ThEl_interaction,detKinBeamRot_E_clus3x3,wgt_full);}
    
  /*if (photon_n_cell_ph!=0 && photon_n_cell_ph!=1 && photon_n_cell_ph!=2 && photon_n_cell_ph!=3 && photon_n_cell_ph!=4 && photon_n_cell_ph!=5 && photon_n_cell_ph!=10 && photon_n_cell_ph!=15 && photon_n_cell_ph!=20 && photon_n_cell_ph!=25 && photon_n_cell_ph!=24 && photon_n_cell_ph!=23 && photon_n_cell_ph!=22 && photon_n_cell_ph!=21 && photon_n_cell_ph!=16 && photon_n_cell_ph!=11 && photon_n_cell_ph!=6)
  { Ephout->Fill(photon_energy,wgt_full);
    Thph->Fill(photon_theta,wgt_full);
    The->Fill(detKinBeamRot_the,wgt_full);
    
   Double_t d_e_ph=sqrt( (detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy) )/Rm; 
    double Dtheta=detKinBeamRot_the-photon_theta;

      n_tot_eph+=wgt_full; // e+gamma sul calorimetro
//      hist_E9_eph->Fill(E9,wgt_full);
      hist_dist->Fill(d_e_ph,wgt_full);
    hist_ang->Fill(Dtheta,wgt_full);
        if (photon_n_cell_ph==detKinBeamRot_n_cell_e)
        {same_cell+=wgt_full;
         hist_E9_eph_same->Fill(E9,wgt_full);
         hist_dist_same->Fill(d_e_ph,wgt_full);
        hist_ang_same->Fill(Dtheta,wgt_full);
         
      hist_Eout_9_eph_same->Fill(Eout_9,wgt_full);
        }//stessa cella
        else {different_cell+=wgt_full;
              hist_E9_eph_diff->Fill(E9,wgt_full);
              hist_dist_diff->Fill(d_e_ph,wgt_full);
                hist_Eout_9_eph_diff->Fill(Eout_9,wgt_full);
              hist_ang_diff->Fill(Dtheta,wgt_full);
             }  // cella diversa
  } else {
      n_tot_NOph+=wgt_full;   
      hist_E9_NOph->Fill(E9,wgt_full);
      hist_Eout_9_NOph->Fill(Eout_9,wgt_full);
      //if (detKinBeamRot_E_clus3x3!=0) {E3x3noph->Fill(detKinBeamRot_ThEl_interaction,detKinBeamRot_E_clus3x3);}
        }*/

  }       
/*if (photon_n_cell_ph==0 && detKinBeamRot_n_cell_e!=0)   
{ hist_E9_e->Fill(E9,wgt_full);
if (Eout_9!=0) hist_Eout_9_e->Fill(Eout_9,wgt_full);
}*/
   
       
/*if (detKinBeamRot_n_cell_e!=0)
    {n_tot+=wgt_full;

// SE IL FOTONE E' PRODOTTO DENTRO AL CALORIMETRO
           if (photon_n_cell_ph!=0)    
           {
// SE IL FOTONE E' NEL CALORIMETRO AD UNA d=N*RM DALL'ELETTRONE             
                if (d_e_ph>3*Rm)
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
                if (d_e_ph>3*Rm)
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
                if (d_e_ph>3*Rm)
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
                if (d_e_ph>3*Rm)
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
                if (d_e_ph>3*Rm)
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
                if (d_e_ph>3*Rm)
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
ratio_cut1=n_two_cut1/n_tot_cut1;    */    
 
}
 cout << " Numero elettroni totali " << n_tot_e << " su un totale di " << n_tot << " eventi " << endl;
cout << "Numero el+fotoni sul calorimetro = " << n_tot_eph << " dove nella stessa cella ce ne sono: " << same_cell << " in una diversa cella: " << different_cell << endl;   
cout << "Numero elettorni senza fotoni sul calorimetro = " << n_tot_NOph << endl;
cout << "-------------------------------------------"<<endl;

    
    
/*cout << "Elettroni totali nel calorimetro: " << n_tot << endl;
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
cout << "Frazione di eventi scartabili CON TAGLIO TAR 1: " << ratio_cut1 <<endl;*/
    
TCanvas * c1= new TCanvas("c1","c1",1000,100,2500,2000);

hist_E9_eBIG->GetXaxis()->SetTitle("Ecentral/E3x3");
hist_E9_eBIG->SetLineWidth(3);
hist_E9_eBIG->SetLineColor(kRed);
    
hist_E9_eBIG->Draw("HIST");     
    
hist_E9_e->GetXaxis()->SetTitle("Ecentral/E3x3");
hist_E9_e->SetLineWidth(3);
hist_E9_e->Draw("HIST same"); 
gPad->BuildLegend(0.25,0.15,0.25,0.15);
/*hist_E9_eph_same->SetLineColor(kRed);
hist_E9_eph_same->SetLineWidth(3);
hist_E9_eph_same->Draw("HIST same"); 
    
hist_E9_eph_diff->SetLineColor(kGreen);
hist_E9_eph_diff->SetLineWidth(3);
hist_E9_eph_diff->Draw("HIST same"); 
    
hist_E9_NOph->SetLineColor(kOrange);
hist_E9_NOph->SetLineWidth(3);
hist_E9_NOph->Draw("HIST same");
gPad->BuildLegend(0.25,0.15,0.25,0.15);*/

c1->SaveAs("/home/LHCB-T3/espedicato/tesi/E9.png");
    
/*TCanvas * c2= new TCanvas("c2","c2",1000,100,2500,2000);

hist_dist->GetXaxis()->SetTitle("r[Rm]");
hist_dist->SetLineWidth(3);
hist_dist->Draw("HIST"); 

hist_dist_same->SetLineColor(kRed);
hist_dist_same->SetLineWidth(3);
hist_dist_same->Draw("HIST same"); 
    
hist_dist_diff->SetLineColor(kGreen);
hist_dist_diff->SetLineWidth(3);
hist_dist_diff->Draw("HIST same");
gPad->BuildLegend(0.25,0.15,0.25,0.15);

c2->SaveAs("/home/LHCB-T3/espedicato/tesi/dist.png");
    
TCanvas * c2a= new TCanvas("c2a","c2a",1000,100,2500,2000);

hist_ang->GetXaxis()->SetTitle("Dth [mrad]");
hist_ang->SetLineWidth(3);
hist_ang->Draw("HIST"); 

hist_ang_same->SetLineColor(kRed);
hist_ang_same->SetLineWidth(3);
hist_ang_same->Draw("HIST same"); 
    
hist_ang_diff->SetLineColor(kGreen);
hist_ang_diff->SetLineWidth(3);
hist_ang_diff->Draw("HIST same");
gPad->BuildLegend(0.25,0.15,0.25,0.15);

c2a->SaveAs("/home/LHCB-T3/espedicato/tesi/Dtheta.png");*/
    
/*TCanvas * c3= new TCanvas("c3","c3",1000,100,2500,2000);

hist_Eout_9_e->GetXaxis()->SetTitle("E3x3/E5x5");
hist_Eout_9_e->SetLineWidth(3);
hist_Eout_9_e->Draw("HIST"); 

hist_Eout_9_eph_same->SetLineColor(kRed);
hist_Eout_9_eph_same->SetLineWidth(3);
hist_Eout_9_eph_same->Draw("HIST same"); 
    
hist_Eout_9_eph_diff->SetLineColor(kGreen);
hist_Eout_9_eph_diff->SetLineWidth(3);
hist_Eout_9_eph_diff->Draw("HIST same"); 
    
hist_Eout_9_NOph->SetLineColor(kOrange);
hist_Eout_9_NOph->SetLineWidth(3);
hist_Eout_9_NOph->Draw("HIST same");
gPad->BuildLegend(0.25,0.15,0.25,0.15);

c3->SaveAs("/home/LHCB-T3/espedicato/tesi/out+3x3.png");*/

Int_t nx = E3x3->GetNbinsX();
Int_t ny = E3x3->GetNbinsY();
for (Int_t i=1; i<nx+1; i++) {
for (Int_t j=1; j<ny+1; j++) {
if (E3x3->GetBinContent(i,j)<1) E3x3->SetBinContent(i,j,0);}}
    
TCanvas * c4= new TCanvas("c4","c4",1000,100,2500,2000);
gStyle->SetPalette(kLake);
TColor::InvertPalette(); 
E3x3->GetXaxis()->SetTitle("Theta_el");
E3x3->GetYaxis()->SetTitle("Ereco3x3");
E3x3->Draw("COLZ");

c4->SaveAs("/home/LHCB-T3/espedicato/tesi/thE.png");   

/*TCanvas * c5= new TCanvas("c5","c5",1000,100,2500,2000);
c5->Divide(1,2);
c5->cd(1);
Ephout->GetXaxis()->SetTitle("E_ph[GeV]");
Ephout->SetLineWidth(3);
Ephout->SetMinimum(1);
gPad->SetLogy();
Ephout->Draw("HIST"); 
c5->cd(2);
The->SetLineColor(kRed);
The->Draw("HIST"); 
Thph->Draw("HIST same"); 
    
c5->SaveAs("/home/LHCB-T3/espedicato/tesi/ph_energy.png");     */
    
TCanvas * c5= new TCanvas("c5","c5",1000,100,2500,2000);
c5->Divide(1,3);
c5->cd(1);
TheBIG->SetLineColor(kRed);
TheBIG->Draw("HIST"); 
The->Draw("HIST same"); 
c5->cd(2);
TheBIG1->SetLineColor(kRed);
TheBIG1->Draw("HIST"); 
The1->Draw("HIST same");     
c5->cd(3);
TheBIG2->SetLineColor(kRed);
TheBIG2->Draw("HIST"); 
The2->Draw("HIST same"); 
    
c5->SaveAs("/home/LHCB-T3/espedicato/tesi/th_el.png"); 
  
    TCanvas * theC= new TCanvas("tar","tar",1500,1000,3500,2000);
    theC->Divide(2,2);
    theC->cd(1);
    hist_thxz_e->SetLineColor(46);
    hist_thxz_e->SetLineWidth(3);
    hist_thxz_e->Draw("HIST");
    hist_thxz_e1->SetLineColor(8);
    hist_thxz_e1->SetLineWidth(3);
    hist_thxz_e1->Draw("HIST SAME");
    hist_thxz_e2->SetLineColor(kBlack);
    hist_thxz_e2->SetLineWidth(3);
    hist_thxz_e2->Draw("HIST SAME");
    hist_thxz_e->GetXaxis()->SetTitle("Theta XZ [rad]");
    theC->cd(2);
    hist_thyz_e->SetLineColor(46);
    hist_thyz_e->SetLineWidth(3);
    hist_thyz_e->Draw("HIST");
    hist_thyz_e1->SetLineColor(8);
    hist_thyz_e1->SetLineWidth(3);
    hist_thyz_e1->Draw("HIST SAME");
    hist_thyz_e2->SetLineColor(kBlack);
    hist_thyz_e2->SetLineWidth(3);
    hist_thyz_e2->Draw("HIST SAME");
    hist_thyz_e->GetXaxis()->SetTitle("Theta YZ [rad]");
    

    theC->cd(3);
    hist_thxz_eBIG->SetLineColor(46);
    hist_thxz_eBIG->SetLineWidth(3);
    hist_thxz_eBIG->Draw("HIST");
    hist_thxz_e1BIG->SetLineColor(8);
    hist_thxz_e1BIG->SetLineWidth(3);
    hist_thxz_e1BIG->Draw("HIST SAME");
    hist_thxz_e2BIG->SetLineColor(kBlack);
    hist_thxz_e2BIG->SetLineWidth(3);
    hist_thxz_e2BIG->Draw("HIST SAME");
    hist_thxz_eBIG->GetXaxis()->SetTitle("Theta XZ [rad]");
    theC->cd(4);
    hist_thyz_eBIG->SetLineColor(46);
    hist_thyz_eBIG->SetLineWidth(3);
    hist_thyz_eBIG->Draw("HIST");
    hist_thyz_e1BIG->SetLineColor(8);
    hist_thyz_e1BIG->SetLineWidth(3);
    hist_thyz_e1BIG->Draw("HIST SAME");
    hist_thyz_e2BIG->SetLineColor(kBlack);
    hist_thyz_e2BIG->SetLineWidth(3);
    hist_thyz_e2BIG->Draw("HIST SAME");
    hist_thyz_eBIG->GetXaxis()->SetTitle("Theta YZ [rad]");
    
    theC->SaveAs("/home/LHCB-T3/espedicato/tesi/Th_el_XZYZ.png");     
    
    
}