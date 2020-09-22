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
Double_t nElNOph=0; // numero di elettroni senza fotone
    
    
Double_t nEl0=0; // numero di elettroni con un fotone dal target 0
Double_t nElPh0=0; // numero elettroni con fotone in generale dal target 0
Double_t nElNO0=0; // numero di elettroni nel cal dal target 0
Double_t nElNORm0=0; // numero di elettroni con un fotone ma senza richesta Rm dal target 0
Double_t nElNOph0=0; // numero di elettroni senza fotone dal target 0
    
    
Double_t nEl1=0; // numero di elettroni con un fotone dal target 1
Double_t nElPh1=0; // numero elettroni con fotone in generale dal target 1
Double_t nElNO1=0; // numero di elettroni nel cal dal target 1
Double_t nElNORm1=0; // numero di elettroni con un fotone ma senza richesta Rm dal target 1
Double_t nElNOph1=0; // numero di elettroni senza fotone dal target 1


TH1F* EnCalNORm=new TH1F("h2aN", "Energy e not 2 Rm distant from photons", 200,0,160);
TH1F* EnCalNORm0=new TH1F("h2aN", "Energy e not 2 Rm distant from photons TAR 0", 200,0,160);
TH1F* EnCalNORm1=new TH1F("h2aN", "Energy e not 2 Rm distant from photons TAR 1", 200,0,160);
TH1F* EnCalNoPh=new TH1F("h2aN", "Energy e wh/out ph", 200,0,160);
TH1F* EnCalNoPh0=new TH1F("h2aN", "Energy e wh/out ph TAR 0", 200,0,160);
TH1F* EnCalNoPh1=new TH1F("h2aN", "Energy e wh/out ph TAR 1", 200,0,160);   
TH1F* EnPhNoCal=new TH1F("h2aN", "Energy e with ph out of calorimeter", 200,0,160);
TH1F* EnPhNoCal0=new TH1F("h2aN", "Energy e with ph out of calorimeter TAR 0", 200,0,160);
TH1F* EnPhNoCal1=new TH1F("h2aN", "Energy e with ph out of calorimeter TAR 1", 200,0,160);   
    
TH1F* En_tot=new TH1F("h2aN", "Energy e + energy ph inside 2 RM", 200,0,160);
    
TH1F* ThCalNORm=new TH1F("h2aN", "Theta e not 2 Rm distant from photons",180,0,100);
TH1F* ThCalNORm0=new TH1F("h2aN", "Theta e not 2 Rm distant from photons TAR 0", 180,0,100);
TH1F* ThCalNORm1=new TH1F("h2aN", "Theta e not 2 Rm distant from photons TAR 1", 180,0,100);
TH1F* ThCalNoPh=new TH1F("h2aN", "Theta e wh/out ph", 180,0,100);
TH1F* ThCalNoPh0=new TH1F("h2aN", "Theta e wh/out ph TAR 0", 180,0,100);
TH1F* ThCalNoPh1=new TH1F("h2aN", "Theta e wh/out ph TAR 1", 180,0,100);   
TH1F* ThPhNoCal=new TH1F("h2aN", "Theta e with ph out of calorimete", 180,0,100);
TH1F* ThPhNoCal0=new TH1F("h2aN", "Theta e with ph out of calorimete TAR 0", 180,0,100);
TH1F* ThPhNoCal1=new TH1F("h2aN", "Theta e with ph out of calorimete TAR 1", 180,0,100); 

TH1F* ThCalNORmMU=new TH1F("h2aN", "Theta MU not 2 Rm distant from photons",180,0,5);
TH1F* ThCalNORm0MU=new TH1F("h2aN", "Theta MU not 2 Rm distant from photons TAR 0", 180,0,5);
TH1F* ThCalNORm1MU=new TH1F("h2aN", "Theta MU not 2 Rm distant from photons TAR 1", 180,0,5);
TH1F* ThCalNoPhMU=new TH1F("h2aN", "Theta MU wh/out ph", 180,0,5);
TH1F* ThCalNoPh0MU=new TH1F("h2aN", "Theta MU wh/out ph TAR 0", 180,0,5);
TH1F* ThCalNoPh1MU=new TH1F("h2aN", "Theta MU wh/out ph TAR 1", 180,0,5);   
TH1F* ThPhNoCalMU=new TH1F("h2aN", "Theta MU with ph out of calorimete", 180,0,5);
TH1F* ThPhNoCal0MU=new TH1F("h2aN", "Theta MU with ph out of calorimete TAR 0", 180,0,5);
TH1F* ThPhNoCal1MU=new TH1F("h2aN", "Theta MU with ph out of calorimete TAR 1", 180,0,5); 
    
TH2F  *Th_E_noph  = new TH2F("h2da" , " Th e Vs. Th  of the electrons whitout photons (LO)",140,0,100,140,0,160);
TH2F  *Th_E_noRm = new TH2F("h2da" , " Th  Vs. E of the electrons with photons <2Rm",140,0,100,140,0,160);
TH2F  *Th_E_PhNoCal = new TH2F("h2da" , " Th  Vs. E of the electrons with photons out of cal",140,0,100,140,0,160);

TH2F  *Th_E_nophMU  = new TH2F("h2da" , " Th MU Vs. Ee whitout photons (LO)",140,0,5,140,0,160);
TH2F  *Th_E_noRmMU = new TH2F("h2da" , " Th MU Vs. Ee with photons <2Rm",140,0,5,140,0,160);
TH2F  *Th_E_PhNoCalMU = new TH2F("h2da" , " Th MU Vs. Ee with photons out of cal",140,0,5,140,0,160);

    
TH2F  *Th_E_noph0  = new TH2F("h2da0" , " Th  Vs. E of the electrons whitout photons (LO) tar 0 ",140,0,100,140,0,160);
TH2F  *Th_E_noRm0 = new TH2F("h2da0" , " Th  Vs. E of the electrons with photons <2Rm tar 0",140,0,100,140,0,160);
TH2F  *Th_E_PhNoCal0 = new TH2F("h2da" , " Th  Vs. E of the electrons with photons out of cal tar 0",140,0,100,140,0,160);
    
TH2F  *Th_E_noph0MU  = new TH2F("h2da" , " Th MU Vs. Ee whitout photons (LO) tar 0",140,0,5,140,0,160);
TH2F  *Th_E_noRm0MU = new TH2F("h2da" , " Th MU Vs. Ee with photons <2Rm tar 0",140,0,5,140,0,160);
TH2F  *Th_E_PhNoCal0MU = new TH2F("h2da" , " Th MU Vs. Ee with photons out of cal tar 0",140,0,5,140,0,160);
    
    
TH2F  *Th_E_noph1  = new TH2F("h2da1" , " Th  Vs. E of the electrons whitout photons (LO) tar 1",140,0,100,140,0,160);
TH2F  *Th_E_noRm1 = new TH2F("h2da1" , " Th  Vs. E of the electrons with photons <2Rm tar 1",140,0,100,140,0,160);
TH2F  *Th_E_PhNoCal1 = new TH2F("h2da1" , " Th  Vs. E of the electrons with photons out of cal tar 1",140,0,100,140,0,160);
    
TH2F  *Th_E_noph1MU  = new TH2F("h2da1" , " Th MU Vs. Ee whitout photons (LO) tar 1",140,0,5,140,0,160);
TH2F  *Th_E_noRm1MU = new TH2F("h2da1" , " Th MU Vs. Ee with photons <2Rm tar 1",140,0,5,140,0,160);
TH2F  *Th_E_PhNoCal1MU = new TH2F("h2da1" , " Th MU Vs. Ee with photons out of cal tar 1",140,0,5,140,0,160);
    
TH2F  *Th_E_mu = new TH2F("h2da1" , " Th mu Vs. Th e when the photons are produced but out of CAL ",140,0,100,140,0,5);
TH2F  *Th_E_eph = new TH2F("h2da1" , " Th mu Vs. Th e when photons <2RM (NLO)",140,0,100,140,0,5);
TH2F  *Th_E_eNoph = new TH2F("h2da1" , " Th mu Vs. Th e when the electrons are without photons",140,0,100,140,0,5);

    
TH2F  *Th_E_mu0 = new TH2F("h2da1" , " Th mu Vs. Th e when e and mu when the photons are produced but out of CAL tar 0",140,0,100,140,0,5);
TH2F  *Th_E_eph0 = new TH2F("h2da1" , " Th mu Vs. Th e with photons <2RM(NLO) tar 0",140,0,100,140,0,5);
TH2F  *Th_E_eNoph0 = new TH2F("h2da1" , " Th mu Vs. Th e when electrons without photons tar 0",140,0,100,140,0,5);

    
TH2F  *Th_E_mu1 = new TH2F("h2da1" , " Th mu Vs. Th e when e and mu when the photons are produced but out of CAL tar 1",140,0,100,140,0,5);
TH2F  *Th_E_eph1 = new TH2F("h2da1" , " Th mu Vs. Th e with photons <2RM(NLO) tar 1",140,0,100,140,0,5);
TH2F  *Th_E_eNoph1 = new TH2F("h2da1" , " Th mu Vs. Th e when electrons without photons tar 1",140,0,100,140,0,5);
    
TH1F* E_ph=new TH1F("h2aN1", "E ph d>2Rm",180,0,100);
TH1F* E_ph0=new TH1F("h2aN2", "E ph d>2Rm tar 0",180,0,100);
TH1F* E_ph1=new TH1F("h2aN3", "E ph d>2Rm tar 1",180,0,100);
    
TH2F  *Ee_Eph = new TH2F("h2da1" , " E e Vs. E oh of the photons with d<2Rm",140,0,100,140,0,100);

    
     if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
    
   //TGraph *energyThEl= new TGraph(nentries); 

    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
   
       //energyThEl->SetPoint(jentry,detKinBeamRot_the,detKinBeamRot_Ee);
       
    //if ((detKinBeamRot_Ee>3.5 || detKinBeamRot_the>12) && (detKinBeamRot_Ee>40 || detKinBeamRot_the>4) && (detKinBeamRot_Ee>(-168/3.4)+((105/3.4)*detKinBeamRot_thmu)) && ( detKinBeamRot_Ee>125 || detKinBeamRot_Ee<25+(100/3.5)*detKinBeamRot_thmu ) ){   
       
       if (detKinBeamRot_cooXe < 0.07 && detKinBeamRot_cooYe < 0.07 && detKinBeamRot_cooXe > -0.07 && detKinBeamRot_cooYe > -0.07)
        {   
// SE VIENE DAL TARGET ZERO E HA ENERGIA MAGGIORE DI 1 GEV
           if (detKinBeamRot_tar==0 && detKinBeamRot_Ee>1)
       
            
        {   nElNO++;
      
         
//SE CI SONO FOTONI
       if (photon_coox != -1)
       {

           d=sqrt((detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy));   
  
// SE QUESTI FOTONI SONO NEL TARGET 
      if (photon_coox < 0.07 && photon_cooy < 0.07 && photon_coox > -0.07 && photon_cooy > -0.07)
       {
           
// SE HANNO UNA DISTANZA DALL'ELETTRONE MAGGIORE DI 2 RAGGI DI MOLIERE
          if (abs(d)>2*Rm)
           {
               nEl++;
           }
// SE INVECE IL PUNTO D'IMPATTO è AD UNA DSTANZA MINORE
          else {
              nElNORm++; 
              EnCalNORm->Fill(detKinBeamRot_Ee,wgt_full); 
              ThCalNORm->Fill(detKinBeamRot_the,wgt_full); 
              Th_E_noRm->Fill(detKinBeamRot_the,detKinBeamRot_Ee,wgt_full); 
              E_ph->Fill(photon_energy,wgt_full);
              
// MUONI ASSOCIATI AD ELETTRONI CON FOTONI CON PUNTO D'IMPATTO <2*RM          
                if (abs(detKinBeamRot_cooXmu)<0.07 && abs(detKinBeamRot_cooYmu)<0.07)
                {
              ThCalNORmMU->Fill(detKinBeamRot_thmu,wgt_full);
              Th_E_noRmMU->Fill(detKinBeamRot_thmu,detKinBeamRot_Ee,wgt_full);
                   Th_E_eph->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full); 
                }}
       }
// SE I FOTONI NON CI SONO NEL CALROIMETRO MA SONO STATI PRODOTTI, QUINDI SONO EVENTI NLO
      else {nElPh++; 
            ThPhNoCal->Fill(detKinBeamRot_the,wgt_full);
            EnPhNoCal->Fill(detKinBeamRot_Ee,wgt_full);
            Th_E_PhNoCal->Fill(detKinBeamRot_the, detKinBeamRot_Ee,wgt_full);
        if (abs(detKinBeamRot_cooXmu)<0.07 && abs(detKinBeamRot_cooYmu)<0.07)
                { 
            ThPhNoCalMU->Fill(detKinBeamRot_thmu,wgt_full);
            Th_E_PhNoCalMU->Fill(detKinBeamRot_thmu, detKinBeamRot_Ee,wgt_full);
                   Th_E_mu->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
                }
           }     
           }     
// SE NON CI SONO FOTONI QUINDI EVENTI LO
else {nElNOph++;
      EnCalNoPh->Fill(detKinBeamRot_Ee,wgt_full); 
      ThCalNoPh->Fill(detKinBeamRot_the,wgt_full);
      Th_E_noph->Fill(detKinBeamRot_the, detKinBeamRot_Ee,wgt_full);
     
     if (abs(detKinBeamRot_cooXmu)<0.07 && abs(detKinBeamRot_cooYmu)<0.07)
     {
      ThCalNoPhMU->Fill(detKinBeamRot_thmu,wgt_full);
      Th_E_nophMU->Fill(detKinBeamRot_thmu, detKinBeamRot_Ee,wgt_full);
         Th_E_eNoph->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
     }
     }}
      
         
// SE VIENE DAL TARGET UNO
         if (detKinBeamRot_tar==1)
         { nElNO++;
      
         
//SE CI SONO FOTONI
       if (photon_coox != -1)
       {

           d=sqrt((detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy));   
  
// SE QUESTI FOTONI SONO NEL TARGET 
      if (photon_coox < 0.07 && photon_cooy < 0.07 && photon_coox > -0.07 && photon_cooy > -0.07)
       {
           
// SE HANNO UNA DISTANZA DALL'ELETTRONE MAGGIORE DI 2 RAGGI DI MOLIERE
          if (abs(d)>2*Rm)
           {
               nEl++;
           }
// SE INVECE IL PUNTO D'IMPATTO è AD UNA DSTANZA MINORE
          else {
              nElNORm++; 
              EnCalNORm->Fill(detKinBeamRot_Ee,wgt_full); 
              ThCalNORm->Fill(detKinBeamRot_the,wgt_full); 
              Th_E_noRm->Fill(detKinBeamRot_the,detKinBeamRot_Ee,wgt_full); 
              E_ph->Fill(photon_energy,wgt_full);
              
// MUONI ASSOCIATI AD ELETTRONI CON FOTONI CON PUNTO D'IMPATTO <2*RM          
                if (abs(detKinBeamRot_cooXmu)<0.07 && abs(detKinBeamRot_cooYmu)<0.07)
                {
              ThCalNORmMU->Fill(detKinBeamRot_thmu,wgt_full);
              Th_E_noRmMU->Fill(detKinBeamRot_thmu,detKinBeamRot_Ee,wgt_full);
                   Th_E_eph->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full); 
                }}
       }
// SE I FOTONI NON CI SONO NEL CALROIMETRO MA SONO STATI PRODOTTI, QUINDI SONO EVENTI NLO
      else {nElPh++; 
            ThPhNoCal->Fill(detKinBeamRot_the,wgt_full); 
            EnPhNoCal->Fill(detKinBeamRot_Ee,wgt_full);
            Th_E_PhNoCal->Fill(detKinBeamRot_the, detKinBeamRot_Ee,wgt_full);
        if (abs(detKinBeamRot_cooXmu)<0.07 && abs(detKinBeamRot_cooYmu)<0.07)
                {
            ThPhNoCalMU->Fill(detKinBeamRot_thmu,wgt_full);
            Th_E_PhNoCalMU->Fill(detKinBeamRot_thmu, detKinBeamRot_Ee,wgt_full);
            Th_E_mu->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
                }
           }     
           }     
// SE NON CI SONO FOTONI QUINDI EVENTI LO
else {nElNOph++;
      EnCalNoPh->Fill(detKinBeamRot_Ee,wgt_full); 
      ThCalNoPh->Fill(detKinBeamRot_the,wgt_full);
      Th_E_noph->Fill(detKinBeamRot_the, detKinBeamRot_Ee,wgt_full);
     
     if (abs(detKinBeamRot_cooXmu)<0.07 && abs(detKinBeamRot_cooYmu)<0.07)
     {
         
      ThCalNoPhMU->Fill(detKinBeamRot_thmu,wgt_full);
      Th_E_nophMU->Fill(detKinBeamRot_thmu, detKinBeamRot_Ee,wgt_full);
      Th_E_eNoph->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
     }
     }}
           
       }
              
  
        
        
if (detKinBeamRot_cooXe < 0.07 && detKinBeamRot_cooYe < 0.07 && detKinBeamRot_cooXe > -0.07 && detKinBeamRot_cooYe > -0.07 && detKinBeamRot_tar==1)
        {
       nElNO1++;
              if (photon_coox != -1)
              {
          
d=sqrt((detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy)); 
  
       if (photon_coox < 0.07 && photon_cooy < 0.07 && photon_coox > -0.07 && photon_cooy > -0.07)
       {
           if (abs(d)>2*Rm)
           {
               nEl1++;
           }
           else {nElNORm1++; 
                 EnCalNORm1->Fill(detKinBeamRot_Ee,wgt_full); 
                 ThCalNORm1->Fill(detKinBeamRot_the,wgt_full); 
                 Th_E_noRm1->Fill(detKinBeamRot_the, detKinBeamRot_Ee,wgt_full);
                 E_ph1->Fill(photon_energy,wgt_full);
                 Ee_Eph->Fill(detKinBeamRot_Ee,photon_energy,wgt_full);
                
                 Double_t Etot=photon_energy+detKinBeamRot_Ee;
                En_tot->Fill(Etot,wgt_full);
// MUONI ASSOCIATI AD ELETTRONI CON FOTONI CON PUNTO D'IMPATTO <2*RM          
                if (abs(detKinBeamRot_cooXmu)<0.07 && abs(detKinBeamRot_cooYmu)<0.07)
                {
                 ThCalNORm1MU->Fill(detKinBeamRot_thmu,wgt_full); 
                 Th_E_noRm1MU->Fill(detKinBeamRot_thmu, detKinBeamRot_Ee,wgt_full); 
                   Th_E_eph1->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full); 
                }
                }
       }
   else  {nElPh1++; 
          ThPhNoCal1->Fill(detKinBeamRot_the,wgt_full); 
          EnPhNoCal1->Fill(detKinBeamRot_Ee,wgt_full); 
         Th_E_PhNoCal1->Fill(detKinBeamRot_the, detKinBeamRot_Ee,wgt_full);
      
        if (abs(detKinBeamRot_cooXmu)<0.07 && abs(detKinBeamRot_cooYmu)<0.07)
                {   
                   Th_E_PhNoCal1MU->Fill(detKinBeamRot_thmu, detKinBeamRot_Ee,wgt_full);
                   Th_E_mu1->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
                   ThPhNoCal1MU->Fill(detKinBeamRot_thmu,wgt_full); 
                }
         }}
    else {nElNOph1++; 
          EnCalNoPh1->Fill(detKinBeamRot_Ee,wgt_full); 
          ThCalNoPh1->Fill(detKinBeamRot_the,wgt_full); 
          Th_E_noph1->Fill(detKinBeamRot_the, detKinBeamRot_Ee,wgt_full);
              if (abs(detKinBeamRot_cooXmu)<0.07 && abs(detKinBeamRot_cooYmu)<0.07)
     {
          ThCalNoPh1MU->Fill(detKinBeamRot_thmu,wgt_full); 
          Th_E_noph1MU->Fill(detKinBeamRot_thmu, detKinBeamRot_Ee,wgt_full);
         Th_E_eNoph1->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
     }}

           }
        
        
    
   if (detKinBeamRot_cooXe < 0.07 && detKinBeamRot_cooYe < 0.07 && detKinBeamRot_cooXe > -0.07 && detKinBeamRot_cooYe > -0.07 && detKinBeamRot_tar==0 && detKinBeamRot_Ee>1)
        {
       nElNO0++;
              if (photon_coox != -1)
              {
          
d=sqrt((detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy)); 
  
       if (photon_coox < 0.07 && photon_cooy < 0.07 && photon_coox > -0.07 && photon_cooy > -0.07)
       {
           if (abs(d)>2*Rm)
           {
               nEl0++;
           }
           else {nElNORm0++; 
                 EnCalNORm0->Fill(detKinBeamRot_Ee,wgt_full); 
                 ThCalNORm0->Fill(detKinBeamRot_the,wgt_full); 
                 Th_E_noRm0->Fill(detKinBeamRot_the, detKinBeamRot_Ee,wgt_full);
                 E_ph0->Fill(photon_energy,wgt_full);
                 Ee_Eph->Fill(detKinBeamRot_Ee,photon_energy,wgt_full);
                
                 Double_t Etot=photon_energy+detKinBeamRot_Ee;
                En_tot->Fill(Etot,wgt_full);
// MUONI ASSOCIATI AD ELETTRONI CON FOTONI CON PUNTO D'IMPATTO <2*RM          
                if (abs(detKinBeamRot_cooXmu)<0.07 && abs(detKinBeamRot_cooYmu)<0.07)
                {
                    ThCalNORm0MU->Fill(detKinBeamRot_thmu,wgt_full); 
                    Th_E_noRm0MU->Fill(detKinBeamRot_thmu, detKinBeamRot_Ee,wgt_full); 
                    Th_E_eph0->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full); 
                }
                }
       }
   else  {nElPh0++; 
          ThPhNoCal0->Fill(detKinBeamRot_the,wgt_full);
          EnPhNoCal0->Fill(detKinBeamRot_Ee,wgt_full); 
         Th_E_PhNoCal0->Fill(detKinBeamRot_the, detKinBeamRot_Ee,wgt_full);

        if (abs(detKinBeamRot_cooXmu)<0.07 && abs(detKinBeamRot_cooYmu)<0.07)
                { 
            ThPhNoCal0MU->Fill(detKinBeamRot_thmu,wgt_full); 
            Th_E_PhNoCal0MU->Fill(detKinBeamRot_thmu, detKinBeamRot_Ee,wgt_full);
            Th_E_mu0->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
                }
         }}
    else {nElNOph0++; 
          EnCalNoPh0->Fill(detKinBeamRot_Ee,wgt_full); 
          ThCalNoPh0->Fill(detKinBeamRot_the,wgt_full); 
          Th_E_noph0->Fill(detKinBeamRot_the, detKinBeamRot_Ee,wgt_full);
            if (abs(detKinBeamRot_cooXmu)<0.07 && abs(detKinBeamRot_cooYmu)<0.07)
     {  
          ThCalNoPh0MU->Fill(detKinBeamRot_thmu,wgt_full); 
          Th_E_noph0MU->Fill(detKinBeamRot_thmu, detKinBeamRot_Ee,wgt_full);
         Th_E_eNoph0->Fill(detKinBeamRot_the,detKinBeamRot_thmu,wgt_full);
     }}
  
           }
            
        
        
       
   }
   //}

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
    cout << "numero di elettroni senza fotoni: " << nElNOph << endl;
    cout<<endl;
    
    cout << "numero di elettroni nel calorimetro dal tar 0: " << nElNO0 << endl;
    cout << "numero di elettroni con fotoni fuori dal calorimetro dal tar 0: " << nElPh0 << endl;
    cout << "numero di elettroni con fotoni nel calorimetro dal tar 0: " << nEl0 << endl;
    cout << "numero di elettroni con fotoni senza richiesta raggio Moliere dal tar 0: " << nElNORm0 << endl;
    cout << "numero di elettroni senza fotoni dal tar 0: " << nElNOph0 << endl;
    
    cout<<endl;
    cout << "numero di elettroni nel calorimetro dal tar 1: " << nElNO1 << endl;
    cout << "numero di elettroni con fotoni fuori dal calorimetro dal tar 1: " << nElPh1 << endl;
    cout << "numero di elettroni con fotoni nel calorimetro dal tar 1: " << nEl1 << endl;
    cout << "numero di elettroni con fotoni senza richiesta raggio Moliere dal tar 1: " << nElNORm1 << endl;
    cout << "numero di elettroni senza fotoni: " << nElNOph1 << endl;
    
    
    
    TCanvas * e= new TCanvas("e","e",200,10,1000,900);
    e->Divide(1,3);
    e->cd(1);
    EnCalNoPh->SetLineColor(kRed);
    EnCalNoPh->Draw("HIST");
    EnCalNORm->Draw("HIST same");
    EnPhNoCal->SetLineColor(kOrange);
    EnPhNoCal->Draw("HIST same");
    gPad->SetLogx();
    
    e->cd(2);
    EnCalNoPh0->SetLineColor(kRed);
    EnCalNoPh0->Draw("HIST");
    EnCalNORm0->Draw("HIST same");
    EnPhNoCal0->SetLineColor(kOrange);
    EnPhNoCal0->Draw("HIST same");
    gPad->SetLogx();
    
    e->cd(3);
    EnCalNoPh1->SetLineColor(kRed);
    EnCalNoPh1->Draw("HIST");
    EnCalNORm1->Draw("HIST same");
    EnPhNoCal1->SetLineColor(kOrange);
    EnPhNoCal1->Draw("HIST same");
    gPad->SetLogx();
    
    /*e->cd(4);
    E_ph->SetLineColor(31);
    E_ph->Draw("HIST");
    E_ph0->Draw("HIST same");
    E_ph1->SetLineColor(49);
    E_ph1->Draw("HIST same");
    gPad->SetLogx();
 
    e->cd(5);
    En_tot->Draw("HIST");
    gPad->SetLogx();*/
    
    e->SaveAs("EnergyElnoRm.png");
    
    TCanvas * te= new TCanvas("te","te",200,10,1000,1000);
    te->Divide(1,3);
    te->cd(1);
    ThCalNoPh->SetLineColor(kRed);
    ThCalNoPh->Draw("HIST");
    ThCalNORm->Draw("HIST same");
    ThPhNoCal->SetLineColor(kOrange);
    ThPhNoCal->Draw("HIST same");

    
    te->cd(2);
    ThCalNoPh0->SetLineColor(kRed);
    ThCalNoPh0->Draw("HIST");
    ThCalNORm0->Draw("HIST same");
    ThPhNoCal0->SetLineColor(kOrange);
    ThPhNoCal0->Draw("HIST same");
   
    
    te->cd(3);
    ThCalNoPh1->SetLineColor(kRed);
    ThCalNoPh1->Draw("HIST");
    ThCalNORm1->Draw("HIST same");
    ThPhNoCal1->SetLineColor(kOrange);
    ThPhNoCal1->Draw("HIST same");
    
    te->SaveAs("ThElnoRm.png");

    TCanvas * teMU= new TCanvas("teMU","teMU",200,10,1000,1000);
    teMU->Divide(1,3);
    teMU->cd(1);
    ThCalNoPhMU->SetLineColor(kRed);
    ThCalNoPhMU->Draw("HIST");
    ThCalNORmMU->Draw("HIST same");
    ThPhNoCalMU->SetLineColor(kOrange);
    ThPhNoCalMU->Draw("HIST same");

    
    teMU->cd(2);
    ThCalNoPh0MU->SetLineColor(kRed);
    ThCalNoPh0MU->Draw("HIST");
    ThCalNORm0MU->Draw("HIST same");
    ThPhNoCal0MU->SetLineColor(kOrange);
    ThPhNoCal0MU->Draw("HIST same");
   
    
    teMU->cd(3);
    ThCalNoPh1MU->SetLineColor(kRed);
    ThCalNoPh1MU->Draw("HIST");
    ThCalNORm1MU->Draw("HIST same");
    ThPhNoCal1MU->SetLineColor(kOrange);
    ThPhNoCal1MU->Draw("HIST same");
    
    teMU->SaveAs("ThElnoRmMU.png");
    
    TCanvas * dued1= new TCanvas("dued1","dued1",1000,100,2500,2000);
    Th_E_noph1->SetMarkerColor(kRed);
    Th_E_noph1->Draw("HIST");
    Th_E_noph1->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noph1->GetYaxis()->SetTitle("E [GeV]");

    Th_E_noRm1->Draw("HIST same");
    Th_E_noRm1->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noRm1->GetYaxis()->SetTitle("E [GeV]");
        
    Th_E_PhNoCal1->SetMarkerColor(kOrange);
    Th_E_PhNoCal1->Draw("HIST same");
    Th_E_PhNoCal1->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_PhNoCal1->GetYaxis()->SetTitle("E [GeV]");
  dued1->SaveAs("Eth1.png");
    
    TCanvas * dued= new TCanvas("dued","dued",1000,100,2500,2000);

    Th_E_noph->SetMarkerColor(kRed);
    Th_E_noph->Draw("HIST");
    Th_E_noph->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noph->GetYaxis()->SetTitle("E [GeV]");

    Th_E_noRm->Draw("HIST same");
    Th_E_noRm->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noRm->GetYaxis()->SetTitle("E [GeV]");
        
    Th_E_PhNoCal->SetMarkerColor(kOrange);
    Th_E_PhNoCal->Draw("HIST same");
    Th_E_PhNoCal->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_PhNoCal->GetYaxis()->SetTitle("E [GeV]");
  dued->SaveAs("Eth.png");
    
    TCanvas * dued0= new TCanvas("dued0","dued0",1000,100,2500,2000);

    Th_E_noph0->SetMarkerColor(kRed);
    Th_E_noph0->Draw("HIST");
    Th_E_noph0->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noph0->GetYaxis()->SetTitle("E [GeV]");

    Th_E_noRm0->Draw("HIST same");
    Th_E_noRm0->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noRm0->GetYaxis()->SetTitle("E [GeV]");
        
    Th_E_PhNoCal0->SetMarkerColor(kOrange);
    Th_E_PhNoCal0->Draw("HIST same");
    Th_E_PhNoCal0->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_PhNoCal0->GetYaxis()->SetTitle("E [GeV]");
  dued0->SaveAs("Eth0.png");
    
    
    
        TCanvas * separati= new TCanvas("dued1","dued1",1000,100,2500,2000);
    separati->Divide(3,3);
    
    separati->cd(1);
    Th_E_noph1->SetMarkerColor(kRed);
    Th_E_noph1->Draw("HIST");
    Th_E_noph1->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noph1->GetYaxis()->SetTitle("E [GeV]");
    
    separati->cd(2);
    Th_E_noRm1->Draw("HIST");
    Th_E_noRm1->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noRm1->GetYaxis()->SetTitle("E [GeV]");
    
    separati->cd(3); 
    Th_E_PhNoCal1->SetMarkerColor(kOrange);
    Th_E_PhNoCal1->Draw("HIST ");
    Th_E_PhNoCal1->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_PhNoCal1->GetYaxis()->SetTitle("E [GeV]");

    separati->cd(4);
    Th_E_noph->SetMarkerColor(kRed);
    Th_E_noph->Draw("HIST");
    Th_E_noph->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noph->GetYaxis()->SetTitle("E [GeV]");

    separati->cd(5);
    Th_E_noRm->Draw("HIST");
    Th_E_noRm->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noRm->GetYaxis()->SetTitle("E [GeV]");
        
    separati->cd(6);
    Th_E_PhNoCal->SetMarkerColor(kOrange);
    Th_E_PhNoCal->Draw("HIST");
    Th_E_PhNoCal->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_PhNoCal->GetYaxis()->SetTitle("E [GeV]");
    
  
    separati->cd(7);
    Th_E_noph0->SetMarkerColor(kRed);
    Th_E_noph0->Draw("HIST");
    Th_E_noph0->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noph0->GetYaxis()->SetTitle("E [GeV]");

    separati->cd(8);
    Th_E_noRm0->Draw("HIST");
    Th_E_noRm0->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noRm0->GetYaxis()->SetTitle("E [GeV]");
        
    separati->cd(9);
    Th_E_PhNoCal0->SetMarkerColor(kOrange);
    Th_E_PhNoCal0->Draw("HIST same");
    Th_E_PhNoCal0->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_PhNoCal0->GetYaxis()->SetTitle("E [GeV]");
  separati->SaveAs("Energy-theta.png");
       
    TCanvas * separatiA= new TCanvas("A","A",1000,100,2500,2000);
    separatiA->Divide(3,3);
    
    separatiA->cd(1);
    Th_E_eNoph1->SetMarkerColor(46);
    Th_E_eNoph1->Draw("HIST");
    Th_E_eNoph1->GetXaxis()->SetTitle("Th el[mrad]");
    Th_E_eNoph1->GetYaxis()->SetTitle("Th mu [mrad]");
    
    separatiA->cd(2);
    Th_E_eph1->Draw("HIST");
    Th_E_eph1->GetXaxis()->SetTitle("Th el[mrad]");
    Th_E_eph1->GetYaxis()->SetTitle("Th mu [mrad]");
    
    separatiA->cd(3); 
    Th_E_mu1->SetMarkerColor(29);
    Th_E_mu1->Draw("HIST ");
    Th_E_mu1->GetXaxis()->SetTitle("Th el[mrad]");
    Th_E_mu1->GetYaxis()->SetTitle("Th mu [mrad]");

    separatiA->cd(4);
    Th_E_eNoph->SetMarkerColor(46);
    Th_E_eNoph->Draw("HIST");
    Th_E_eNoph->GetXaxis()->SetTitle("Th el[mrad]");
    Th_E_eNoph->GetYaxis()->SetTitle("Th mu [mrad]");

    separatiA->cd(5);
    Th_E_eph->Draw("HIST");
    Th_E_eph->GetXaxis()->SetTitle("Th el[mrad]");
    Th_E_eph->GetYaxis()->SetTitle("Th mu [mrad]");
        
    separatiA->cd(6);
    Th_E_mu->SetMarkerColor(29);
    Th_E_mu->Draw("HIST");
    Th_E_mu->GetXaxis()->SetTitle("Th el[mrad]");
    Th_E_mu->GetYaxis()->SetTitle("Th mu [mrad]");

    
  
    separatiA->cd(7);
    Th_E_eNoph0->SetMarkerColor(46);
    Th_E_eNoph0->Draw("HIST");
    Th_E_eNoph0->GetXaxis()->SetTitle("Th el [mrad]");
    Th_E_eNoph0->GetYaxis()->SetTitle("Th mu [mrad]");

    separatiA->cd(8);
    Th_E_eph0->Draw("HIST");
    Th_E_eph0->GetXaxis()->SetTitle("Th el [mrad]");
    Th_E_eph0->GetYaxis()->SetTitle("Th mu [mrad]");
        
    separatiA->cd(9);
    Th_E_mu0->SetMarkerColor(29);
    Th_E_mu0->Draw("HIST same");
    Th_E_mu0->GetXaxis()->SetTitle("Th el [mrad]");
    Th_E_mu0->GetYaxis()->SetTitle("Th mu [mrad]");
  separatiA->SaveAs("thetae-thetamu.png");   
    
TCanvas * E= new TCanvas("A","A",1000,100,2500,2000);
    
    Ee_Eph->Draw("HIST");
  E->SaveAs("Ee-Eph.png"); 



        TCanvas * separatiMU= new TCanvas("dued1MU","dued1MU",1000,100,2500,2000);
    separatiMU->Divide(3,3);
    
    separatiMU->cd(1);
    Th_E_noph1MU->SetMarkerColor(kRed);
    Th_E_noph1MU->Draw("HIST");
    Th_E_noph1MU->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noph1MU->GetYaxis()->SetTitle("E [GeV]");
    
    separatiMU->cd(2);
    Th_E_noRm1MU->Draw("HIST");
    Th_E_noRm1MU->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noRm1MU->GetYaxis()->SetTitle("E [GeV]");
    
    separatiMU->cd(3); 
    Th_E_PhNoCal1MU->SetMarkerColor(kOrange);
    Th_E_PhNoCal1MU->Draw("HIST ");
    Th_E_PhNoCal1MU->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_PhNoCal1MU->GetYaxis()->SetTitle("E [GeV]");

    separatiMU->cd(4);
    Th_E_nophMU->SetMarkerColor(kRed);
    Th_E_nophMU->Draw("HIST");
    Th_E_nophMU->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_nophMU->GetYaxis()->SetTitle("E [GeV]");

    separatiMU->cd(5);
    Th_E_noRmMU->Draw("HIST");
    Th_E_noRmMU->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noRmMU->GetYaxis()->SetTitle("E [GeV]");
        
    separatiMU->cd(6);
    Th_E_PhNoCalMU->SetMarkerColor(kOrange);
    Th_E_PhNoCalMU->Draw("HIST");
    Th_E_PhNoCalMU->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_PhNoCalMU->GetYaxis()->SetTitle("E [GeV]");
    
  
    separatiMU->cd(7);
    Th_E_noph0MU->SetMarkerColor(kRed);
    Th_E_noph0MU->Draw("HIST");
    Th_E_noph0MU->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noph0MU->GetYaxis()->SetTitle("E [GeV]");

    separatiMU->cd(8);
    Th_E_noRm0MU->Draw("HIST");
    Th_E_noRm0MU->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_noRm0MU->GetYaxis()->SetTitle("E [GeV]");
        
    separatiMU->cd(9);
    Th_E_PhNoCal0MU->SetMarkerColor(kOrange);
    Th_E_PhNoCal0MU->Draw("HIST same");
    Th_E_PhNoCal0MU->GetXaxis()->SetTitle("Th [mrad]");
    Th_E_PhNoCal0MU->GetYaxis()->SetTitle("E [GeV]");
  separatiMU->SaveAs("Energye-thetaMU.png");
    
       
    TCanvas * incrocio= new TCanvas("incrocio","incrocio",1000,100,2500,2000);
    Th_E_noRm->Draw("HIST");
    Th_E_PhNoCal->Draw("HIST same");
    Th_E_noRmMU->SetMarkerColor(kRed);
    Th_E_PhNoCalMU->SetMarkerColor(kRed);
    Th_E_noRmMU->Draw("HIST same");
    Th_E_PhNoCalMU->Draw("HIST same");
  incrocio->SaveAs("incrocioEnergye-thetaNLO.png");
    
        TCanvas * incrocioLO= new TCanvas("incrocioLO","incrocioLO",1000,100,2500,2000);
    //Th_E_noph->SetXaxis
    //Th_E_noph->SetXaxis
    Th_E_noph->Draw("HIST");
    Th_E_nophMU->SetMarkerColor(kBlack);
    Th_E_nophMU->Draw("HIST same");
  incrocioLO->SaveAs("incrocioEnergye-thetaLO.png");
    
}