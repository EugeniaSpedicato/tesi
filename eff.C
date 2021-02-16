#define atree_cxx
#include "atree.h"
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>

#include <TGraph.h>
#include <TTree.h>
#include <cmath>
using namespace std;

#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{
    TH1::SetDefaultSumw2();


Double_t E_CAL=0.;

TH1F* The_trueCUT=new TH1F("th", "th El true CUT", 120,0,30);
TH1F* The_trueCUT1=new TH1F("th", "th El true TAR 1 CUT", 120,0,30);
TH1F* The_trueCUT2=new TH1F("th", "th El true TAR 2 CUT", 120,0,30);

TH1F* The_trueCUTmu=new TH1F("th", "th El  true CUT th_mu>0.2mrad", 120,0,30);
TH1F* The_trueCUT1mu=new TH1F("th", "th El true TAR 1 CUT th_mu>0.2mrad", 120,0,30);
TH1F* The_trueCUT2mu=new TH1F("th", "th El true TAR 2 CUT th_mu>0.2mrad", 120,0,30);
    
TH1F* The_trueCUTEe=new TH1F("th", "th El true CUT on E_e", 120,0,30);
TH1F* The_trueCUT1Ee=new TH1F("th", "th El true TAR 1 CUT on E_e", 120,0,30);
TH1F* The_trueCUT2Ee=new TH1F("th", "th El true TAR 2 CUT on E_e", 120,0,30);
    
TH1F* The_trueCUTtot=new TH1F("th", "th El true CUT th_mu+Ee", 120,0,30);
TH1F* The_trueCUT1tot=new TH1F("th", "th El true TAR 1 CUT th_mu+Ee", 120,0,30);
TH1F* The_trueCUT2tot=new TH1F("th", "th El true TAR 2 CUT th_mu+Ee", 120,0,30);

TH1F* The_true=new TH1F("th", "th El true", 120,0,30);    
TH1F* The=new TH1F("th", "th El core", 120,0,30); 
/*TH1F* TheBIG=new TH1F("th", "th El 5X5", 120,0,30); 
TH1F* The2P=new TH1F("th", "th El crown", 120,0,30); */

TH1F* The_true1=new TH1F("th", "th El true TAR 1", 120,0,30);    
TH1F* The1=new TH1F("th", "th El core TAR 1", 120,0,30); 
/*TH1F* TheBIG1=new TH1F("th", "th El 5X5 TAR 1", 120,0,30); 
TH1F* The2P1=new TH1F("th", "th El crown TAR 1", 120,0,30); */

    
TH1F* The_true2=new TH1F("th", "th El true TAR 2", 120,0,30);     
TH1F* The2=new TH1F("th", "th El core TAR 2", 120,0,30); 
/*TH1F* TheBIG2=new TH1F("th", "th El 5X5 TAR 2", 120,0,30); 
TH1F* The2P2=new TH1F("th", "th El crown TAR 2", 120,0,30);*/ 
    
    
TH1F* TheCUT=new TH1F("th", "th El core CUT", 120,0,30); 
TH1F* The1CUT=new TH1F("th", "th El TAR 1 core CUT", 120,0,30); 
TH1F* The2CUT=new TH1F("th", "th El TAR 2 core CUT", 120,0,30); 

/*TH1F* TheBIGCUT=new TH1F("th", "th El 5X5 CUT", 120,0,30); 
TH1F* The2PCUT=new TH1F("th", "th El crown CUT", 120,0,30);*/ 
    
    
TH1F* TheCUTmu=new TH1F("th", "th El core CUT th_mu>0.2mrad", 120,0,30); 
TH1F* The1CUTmu=new TH1F("th", "th El TAR 1 core CUT th_mu>0.2mrad", 120,0,30); 
TH1F* The2CUTmu=new TH1F("th", "th El TAR 2 core CUT th_mu>0.2mrad", 120,0,30); 

/*TH1F* TheBIG1CUT=new TH1F("th", "th El 5X5 TAR 1 CUT", 120,0,30); 
TH1F* The2P1CUT=new TH1F("th", "th El crown TAR 1 CUT", 120,0,30); */
    
    
TH1F* TheCUTEe=new TH1F("th", "th El core CUT on E_e", 120,0,30); 
TH1F* The1CUTEe=new TH1F("th", "th El core TAR 1 CUT on E_e", 120,0,30); 
TH1F* The2CUTEe=new TH1F("th", "th El core TAR 2 CUT on E_e", 120,0,30); 
    
/*TH1F* TheBIG2CUT=new TH1F("th", "th El 5X5 TAR 2 CUT", 120,0,30);
TH1F* The2P2CUT=new TH1F("th", "th El crown TAR 2 CUT", 120,0,30); */

TH1F* TheCUTtot=new TH1F("th", "th El core CUT th_mu+Ee", 120,0,30); 
TH1F* The1CUTtot=new TH1F("th", "th El core TAR 1 CUT th_mu+Ee", 120,0,30); 
TH1F* The2CUTtot=new TH1F("th", "th El core TAR 2 CUT th_mu+Ee", 120,0,30); 
    
TH1F* rmu=new TH1F("rmu", "impact point", 100,0,10); 

TH2F  *E3x31  = new TH2F("ThEel1" , " Th_el Vs. E_E3x3 core TAR1",100,0,50,280,0,140);
TH2F  *E3x32  = new TH2F("ThEel2" , " Th_el Vs. E_E3x3 core TAR2",100,0,50,280,0,140);
    

TH2F  *E3x31CUT  = new TH2F("ThEel1" , " Th_el Vs. E_E3x3 core TAR1 CUT",100,0,50,280,0,140);
TH2F  *E3x32CUT  = new TH2F("ThEel2" , " Th_el Vs. E_E3x3 core TAR2 CUT",100,0,50,280,0,140);
    
TH2F  *E3x31CUTmu  = new TH2F("ThEel1" , " Th_el Vs. E_E3x3 core TAR1 CUT th_mu>0.2mrad",100,0,50,280,0,140);
TH2F  *E3x32CUTmu  = new TH2F("ThEel2" , " Th_el Vs. E_E3x3 core TAR2 CUT th_mu>0.2mrad",100,0,50,280,0,140);
    
TH2F  *E3x31CUTEe  = new TH2F("ThEel1" , " Th_el Vs. E_E3x3 core TAR1 CUT on E_e",100,0,50,280,0,140);
TH2F  *E3x32CUTEe  = new TH2F("ThEel2" , " Th_el Vs. E_E3x3 core TAR2 CUT on E_e",100,0,50,280,0,140);
    
TH2F  *E3x31CUTtot  = new TH2F("ThEel1" , " Th_el Vs. E_E3x3 core TAR1 CUT th_mu+Ee",100,0,50,280,0,140);
TH2F  *E3x32CUTtot  = new TH2F("ThEel2" , " Th_el Vs. E_E3x3 core TAR2 CUT th_mu+Ee",100,0,50,280,0,140);

    
    
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
    {E_CAL=detKinBeamRot_Ee;}

        detKinBeamRot_cooXe=detKinBeamRot_cooXe*100; // cm
        detKinBeamRot_cooYe=detKinBeamRot_cooYe*100; // cm
        detKinBeamRot_x_in=detKinBeamRot_x_in*100; // cm coo muone entrante
        detKinBeamRot_y_in=detKinBeamRot_y_in*100; // cm coo muone entrante
        photon_coox=photon_coox*100; // cm
        photon_cooy=photon_cooy*100; // cm
               
       double r_mu=sqrt((detKinBeamRot_x_in*detKinBeamRot_x_in)+(detKinBeamRot_y_in*detKinBeamRot_y_in));
              

       
The_true->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3) The_trueCUT->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3) The_trueCUTmu->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_E_clus3x3>1) The_trueCUTEe->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3 && detKinBeamRot_E_clus3x3>1) The_trueCUTtot->Fill(detKinBeamRot_def_angle_e,wgt_full);
    

    if (detKinBeamRot_tar==0)
    {The_true1->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3) The_trueCUT1->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3) The_trueCUT1mu->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_E_clus3x3>1) The_trueCUT1Ee->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3 && detKinBeamRot_E_clus3x3>1) The_trueCUT1tot->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    if (detKinBeamRot_tar==1)
    {The_true2->Fill(detKinBeamRot_def_angle_e,wgt_full);
     
if(r_mu<3) The_trueCUT2->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3) The_trueCUT2mu->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_E_clus3x3>1) The_trueCUT2Ee->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3 && detKinBeamRot_E_clus3x3>1) The_trueCUT2tot->Fill(detKinBeamRot_def_angle_e,wgt_full);}
       
    if (detKinBeamRot_tar==0) rmu->Fill(r_mu,wgt_full);
       
if (detKinBeamRot_n_cell_e!=0 && E_CAL>0)
{  

    /*TheBIG->Fill(detKinBeamRot_def_angle_e,wgt_full);
    TheBIGCUT->Fill(detKinBeamRot_def_angle_e,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {TheBIG1->Fill(detKinBeamRot_def_angle_e,wgt_full);
    TheBIG1CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    if (detKinBeamRot_tar==1)
    {TheBIG2->Fill(detKinBeamRot_def_angle_e,wgt_full);
    TheBIG2CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);}*/

    
if(abs(detKinBeamRot_cooXe)<4.275 && abs(detKinBeamRot_cooYe)<4.275) {

The->Fill(detKinBeamRot_def_angle_e,wgt_full);
    
if(r_mu<3) TheCUT->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3) TheCUTmu->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_E_clus3x3>1) TheCUTEe->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3 && detKinBeamRot_E_clus3x3>1)TheCUTtot->Fill(detKinBeamRot_def_angle_e,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {The1->Fill(detKinBeamRot_def_angle_e,wgt_full);
     
if(r_mu<3) The1CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3) The1CUTmu->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_E_clus3x3>1) The1CUTEe->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3 && detKinBeamRot_E_clus3x3>1)The1CUTtot->Fill(detKinBeamRot_def_angle_e,wgt_full);

     if (detKinBeamRot_E_clus3x3!=0) 
     {
E3x31->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);
if(r_mu<3) E3x31CUT->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3) E3x31CUTmu->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);
if(r_mu<3 && detKinBeamRot_E_clus3x3>1) E3x31CUTEe->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3 && detKinBeamRot_E_clus3x3>1) E3x31CUTtot->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);
} }
    if (detKinBeamRot_tar==1)
    {The2->Fill(detKinBeamRot_def_angle_e,wgt_full);
     
if(r_mu<3) The2CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3) The2CUTmu->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_E_clus3x3>1) The2CUTEe->Fill(detKinBeamRot_def_angle_e,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3 && detKinBeamRot_E_clus3x3>1)The2CUTtot->Fill(detKinBeamRot_def_angle_e,wgt_full);
     if (detKinBeamRot_E_clus3x3!=0) 
     {
E3x32->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);
if(r_mu<3) E3x32CUT->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3) E3x32CUTmu->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);
if(r_mu<3 && detKinBeamRot_E_clus3x3>1) E3x32CUTEe->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);
if(r_mu<3 && detKinBeamRot_def_angle_mu>0.3 && detKinBeamRot_E_clus3x3>1) E3x32CUTtot->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);
} 
    
    } /*else { 

    The2P->Fill(detKinBeamRot_def_angle_e,wgt_full);
    The2PCUT->Fill(detKinBeamRot_def_angle_e,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {The2P1->Fill(detKinBeamRot_def_angle_e,wgt_full);
    The2P1CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    if (detKinBeamRot_tar==1)
    {The2P2->Fill(detKinBeamRot_def_angle_e,wgt_full);
    The2P2CUT->Fill(detKinBeamRot_def_angle_e,wgt_full);}
}*/
     
}}}
TH1F *Eff1 = new TH1F("ef1", "Eff Tar1", 120,0,30);
TH1F *Eff2 = new TH1F("ef2", "Eff Tar2", 120,0,30);
Eff1->Divide(The1,The_true1,1,1,"B");
Eff2->Divide(The2,The_true2,1,1,"B");
    
TH1F *Eff1CUT = new TH1F("ef1cut", "Eff Tar1 cut", 120,0,30);
TH1F *Eff2CUT = new TH1F("ef2cut", "Eff Tar2 cut", 120,0,30);
Eff1CUT->Divide(The1CUT,The_trueCUT1,1,1,"B");
Eff2CUT->Divide(The2CUT,The_trueCUT2,1,1,"B");

TH1F *Eff1CUTmu = new TH1F("ef1cutmu", "Eff Tar1 cut th_mu>0.2mrad", 120,0,30);
TH1F *Eff2CUTmu = new TH1F("ef2cutmu", "Eff Tar2 cut th_mu>0.2mrad", 120,0,30);
Eff1CUTmu->Divide(The1CUTmu,The_trueCUT1mu,1,1,"B");
Eff2CUTmu->Divide(The2CUTmu,The_trueCUT2mu,1,1,"B");
    
TH1F *Eff1CUTEe = new TH1F("ef1cutEe", "Eff Tar1 cut E_e", 120,0,30);
TH1F *Eff2CUTEe = new TH1F("ef2cutEe", "Eff Tar2 cut E_e", 120,0,30);
Eff1CUTEe->Divide(The1CUTEe,The_trueCUT1Ee,1,1,"B");
Eff2CUTEe->Divide(The2CUTEe,The_trueCUT2Ee,1,1,"B");

TH1F *Eff1CUTtot = new TH1F("ef1cuttot", "Eff Tar1 cut th_mu+Ee", 120,0,30);
TH1F *Eff2CUTtot = new TH1F("ef2cuttot", "Eff Tar2 cut th_mu+Ee", 120,0,30);
Eff1CUTtot->Divide(The1CUTtot,The_trueCUT1tot,1,1,"B");
Eff2CUTtot->Divide(The2CUTtot,The_trueCUT2tot,1,1,"B");

       

TCanvas * ef= new TCanvas("ef","ef",1000,100,2500,2000);
ef->Divide(1,2);
ef->cd(1);
Eff1->GetXaxis()->SetTitle("Theta el[mrad]");
Eff1->GetYaxis()->SetTitle("Efficency");
Eff1->SetLineWidth(3);
Eff1->SetLineColor(kBlue);
Eff1->SetMaximum(1.2);
//Eff1->SetMinimum(0);
Eff1->Draw();  
Eff1CUT->GetXaxis()->SetTitle("Theta el[mrad]");
Eff1CUT->GetYaxis()->SetTitle("Efficency");
Eff1CUT->SetLineWidth(3);
Eff1CUT->SetLineColor(kRed);
Eff1CUT->SetMaximum(1.2);
//Eff1CUT->SetMinimum(0);
Eff1CUT->Draw("same"); 

Eff1CUTmu->GetXaxis()->SetTitle("Theta el[mrad]");
Eff1CUTmu->GetYaxis()->SetTitle("Efficency");
Eff1CUTmu->SetLineWidth(3);
Eff1CUTmu->SetLineColor(kOrange);
Eff1CUTmu->SetMaximum(1.2);
//Eff1CUTmu->SetMinimum(0);
Eff1CUTmu->Draw("same"); 
       
Eff1CUTEe->GetXaxis()->SetTitle("Theta el[mrad]");
Eff1CUTEe->GetYaxis()->SetTitle("Efficency");
Eff1CUTEe->SetLineWidth(3);
Eff1CUTEe->SetLineColor(kGreen);
Eff1CUTEe->SetMaximum(1.2);
//Eff1CUTEe->SetMinimum(0);
Eff1CUTEe->Draw("same"); 

Eff1CUTtot->GetXaxis()->SetTitle("Theta el[mrad]");
Eff1CUTtot->GetYaxis()->SetTitle("Efficency");
Eff1CUTtot->SetLineWidth(3);
Eff1CUTtot->SetLineColor(kBlack);
Eff1CUTtot->SetMaximum(1.2);
//Eff1CUTtot->SetMinimum(0);
Eff1CUTtot->Draw("same"); 
gPad->BuildLegend(0.25,0.15,0.25,0.15);
ef->cd(2);
Eff2->GetXaxis()->SetTitle("Theta el[mrad]");
Eff2->GetYaxis()->SetTitle("Efficency");
Eff2->SetLineWidth(3);
Eff2->SetLineColor(kBlue);
Eff2->SetMaximum(1.2);
//Eff2->SetMinimum(0);
Eff2->Draw();
Eff2CUT->GetXaxis()->SetTitle("Theta el[mrad]");
Eff2CUT->GetYaxis()->SetTitle("Efficency");
Eff2CUT->SetLineWidth(3);
Eff2CUT->SetLineColor(kRed);
Eff2CUT->SetMaximum(1.2);
//Eff2CUT->SetMinimum(0);
Eff2CUT->Draw("same");   
       
Eff2CUTmu->GetXaxis()->SetTitle("Theta el[mrad]");
Eff2CUTmu->GetYaxis()->SetTitle("Efficency");
Eff2CUTmu->SetLineWidth(3);
Eff2CUTmu->SetLineColor(kOrange);
Eff2CUTmu->SetMaximum(1.2);
//Eff2CUTmu->SetMinimum(0);
Eff2CUTmu->Draw("same"); 
       
Eff2CUTEe->GetXaxis()->SetTitle("Theta el[mrad]");
Eff2CUTEe->GetYaxis()->SetTitle("Efficency");
Eff2CUTEe->SetLineWidth(3);
Eff2CUTEe->SetLineColor(kGreen);
Eff2CUTEe->SetMaximum(1.2);
//Eff2CUTEe->SetMinimum(0);
Eff2CUTEe->Draw("same"); 

Eff2CUTtot->GetXaxis()->SetTitle("Theta el[mrad]");
Eff2CUTtot->GetYaxis()->SetTitle("Efficency");
Eff2CUTtot->SetLineWidth(3);
Eff2CUTtot->SetLineColor(kBlack);
Eff2CUTtot->SetMaximum(1.2);
//Eff2CUTtot->SetMinimum(0);
Eff2CUTtot->Draw("same"); 
gPad->BuildLegend(0.25,0.15,0.25,0.15);

ef->SaveAs("/home/LHCB-T3/espedicato/tesi/eff/Eff.png");    


    
/*TCanvas * efCUT= new TCanvas("efcut","efcut",1000,100,2500,2000);
efCUT->Divide(1,2);
efCUT->cd(1);
Eff1CUT->GetXaxis()->SetTitle("Theta el[mrad]");
Eff1CUT->GetYaxis()->SetTitle("Efficency");
Eff1CUT->SetLineWidth(3);
Eff1CUT->SetLineColor(kBlack);
Eff1CUT->SetMaximum(1);
Eff1CUT->SetMinimum(0);
Eff1CUT->Draw();  
gPad->BuildLegend(0.25,0.15,0.25,0.15);
efCUT->cd(2);
Eff2CUT->GetXaxis()->SetTitle("Theta el[mrad]");
Eff2CUT->GetYaxis()->SetTitle("Efficency");
Eff2CUT->SetLineWidth(3);
Eff2CUT->SetLineColor(kRed);
Eff2CUT->SetMaximum(1);
Eff2CUT->SetMinimum(0);
Eff2CUT->Draw();   
gPad->BuildLegend(0.25,0.15,0.25,0.15);
efCUT->SaveAs("/home/LHCB-T3/espedicato/tesi/eff/EffCUT.png");*/
    
    
TCanvas * c5= new TCanvas("c5","c5",1000,1000,2500,2000);
c5->Divide(2,1);
/*c5->cd(1);
The_true->SetLineColor(kBlack);
The_true->SetLineWidth(3);
The_true->Draw("HIST"); 
The->SetLineWidth(3);
The->Draw("HIST same");
gPad->BuildLegend(0.3,0.21,0.3,0.21);*/
c5->cd(1);
The_true1->SetLineColor(kBlack);
The_true1->SetLineWidth(3);
The_true1->Draw("HIST");
The1->SetLineWidth(3);
The1->Draw("HIST same");  

gPad->BuildLegend(0.3,0.21,0.3,0.21);
c5->cd(2);
The_true2->SetLineColor(kBlack);
The_true2->SetLineWidth(3);
The_true2->Draw("HIST");   
The2->SetLineWidth(3);
The2->Draw("HIST same"); 

gPad->BuildLegend(0.3,0.21,0.3,0.21);
    
c5->SaveAs("/home/LHCB-T3/espedicato/tesi/eff/th_el.png"); 
    
TCanvas * c5MCS= new TCanvas("c5MCS","c5MCS",1000,100,2500,2000);
c5MCS->Divide(2,4);
c5MCS->cd(1);
The_trueCUT1->SetLineColor(kRed);
The_trueCUT1->SetLineWidth(3);
The_trueCUT1->Draw("HIST");
The1CUT->SetLineWidth(3);
The1CUT->Draw("HIST same");  

gPad->BuildLegend(0.3,0.21,0.3,0.21);
c5MCS->cd(2);
The_trueCUT2->SetLineColor(kRed);
The_trueCUT2->SetLineWidth(3);
The_trueCUT2->Draw("HIST");     

The2CUT->SetLineWidth(3);
The2CUT->Draw("HIST same"); 

gPad->BuildLegend(0.3,0.21,0.3,0.21);
    
    
c5MCS->cd(3);
The_trueCUT1mu->SetLineColor(kOrange);
The_trueCUT1mu->SetLineWidth(3);
The_trueCUT1mu->Draw("HIST");
The1CUTmu->SetLineWidth(3);
The1CUTmu->Draw("HIST same");  

gPad->BuildLegend(0.3,0.21,0.3,0.21);
c5MCS->cd(4);
The_trueCUT2mu->SetLineColor(kOrange);
The_trueCUT2mu->SetLineWidth(3);
The_trueCUT2mu->Draw("HIST");     

The2CUTmu->SetLineWidth(3);
The2CUTmu->Draw("HIST same"); 

gPad->BuildLegend(0.3,0.21,0.3,0.21);
    
    
c5MCS->cd(5);
The_trueCUT1Ee->SetLineColor(kGreen);
The_trueCUT1Ee->SetLineWidth(3);
The_trueCUT1Ee->Draw("HIST");
The1CUTEe->SetLineWidth(3);
The1CUTEe->Draw("HIST same");  

gPad->BuildLegend(0.3,0.21,0.3,0.21);
c5MCS->cd(6);
The_trueCUT2Ee->SetLineColor(kGreen);
The_trueCUT2Ee->SetLineWidth(3);
The_trueCUT2Ee->Draw("HIST");     

The2CUTEe->SetLineWidth(3);
The2CUTEe->Draw("HIST same"); 

gPad->BuildLegend(0.3,0.21,0.3,0.21);
    
c5MCS->cd(7);
The_trueCUT1tot->SetLineColor(kBlack);
The_trueCUT1tot->SetLineWidth(3);
The_trueCUT1tot->Draw("HIST");
The1CUTtot->SetLineWidth(3);
The1CUTtot->Draw("HIST same");  

gPad->BuildLegend(0.3,0.21,0.3,0.21);
c5MCS->cd(8);
The_trueCUT2tot->SetLineColor(kBlack);
The_trueCUT2tot->SetLineWidth(3);
The_trueCUT2tot->Draw("HIST");     

The2CUTtot->SetLineWidth(3);
The2CUTtot->Draw("HIST same"); 

gPad->BuildLegend(0.3,0.21,0.3,0.21);
        
    
c5MCS->SaveAs("/home/LHCB-T3/espedicato/tesi/eff/th_elCUT.png"); 
  
TCanvas * c= new TCanvas("c5MCS","c5MCS",1000,100,2500,2000);
rmu->Draw();
c->SaveAs("/home/LHCB-T3/espedicato/tesi/eff/rmu.png"); 
    

Int_t nx1 = E3x31->GetNbinsX();
Int_t ny1 = E3x31->GetNbinsY();
for (Int_t i=1; i<nx1+1; i++) {
for (Int_t j=1; j<ny1+1; j++) {
if (E3x31->GetBinContent(i,j)<1) E3x31->SetBinContent(i,j,0);}}
    
Int_t nx2 = E3x32->GetNbinsX();
Int_t ny2 = E3x32->GetNbinsY();
for (Int_t i=1; i<nx2+1; i++) {
for (Int_t j=1; j<ny2+1; j++) {
if (E3x32->GetBinContent(i,j)<1) E3x32->SetBinContent(i,j,0);}}


    
Int_t nx1CUT = E3x31CUT->GetNbinsX();
Int_t ny1CUT = E3x31CUT->GetNbinsY();
for (Int_t i=1; i<nx1CUT+1; i++) {
for (Int_t j=1; j<ny1CUT+1; j++) {
if (E3x31CUT->GetBinContent(i,j)<1) E3x31CUT->SetBinContent(i,j,0);}}
    
Int_t nx2CUT = E3x32CUT->GetNbinsX();
Int_t ny2CUT = E3x32CUT->GetNbinsY();
for (Int_t i=1; i<nx2CUT+1; i++) {
for (Int_t j=1; j<ny2CUT+1; j++) {
if (E3x32CUT->GetBinContent(i,j)<1) E3x32CUT->SetBinContent(i,j,0);}}
    
    
    
Int_t nx1CUTmu = E3x31CUTmu->GetNbinsX();
Int_t ny1CUTmu = E3x31CUTmu->GetNbinsY();
for (Int_t i=1; i<nx1CUTmu+1; i++) {
for (Int_t j=1; j<ny1CUTmu+1; j++) {
if (E3x31CUTmu->GetBinContent(i,j)<1) E3x31CUTmu->SetBinContent(i,j,0);}}
    
Int_t nx2CUTmu = E3x32CUTmu->GetNbinsX();
Int_t ny2CUTmu = E3x32CUTmu->GetNbinsY();
for (Int_t i=1; i<nx2CUTmu+1; i++) {
for (Int_t j=1; j<ny2CUTmu+1; j++) {
if (E3x32CUTmu->GetBinContent(i,j)<1) E3x32CUTmu->SetBinContent(i,j,0);}} 
  

Int_t nx1CUTEe = E3x31CUTEe->GetNbinsX();
Int_t ny1CUTEe = E3x31CUTEe->GetNbinsY();
for (Int_t i=1; i<nx1CUTEe+1; i++) {
for (Int_t j=1; j<ny1CUTEe+1; j++) {
if (E3x31CUTEe->GetBinContent(i,j)<1) E3x31CUTEe->SetBinContent(i,j,0);}}
    
Int_t nx2CUTEe = E3x32CUTEe->GetNbinsX();
Int_t ny2CUTEe = E3x32CUTEe->GetNbinsY();
for (Int_t i=1; i<nx2CUTEe+1; i++) {
for (Int_t j=1; j<ny2CUTEe+1; j++) {
if (E3x32CUTEe->GetBinContent(i,j)<1) E3x32CUTEe->SetBinContent(i,j,0);}}
    
    
Int_t nx1CUTtot = E3x31CUTtot->GetNbinsX();
Int_t ny1CUTtot = E3x31CUTtot->GetNbinsY();
for (Int_t i=1; i<nx1CUT+1; i++) {
for (Int_t j=1; j<ny1CUT+1; j++) {
if (E3x31CUTtot->GetBinContent(i,j)<1) E3x31CUTtot->SetBinContent(i,j,0);}}
    
Int_t nx2CUTtot = E3x32CUTtot->GetNbinsX();
Int_t ny2CUTtot = E3x32CUTtot->GetNbinsY();
for (Int_t i=1; i<nx2CUTtot+1; i++) {
for (Int_t j=1; j<ny2CUTtot+1; j++) {
if (E3x32CUTtot->GetBinContent(i,j)<1) E3x32CUTtot->SetBinContent(i,j,0);}}
    
TCanvas * c4a= new TCanvas("c4a","c4a",100,100,2500,2000);
c4a->Divide(2,5);
c4a->cd(1);
gStyle->SetPalette(kLake);
TColor::InvertPalette(); 
E3x31->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x31->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x31->Draw("COLZ");
c4a->cd(2);   
E3x32->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x32->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x32->Draw("COLZ");
c4a->cd(3);   
E3x31CUT->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x31CUT->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x31CUT->Draw("COLZ");
c4a->cd(4);   
E3x32CUT->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x32CUT->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x32CUT->Draw("COLZ");
c4a->cd(5);   
E3x31CUTmu->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x31CUTmu->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x31CUTmu->Draw("COLZ");
c4a->cd(6);   
E3x32CUTmu->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x32CUTmu->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x32CUTmu->Draw("COLZ");
c4a->cd(7);   
E3x31CUTEe->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x31CUTEe->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x31CUTEe->Draw("COLZ");
c4a->cd(8);   
E3x32CUTEe->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x32CUTEe->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x32CUTEe->Draw("COLZ");
c4a->cd(9);   
E3x31CUTtot->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x31CUTtot->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x31CUTtot->Draw("COLZ");
c4a->cd(10);   
E3x32CUTtot->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x32CUTtot->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x32CUTtot->Draw("COLZ");


c4a->SaveAs("/home/LHCB-T3/espedicato/tesi/eff/grthE_cut.png");    
    

}