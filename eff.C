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

TH1F* The_trueCUT=new TH1F("th", "th El out true cut", 250,0,50);
TH1F* The_trueCUT1=new TH1F("th", "th El out true TAR 1 cut", 250,0,50);
TH1F* The_trueCUT2=new TH1F("th", "th El out true TAR 2 cut", 250,0,50);

TH1F* The_true=new TH1F("th", "th El out true", 500,0,50);    
TH1F* The=new TH1F("th", "th El out core", 250,0,50); 
/*TH1F* TheBIG=new TH1F("th", "th El out 5X5", 250,0,50); 
TH1F* The2P=new TH1F("th", "th El out crown", 250,0,50); */

TH1F* The_true1=new TH1F("th", "th El out true TAR 1", 250,0,50);    
TH1F* The1=new TH1F("th", "th El out core TAR 1", 250,0,50); 
/*TH1F* TheBIG1=new TH1F("th", "th El out 5X5 TAR 1", 250,0,50); 
TH1F* The2P1=new TH1F("th", "th El out crown TAR 1", 250,0,50); */

    
TH1F* The_true2=new TH1F("th", "th El out true TAR 2", 250,0,50);     
TH1F* The2=new TH1F("th", "th El out core TAR 2", 250,0,50); 
/*TH1F* TheBIG2=new TH1F("th", "th El out 5X5 TAR 2", 250,0,50); 
TH1F* The2P2=new TH1F("th", "th El out crown TAR 2", 250,0,50);*/ 
    
    
TH1F* TheMCS=new TH1F("th", "th El out core MCS", 250,0,50); 
/*TH1F* TheBIGMCS=new TH1F("th", "th El out 5X5 MCS", 250,0,50); 
TH1F* The2PMCS=new TH1F("th", "th El out crown MCS", 250,0,50);*/ 
    
    
TH1F* The1MCS=new TH1F("th", "th El out TAR 1 core MCS", 250,0,50); 
/*TH1F* TheBIG1MCS=new TH1F("th", "th El out 5X5 TAR 1 MCS", 250,0,50); 
TH1F* The2P1MCS=new TH1F("th", "th El out crown TAR 1 MCS", 250,0,50); */
    
    
TH1F* The2MCS=new TH1F("th", "th El out core TAR 2 MCS", 250,0,50); 
/*TH1F* TheBIG2MCS=new TH1F("th", "th El out 5X5 TAR 2 MCS", 250,0,50);
TH1F* The2P2MCS=new TH1F("th", "th El out crown TAR 2 MCS", 250,0,50); */

    

    
if (fChain == 0) return;

Long64_t nentries = fChain->GetEntriesFast();



    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
The_true->Fill(detKinBeamRot_def_angle_e,wgt_full);
The_trueCUT->Fill(detKinBeamRot_def_angle_e,wgt_full);

    if (photon_coox!=-1 && photon_cooy!=-1)
    {  
     E_CAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {E_CAL=detKinBeamRot_Ee;}

        detKinBeamRot_cooXe=detKinBeamRot_cooXe*100; // cm
        detKinBeamRot_cooYe=detKinBeamRot_cooYe*100; // cm
        photon_coox=photon_coox*100; // cm
        photon_cooy=photon_cooy*100; // cm
     
    if (detKinBeamRot_tar==0)
    {The_true1->Fill(detKinBeamRot_def_angle_e,wgt_full);
    if(detKinBeamRot_def_angle_mu>0.2) The_trueCUT1->Fill(detKinBeamRot_def_angle_e,wgt_full);
    }
    if (detKinBeamRot_tar==1)
    {The_true2->Fill(detKinBeamRot_def_angle_e,wgt_full);
    if(detKinBeamRot_def_angle_mu>0.2) The_trueCUT2->Fill(detKinBeamRot_def_angle_e,wgt_full);
    }
       
       
if (detKinBeamRot_n_cell_e!=0 && E_CAL>0)
{  

    /*TheBIG->Fill(detKinBeamRot_def_angle_e,wgt_full);
    TheBIGMCS->Fill(detKinBeamRot_def_angle_e,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {TheBIG1->Fill(detKinBeamRot_def_angle_e,wgt_full);
    TheBIG1MCS->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    if (detKinBeamRot_tar==1)
    {TheBIG2->Fill(detKinBeamRot_def_angle_e,wgt_full);
    TheBIG2MCS->Fill(detKinBeamRot_def_angle_e,wgt_full);}*/

    
if(abs(detKinBeamRot_cooXe)<4.275 && abs(detKinBeamRot_cooYe)<4.275) {

    The->Fill(detKinBeamRot_def_angle_e,wgt_full);
    if(detKinBeamRot_def_angle_mu>0.2) TheMCS->Fill(detKinBeamRot_def_angle_e,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {The1->Fill(detKinBeamRot_def_angle_e,wgt_full);
    if(detKinBeamRot_def_angle_mu>0.2) The1MCS->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    if (detKinBeamRot_tar==1)
    {The2->Fill(detKinBeamRot_def_angle_e,wgt_full);
    if(detKinBeamRot_def_angle_mu>0.2) The2MCS->Fill(detKinBeamRot_def_angle_e,wgt_full);}

} /*else { 

    The2P->Fill(detKinBeamRot_def_angle_e,wgt_full);
    The2PMCS->Fill(detKinBeamRot_def_angle_e,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {The2P1->Fill(detKinBeamRot_def_angle_e,wgt_full);
    The2P1MCS->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    if (detKinBeamRot_tar==1)
    {The2P2->Fill(detKinBeamRot_def_angle_e,wgt_full);
    The2P2MCS->Fill(detKinBeamRot_def_angle_e,wgt_full);}
}*/
     
}}
TH1F *Eff1 = new TH1F("ef1", "Eff Tar1", 250,0,50);
TH1F *Eff2 = new TH1F("ef2", "Eff Tar2", 250,0,50);
Eff1->Divide(The1,The_true1,1,1,"B");
Eff2->Divide(The2,The_true2,1,1,"B");

TCanvas * ef= new TCanvas("ef","ef",1000,100,2500,2000);
ef->Divide(1,2);
ef->cd(1);
Eff1->GetXaxis()->SetTitle("Theta el[mrad]");
Eff1->GetYaxis()->SetTitle("Efficency");
Eff1->SetLineWidth(3);
Eff1->SetLineColor(kBlack);
Eff1->SetMaximum(1);
Eff1->SetMinimum(0);
Eff1->Draw();  
gPad->BuildLegend(0.25,0.15,0.25,0.15);
ef->cd(2);
Eff2->GetXaxis()->SetTitle("Theta el[mrad]");
Eff2->GetYaxis()->SetTitle("Efficency");
Eff2->SetLineWidth(3);
Eff2->SetLineColor(kRed);
Eff2->SetMaximum(1);
Eff2->SetMinimum(0);
Eff2->Draw();   
gPad->BuildLegend(0.25,0.15,0.25,0.15);

ef->SaveAs("/home/LHCB-T3/espedicato/tesi/Eff.png");    
    
TH1F *Eff1CUT = new TH1F("ef1cut", "Eff Tar1 cut th_mu>0.2 mrad", 250,0,50);
TH1F *Eff2CUT = new TH1F("ef2cut", "Eff Tar2 cut th_mu>0.2 mrad", 250,0,50);
Eff1CUT->Divide(The1MCS,The_trueCUT1,1,1,"B");
Eff2CUT->Divide(The2MCS,The_trueCUT2,1,1,"B");

    
TCanvas * efCUT= new TCanvas("ef","ef",1000,100,2500,2000);
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
ef->cd(2);
Eff2CUT->GetXaxis()->SetTitle("Theta el[mrad]");
Eff2CUT->GetYaxis()->SetTitle("Efficency");
Eff2CUT->SetLineWidth(3);
Eff2CUT->SetLineColor(kRed);
Eff2CUT->SetMaximum(1);
Eff2CUT->SetMinimum(0);
Eff2CUT->Draw();   
gPad->BuildLegend(0.25,0.15,0.25,0.15);
efCUT->SaveAs("/home/LHCB-T3/espedicato/tesi/EffCUT.png");
 
TCanvas * t= new TCanvas("ef","ef",1000,100,2500,2000);
The_true->Draw();
The_true->SetLineColor(kRed);
The_trueCUT->Draw("same");

t->SaveAs("/home/LHCB-T3/espedicato/tesi/t.png");


}