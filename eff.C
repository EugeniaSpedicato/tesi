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

TH1F* The_true=new TH1F("th", "th El out true", 500,0,50);
TH1F* The_trueCUT=new TH1F("th", "th El out true", 300,0,30);
TH1F* The_trueCUT1=new TH1F("th", "th El out true TAR 1", 300,0,30);
TH1F* The_trueCUT2=new TH1F("th", "th El out true TAR 2", 300,0,30);

TH1F* The=new TH1F("th", "th El out core", 250,0,50); 
TH1F* TheBIG=new TH1F("th", "th El out 5X5", 250,0,50); 
TH1F* The2P=new TH1F("th", "th El out crown", 250,0,50); 

TH1F* The_true1=new TH1F("th", "th El out true TAR 1", 250,0,50);    
TH1F* The1=new TH1F("th", "th El out core TAR 1", 250,0,50); 
TH1F* TheBIG1=new TH1F("th", "th El out 5X5 TAR 1", 250,0,50); 
TH1F* The2P1=new TH1F("th", "th El out crown TAR 1", 250,0,50); 

    
TH1F* The_true2=new TH1F("th", "th El out true TAR 2", 250,0,50);     
TH1F* The2=new TH1F("th", "th El out core TAR 2", 250,0,50); 
TH1F* TheBIG2=new TH1F("th", "th El out 5X5 TAR 2", 250,0,50); 
TH1F* The2P2=new TH1F("th", "th El out crown TAR 2", 250,0,50); 
    
    
TH1F* TheMCS=new TH1F("th", "th El out core MCS", 300,0,30); 
TH1F* TheBIGMCS=new TH1F("th", "th El out 5X5 MCS", 300,0,30); 
TH1F* The2PMCS=new TH1F("th", "th El out crown MCS", 300,0,30); 
    
    
TH1F* The1MCS=new TH1F("th", "th El out TAR 1 core MCS", 300,0,30); 
TH1F* TheBIG1MCS=new TH1F("th", "th El out 5X5 TAR 1 MCS", 300,0,30); 
TH1F* The2P1MCS=new TH1F("th", "th El out crown TAR 1 MCS", 300,0,30); 
    
    
TH1F* The2MCS=new TH1F("th", "th El out core TAR 2 MCS", 300,0,30); 
TH1F* TheBIG2MCS=new TH1F("th", "th El out 5X5 TAR 2 MCS", 300,0,30);
TH1F* The2P2MCS=new TH1F("th", "th El out crown TAR 2 MCS", 300,0,30); 

    
//TH2F  *E3x3noph  = new TH2F("ThEel" , " Theta el Vs. E_ECAL no ph",100,0,70,70,0.2,140);

    
if (fChain == 0) return;

Long64_t nentries = fChain->GetEntriesFast();



    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
       n_tot+=wgt_full;
       ++n_tot_true;
       

    if (photon_coox!=-1 && photon_cooy!=-1)
    {  
     E_CAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {E_CAL=detKinBeamRot_Ee;}

        detKinBeamRot_cooXe=detKinBeamRot_cooXe*100; // cm
        detKinBeamRot_cooYe=detKinBeamRot_cooYe*100; // cm
        photon_coox=photon_coox*100; // cm
        photon_cooy=photon_cooy*100; // cm
       
       E9=detKinBeamRot_E_1/detKinBeamRot_E_clus3x3;
        Double_t anglex_e = atan2(detKinBeamRot_pXe_out, detKinBeamRot_pZe_out);
        Double_t angley_e = atan2(detKinBeamRot_pYe_out, detKinBeamRot_pZe_out);


The_true->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
    if (detKinBeamRot_tar==0){The_true1->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
                              The_trueCUT1->Fill(detKinBeamRot_ThEl_interaction,wgt_full);}
    if (detKinBeamRot_tar==1){The_true2->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
                              The_trueCUT2->Fill(detKinBeamRot_ThEl_interaction,wgt_full);}
       
The_trueCUT->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
       
       
if (detKinBeamRot_n_cell_e!=0 && E_CAL>0)
{
   
        
    n_tot_eBIG+=wgt_full;
    ++n_tot_eBIG_true;
    hist_E9_eBIG->Fill(E9,wgt_full);
    hist_thxz_eBIG->Fill(anglex_e,wgt_full);
    hist_thyz_eBIG->Fill(angley_e,wgt_full);
    TheBIG->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
    TheBIGMCS->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {hist_thxz_e1BIG->Fill(anglex_e,wgt_full);
    hist_thyz_e1BIG->Fill(angley_e,wgt_full);
    TheBIG1->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
    TheBIG1MCS->Fill(detKinBeamRot_ThEl_interaction,wgt_full);}
    if (detKinBeamRot_tar==1)
    {hist_thxz_e2BIG->Fill(anglex_e,wgt_full);
    hist_thyz_e2BIG->Fill(angley_e,wgt_full);
    TheBIG2->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
    TheBIG2MCS->Fill(detKinBeamRot_ThEl_interaction,wgt_full);}
    
    //hist_Eout_9_e->Fill(Eout_9,wgt_full);
    
if (detKinBeamRot_E_clus3x3!=0 ) {E3x3BIG->Fill(detKinBeamRot_ThEl_interaction,detKinBeamRot_E_clus3x3,wgt_full);
                                 E3x3true->Fill(detKinBeamRot_ThEl_interaction,detKinBeamRot_Ee,wgt_full);}
    
if(abs(detKinBeamRot_cooXe)<4.275 && abs(detKinBeamRot_cooYe)<4.275) {
    n_tot_e+=wgt_full;
    ++n_tot_e_true;
    hist_E9_e->Fill(E9,wgt_full);
    hist_thxz_e->Fill(anglex_e,wgt_full);
    hist_thyz_e->Fill(angley_e,wgt_full);
    The->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
    TheMCS->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {hist_thxz_e1->Fill(anglex_e,wgt_full);
    hist_thyz_e1->Fill(angley_e,wgt_full);
    The1->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
    The1MCS->Fill(detKinBeamRot_ThEl_interaction,wgt_full);}
    if (detKinBeamRot_tar==1)
    {hist_thxz_e2->Fill(anglex_e,wgt_full);
    hist_thyz_e2->Fill(angley_e,wgt_full);
    The2->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
    The2MCS->Fill(detKinBeamRot_ThEl_interaction,wgt_full);}
    
    //hist_Eout_9_e->Fill(Eout_9,wgt_full);
    
    if (detKinBeamRot_E_clus3x3!=0) {
        E3x3->Fill(detKinBeamRot_ThEl_interaction,detKinBeamRot_E_clus3x3,wgt_full);
        //E3x3true->Fill(detKinBeamRot_ThEl_interaction,detKinBeamRot_Ee,wgt_full);
    } //(detKinBeamRot_E_clus3x3!=0 && E9<0.95 && E9>0.5) 
} else { 
    n_tot_e2P+=wgt_full;
    ++n_tot_e2P_true;
    hist_E9_e2P->Fill(E9,wgt_full);
    hist_thxz_e2P->Fill(anglex_e,wgt_full);
    hist_thyz_e2P->Fill(angley_e,wgt_full);
    The2P->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
    The2PMCS->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {hist_thxz_e12P->Fill(anglex_e,wgt_full);
    hist_thyz_e12P->Fill(angley_e,wgt_full);
    The2P1->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
    The2P1MCS->Fill(detKinBeamRot_ThEl_interaction,wgt_full);}
    if (detKinBeamRot_tar==1)
    {hist_thxz_e22P->Fill(anglex_e,wgt_full);
    hist_thyz_e22P->Fill(angley_e,wgt_full);
    The2P2->Fill(detKinBeamRot_ThEl_interaction,wgt_full);
    The2P2MCS->Fill(detKinBeamRot_ThEl_interaction,wgt_full);}
    
    //hist_Eout_9_e->Fill(Eout_9,wgt_full);
    
if (detKinBeamRot_E_clus3x3!=0) {E3x32P->Fill(detKinBeamRot_ThEl_interaction,detKinBeamRot_E_clus3x3,wgt_full);}}}
     
}


TH1F *Eff1 = new TH1F("ef1", "Efficency Tar1", 250,0,50);
TH1F *Eff2 = new TH1F("ef2", "Efficency Tar2", 250,0,50);
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
Eff1->Draw("HIST");  
gPad->BuildLegend(0.25,0.15,0.25,0.15);
ef->cd(2);
Eff2->GetXaxis()->SetTitle("Theta el[mrad]");
Eff2->GetYaxis()->SetTitle("Efficency");
Eff2->SetLineWidth(3);
Eff2->SetLineColor(kRed);
Eff2->SetMaximum(1);
Eff2->SetMinimum(0);
Eff2->Draw("HIST");   
gPad->BuildLegend(0.25,0.15,0.25,0.15);

ef->SaveAs("/home/LHCB-T3/espedicato/tesi/Eff.png");
    

 
}