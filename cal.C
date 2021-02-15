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
    
   // typedef map<int, double>  energy_cell; 
   // energy_cell en_c;    

Double_t n=0.; Double_t n_true=0.;
Double_t n_tot=0.; Double_t n_tot_true=0.; 
Double_t n_tot_e=0.; Double_t n_tot_e_true=0.;
Double_t n_tot_eBIG=0.; Double_t n_tot_eBIG_true=0.;
Double_t n_tot_e2P=0.; Double_t n_tot_e2P_true=0.;
 
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

TH1F* hist_E9_e=new TH1F("E9e", "E9 e- tot", 500,0,1);
TH1F* hist_thxz_e=new TH1F("thXZ", "th XZ e- CORE", 150,-0.1,0.1);
TH1F* hist_thyz_e=new TH1F("thYZ", "th YZ e- CORE", 150,-0.1,0.1);
TH1F* hist_thxz_e1=new TH1F("thXZ", "th XZ TAR 1 e- CORE", 150,-0.1,0.1);
TH1F* hist_thyz_e1=new TH1F("thYZ", "th YZ e- TAR 1 CORE", 150,-0.1,0.1);
TH1F* hist_thxz_e2=new TH1F("thXZ", "th XZ e- TAR 2 CORE", 150,-0.1,0.1);
TH1F* hist_thyz_e2=new TH1F("thYZ", "th YZ e- TAR 2 CORE", 150,-0.1,0.1);
    
TH1F* hist_E9_eBIG=new TH1F("E9e5X5", "E9 e- tot 5X5", 500,0,1);
TH1F* hist_thxz_eBIG=new TH1F("thXZ5X5", "th XZ e- 5X5", 150,-0.1,0.1);
TH1F* hist_thyz_eBIG=new TH1F("thYZ5X5", "th YZ e- 5X5", 150,-0.1,0.1);
TH1F* hist_thxz_e1BIG=new TH1F("thXZ5X5", "th XZ e- TAR 1 5X5", 150,-0.1,0.1);
TH1F* hist_thyz_e1BIG=new TH1F("thYZ5X5", "th YZ e- TAR 1 5X5", 150,-0.1,0.1);
TH1F* hist_thxz_e2BIG=new TH1F("thXZ5X5", "th XZ e- TAR 2 5X5", 150,-0.1,0.1);
TH1F* hist_thyz_e2BIG=new TH1F("thYZ5X5", "th YZ e- TAR 2 5X5", 150,-0.1,0.1);

TH1F* hist_E9_e2P=new TH1F("E9ecrown", "E9 e- tot crown", 500,0,1);
TH1F* hist_thxz_e2P=new TH1F("thXZcrown", "th XZ e- crown", 150,-0.1,0.1);
TH1F* hist_thyz_e2P=new TH1F("thYZcrown", "th YZ e- crown", 150,-0.1,0.1);
TH1F* hist_thxz_e12P=new TH1F("thXZcrown", "th XZ e- TAR 1 crown", 150,-0.1,0.1);
TH1F* hist_thyz_e12P=new TH1F("thYZcrown", "th YZ e- TAR 1 crown", 150,-0.1,0.1);
TH1F* hist_thxz_e22P=new TH1F("thXZcrown", "th XZ e- TAR 2 crown", 150,-0.1,0.1);
TH1F* hist_thyz_e22P=new TH1F("thYZcrown", "th YZ e- TAR 2 crown", 150,-0.1,0.1);

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
TH1F* The_true=new TH1F("th", "th El out true", 1000,0,100);
TH1F* The_trueCUT=new TH1F("th", "th El out true", 300,0,30);
TH1F* The_trueCUT1=new TH1F("th", "th El out true TAR 1", 300,0,30);
TH1F* The_trueCUT2=new TH1F("th", "th El out true TAR 2", 300,0,30);

    
    
TH1F* The=new TH1F("th", "th El out core", 1000,0,100); 
TH1F* TheBIG=new TH1F("th", "th El out 5X5", 1000,0,100); 
TH1F* The2P=new TH1F("th", "th El out crown", 1000,0,100); 

TH1F* The_true1=new TH1F("th", "th El out true TAR 1", 1000,0,100);    
TH1F* The1=new TH1F("th", "th El out core TAR 1", 1000,0,100); 
TH1F* TheBIG1=new TH1F("th", "th El out 5X5 TAR 1", 1000,0,100); 
TH1F* The2P1=new TH1F("th", "th El out crown TAR 1", 1000,0,100); 

    
TH1F* The_true2=new TH1F("th", "th El out true TAR 2", 1000,0,100);     
TH1F* The2=new TH1F("th", "th El out core TAR 2", 1000,0,100); 
TH1F* TheBIG2=new TH1F("th", "th El out 5X5 TAR 2", 1000,0,100); 
TH1F* The2P2=new TH1F("th", "th El out crown TAR 2", 1000,0,100); 
    
    
TH1F* TheMCS=new TH1F("th", "th El out core MCS", 300,0,30); 
TH1F* TheBIGMCS=new TH1F("th", "th El out 5X5 MCS", 300,0,30); 
TH1F* The2PMCS=new TH1F("th", "th El out crown MCS", 300,0,30); 
    
    
TH1F* The1MCS=new TH1F("th", "th El out TAR 1 core MCS", 300,0,30); 
TH1F* TheBIG1MCS=new TH1F("th", "th El out 5X5 TAR 1 MCS", 300,0,30); 
TH1F* The2P1MCS=new TH1F("th", "th El out crown TAR 1 MCS", 300,0,30); 
    
    
TH1F* The2MCS=new TH1F("th", "th El out core TAR 2 MCS", 300,0,30); 
TH1F* TheBIG2MCS=new TH1F("th", "th El out 5X5 TAR 2 MCS", 300,0,30);
TH1F* The2P2MCS=new TH1F("th", "th El out crown TAR 2 MCS", 300,0,30); 
    
    
TH2F  *E3x3  = new TH2F("ThEel" , " Theta el Vs. E_ECAL core",100,0,50,280,0,140);
TH2F  *E3x3true  = new TH2F("ThEel" , " Theta el Vs. E_ECAL true",100,0,50,280,0,140);

TH2F  *E3x3BIG  = new TH2F("ThEelbig" , " Theta el Vs. E_ECAL 5x5",100,0,50,280,0,140);
TH2F  *E3x32P  = new TH2F("ThEel2p" , " Theta el Vs. E_ECAL crown",100,0,50,280,0,140);
    

    
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
       

        /*en_c[1]=detKinBeamRot_Ecell1; en_c[2]=detKinBeamRot_Ecell2; en_c[3]=detKinBeamRot_Ecell3; en_c[4]=detKinBeamRot_Ecell4; en_c[5]=detKinBeamRot_Ecell5;
        en_c[6]=detKinBeamRot_Ecell6; en_c[7]=detKinBeamRot_Ecell7; en_c[8]=detKinBeamRot_Ecell8; en_c[9]=detKinBeamRot_Ecell9; en_c[10]=detKinBeamRot_Ecell10;
        en_c[11]=detKinBeamRot_Ecell11; en_c[12]=detKinBeamRot_Ecell12; en_c[13]=detKinBeamRot_Ecell13; en_c[14]=detKinBeamRot_Ecell14; en_c[15]=detKinBeamRot_Ecell15;
        en_c[16]=detKinBeamRot_Ecell16; en_c[17]=detKinBeamRot_Ecell17; en_c[18]=detKinBeamRot_Ecell18; en_c[19]=detKinBeamRot_Ecell19; en_c[20]=detKinBeamRot_Ecell20;
        en_c[21]=detKinBeamRot_Ecell21; en_c[22]=detKinBeamRot_Ecell22; en_c[23]=detKinBeamRot_Ecell23; en_c[24]=detKinBeamRot_Ecell24; en_c[25]=detKinBeamRot_Ecell25;*/
       
    if (photon_coox!=-1 && photon_cooy!=-1)
    {  
     E_CAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {E_CAL=detKinBeamRot_Ee;}

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

//if (detKinBeamRot_n_cell_e!=0 && abs(detKinBeamRot_cooXe)<4.275 && abs(detKinBeamRot_cooYe)<4.275)  

The_true->Fill(detKinBeamRot_def_angle_e,wgt_full);
    if (detKinBeamRot_tar==0){The_true1->Fill(detKinBeamRot_def_angle_e,wgt_full);
                              The_trueCUT1->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    if (detKinBeamRot_tar==1){The_true2->Fill(detKinBeamRot_def_angle_e,wgt_full);
                              The_trueCUT2->Fill(detKinBeamRot_def_angle_e,wgt_full);}
       
The_trueCUT->Fill(detKinBeamRot_def_angle_e,wgt_full);
       
       
if (detKinBeamRot_n_cell_e!=0 && E_CAL>0)
{
   
        
    n_tot_eBIG+=wgt_full;
    ++n_tot_eBIG_true;
    hist_E9_eBIG->Fill(E9,wgt_full);
    hist_thxz_eBIG->Fill(anglex_e,wgt_full);
    hist_thyz_eBIG->Fill(angley_e,wgt_full);
    TheBIG->Fill(detKinBeamRot_def_angle_e,wgt_full);
    TheBIGMCS->Fill(detKinBeamRot_def_angle_e,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {hist_thxz_e1BIG->Fill(anglex_e,wgt_full);
    hist_thyz_e1BIG->Fill(angley_e,wgt_full);
    TheBIG1->Fill(detKinBeamRot_def_angle_e,wgt_full);
    TheBIG1MCS->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    if (detKinBeamRot_tar==1)
    {hist_thxz_e2BIG->Fill(anglex_e,wgt_full);
    hist_thyz_e2BIG->Fill(angley_e,wgt_full);
    TheBIG2->Fill(detKinBeamRot_def_angle_e,wgt_full);
    TheBIG2MCS->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    
    //hist_Eout_9_e->Fill(Eout_9,wgt_full);
    
if (detKinBeamRot_E_clus3x3!=0 ) {E3x3BIG->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);
                                 E3x3true->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_Ee,wgt_full);}
    
if(abs(detKinBeamRot_cooXe)<4.275 && abs(detKinBeamRot_cooYe)<4.275) {
    n_tot_e+=wgt_full;
    ++n_tot_e_true;
    hist_E9_e->Fill(E9,wgt_full);
    hist_thxz_e->Fill(anglex_e,wgt_full);
    hist_thyz_e->Fill(angley_e,wgt_full);
    The->Fill(detKinBeamRot_def_angle_e,wgt_full);
    TheMCS->Fill(detKinBeamRot_def_angle_e,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {hist_thxz_e1->Fill(anglex_e,wgt_full);
    hist_thyz_e1->Fill(angley_e,wgt_full);
    The1->Fill(detKinBeamRot_def_angle_e,wgt_full);
    The1MCS->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    if (detKinBeamRot_tar==1)
    {hist_thxz_e2->Fill(anglex_e,wgt_full);
    hist_thyz_e2->Fill(angley_e,wgt_full);
    The2->Fill(detKinBeamRot_def_angle_e,wgt_full);
    The2MCS->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    
    //hist_Eout_9_e->Fill(Eout_9,wgt_full);
    
    if (detKinBeamRot_E_clus3x3!=0) {
        E3x3->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);
        //E3x3true->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_Ee,wgt_full);
    } //(detKinBeamRot_E_clus3x3!=0 && E9<0.95 && E9>0.5) 
} else { 
    n_tot_e2P+=wgt_full;
    ++n_tot_e2P_true;
    hist_E9_e2P->Fill(E9,wgt_full);
    hist_thxz_e2P->Fill(anglex_e,wgt_full);
    hist_thyz_e2P->Fill(angley_e,wgt_full);
    The2P->Fill(detKinBeamRot_def_angle_e,wgt_full);
    The2PMCS->Fill(detKinBeamRot_def_angle_e,wgt_full);
    
    if (detKinBeamRot_tar==0)
    {hist_thxz_e12P->Fill(anglex_e,wgt_full);
    hist_thyz_e12P->Fill(angley_e,wgt_full);
    The2P1->Fill(detKinBeamRot_def_angle_e,wgt_full);
    The2P1MCS->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    if (detKinBeamRot_tar==1)
    {hist_thxz_e22P->Fill(anglex_e,wgt_full);
    hist_thyz_e22P->Fill(angley_e,wgt_full);
    The2P2->Fill(detKinBeamRot_def_angle_e,wgt_full);
    The2P2MCS->Fill(detKinBeamRot_def_angle_e,wgt_full);}
    
    //hist_Eout_9_e->Fill(Eout_9,wgt_full);
    
if (detKinBeamRot_E_clus3x3!=0) {E3x32P->Fill(detKinBeamRot_def_angle_e,detKinBeamRot_E_clus3x3,wgt_full);}}}
    

    
 
}
 cout << " Numero elettroni totali IN 5X5 " << n_tot_eBIG << " su un totale di " << n_tot << " eventi " << endl;
 cout << " Numero elettroni totali IN CORE " << n_tot_e << " su un totale di " << n_tot << " eventi " << endl;
 cout << " Numero elettroni totali su CORONA " << n_tot_e2P << " su un totale di " << n_tot << " eventi " << endl;

 cout << " Numero elettroni totali IN 5X5 true " << n_tot_eBIG_true << " su un totale di " << n_tot_true << " eventi " << endl;
 cout << " Numero elettroni totali IN CORE true " << n_tot_e_true << " su un totale di " << n_tot_true << " eventi " << endl;
 cout << " Numero elettroni totali su CORONA true " << n_tot_e2P_true << " su un totale di " << n_tot_true << " eventi " << endl;


    
/*TCanvas * c1= new TCanvas("c1","c1",1000,100,2500,2000);

hist_E9_eBIG->GetXaxis()->SetTitle("Ecentral/E3x3");
hist_E9_eBIG->SetLineWidth(3);
hist_E9_eBIG->SetLineColor(kRed);
    
hist_E9_eBIG->Draw("HIST");     
    
hist_E9_e->GetXaxis()->SetTitle("Ecentral/E3x3");
hist_E9_e->SetLineWidth(3);
hist_E9_e->Draw("HIST same"); 
    
hist_E9_e2P->GetXaxis()->SetTitle("Ecentral/E3x3");
hist_E9_e2P->SetLineWidth(3);
hist_E9_e2P->SetLineColor(kViolet);
hist_E9_e2P->Draw("HIST same"); 
gPad->BuildLegend(0.25,0.15,0.25,0.15);

c1->SaveAs("/home/LHCB-T3/espedicato/tesi/cal/E9.png");*/
    
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

c2->SaveAs("/home/LHCB-T3/espedicato/tesi/cal/dist.png");
    
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

c2a->SaveAs("/home/LHCB-T3/espedicato/tesi/cal/Dtheta.png");*/


Int_t nx = E3x3->GetNbinsX();
Int_t ny = E3x3->GetNbinsY();
for (Int_t i=1; i<nx+1; i++) {
for (Int_t j=1; j<ny+1; j++) {
if (E3x3->GetBinContent(i,j)<1) E3x3->SetBinContent(i,j,0);}}

Int_t nxBIG = E3x3BIG->GetNbinsX();
Int_t nyBIG = E3x3BIG->GetNbinsY();
for (Int_t i=1; i<nxBIG+1; i++) {
for (Int_t j=1; j<nyBIG+1; j++) {
if (E3x3BIG->GetBinContent(i,j)<1) E3x3BIG->SetBinContent(i,j,0);}}
    
Int_t nx2P = E3x32P->GetNbinsX();
Int_t ny2P = E3x32P->GetNbinsY();
for (Int_t i=1; i<nx2P+1; i++) {
for (Int_t j=1; j<ny2P+1; j++) {
if (E3x32P->GetBinContent(i,j)<1) E3x32P->SetBinContent(i,j,0);}}
    
Int_t nxtrue = E3x3true->GetNbinsX();
Int_t nytrue = E3x3true->GetNbinsY();
for (Int_t i=1; i<nxtrue+1; i++) {
for (Int_t j=1; j<nytrue+1; j++) {
if (E3x3true->GetBinContent(i,j)<1) E3x3true->SetBinContent(i,j,0);}}
    
TProfile *px = E3x3BIG->ProfileX("px");
   // px->SetErrorOption("S");
TProfile *px1 = E3x3->ProfileX("px1");  
   // px1->SetErrorOption("S");
TProfile *px2 = E3x32P->ProfileX("px2");   
   // px2->SetErrorOption("S");
TProfile *pxtrue = E3x3true->ProfileX("pxtrue");   

    
TCanvas * c4= new TCanvas("c4","c4",100,100,2500,2000);

px->GetXaxis()->SetTitle("Theta_el[mrad]");
px->GetYaxis()->SetTitle("Ereco3x3[GeV]");
px->SetLineColor(kRed);
px->SetLineWidth(3);
px->Draw();
 
px1->GetXaxis()->SetTitle("Theta_el[mrad]");
px1->GetYaxis()->SetTitle("Ereco3x3[GeV]");
px1->SetLineWidth(3);
px1->Draw("same");
  
px2->GetXaxis()->SetTitle("Theta_el[mrad]");
px2->GetYaxis()->SetTitle("Ereco3x3[GeV]");
px2->SetLineColor(kViolet);
px2->SetLineWidth(3);
px2->Draw("same");
    
pxtrue->GetXaxis()->SetTitle("Theta_el[mrad]");
pxtrue->GetYaxis()->SetTitle("E[GeV]");
pxtrue->SetLineColor(kGreen);
pxtrue->SetLineWidth(3);
pxtrue->Draw("same");    
    
gPad->BuildLegend(0.25,0.15,0.25,0.15);
c4->SaveAs("/home/LHCB-T3/espedicato/tesi/cal/thEPROFILE.png");   
    
TCanvas * cp= new TCanvas("c4p","c4p",100,100,2500,2000);

pxtrue->GetXaxis()->SetTitle("Theta_el[mrad]");
pxtrue->GetYaxis()->SetTitle("E[GeV]");
pxtrue->SetLineColor(kGreen);
pxtrue->SetLineWidth(3);
pxtrue->Draw();
 
px->GetXaxis()->SetTitle("Theta_el[mrad]");
px->GetYaxis()->SetTitle("Ereco3x3[GeV]");
px->SetLineWidth(3);
px->Draw("same");
    
px1->GetXaxis()->SetTitle("Theta_el[mrad]");
px1->GetYaxis()->SetTitle("Ereco3x3[GeV]");
px1->SetLineWidth(3);
px1->Draw("same");
gPad->BuildLegend(0.25,0.15,0.25,0.15);
cp->SaveAs("/home/LHCB-T3/espedicato/tesi/cal/thEprofTRUE.png");        
    
TCanvas * c4a= new TCanvas("c4","c4",100,100,2500,2000);
c4a->Divide(3,1);
c4a->cd(1);
gStyle->SetPalette(kLake);
TColor::InvertPalette(); 
E3x3BIG->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x3BIG->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x3BIG->Draw("COLZ");
c4a->cd(2);   
E3x3->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x3->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x3->Draw("COLZ");
c4a->cd(3);   
E3x32P->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x32P->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x32P->Draw("COLZ");

c4a->SaveAs("/home/LHCB-T3/espedicato/tesi/cal/thE.png");
    
TCanvas * cc= new TCanvas("c4","c4",100,100,2500,2000); 
cc->Divide(1,2);
cc->cd(1);
E3x3->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x3->GetYaxis()->SetTitle("Ereco3x3[GeV]");
E3x3->Draw("COLZ");
cc->cd(2);
E3x3true->GetXaxis()->SetTitle("Theta_el[mrad]");
E3x3true->GetYaxis()->SetTitle("Etrue");
E3x3true->Draw("COLZ");
cc->SaveAs("/home/LHCB-T3/espedicato/tesi/cal/thEcore.png");

    
TCanvas * c5= new TCanvas("c5","c5",1000,100,2500,2000);
c5->Divide(1,3);
c5->cd(1);
The_true->SetLineColor(kBlack);
The_true->SetLineWidth(3);
The_true->Draw("HIST"); 
TheBIG->SetLineColor(kRed);
TheBIG->SetLineWidth(3);
TheBIG->Draw("HIST same"); 
The->SetLineWidth(3);
The->Draw("HIST same");
The2P->SetLineColor(kViolet);
The2P->SetLineWidth(3);
The2P->Draw("HIST same"); 
gPad->BuildLegend(0.3,0.21,0.3,0.21);
c5->cd(2);
The_true1->SetLineColor(kBlack);
The_true1->SetLineWidth(3);
The_true1->Draw("HIST");
TheBIG1->SetLineColor(kRed);
TheBIG1->SetLineWidth(3);
TheBIG1->Draw("HIST same");
The1->SetLineWidth(3);
The1->Draw("HIST same");  
The2P1->SetLineColor(kViolet);
The2P1->SetLineWidth(3);
The2P1->Draw("HIST same");
gPad->BuildLegend(0.3,0.21,0.3,0.21);
c5->cd(3);
The_true2->SetLineColor(kBlack);
The_true2->SetLineWidth(3);
The_true2->Draw("HIST");     
TheBIG2->SetLineColor(kRed);
TheBIG2->SetLineWidth(3);
TheBIG2->Draw("HIST same"); 
The2->SetLineWidth(3);
The2->Draw("HIST same"); 
The2P2->SetLineColor(kViolet);
The2P2->SetLineWidth(3);
The2P2->Draw("HIST same"); 
gPad->BuildLegend(0.3,0.21,0.3,0.21);
    
c5->SaveAs("/home/LHCB-T3/espedicato/tesi/cal/th_el.png"); 
    
TCanvas * c5MCS= new TCanvas("c5MCS","c5MCS",1000,100,2500,2000);
c5MCS->Divide(1,3);
c5MCS->cd(1);
The_trueCUT->SetLineColor(kBlack);
The_trueCUT->SetLineWidth(3);
The_trueCUT->Draw("HIST"); 
TheBIG->SetLineColor(kRed);
TheBIG->SetLineWidth(3);
TheBIG->Draw("HIST same"); 
The->SetLineWidth(3);
The->Draw("HIST same");
The2P->SetLineColor(kViolet);
The2P->SetLineWidth(3);
The2P->Draw("HIST same"); 
gPad->BuildLegend(0.3,0.21,0.3,0.21);
c5MCS->cd(2);
The_trueCUT1->SetLineColor(kBlack);
The_trueCUT1->SetLineWidth(3);
The_trueCUT1->Draw("HIST");
TheBIG1->SetLineColor(kRed);
TheBIG1->SetLineWidth(3);
TheBIG1->Draw("HIST same");
The1->SetLineWidth(3);
The1->Draw("HIST same");  
The2P1->SetLineColor(kViolet);
The2P1->SetLineWidth(3);
The2P1->Draw("HIST same");
gPad->BuildLegend(0.3,0.21,0.3,0.21);
c5MCS->cd(3);
The_trueCUT2->SetLineColor(kBlack);
The_trueCUT2->SetLineWidth(3);
The_trueCUT2->Draw("HIST");     
TheBIG2->SetLineColor(kRed);
TheBIG2->SetLineWidth(3);
TheBIG2->Draw("HIST same"); 
The2->SetLineWidth(3);
The2->Draw("HIST same"); 
The2P2->SetLineColor(kViolet);
The2P2->SetLineWidth(3);
The2P2->Draw("HIST same"); 
gPad->BuildLegend(0.3,0.21,0.3,0.21);
    
    
c5MCS->SaveAs("/home/LHCB-T3/espedicato/tesi/cal/th_elMCS.png"); 
  
  /*  TCanvas * theC= new TCanvas("tar","tar",1500,1000,3500,2000);
    theC->Divide(2,3);
    theC->cd(1);
    hist_thxz_e->SetLineColor(kBlue);
    hist_thxz_e->SetLineWidth(3);
    hist_thxz_e->Draw("HIST");
    hist_thxz_e1->SetLineColor(8);
    hist_thxz_e1->SetLineWidth(3);
    hist_thxz_e1->Draw("HIST SAME");
    hist_thxz_e2->SetLineColor(kBlack);
    hist_thxz_e2->SetLineWidth(3);
    hist_thxz_e2->Draw("HIST SAME");
    hist_thxz_e->GetXaxis()->SetTitle("Theta XZ [rad]");
    gPad->BuildLegend(0.3,0.21,0.3,0.21);
    theC->cd(2);
    hist_thyz_e->SetLineColor(kBlue);
    hist_thyz_e->SetLineWidth(3);
    hist_thyz_e->Draw("HIST");
    hist_thyz_e1->SetLineColor(8);
    hist_thyz_e1->SetLineWidth(3);
    hist_thyz_e1->Draw("HIST SAME");
    hist_thyz_e2->SetLineColor(kBlack);
    hist_thyz_e2->SetLineWidth(3);
    hist_thyz_e2->Draw("HIST SAME");
    hist_thyz_e->GetXaxis()->SetTitle("Theta YZ [rad]");
    gPad->BuildLegend(0.3,0.21,0.3,0.21);
    

    theC->cd(3);
    hist_thxz_eBIG->SetLineColor(kRed);
    hist_thxz_eBIG->SetLineWidth(3);
    hist_thxz_eBIG->Draw("HIST");
    hist_thxz_e1BIG->SetLineColor(8);
    hist_thxz_e1BIG->SetLineWidth(3);
    hist_thxz_e1BIG->Draw("HIST SAME");
    hist_thxz_e2BIG->SetLineColor(kBlack);
    hist_thxz_e2BIG->SetLineWidth(3);
    hist_thxz_e2BIG->Draw("HIST SAME");
    hist_thxz_eBIG->GetXaxis()->SetTitle("Theta XZ [rad]");
    gPad->BuildLegend(0.3,0.21,0.3,0.21);
    theC->cd(4);
    hist_thyz_eBIG->SetLineColor(kRed);
    hist_thyz_eBIG->SetLineWidth(3);
    hist_thyz_eBIG->Draw("HIST");
    hist_thyz_e1BIG->SetLineColor(8);
    hist_thyz_e1BIG->SetLineWidth(3);
    hist_thyz_e1BIG->Draw("HIST SAME");
    hist_thyz_e2BIG->SetLineColor(kBlack);
    hist_thyz_e2BIG->SetLineWidth(3);
    hist_thyz_e2BIG->Draw("HIST SAME");
    hist_thyz_eBIG->GetXaxis()->SetTitle("Theta YZ [rad]");
    gPad->BuildLegend(0.3,0.21,0.3,0.21);
    
    theC->cd(5);
    hist_thxz_e2P->SetLineColor(kViolet);
    hist_thxz_e2P->SetLineWidth(3);
    hist_thxz_e2P->Draw("HIST");
    hist_thxz_e12P->SetLineColor(8);
    hist_thxz_e12P->SetLineWidth(3);
    hist_thxz_e12P->Draw("HIST SAME");
    hist_thxz_e22P->SetLineColor(kBlack);
    hist_thxz_e22P->SetLineWidth(3);
    hist_thxz_e22P->Draw("HIST SAME");
    hist_thxz_e2P->GetXaxis()->SetTitle("Theta XZ [rad]");
    gPad->BuildLegend(0.3,0.21,0.3,0.21);
    theC->cd(6);
    hist_thyz_e2P->SetLineColor(kViolet);
    hist_thyz_e2P->SetLineWidth(3);
    hist_thyz_e2P->Draw("HIST");
    hist_thyz_e12P->SetLineColor(8);
    hist_thyz_e12P->SetLineWidth(3);
    hist_thyz_e12P->Draw("HIST SAME");
    hist_thyz_e22P->SetLineColor(kBlack);
    hist_thyz_e22P->SetLineWidth(3);
    hist_thyz_e22P->Draw("HIST SAME");
    hist_thyz_e2P->GetXaxis()->SetTitle("Theta YZ [rad]");
    gPad->BuildLegend(0.3,0.21,0.3,0.21);
    
    theC->SaveAs("/home/LHCB-T3/espedicato/tesi/cal/Th_el_XZYZ.png");     */
    
    TCanvas * the2= new TCanvas("tar2","tar2",1500,1000,3500,2000);
    the2->Divide(2,3);
    the2->cd(1);
    hist_thxz_eBIG->SetLineColor(kRed);
    hist_thxz_eBIG->SetLineWidth(3);
    hist_thxz_eBIG->Draw("HIST");
    hist_thxz_e->SetLineColor(kBlue);
    hist_thxz_e->SetLineWidth(3);
    hist_thxz_e->Draw("HIST same");
    hist_thxz_e2P->SetLineColor(kViolet);
    hist_thxz_e2P->SetLineWidth(3);
    hist_thxz_e2P->Draw("HIST same");
    gPad->BuildLegend(0.3,0.21,0.3,0.21);
    the2->cd(2);
    hist_thyz_eBIG->SetLineColor(kRed);
    hist_thyz_eBIG->SetLineWidth(3);
    hist_thyz_eBIG->Draw("HIST");
    hist_thyz_e->SetLineColor(kBlue);
    hist_thyz_e->SetLineWidth(3);
    hist_thyz_e->Draw("HIST same");
    hist_thyz_e2P->SetLineColor(kViolet);
    hist_thyz_e2P->SetLineWidth(3);
    hist_thyz_e2P->Draw("HIST same");
    gPad->BuildLegend(0.3,0.21,0.3,0.21);
    
    the2->cd(3);
    hist_thxz_e1BIG->SetLineColor(kRed);
    hist_thxz_e1BIG->SetLineWidth(3);
    hist_thxz_e1BIG->Draw("HIST");
    hist_thxz_e1->SetLineColor(kBlue);
    hist_thxz_e1->SetLineWidth(3);
    hist_thxz_e1->Draw("HIST same");
    hist_thxz_e12P->SetLineColor(kViolet);
    hist_thxz_e12P->SetLineWidth(3);
    hist_thxz_e12P->Draw("HIST same");
    gPad->BuildLegend(0.3,0.21,0.3,0.21);
    the2->cd(4);
    hist_thyz_e1BIG->SetLineColor(kRed);
    hist_thyz_e1BIG->SetLineWidth(3);
    hist_thyz_e1BIG->Draw("HIST");
    hist_thyz_e1->SetLineColor(kBlue);
    hist_thyz_e1->SetLineWidth(3);
    hist_thyz_e1->Draw("HIST same");
    hist_thyz_e12P->SetLineColor(kViolet);
    hist_thyz_e12P->SetLineWidth(3);
    hist_thyz_e12P->Draw("HIST same");
    gPad->BuildLegend(0.3,0.21,0.3,0.21);
    
    the2->cd(5);
    hist_thxz_e2BIG->SetLineColor(kRed);
    hist_thxz_e2BIG->SetLineWidth(3);
    hist_thxz_e2BIG->Draw("HIST");
    hist_thxz_e2->SetLineColor(kBlue);
    hist_thxz_e2->SetLineWidth(3);
    hist_thxz_e2->Draw("HIST same");
    hist_thxz_e22P->SetLineColor(kViolet);
    hist_thxz_e22P->SetLineWidth(3);
    hist_thxz_e22P->Draw("HIST same");
    gPad->BuildLegend(0.3,0.21,0.3,0.21);
    the2->cd(6);
    hist_thyz_e2BIG->SetLineColor(kRed);
    hist_thyz_e2BIG->SetLineWidth(3);
    hist_thyz_e2BIG->Draw("HIST");
    hist_thyz_e2->SetLineColor(kBlue);
    hist_thyz_e2->SetLineWidth(3);
    hist_thyz_e2->Draw("HIST same");
    hist_thyz_e22P->SetLineColor(kViolet);
    hist_thyz_e22P->SetLineWidth(3);
    hist_thyz_e22P->Draw("HIST same");
    gPad->BuildLegend(0.3,0.21,0.3,0.21);
        
    the2->SaveAs("/home/LHCB-T3/espedicato/tesi/cal/thxzyz.png");   
}