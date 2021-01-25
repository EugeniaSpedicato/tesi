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
    
    typedef map<int, double>  energy_cell; 
    energy_cell en_c;    
    
Double_t n_tot=0.;
Double_t n_tot_e=0.;
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
TH1F* hist_E9_eph_same=new TH1F("E9eph", "E9 e+ph same cell", 500,0,1);
TH1F* hist_E9_eph_diff=new TH1F("E9eph", "E9 e+ph different cell", 500,0,1);
TH1F* hist_E9_e=new TH1F("E9e", "E9 e tot", 500,0,1);
TH1F* hist_E9_NOph=new TH1F("E9noph", "E9 NO photons", 500,0,1);

    
    
TH1F* hist_Eout_9_eph_same=new TH1F("E9outeph", "Eout e+ph same cell", 300,0,0.3);
TH1F* hist_Eout_9_eph_diff=new TH1F("E9outeph", "Eout9 e+ph diff cell", 300,0,0.3);
TH1F* hist_Eout_9_e=new TH1F("E9oute", "Eout9 e",300,0,0.3);
TH1F* hist_Eout_9_NOph=new TH1F("E9outnoph", "Eout9 NO photons", 300,0,0.3);

TH1F* hist_dist=new TH1F("dist", "Dist e-gamma", 400,0,4);
TH1F* hist_dist_same=new TH1F("dist", "Dist e-gamma same cell", 400,0,4);
TH1F* hist_dist_diff=new TH1F("dist", "Dist e-gamma diff cel", 400,0,4);
    


    
if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
TGraph* E3x3 = new TGraph(nentries);
TGraph* E3x3noph = new TGraph(nentries);
    


    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
       
        en_c[1]=detKinBeamRot_Ecell1; en_c[2]=detKinBeamRot_Ecell2; en_c[3]=detKinBeamRot_Ecell3; en_c[4]=detKinBeamRot_Ecell4; en_c[5]=detKinBeamRot_Ecell5;
        en_c[6]=detKinBeamRot_Ecell6; en_c[7]=detKinBeamRot_Ecell7; en_c[8]=detKinBeamRot_Ecell8; en_c[9]=detKinBeamRot_Ecell9; en_c[10]=detKinBeamRot_Ecell10;
        en_c[11]=detKinBeamRot_Ecell11; en_c[12]=detKinBeamRot_Ecell12; en_c[13]=detKinBeamRot_Ecell13; en_c[14]=detKinBeamRot_Ecell14; en_c[15]=detKinBeamRot_Ecell15;
        en_c[16]=detKinBeamRot_Ecell16; en_c[17]=detKinBeamRot_Ecell17; en_c[18]=detKinBeamRot_Ecell18; en_c[19]=detKinBeamRot_Ecell19; en_c[20]=detKinBeamRot_Ecell20;
        en_c[21]=detKinBeamRot_Ecell21; en_c[22]=detKinBeamRot_Ecell22; en_c[23]=detKinBeamRot_Ecell23; en_c[24]=detKinBeamRot_Ecell24; en_c[25]=detKinBeamRot_Ecell25;
       
    if (photon_coox!=-1 && photon_cooy!=-1)
    {  
     E_CAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {  E_CAL=detKinBeamRot_Ee;}
       
        detKinBeamRot_cooXe=detKinBeamRot_cooXe*100; // cm
        detKinBeamRot_cooYe=detKinBeamRot_cooYe*100; // cm
        photon_coox=photon_coox*100; // cm
        photon_cooy=photon_cooy*100; // cm
       
       double Etotcal =0.;
       for(int i=1;i<26;++i){Etotcal+=en_c[i];}
       double Eout_9=(Etotcal-detKinBeamRot_E_clus3x3)/detKinBeamRot_E_clus3x3;
       E9=en_c[detKinBeamRot_n_max_Cell]/detKinBeamRot_E_clus3x3;
       /*cout << detKinBeamRot_n_max_Cell << " cella impatto elettrone " << detKinBeamRot_n_cell_e << "con energia " <<detKinBeamRot_Ee << " cella impatto fotone " << photon_n_cell_ph<< "con energia " <<photon_energy <<endl;*/
       

       
       
    Double_t d_e_ph=sqrt( (detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy) )/Rm; 
       
    /*Double_t d_e_mu=sqrt( (detKinBeamRot_cooXe-detKinBeamRot_cooXmu)*(detKinBeamRot_cooXe-detKinBeamRot_cooXmu)+(detKinBeamRot_cooYe-detKinBeamRot_cooYmu)*(detKinBeamRot_cooYe-detKinBeamRot_cooYmu) ); */
    n_tot+=wgt_full;
if (detKinBeamRot_n_cell_e!=0)  {     
    
    n_tot_e+=wgt_full;
    hist_E9_e->Fill(E9,wgt_full);
    hist_Eout_9_e->Fill(Eout_9,wgt_full);
    
    
  if (photon_n_cell_ph!=0)
  {   
      n_tot_eph+=wgt_full; // e+gamma sul calorimetro
//      hist_E9_eph->Fill(E9,wgt_full);
      hist_dist->Fill(d_e_ph,wgt_full);
    E3x3->SetPoint(j,detKinBeamRot_ThEl_interaction,detKinBeamRot_E_clus3x3); ++j;
        if (photon_n_cell_ph==detKinBeamRot_n_cell_e)
        {same_cell+=wgt_full;
         hist_E9_eph_same->Fill(E9,wgt_full);
         hist_dist_same->Fill(d_e_ph,wgt_full);
      hist_Eout_9_eph_same->Fill(Eout_9,wgt_full);
        }//stessa cella
        else {different_cell+=wgt_full;
              hist_E9_eph_diff->Fill(E9,wgt_full);
              hist_dist_diff->Fill(d_e_ph,wgt_full);
                hist_Eout_9_eph_diff->Fill(Eout_9,wgt_full);
             }  // cella diversa
  } else {
      n_tot_NOph+=wgt_full;   
      hist_E9_NOph->Fill(E9,wgt_full);
      hist_Eout_9_NOph->Fill(Eout_9,wgt_full);
      E3x3noph->SetPoint(i,detKinBeamRot_ThEl_interaction,detKinBeamRot_E_clus3x3); ++i;
        }

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

hist_E9_e->GetXaxis()->SetTitle("Ecentral/E3x3");
hist_E9_e->SetLineWidth(3);
hist_E9_e->Draw("HIST"); 

hist_E9_eph_same->SetLineColor(kRed);
hist_E9_eph_same->SetLineWidth(3);
hist_E9_eph_same->Draw("HIST same"); 
    
hist_E9_eph_diff->SetLineColor(kGreen);
hist_E9_eph_diff->SetLineWidth(3);
hist_E9_eph_diff->Draw("HIST same"); 
    
hist_E9_NOph->SetLineColor(kOrange);
hist_E9_NOph->SetLineWidth(3);
hist_E9_NOph->Draw("HIST same");
gPad->BuildLegend(0.25,0.15,0.25,0.15);

c1->SaveAs("/home/LHCB-T3/espedicato/tesi/E9.png");
    
TCanvas * c2= new TCanvas("c2","c2",1000,100,2500,2000);

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
    
TCanvas * c3= new TCanvas("c3","c3",1000,100,2500,2000);

hist_Eout_9_e->GetXaxis()->SetTitle("Ecentral/E3x3");
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

c3->SaveAs("/home/LHCB-T3/espedicato/tesi/out+3x3.png");

 
TCanvas * c4= new TCanvas("c4","c4",1000,100,2500,2000);
E3x3noph->SetMarkerColor(kBlack);
E3x3noph->Draw("AP");
E3x3->SetMarkerColor(kRed);
E3x3->Draw("AP same");
c4->SaveAs("/home/LHCB-T3/espedicato/tesi/thE.png");   
    
}