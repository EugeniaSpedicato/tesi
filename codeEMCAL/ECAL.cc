#include <TH2F.h>
#include <map>
#include <cmath>

#include <iostream>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"


#include "ECAL.h"


using namespace std;
using namespace RooFit;


ECAL::ECAL(double nbinsx, 
    double xlow, 
    double xup, 
    double nbinsy, 
    double ylow, 
    double yup)
    :
    nbinX(nbinsx),nbinY(nbinsy),Xlow(xlow),Xup(xup),Ylow(ylow),Yup(yup) 
    {
        //Queste mappe servono a mappare il numero di bin nel numero vero della cella e viceversa
        //perchè i numeri dei bin sono sballati a causa degli overflow e underflow bins
        number[36]=1; number[37]=2; number[38]=3; number[39]=4; number[40]=5;
        number[29]=6; number[30]=7; number[31]=8; number[32]=9; number[33]=10;
        number[22]=11; number[23]=12; number[24]=13; number[25]=14; number[26]=15;
        number[15]=16; number[16]=17; number[17]=18; number[18]=19; number[19]=20;
        number[8]=21; number[9]=22; number[10]=23; number[11]=24; number[12]=25;
        
        Rev_number[1]=36; Rev_number[2]=37; Rev_number[3]=38; Rev_number[4]=39; Rev_number[5]=40;
        Rev_number[6]=29; Rev_number[7]=30; Rev_number[8]=31; Rev_number[9]=32; Rev_number[10]=33;
        Rev_number[11]=22; Rev_number[12]=23; Rev_number[13]=24; Rev_number[14]=25; Rev_number[15]=26;
        Rev_number[16]=15; Rev_number[17]=16; Rev_number[18]=17; Rev_number[19]=18; Rev_number[20]=19;
        Rev_number[21]=8; Rev_number[22]=9; Rev_number[23]=10; Rev_number[24]=11; Rev_number[25]=12;
      

        
    EnRad_3 = new TProfile("Step3", "Radial Profile 2-3 X0", 20, 0, 4, 0, 5);
    EnRad_3->SetErrorOption("S");
    EnRad_6 = new TProfile("Step6", "Radial Profile 6-7 X0", 20, 0, 4, 0, 5);
    EnRad_6->SetErrorOption("S");
    EnRad_13 = new TProfile("Step13", "Radial Profile 19-20 X0", 20, 0, 4, 0, 5);
    EnRad_13->SetErrorOption("S");
    EnRad_20 = new TProfile("Step20", "Radial Profile 22-23 X0", 20, 0, 4, 0, 5);
    EnRad_20->SetErrorOption("S");
    EnRad_tot = new TProfile("Step20", "Radial Profile Total", 20, 0, 4, 0, 2);
    EnRad_tot->SetErrorOption("S");
        
    EnRad_3ERR = new TProfile("Step3", "Radial Profile 1-2 X0 with RMS", 20, 0, 4, 0, 5);
    EnRad_3ERR->SetErrorOption("S");
    EnRad_6ERR = new TProfile("Step6", "Radial Profile 5-6 X0 with RMS", 20, 0, 4, 0, 5);
    EnRad_6ERR->SetErrorOption("S");
    EnRad_13ERR = new TProfile("Step13", "Radial Profile 13-14 X0 with RMS", 20, 0, 4, 0, 5);
    EnRad_13ERR->SetErrorOption("S");
    EnRad_20ERR = new TProfile("Step20", "Radial Profile 22-23 X0 with RMS", 20, 0, 4, 0, 5);
    EnRad_20ERR->SetErrorOption("S");
        
    EnLong = new TProfile("Long", "Longitudinal Profile", 25, 0, 25);
    EnLong->SetErrorOption("S");
    EnLongERR = new TProfile("Long", "Longitudinal Profile with RMS", 25, 0, 25);
    EnLongERR->SetErrorOption("S");
    
    Er= new TProfile("<E(r)>/E", "Mean Energy Fraction in r <E(r)>/E ", 20, 0, 4, 0, 5);
    Er->SetErrorOption("S");
        
    Er2= new TProfile("<E(r)>", "Mean Energy Fraction in r <E(r)> ", 20, 0, 4);
    Er2->SetErrorOption("S");
    
    Energy_dist =new TH1F("Energy", "Energy",100,90,100);
    Energy_dist1 =new TH1F("Energy", "Energy 1 cell",150,90,105);
    Energy_dist3x3 =new TH1F("Energy", "Energy 3x3 cells",150,80,100);
    

    sigma =  new TProfile("Res", "Stochastic term",20, 0, 4, 0, 5);
    sigma->SetErrorOption("S");
    
    Array9=0;
        
    }



// metodo che crea l'istogramma rappresentante il calorimetro
TH2F* ECAL::CreateGrid(double nbinsx,double xlow,double xup,double nbinsy,double ylow,double yup)
{
    TH2F* EcalGrid = new TH2F("EcalGrid" , "EM Calorimeter with E in GeV",nbinsx,xlow,xup,nbinsy,ylow,yup);
    return EcalGrid;
};

void ECAL::SetEnergy(double energy)
{
    energy_IN=energy;
}
// metodo che assegna il numero della cella che viene colpita dalla particella 
double ECAL::GiveCentralCell(double coox,double cooy,TH2F* a)
{   
    int binx = a->GetXaxis()->FindBin(coox);
    int biny = a->GetYaxis()->FindBin(cooy);
    int nbin = a->GetBin(binx,biny);

    cout <<"Number of the cell:" << number[nbin] << endl;

    return number[nbin];
};

int* ECAL::GiveArray3x3(int n)
{
    if (n==1) {Array9= new int[9]{1,6,7,2,0,0,0,0,0};}
    if (n==2) {Array9= new int[9]{1,2,6,7,8,3,0,0,0};}
    if (n==3) {Array9= new int[9]{2,3,7,8,9,4,0,0,0};}
    if (n==4) {Array9= new int[9]{3,4,8,9,10,0,0,0};}
    if (n==5) {Array9= new int[9]{4,5,9,10,0,0,0};}
    if (n==6) {Array9= new int[9]{1,2,6,7,12,11,0,0,0};}
    if (n==7) {Array9= new int[9]{1,2,3,6,6,8,11,12,13};}
    if (n==8) {Array9= new int[9]{2,3,4,7,8,9,12,13,14};}
    if (n==9) {Array9= new int[9]{3,4,5,8,9,10,13,14,15};}
    if (n==10) {Array9= new int[9]{4,5,9,10,14,15,0,0,0};}
    if (n==11) {Array9= new int[9]{6,7,11,12,16,17,0,0,0};}
    if (n==12) {Array9= new int[9]{6,7,8,11,12,13,16,17,18};}
    if (n==13) {Array9= new int[9]{7,8,9,12,13,14,17,18,19};}
    if (n==14) {Array9= new int[9]{8,9,10,13,14,15,18,19,20};}
    if (n==15) {Array9= new int[9]{9,10,14,15,19,20,0,0,0};}
    if (n==16) {Array9= new int[9]{11,12,16,17,21,22,0,0,0};}
   if (n==17) {Array9= new int[9]{11,12,13,16,17,18,21,22,23};}
   if (n==18) {Array9= new int[9]{12,13,14,17,18,19,22,23,24};}
   if (n==19) {Array9= new int[9]{13,14,15,18,19,20,23,24,25};}
    if (n==20) {Array9= new int[9]{14,15,19,20,24,25,0,0,0};}
    if (n==21) {Array9= new int[9]{16,17,21,22,0,0,0,0,0};}
    if (n==22) {Array9= new int[9]{16,17,18,21,22,23,0,0,0};}
    if (n==23) {Array9= new int[9]{17,18,19,22,23,24,0,0,0};}
    if (n==24) {Array9= new int[9]{18,19,20,23,24,25,0,0,0};}
    if (n==25) {Array9= new int[9]{19,20,24,25,0,0,0,0,0};}
    
    return 0;
}

// metodo che aggiunge il punto di coo(x,y) all'istogramma, quindi al calorimetro e dà numero cella
double ECAL::AddHitCoo(double r, double phi,double xi, double yi, double w, TH2F* a)
{   r *= 2.19;
    double x=r*cos(phi)+xi; // coo x in cm
    double y=r*sin(phi)+yi; // coo y in cm
    a->Fill(x,y,w);   
 
double number=ECAL::GiveCentralCell(x,y,a);
return number;
};

void ECAL::AddHitCooDepth(double r, double phi,double xi, double yi, double w, double depth, double X0depth, TH2F* a)
{   depth += X0depth;
    r *= 2.19;
    double x=r*cos(phi)+xi; // coo x in cm
    double y=r*sin(phi)+yi; // coo y in cm
 if (24.7-X0depth>depth) 
 {a->Fill(x,y,w);   
double number=ECAL::GiveCentralCell(x,y,a); cout <<"è giusto"<< endl;}

};

// metodo che disegna l'evento nel calorimetro e le celle che vengono colpite
void ECAL::Draw_ECAL(TH2F* a){

TCanvas * Ecal_= new TCanvas("Ecal_","Ecal_",1500,100,3500,2000);
Ecal_->Divide(2,1);
Ecal_->cd(1);
gStyle->SetPalette(kAquamarine);
//TColor::InvertPalette();
a->SetXTitle("x (cm)");
a->SetYTitle("y (cm)");
a->Draw("COL");
a->Draw("TEXT SAME");
Ecal_->cd(2);
a->Draw("LEGO");
Ecal_->SaveAs("/home/LHCB-T3/espedicato/tesi/Ecal.png");


// riempi celle    
int binMax=a->GetMaximumBin();  
int CentralCell=number[binMax];
cout << "cella centrale rev " << Rev_number[CentralCell] <<" and vera " << CentralCell << endl;
Energy_dist1->Fill(a->GetBinContent(binMax));

double energy3x3=0.;    
ECAL::GiveArray3x3(CentralCell);
for (int i=0; i<9; ++i)
{
    if (Array9[i]>0 & Array9[i]<25 & Array9[i]!=0) energy3x3+=a->GetBinContent(Rev_number[Array9[i]]);
    cout << Rev_number[Array9[i]] << " and vera " << Array9[i]<< " c'è energia " << energy3x3 << endl;
}
Energy_dist3x3->Fill((energy3x3/energy_IN)*100);  
};


}

void ECAL::Print_()
{
   
TCanvas * encell= new TCanvas("Energy cells","Energy cells",1000,100,2500,2000);
Energy_dist3x3->Fit("gaus");  
Energy_dist3x3->SetLineWidth(2);
Energy_dist3x3->Draw("same");
encell->SaveAs("/home/LHCB-T3/espedicato/tesi/EnCell.png");

    
// Observable
RooRealVar energy3("energy3","energy3",80,100) ;
RooRealVar mean("mean","mean",94,96.5) ;
RooRealVar sigma("sigma","sigma",0.2,1.6) ;
RooRealVar alpha("alpha","alpha",1,0,20) ;
RooRealVar n("n","n",4,1,8) ;
RooCBShape CrystallBall("CrystallBall", "CrystallBall", energy3, mean, sigma, alpha, n);


RooDataHist en3("en3x3","en3x3",energy3,Import(*Energy_dist3x3));

RooPlot *frame = energy3.frame(Title("energy 3x3 cells"));
en3.plotOn(frame,MarkerStyle(kFullDotMedium));
CrystallBall.fitTo(en3);
CrystallBall.plotOn(frame);
CrystallBall.paramOn(frame,Layout(0.12,0.50));
frame->Draw();
    
TCanvas* cROO= new TCanvas("cROO","cROO",400,10,1100,800);
frame->GetXaxis()->SetTitle("Energy [GeV]");
frame->Draw();
cROO->SaveAs("/home/LHCB-T3/espedicato/tesi/frame.png");

    
}