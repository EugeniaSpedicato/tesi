#include <TH2.h>
#include <TH1.h>
#include <TGraph.h>
#include <cmath>
#include <TMath.h>
using namespace std;

#include <TStyle.h>
#include <TCanvas.h>

void ECAL()
{
    
Double_t Rm = 1.959 ; //raggio di Moliere in centimetri 
    
TF2 *fxy = new TF2("fxy","(1/(2*3.14159265358979323846*[0]*[1]))*(exp(-((x-[2])*(x-[2])+(y-[3])*(y-[3]))/(2*[0]*[1])))",-7.125,7.125,-7.125,7.125 );
 
TCanvas * ER= new TCanvas("ER","ER",1000,1000,2000,2500);

    fxy->SetParameter(0,Rm);
    fxy->SetParameter(1,Rm);
    fxy->SetParameter(2,0);
    fxy->SetParameter(3,0);
    fxy->Draw("surf1z"); 
    
    ER->SaveAs("functionE.png");
}