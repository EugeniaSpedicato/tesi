#include <cmath>
#include <iostream>
#include "TRandom.h"
#include "TMath.h"
#include "ElasticState.h"
#include "ResolutionModels.h"
#include "FastSim.h"
#include "TMatrixF.h" 
#include "TF2.h" 
#include "TMatrixD.h"
#include "TMatrixFBase.h"
#include <TMatrixFSym.h>
#include "TString.h"
#include "TVector3.h"
#include <TApplication.h>

using namespace MuE;
using namespace std;
using namespace ROOT::Math;

// from PDG Book 2018:
const Double_t FastSim::mm_PDG = 105.6583745 *0.001;
const Double_t FastSim::me_PDG = 0.5109989461 *0.001;

FastSim::FastSim(const MuE::MCpara & pargen, const MuE::FS_Input & fsi, bool _debug_):
  mm(pargen.mass_mu), me(pargen.mass_e), Ebeam(pargen.Ebeam), EbeamRMS(pargen.EbeamRMS),
  model(fsi.model), MSopt(fsi.MSopt), thickness(fsi.thickness), intrinsic_resolution(fsi.resolution), debug(_debug_), Minv(0)
  
{
  if (std::abs(mm - mm_PDG)/mm_PDG > 1e-7) {
    cout<<"\n"<< "***WARNING: muon mass = "<<mm<<" is different from the PDG mass: "<<mm_PDG<<endl;
  }
  if (std::abs(me - me_PDG)/me_PDG > 1e-7) {
    cout<<"\n"<< "***WARNING: electron mass = "<<me<<" is different from the PDG mass: "<<me_PDG<<endl;
  }
  
  if (model == 0) {
    cout<<"\n Simple detector model, DetModel = "<<model<<", option = "<<MSopt<<endl;
    cout<<"material thickness = "<<thickness<<" X0, intrinsic angular resolution = "<<intrinsic_resolution<<" mrad "<<endl;
  }
  else if (model == 1) {
    cout<<"\n Antonios detector model, 3 parameters"<< endl;
  }
  else cout<<"\n"<<"*** ERROR : FastSim, undefined detector model with DetModel = "<<model<<endl; 
}

// process an event
//
void FastSim::Process(const MuE::Event & event) {
  
  if (debug) cout<<"\n Process:  Run = "<<event.RunNr << " , Event = "<< event.EventNr << endl;
  
  Double_t emuin = event.E_mu_in;
  Double_t pmuin = sqrt(emuin*emuin-mm*mm);
  Double_t s = mm*mm + me*me + 2*me*emuin;
  
  PxPyPzEVector p_mu_in(0.,0.,pmuin,emuin);
  PxPyPzEVector p_e_in(0.,0.,0.,me);

  p_system = p_mu_in + p_e_in;
  Minv = sqrt(s);

  Double_t pcm = P_2bodies_CoM(Minv, mm, me);
  
  if (debug) cout<<"\n"<<"Incoming muon energy = "<< emuin 
		       << " GeV, s = "<<s<<" GeV^2, sqrt(s) = "<<Minv<<" GeV, pcm = "<<pcm<<" GeV"<<endl;
  
  PxPyPzEVector p_mu_out(event.P_mu_out.px, event.P_mu_out.py, event.P_mu_out.pz, event.P_mu_out.E);
  PxPyPzEVector p_e_out(event.P_e_out.px, event.P_e_out.py, event.P_e_out.pz, event.P_e_out.E);
   
sS    = 6*0.00064; //m spessore silicio !!!EHIIII CAMBIAMIIIIII
x0S = 0.094; // m
sB    = 0.015; //m spessore berillio
x0B = 0.353; // m
tar=gRandom->Integer(2);  // target where the interaction happens
vertex=gRandom->Uniform(); // where inc. muon interacts in the Beryllium tar.
    
sigSI=(13.6/(p_mu_in.E()*1000))*sqrt(sS/x0S)*(1+0.038*log(sS/x0S)); //rad e energy in MeV
sigBE=(13.6/(p_mu_in.E()*1000))*sqrt(sB/x0B)*(1+0.038*log(sB/x0B)); //rad  
sigBE2in=(13.6/(p_mu_in.E()*1000))*sqrt((sB*vertex)/x0B)*(1+0.038*log((sB*vertex)/x0B)); //rad      

// Appply beam divergence on incoming muon and MCS of station0 + vertex*tar0 (tar 0) 
// or station0 + tar0 + station1 + vertex*tar1
TMatrixD muin=RotDivIN(p_mu_in);
PxPyPzEVector p_mu_in_div(muin[0][0],muin[1][0],muin[2][0],muin[3][0]); // divergent P incoming mu   
PxPyPzEVector p_e_in_div=p_e_in; //it's the same as =0
    
//Applay rotation on P outgoing mu and e- in the reference system of the div. beam with also initial MCS
PxPyPzEVector p_mu_out_div=RotDiv(p_mu_in_div,p_mu_out);
PxPyPzEVector p_e_out_div=RotDiv(p_mu_in_div,p_e_out);
    
//MCS effect on outgoing particles, gives matrix with coo, angles and momenta for e- and mu
TMatrixD b=MCSout(p_mu_in_div,p_mu_out_div,p_e_out_div);
PxPyPzEVector p_mu_out_div_smeared(b[8][3],b[8][4],b[8][5],b[8][6]);
PxPyPzEVector p_e_out_div_smeared(b[17][3],b[17][4],b[17][5],b[17][6]);


TMatrixD coo(5,1);
coo[0][0]=b[8][0];
coo[1][0]=b[8][1];
coo[2][0]=b[17][0];
coo[3][0]=b[17][1];
coo[4][0]=b[8][2];
    
TMatrixD cooIN(5,1);
Double_t xin=b[18][0];
Double_t yin=b[18][1];
Double_t TheINT=b[18][2]; // angolo che ha nel momento dell'interazione, senza effetto di MCS 
    

  /*  if (MSopt ==0) {
    p_mu_out_div_smeared = Smear(p_mu_out_div);
    p_e_out_div_smeared = Smear(p_e_out_div);
  }
  else if (MSopt ==1) {
    p_mu_out_div_smeared = SmearX(p_mu_out_div);
    p_e_out_div_smeared = SmearX(p_e_out_div);
  }
  else if (MSopt ==2) {
    p_mu_out_div_smeared = SmearPolar(p_mu_out_div);
    p_e_out_div_smeared = SmearPolar(p_e_out_div);
  }
  else cout<<"\n"<<"*** ERROR : FastSim, undefined detector MS option = "<<MSopt<<endl;*/
  
  LoadKineVars(p_mu_in_div, p_e_in_div, p_mu_out_div_smeared, p_e_out_div_smeared, coo, TheINT, detKinBeamRot);
  
  LoadPhoton(event, photon, p_mu_in_div,muin[4][0],xin,yin);
    
}


// momentum for 2-body kinematics in the centre-of-mass system 
//
Double_t FastSim::P_2bodies_CoM (Double_t M, Double_t mm, Double_t me) const
{
  Double_t msum = std::abs(mm+me);
  Double_t mdif = std::abs(mm-me); 
  return 0.5/M*sqrt((M+mdif)*(M-mdif)*(M+msum)*(M-msum));
}

// Lorentz transformation to the centre-of-mass system
//
PxPyPzEVector FastSim::Lorentz_ToCoM(const PxPyPzEVector & pLab) const
{ 
  if (p_system.E() == Minv) return pLab;
  
  else {
    Double_t ecm = (pLab.E()*p_system.E()
		    -pLab.Px()*p_system.Px()-pLab.Py()*p_system.Py()-pLab.Pz()*p_system.Pz()) /Minv;
    
    Double_t fn = (ecm+pLab.E()) / (p_system.E()+Minv);
    
    return PxPyPzEVector(pLab.Px() - fn * p_system.Px(),
			  pLab.Py() - fn * p_system.Py(),
			  pLab.Pz() - fn * p_system.Pz(),
			  ecm);
  }
}

// Lorentz transformation to the Laboratory system
//
PxPyPzEVector FastSim::Lorentz_ToLab(const PxPyPzEVector & pCoM) const
{
  if (p_system.E() == Minv) return pCoM;
  
  else {
    Double_t elab = (pCoM.Px()*p_system.Px()+pCoM.Py()*p_system.Py()
		     +pCoM.Pz()*p_system.Pz()+pCoM.E()*p_system.E()) /Minv;
    
    Double_t fn   = (elab+pCoM.E()) / (p_system.E()+Minv);
    
    return PxPyPzEVector(pCoM.Px() + fn * p_system.Px(),
			  pCoM.Py() + fn * p_system.Py(),
			  pCoM.Pz() + fn * p_system.Pz(),
			  elab);
  }
}

TMatrixD FastSim::RotDivIN(const PxPyPzEVector & k) const
 {   

Double_t divthx = gRandom->Gaus(0., 0.00027);
Double_t divthy = gRandom->Gaus(0., 0.00020); 
    
     
if (tar==0){
Double_t Thx=gRandom->Gaus(divthx,sigSI);
Double_t Thy=gRandom->Gaus(divthy,sigSI);
Double_t thetaX1= gRandom->Gaus(Thx,sigBE2in);
Double_t thetaY1= gRandom->Gaus(Thy,sigBE2in); 


//Double_t anglex = atan2(k.Px(), k.Pz());
//Double_t angley = atan2(k.Py(), k.Pz()); 
     
Double_t anglex = thetaX1;
Double_t angley = thetaY1;
     
    //NB questi Px Py Pz sono del nuovo!!
Double_t pmuin=sqrt(k.Px()*k.Px()+k.Py()*k.Py()+k.Pz()*k.Pz());
Double_t pz=pmuin/(1+tan(anglex)*tan(anglex)+tan(angley)*tan(angley));
Double_t py=pz*tan(angley);
Double_t px=pz*tan(anglex);


TMatrixD pnewdiv(4,1);
pnewdiv[0][0]=px;
pnewdiv[1][0]=py;
pnewdiv[2][0]=pz; 
pnewdiv[3][0]=k.E();
return pnewdiv;  }
     
if (tar==1){
Double_t Thx=gRandom->Gaus(divthx,sigSI);
Double_t Thx1=gRandom->Gaus(Thx,sigBE);
Double_t Thx2=gRandom->Gaus(Thx1,sigSI);
    
Double_t Thy=gRandom->Gaus(divthy,sigSI);
Double_t Thy1=gRandom->Gaus(Thy,sigBE);
Double_t Thy2=gRandom->Gaus(Thy1,sigSI); 
    
Double_t thetaX1= gRandom->Gaus(Thx2,sigBE2in);
Double_t thetaY1= gRandom->Gaus(Thy2,sigBE2in); 

//Double_t anglex = atan2(k.Px(), k.Pz());
//Double_t angley = atan2(k.Py(), k.Pz()); 
     
Double_t anglex = thetaX1;
Double_t angley = thetaY1;
     
    //NB questi Px Py Pz sono del nuovo!!
Double_t pmuin=sqrt(k.Px()*k.Px()+k.Py()*k.Py()+k.Pz()*k.Pz());
Double_t pz=sqrt((pmuin*pmuin)/(1+tan(anglex)*tan(anglex)+tan(angley)*tan(angley)));
Double_t py=pz*tan(angley);
Double_t px=pz*tan(anglex);


TMatrixD pnewdiv(5,1);
pnewdiv[0][0]=px;
pnewdiv[1][0]=py;
pnewdiv[2][0]=pz; 
pnewdiv[3][0]=k.E();
pnewdiv[4][0]=tar;
    
return pnewdiv;  }
     
else {TMatrixD p(2,1);
p[0][0]=0;
p[1][0]=0;
return p;}
 }

 PxPyPzEVector FastSim::RotDiv(const PxPyPzEVector & k,const PxPyPzEVector & out) const
 {
Double_t pz=k.Pz();
Double_t py=k.Py();
Double_t px=k.Px();

Double_t ptz=sqrt(px*px+pz*pz);

/* PxPyPzEVector pnewdiv(px, py, pz, k.E());     
return pnewdiv;  */

// costruzione della matrice di rotazione
Double_t psi=atan2(px,pz); 
Double_t phi=atan2(py,ptz);  
     
TMatrixD R(3,3);
R[0][0]=cos(psi);
R[0][1]=0;
R[0][2]=sin(psi);
    R[1][0]=-sin(phi)*sin(psi);
    R[1][1]=cos(phi);
    R[1][2]=sin(phi)*cos(psi);
        R[2][0]=-cos(phi)*sin(psi);
        R[2][1]=-sin(phi);
        R[2][2]=cos(phi)*cos(psi);
     

TMatrixD pO(3,1);
    pO[0][0]=out.Px();
    pO[1][0]=out.Py();
    pO[2][0]=out.Pz();


TMatrixD pN(R, TMatrixD::kMult,pO);

PxPyPzEVector pnewdiv(pN[0][0], pN[1][0], pN[2][0], out.E());     
return pnewdiv;    
 }   


TMatrixD FastSim::MCSout(const PxPyPzEVector & kin, const PxPyPzEVector & k, const PxPyPzEVector & ke) const
{
        
    
Double_t const sSin    = 6*0.00064; //m spessore silicio per fascio entrante

Double_t const sS    = 0.00064; //m spessore silicio
Double_t const x0S = 0.094; // m

Double_t const sB    = 0.015; //m spessore berillio
Double_t const x0B = 0.353; // m



Double_t const dSS = 0.01; // m distanza tra i due 2S


Double_t const d = 0.25-0.005;

// m la distanza tra coppie di silici è 0.25. Però 0.25 è tra i due della coppia, quindi considero -dSS/2=0.005
Double_t const dCAL = 0.10; // m distanza silicio calorimetro
    
    
    
Double_t sigSImu=(13.6/(k.E()*1000))*sqrt(sS/x0S)*(1+0.038*log(sS/x0S)); //rad
// considero sB/2 per quando interagisce a metà 
Double_t sigSIe=(13.6/(ke.E()*1000))*sqrt(sS/x0S)*(1+0.038*log(sS/x0S)); //rad
// considero sB/2 per quando interagisce a metà 
Double_t sigBE2mu=(13.6/(k.E()*1000))*sqrt((sB*(1-vertex))/x0B)*(1+0.038*log((sB*(1-vertex))/x0B)); //rad
// considero sB per quando interagisce a nel target in vertex 
Double_t sigBE2e=(13.6/(ke.E()*1000))*sqrt((sB*(1-vertex))/x0B)*(1+0.038*log((sB*(1-vertex))/x0B)); //rad
// considero sB per quando interagisce a nel target in vertex 
Double_t sigBEmu=(13.6/(k.E()*1000))*sqrt(sB/x0B)*(1+0.038*log(sB/x0B)); //rad   
Double_t sigBEe=(13.6/(ke.E()*1000))*sqrt(sB/x0B)*(1+0.038*log(sB/x0B)); //rad   
    

    
    TMatrixD thetaX(2,7);
    TMatrixD thetaY(2,7);
    //coordinate prima dell'entrata  
    TMatrixD xi(2,7);
    TMatrixD yi(2,7);
    //coordinate dopo uscite con effetto MCS
    TMatrixD x(2,7);
    TMatrixD y(2,7);
    
    TMatrixD thetaXe(2,7);
    TMatrixD thetaYe(2,7);
    //coordinate prima dell'entrata entrata 
    TMatrixD xe(2,7);
    TMatrixD ye(2,7);
    //coordinate dopo uscite con effetto MCS
    TMatrixD xei(2,7);
    TMatrixD yei(2,7);
    
    Double_t anglexin = atan2(kin.Px(), kin.Pz());
    Double_t angleyin = atan2(kin.Py(), kin.Pz()); 
  
    // sono gli angoli dai momenti già ruotati
    Double_t anglex = atan2(k.Px(), k.Pz());
    Double_t angley = atan2(k.Py(), k.Pz()); 
    Double_t anglexe = atan2(ke.Px(), ke.Pz());
    Double_t angleye = atan2(ke.Py(), ke.Pz()); 
    Double_t The = ke.Theta(); // angolo in rad che ha nel momento dell'interazione, senza effetto di MCS 
    
    
    /*Double_t thetaEL = ke.Theta()*1000;//mrad
    && thetaEL<35 && thetaEL<70*/
    
    TMatrixD coo_in(1,3);
    coo_in[0][0]=100;
    coo_in[0][1]=100;
    coo_in[0][2]=100;
    
    Double_t xR = gRandom->Gaus(0,0.026);
    Double_t yR = gRandom->Gaus(0,0.027);

    
    if(tar==0)
    {
       // siamo nelle stazioni con il target di berillio. Ora entra nel berillio, la distanza che percorre nel vuoto il fascio è approssimativamente 1m-spessore_silici=1-sSin che serve per il calcolo della coordinata. Considero interazione in "vertex" nel berillio.
        
                //angoli di uscita dal berillio per e- e mu
                 thetaX[0][0]= gRandom->Gaus(anglex,sigBE2mu);
                 thetaY[0][0]= gRandom->Gaus(angley,sigBE2mu); 
                 thetaXe[0][0]= gRandom->Gaus(anglexe,sigBE2e);
                 thetaYe[0][0]= gRandom->Gaus(angleye,sigBE2e); 
        
                //coo dell'incoming muon dopo l'effetto del materiale (6 silici) della stazione 0
                Double_t x1 = gRandom->Gaus(xR,(1/sqrt(3))*sSin*sigSI);
                Double_t y1 = gRandom->Gaus(yR,(1/sqrt(3))*sSin*sigSI);
        
                //I muoni arrivano sul berillio, uguali ovviamente per elettrone e mu perchè ancora non nati..
                xi[0][0] = x1+(1-sSin)*tan(anglexin);
                yi[0][0] = y1+(1-sSin)*tan(angleyin);
                xei[0][0] = x1+(1-sSin)*tan(anglexin);
                yei[0][0] = y1+(1-sSin)*tan(angleyin);
        
                //coo d'interazione dopo metà berillio
                Double_t x2 = gRandom->Gaus(xi[0][0],(1/sqrt(3))*(sB*vertex)*sigBE2in);
                Double_t y2 = gRandom->Gaus(yi[0][0],(1/sqrt(3))*(sB*vertex)*sigBE2in);

                // dopo l'interazione, coordinate d'uscita dal berillio
                 x[0][0]=gRandom->Gaus(x2,(1/sqrt(3))*(sB*(1-vertex))*sigBE2mu);
                 y[0][0]=gRandom->Gaus(y2,(1/sqrt(3))*(sB*(1-vertex))*sigBE2mu);  
    
                 xe[0][0]=gRandom->Gaus(x2,(1/sqrt(3))*(sB*(1-vertex))*sigBE2e);
                 ye[0][0]=gRandom->Gaus(y2,(1/sqrt(3))*(sB*(1-vertex))*sigBE2e); 
            


        
              
//entro nelle 3 coppie di silici
for (Int_t p=1; p<7; p++)  {
                 if(p==1 || p==3 || p==5){
                 thetaX[0][p]= gRandom->Gaus(thetaX[0][p-1],sigSImu);
                 thetaY[0][p]= gRandom->Gaus(thetaY[0][p-1],sigSImu); 
                 thetaXe[0][p]= gRandom->Gaus(thetaXe[0][p-1],sigSIe);
                 thetaYe[0][p]= gRandom->Gaus(thetaYe[0][p-1],sigSIe); 

                xi[0][p]=x[0][p-1]+d*tan(thetaX[0][p-1]);
                x[0][p]=gRandom->Gaus(xi[0][p],(1/sqrt(3))*sS*sigSImu);
                yi[0][p]=y[0][p-1]+d*tan(thetaY[0][p-1]);
                y[0][p]=gRandom->Gaus(yi[0][p],(1/sqrt(3))*sS*sigSImu);
                xei[0][p]=xe[0][p-1]+d*tan(thetaXe[0][p-1]);
                xe[0][p]=gRandom->Gaus(xei[0][p],(1/sqrt(3))*sS*sigSIe);
                yei[0][p]=ye[0][p-1]+d*tan(thetaYe[0][p-1]);
                ye[0][p]=gRandom->Gaus(yei[0][p],(1/sqrt(3))*sS*sigSIe);
                
                thetaX[0][p+1]= gRandom->Gaus(thetaX[0][p],sigSImu);
                thetaY[0][p+1]= gRandom->Gaus(thetaY[0][p],sigSImu);
                thetaXe[0][p+1]= gRandom->Gaus(thetaXe[0][p],sigSIe);
                thetaYe[0][p+1]= gRandom->Gaus(thetaYe[0][p],sigSIe);

                xi[0][p+1]=x[0][p]+dSS*tan(thetaX[0][p]);
                x[0][p+1]=gRandom->Gaus(xi[0][p+1],(1/sqrt(3))*sS*sigSImu);
                yi[0][p+1]=y[0][p]+dSS*tan(thetaY[0][p]);
                y[0][p+1]=gRandom->Gaus(yi[0][p+1],(1/sqrt(3))*sS*sigSImu);
                xei[0][p+1]=xe[0][p]+dSS*tan(thetaXe[0][p]);
                xe[0][p+1]=gRandom->Gaus(xei[0][p+1],(1/sqrt(3))*sS*sigSIe);
                yei[0][p+1]=ye[0][p]+dSS*tan(thetaYe[0][p]);
                ye[0][p+1]=gRandom->Gaus(yei[0][p+1],(1/sqrt(3))*sS*sigSIe);
                 }}
        
                 thetaX[1][0]= gRandom->Gaus(thetaX[0][6],sigBEmu);
                 thetaY[1][0]= gRandom->Gaus(thetaY[0][6],sigBEmu); 
                 thetaXe[1][0]= gRandom->Gaus(thetaXe[0][6],sigBEe);
                 thetaYe[1][0]= gRandom->Gaus(thetaYe[0][6],sigBEe); 
                
                 xi[1][0]=x[0][6]+d*tan(thetaX[0][6]);
                 x[1][0]=gRandom->Gaus(xi[1][0],(1/sqrt(3))*sB*sigBEmu);
                 yi[1][0]=y[0][6]+d*tan(thetaY[0][6]);
                 y[1][0]=gRandom->Gaus(yi[1][0],(1/sqrt(3))*sB*sigBEmu);  
                 xei[1][0]=xe[0][6]+d*tan(thetaXe[0][6]);
                 xe[1][0]=gRandom->Gaus(xei[1][0],(1/sqrt(3))*sB*sigBEe);
                 yei[1][0]=ye[0][6]+d*tan(thetaYe[0][6]);
                 ye[1][0]=gRandom->Gaus(yei[1][0],(1/sqrt(3))*sB*sigBEe);  
   
        //entro nelle 3 coppie di silici
for (Int_t p=1; p<7; p++)  {
                 if(p==1 || p==3 || p==5){
                 thetaX[1][p]= gRandom->Gaus(thetaX[1][p-1],sigSImu);
                 thetaY[1][p]= gRandom->Gaus(thetaY[1][p-1],sigSImu); 
                 thetaXe[1][p]= gRandom->Gaus(thetaXe[1][p-1],sigSIe);
                 thetaYe[1][p]= gRandom->Gaus(thetaYe[1][p-1],sigSIe); 
                     
                     
                xi[1][p]=x[1][p-1]+d*tan(thetaX[1][p-1]);
                x[1][p]=gRandom->Gaus(xi[1][p],(1/sqrt(3))*sS*sigSImu);
                yi[1][p]=y[1][p-1]+d*tan(thetaY[1][p-1]);
                y[1][p]=gRandom->Gaus(yi[1][p],(1/sqrt(3))*sS*sigSImu);
                xei[1][p]=xe[1][p-1]+d*tan(thetaXe[1][p-1]);
                xe[1][p]=gRandom->Gaus(xei[1][p],(1/sqrt(3))*sS*sigSIe);
                yei[1][p]=ye[1][p-1]+d*tan(thetaYe[1][p-1]);
                ye[1][p]=gRandom->Gaus(yei[1][p],(1/sqrt(3))*sS*sigSIe);
                
                thetaX[1][p+1]= gRandom->Gaus(thetaX[1][p],sigSImu);
                thetaY[1][p+1]= gRandom->Gaus(thetaY[1][p],sigSImu);
                thetaXe[1][p+1]= gRandom->Gaus(thetaXe[1][p],sigSIe);
                thetaYe[1][p+1]= gRandom->Gaus(thetaYe[1][p],sigSIe);
                     
                xi[1][p+1]=x[1][p]+dSS*tan(thetaX[1][p]);
                x[1][p+1]=gRandom->Gaus(xi[1][p+1],(1/sqrt(3))*sS*sigSImu);
                yi[1][p+1]=y[1][p]+dSS*tan(thetaY[1][p]);
                y[1][p+1]=gRandom->Gaus(yi[1][p+1],(1/sqrt(3))*sS*sigSImu);
                xei[1][p+1]=xe[1][p]+dSS*tan(thetaXe[1][p]);
                xe[1][p+1]=gRandom->Gaus(xei[1][p+1],(1/sqrt(3))*sS*sigSIe);
                yei[1][p+1]=ye[1][p]+dSS*tan(thetaYe[1][p]);
                ye[1][p+1]=gRandom->Gaus(yei[1][p+1],(1/sqrt(3))*sS*sigSIe);                 
                 }}
        
    Double_t xf = x[1][6]+dCAL*tan(thetaX[1][6]);
    Double_t yf = y[1][6]+dCAL*tan(thetaY[1][6]);
    Double_t xfe = xe[1][6]+dCAL*tan(thetaXe[1][6]);
    Double_t yfe = ye[1][6]+dCAL*tan(thetaYe[1][6]);
    
    TMatrixD coo_ang_fin(19,7);
        
                for (Int_t i=0; i<7; i++)
    { //coordinate stazione 2 muone
        coo_ang_fin[0][i]=x[0][i]; 
        coo_ang_fin[1][i]=y[0][i]; 
        coo_ang_fin[2][i]=thetaX[0][i]; 
        coo_ang_fin[3][i]=thetaY[0][i]; 
    //coordinate stazione 3 muone         
        coo_ang_fin[4][i]=x[1][i]; 
        coo_ang_fin[5][i]=y[1][i]; 
        coo_ang_fin[6][i]=thetaX[1][i]; 
        coo_ang_fin[7][i]=thetaY[1][i];      
    //coordinate stazione 2 elettrone
        coo_ang_fin[9][i]=xe[0][i]; 
        coo_ang_fin[10][i]=ye[0][i]; 
        coo_ang_fin[11][i]=thetaXe[0][i]; 
        coo_ang_fin[12][i]=thetaYe[0][i]; 
    //coordinate stazione 3 elettrone          
        coo_ang_fin[13][i]=xe[1][i]; 
        coo_ang_fin[14][i]=ye[1][i]; 
        coo_ang_fin[15][i]=thetaXe[1][i]; 
        coo_ang_fin[16][i]=thetaYe[1][i];       
    }  
    
  Double_t dxdz = tan(thetaX[1][6]); // could approx tan ~ angle
  Double_t dydz = tan(thetaY[1][6]);
  
  // assuming z-motion is always forward
Double_t pmu=sqrt(k.Px()*k.Px()+k.Py()*k.Py()+k.Pz()*k.Pz());
  Double_t skz = sqrt((pmu*pmu)/ (dxdz*dxdz + dydz*dydz + 1));
  Double_t skx = skz * dxdz;
  Double_t sky = skz * dydz;
        
  Double_t dxdze = tan(thetaXe[1][6]); // could approx tan ~ angle
  Double_t dydze = tan(thetaYe[1][6]);
  
  // assuming z-motion is always forward
Double_t pe=sqrt(ke.Px()*ke.Px()+ke.Py()*ke.Py()+ke.Pz()*ke.Pz());
  Double_t skze = sqrt((pe*pe)/ (dxdze*dxdze + dydze*dydze + 1));
  Double_t skxe = skze * dxdze;
  Double_t skye = skze * dydze;
        
        // coordinate sul calorimetro e momento out smeared muone
        coo_ang_fin[8][0]=xf;   
        coo_ang_fin[8][1]=yf;  
        coo_ang_fin[8][2]=tar;      
        coo_ang_fin[8][3]=skx; 
        coo_ang_fin[8][4]=sky;      
        coo_ang_fin[8][5]=skz;      
        coo_ang_fin[8][6]=k.E();   
        
                
        // coordinate sul calorimetro e momento out smeared elettrone
        coo_ang_fin[17][0]=xfe;   
        coo_ang_fin[17][1]=yfe;  
        coo_ang_fin[17][2]=tar;      
        coo_ang_fin[17][3]=skxe; 
        coo_ang_fin[17][4]=skye;      
        coo_ang_fin[17][5]=skze;      
        coo_ang_fin[17][6]=ke.E();   
        
        //coordinate entranti beam divergente
        coo_ang_fin[18][0]=x2; 
        coo_ang_fin[18][1]=y2; 
        coo_ang_fin[18][2]=The; // angolo che ha nel momento dell'interazione, senza effetto di MCS 

        
        
        
    return coo_ang_fin;
    
    }
     
    
    if(tar==1) //peso 0 
    {
        // siamo nelle stazioni con il target di berillio. Ora entra nel berillio ad una distanza d=0.25-0.005 m dagli ultimi silici. Qui però non puoi trascurare lo spessore del berillio, cioè dove interagisce? considero a metà (7.5 mm), quindi aggiungo a d il pezzo in cui x non è modificato e poi sommo con il nuovo angolo di scattering


                 thetaX[0][0]=0;
                 thetaY[0][0]=0;
                 thetaXe[0][0]=0;
                 thetaYe[0][0]=0;

                 x[0][0]=0;
                 y[0][0]=0;  
                 xe[0][0]=0;
                 ye[0][0]=0;
       
//entro nelle 3 coppie di silici
for (Int_t p=1; p<7; p++)  {
                 if(p==1 || p==3 || p==5){
                 thetaX[0][p]=0;
                 thetaY[0][p]=0;
                 thetaXe[0][p]= 0;
                 thetaYe[0][p]= 0;

                x[0][p]=0;
                y[0][p]=0;
                xe[0][p]=0;
                ye[0][p]=0;
                

                thetaX[0][p+1]=0;
                thetaY[0][p+1]=0;
                thetaXe[0][p+1]= 0;
                thetaYe[0][p+1]= 0;

                x[0][p+1]=0;
                y[0][p+1]=0;
                xe[0][p+1]=0;
                ye[0][p+1]=0; 
                 }}

                 //coo dopo il silicio
                Double_t x1 = gRandom->Gaus(xR,(1/sqrt(3))*sSin*sigSI);
                Double_t y1 = gRandom->Gaus(yR,(1/sqrt(3))*sSin*sigSI);
        
                //coo dopo berillio
                Double_t x2 = gRandom->Gaus(x1,(1/sqrt(3))*2*sB*sigBEin);
                Double_t y2 = gRandom->Gaus(y1,(1/sqrt(3))*2*sB*sigBEin);
         
                //coo dopo il silicio
                Double_t x3 = gRandom->Gaus(x2,(1/sqrt(3))*sSin*sigSI);
                Double_t y3 = gRandom->Gaus(y2,(1/sqrt(3))*sSin*sigSI);
        
                
                xi[1][0]=(2-2*sSin-sB)*tan(anglexin)+x3;
                yi[1][0]=(2-2*sSin-sB)*tan(angleyin)+y3;
        
                //coo dopo metà berillio, poi avviene interazione
                Double_t x4 = gRandom->Gaus(xi[1][0],(1/sqrt(3))*(sB*vertex)*sigBE2in);
                Double_t y4 = gRandom->Gaus(yi[1][0],(1/sqrt(3))*(sB*vertex)*sigBE2in);
        
                //coo all'uscita del berillio dopo l'interazione
                 x[1][0]=gRandom->Gaus(x4,(1/sqrt(3))*(sB*(1-vertex))*sigBE2mu);
                 y[1][0]=gRandom->Gaus(y4,(1/sqrt(3))*(sB*(1-vertex))*sigBE2mu);  
    
                 xe[1][0]=gRandom->Gaus(x4,(1/sqrt(3))*(sB*(1-vertex))*sigBE2e);
                 ye[1][0]=gRandom->Gaus(y4,(1/sqrt(3))*(sB*(1-vertex))*sigBE2e); 
        
                
                thetaX[1][0]=gRandom->Gaus(anglex,sigBE2mu);
                thetaY[1][0]=gRandom->Gaus(angley,sigBE2mu);
                thetaXe[1][0]=gRandom->Gaus(anglexe,sigBE2e);
                thetaYe[1][0]=gRandom->Gaus(angleye,sigBE2e);

   
        //entro nelle 3 coppie di silici
for (Int_t p=1; p<7; p++)  {
                 if(p==1 || p==3 || p==5){
                 thetaX[1][p]= gRandom->Gaus(thetaX[1][p-1],sigSImu);
                 thetaY[1][p]= gRandom->Gaus(thetaY[1][p-1],sigSImu); 
                 thetaXe[1][p]= gRandom->Gaus(thetaXe[1][p-1],sigSIe);
                 thetaYe[1][p]= gRandom->Gaus(thetaYe[1][p-1],sigSIe); 
                     
                xi[1][p]=x[1][p-1]+d*tan(thetaX[1][p-1]);
                x[1][p]=gRandom->Gaus(xi[1][p],(1/sqrt(3))*sS*sigSImu);
                yi[1][p]=y[1][p-1]+d*tan(thetaY[1][p-1]);
                y[1][p]=gRandom->Gaus(yi[1][p],(1/sqrt(3))*sS*sigSImu);
                xei[1][p]=xe[1][p-1]+d*tan(thetaXe[1][p-1]);
                xe[1][p]=gRandom->Gaus(xei[1][p],(1/sqrt(3))*sS*sigSIe);
                yei[1][p]=ye[1][p-1]+d*tan(thetaYe[1][p-1]);
                ye[1][p]=gRandom->Gaus(yei[1][p],(1/sqrt(3))*sS*sigSIe);
                
                thetaX[1][p+1]= gRandom->Gaus(thetaX[1][p],sigSImu);
                thetaY[1][p+1]= gRandom->Gaus(thetaY[1][p],sigSImu);
                thetaXe[1][p+1]= gRandom->Gaus(thetaXe[1][p],sigSIe);
                thetaYe[1][p+1]= gRandom->Gaus(thetaYe[1][p],sigSIe);

                xi[1][p+1]=x[1][p]+dSS*tan(thetaX[1][p]);
                x[1][p+1]=gRandom->Gaus(xi[1][p+1],(1/sqrt(3))*sS*sigSImu);
                yi[1][p+1]=y[1][p]+dSS*tan(thetaY[1][p]);
                y[1][p+1]=gRandom->Gaus(yi[1][p+1],(1/sqrt(3))*sS*sigSImu);
                xei[1][p+1]=xe[1][p]+dSS*tan(thetaXe[1][p]);
                xe[1][p+1]=gRandom->Gaus(xei[1][p+1],(1/sqrt(3))*sS*sigSIe);
                yei[1][p+1]=ye[1][p]+dSS*tan(thetaYe[1][p]);
                ye[1][p+1]=gRandom->Gaus(yei[1][p+1],(1/sqrt(3))*sS*sigSIe);                 
                 }}
        // sul calorimetro
        
    Double_t xf = x[1][6]+dCAL*tan(thetaX[1][6]);
    Double_t yf = y[1][6]+dCAL*tan(thetaY[1][6]);
    Double_t xfe = xe[1][6]+dCAL*tan(thetaXe[1][6]);
    Double_t yfe = ye[1][6]+dCAL*tan(thetaYe[1][6]);
    
    TMatrixD coo_ang_fin(19,7);
        
                for (Int_t i=0; i<7; i++)
    { //coordinate stazione 2 muone
        coo_ang_fin[0][i]=x[0][i]; 
        coo_ang_fin[1][i]=y[0][i]; 
        coo_ang_fin[2][i]=thetaX[0][i]; 
        coo_ang_fin[3][i]=thetaY[0][i]; 
    //coordinate stazione 3 muone         
        coo_ang_fin[4][i]=x[1][i]; 
        coo_ang_fin[5][i]=y[1][i]; 
        coo_ang_fin[6][i]=thetaX[1][i]; 
        coo_ang_fin[7][i]=thetaY[1][i];      
    //coordinate stazione 2 elettrone
        coo_ang_fin[9][i]=xe[0][i]; 
        coo_ang_fin[10][i]=ye[0][i]; 
        coo_ang_fin[11][i]=thetaXe[0][i]; 
        coo_ang_fin[12][i]=thetaYe[0][i]; 
    //coordinate stazione 3 elettrone          
        coo_ang_fin[13][i]=xe[1][i]; 
        coo_ang_fin[14][i]=ye[1][i]; 
        coo_ang_fin[15][i]=thetaXe[1][i]; 
        coo_ang_fin[16][i]=thetaYe[1][i];       
    }  
    
  Double_t dxdz = tan(thetaX[1][6]); // could approx tan ~ angle
  Double_t dydz = tan(thetaY[1][6]);
  
  // assuming z-motion is always forward
Double_t pmu=sqrt(k.Px()*k.Px()+k.Py()*k.Py()+k.Pz()*k.Pz());
  Double_t skz = sqrt((pmu*pmu)/ (dxdz*dxdz + dydz*dydz + 1));
  Double_t skx = skz * dxdz;
  Double_t sky = skz * dydz;
        
  Double_t dxdze = tan(thetaXe[1][6]); // could approx tan ~ angle
  Double_t dydze = tan(thetaYe[1][6]);
  
  // assuming z-motion is always forward
Double_t pe=sqrt(ke.Px()*ke.Px()+ke.Py()*ke.Py()+ke.Pz()*ke.Pz());
  Double_t skze = sqrt((pe*pe) / (dxdze*dxdze + dydze*dydze + 1));
  Double_t skxe = skze * dxdze;
  Double_t skye = skze * dydze;
        
        // coordinate sul calorimetro e momento out smeared
        coo_ang_fin[8][0]=xf;   
        coo_ang_fin[8][1]=yf;  
        coo_ang_fin[8][2]=tar;      
        coo_ang_fin[8][3]=skx; 
        coo_ang_fin[8][4]=sky;      
        coo_ang_fin[8][5]=skz;      
        coo_ang_fin[8][6]=k.E();   
        
                
        // coordinate sul calorimetro e momento out smeared
        coo_ang_fin[17][0]=xfe;   
        coo_ang_fin[17][1]=yfe;  
        coo_ang_fin[17][2]=tar;      
        coo_ang_fin[17][3]=skxe; 
        coo_ang_fin[17][4]=skye;      
        coo_ang_fin[17][5]=skze;      
        coo_ang_fin[17][6]=ke.E();   
        
        //coordinate entranti beam divergente
        coo_ang_fin[18][0]=x4; 
        coo_ang_fin[18][1]=y4;  
        coo_ang_fin[18][2]=The; // angolo che ha nel momento dell'interazione, senza effetto di MCS 
        
        
    return coo_ang_fin;
        
    }

else return coo_in; 

}



TMatrixD FastSim::MCSphoton(const Double_t & tar, const PxPyPzEVector & kp,const Double_t & xin, const Double_t & yin) const

{ Double_t d0 =2.025; //m
  Double_t d1 =1.025;
// Double_t thetaR = theta*0.01;// rad
         TMatrixD coo (1,2);
        coo[0][0]=0;
        coo[0][1]=0;
 
    Double_t anglex = atan2(kp.Px(), kp.Pz());
    Double_t angley = atan2(kp.Py(), kp.Pz()); 
 
    if(tar==0)
    {
        /*
    Double_t d_xy = d0*tan(thetaR);//vettore nel piano xy
    Double_t xf = xin+d_xy*cos(phi);
    Double_t yf = yin+d_xy*sin(phi);
        TMatrixD coo (1,2);*/
        
        Double_t xf=xin+d0*tan(anglex);
        Double_t yf=yin+d0*tan(angley);
        
        coo[0][0]=xf;
        coo[0][1]=yf;

        return coo;
    }

    if(tar==1)
    {
    /*Double_t d_xy = d1*tan(thetaR);//vettore nel piano xy
    Double_t xf = xin+d_xy*cos(phi);
    Double_t yf = yin+d_xy*sin(phi);*/
        
        Double_t xf=xin+d1*tan(anglex);
        Double_t yf=yin+d1*tan(angley); 
            
        TMatrixD coo (1,2);
        coo[0][0]=xf;
        coo[0][1]=yf;

        return coo;
    }
 else return coo;
}


TMatrixD FastSim::Def_angle(const PxPyPzEVector & p_mu_in_div,const PxPyPzEVector & p_mu_out_div,const PxPyPzEVector & p_e_out_div) const {
    


XYZVector p_mu_in_div3 = p_mu_in_div.Vect().Unit();    
XYZVector p_e_out_div3 = p_e_out_div.Vect().Unit();
XYZVector p_mu_out_div3 = p_mu_out_div.Vect().Unit();


    
Double_t DIR_mu=p_mu_in_div3.Dot(p_mu_out_div3);
Double_t A_DIR_mu=abs(DIR_mu)<1. ? 1000.*acos(DIR_mu) : 100.;
    
Double_t DIR_e=p_mu_in_div3.Dot(p_e_out_div3);
Double_t A_DIR_e=abs(DIR_e)<1. ? 1000.*acos(DIR_e) : 100.;

TMatrixD def_angle(2,1);
  def_angle[0][0]=A_DIR_mu;
  def_angle[1][0]=A_DIR_e;
    
return def_angle;
}


Int_t FastSim::ECALe(const Double_t & cooXe, const Double_t & cooYe) const
{       
    Int_t n_cell; //numero di cella in cui cade l'ELETTRONE
   
    
    Double_t detKinBeamRot_cooXe=cooXe*100;
    Double_t detKinBeamRot_cooYe=cooYe*100;

    
    if (detKinBeamRot_cooXe>-7.125 && detKinBeamRot_cooXe<-4.275 && detKinBeamRot_cooYe<7.125 && detKinBeamRot_cooYe>4.275) {n_cell=1;}
    
    if (detKinBeamRot_cooXe>-4.275 && detKinBeamRot_cooXe<-1.425 && detKinBeamRot_cooYe<7.125 && detKinBeamRot_cooYe>4.275) {n_cell=2;}

    if (detKinBeamRot_cooXe>-1.425 && detKinBeamRot_cooXe<1.425 && detKinBeamRot_cooYe<7.125 && detKinBeamRot_cooYe>4.275) {n_cell=3;}

    if (detKinBeamRot_cooXe>1.425 && detKinBeamRot_cooXe<4.275 && detKinBeamRot_cooYe<7.125 && detKinBeamRot_cooYe>4.275) {n_cell=4;}

    if (detKinBeamRot_cooXe>4.275 && detKinBeamRot_cooXe<7.125 && detKinBeamRot_cooYe<7.125 && detKinBeamRot_cooYe>4.275) {n_cell=5;}
    
    
    if (detKinBeamRot_cooXe>-7.125 && detKinBeamRot_cooXe<-4.275 && detKinBeamRot_cooYe<4.275 && detKinBeamRot_cooYe>1.425) {n_cell=6;}
    
    if (detKinBeamRot_cooXe>-4.275 && detKinBeamRot_cooXe<-1.425 && detKinBeamRot_cooYe<4.275 && detKinBeamRot_cooYe>1.425) {n_cell=7;}

    if (detKinBeamRot_cooXe>-1.425 && detKinBeamRot_cooXe<1.425 && detKinBeamRot_cooYe<4.275 && detKinBeamRot_cooYe>1.425) {n_cell=8;}

    if (detKinBeamRot_cooXe>1.425 && detKinBeamRot_cooXe<4.275 && detKinBeamRot_cooYe<4.275 && detKinBeamRot_cooYe>1.425) {n_cell=9;}

    if (detKinBeamRot_cooXe>4.275 && detKinBeamRot_cooXe<7.125 && detKinBeamRot_cooYe<4.275 && detKinBeamRot_cooYe>1.425) {n_cell=10;}
    

    
    if (detKinBeamRot_cooXe>-7.125 && detKinBeamRot_cooXe<-4.275 && detKinBeamRot_cooYe<1.425 && detKinBeamRot_cooYe>-1.425) {n_cell=11;}
    
    if (detKinBeamRot_cooXe>-4.275 && detKinBeamRot_cooXe<-1.425 && detKinBeamRot_cooYe<1.425 && detKinBeamRot_cooYe>-1.425) {n_cell=12;}

    if (detKinBeamRot_cooXe>-1.425 && detKinBeamRot_cooXe<1.425 && detKinBeamRot_cooYe<1.425 && detKinBeamRot_cooYe>-1.425) {n_cell=13;}

    if (detKinBeamRot_cooXe>1.425 && detKinBeamRot_cooXe<4.275 && detKinBeamRot_cooYe<1.425 && detKinBeamRot_cooYe>-1.425) {n_cell=14;}

    if (detKinBeamRot_cooXe>4.275 && detKinBeamRot_cooXe<7.125 && detKinBeamRot_cooYe<1.425 && detKinBeamRot_cooYe>-1.425) {n_cell=15;}
    
    
    
    
    if (detKinBeamRot_cooXe>-7.125 && detKinBeamRot_cooXe<-4.275 && detKinBeamRot_cooYe<-1.425 && detKinBeamRot_cooYe>-4.275) {n_cell=16;}
    
    if (detKinBeamRot_cooXe>-4.275 && detKinBeamRot_cooXe<-1.425 && detKinBeamRot_cooYe<-1.425 && detKinBeamRot_cooYe>-4.275) {n_cell=17;}

    if (detKinBeamRot_cooXe>-1.425 && detKinBeamRot_cooXe<1.425 && detKinBeamRot_cooYe<-1.425 && detKinBeamRot_cooYe>-4.275) {n_cell=18;}

    if (detKinBeamRot_cooXe>1.425 && detKinBeamRot_cooXe<4.275 && detKinBeamRot_cooYe<-1.425 && detKinBeamRot_cooYe>-4.275) {n_cell=19;}

    if (detKinBeamRot_cooXe>4.275 && detKinBeamRot_cooXe<7.125 && detKinBeamRot_cooYe<-1.425 && detKinBeamRot_cooYe>-4.275) {n_cell=20;}
    
    
    
    if (detKinBeamRot_cooXe>-7.125 && detKinBeamRot_cooXe<-4.275 && detKinBeamRot_cooYe<-4.275 && detKinBeamRot_cooYe>-7.125) {n_cell=21;}
    
    if (detKinBeamRot_cooXe>-4.275 && detKinBeamRot_cooXe<-1.425 && detKinBeamRot_cooYe<-4.275 && detKinBeamRot_cooYe>-7.125) {n_cell=22;}

    if (detKinBeamRot_cooXe>-1.425 && detKinBeamRot_cooXe<1.425 && detKinBeamRot_cooYe<-4.275 && detKinBeamRot_cooYe>-7.125) {n_cell=23;}

    if (detKinBeamRot_cooXe>1.425 && detKinBeamRot_cooXe<4.275 && detKinBeamRot_cooYe<-4.275 && detKinBeamRot_cooYe>-7.125) {n_cell=24;}

    if (detKinBeamRot_cooXe>4.275 && detKinBeamRot_cooXe<7.125 && detKinBeamRot_cooYe<-4.275 && detKinBeamRot_cooYe>-7.125) {n_cell=25;}
        
return n_cell;
}

Int_t FastSim::ECALph(const Double_t & coox, const Double_t & cooy) const
{       
    Int_t n_cell_ph; //numero di cella in cui cade il FOTONE
  
    
    Double_t photon_coox=coox*100;
    Double_t photon_cooy=cooy*100;

    
    if (photon_coox>-7.125 && photon_coox<-4.275 && photon_cooy<7.125 && photon_cooy>4.275) {n_cell_ph=1;}
    
    if (photon_coox>-4.275 && photon_coox<-1.425 && photon_cooy<7.125 && photon_cooy>4.275) {n_cell_ph=2;}

    if (photon_coox>-1.425 && photon_coox<1.425 && photon_cooy<7.125 && photon_cooy>4.275) {n_cell_ph=3;}

    if (photon_coox>1.425 && photon_coox<4.275 && photon_cooy<7.125 && photon_cooy>4.275) {n_cell_ph=4;}

    if (photon_coox>4.275 && photon_coox<7.125 && photon_cooy<7.125 && photon_cooy>4.275) {n_cell_ph=5;}
    
    
    if (photon_coox>-7.125 && photon_coox<-4.275 && photon_cooy<4.275 && photon_cooy>1.425) {n_cell_ph=6;}
    
    if (photon_coox>-4.275 && photon_coox<-1.425 && photon_cooy<4.275 && photon_cooy>1.425) {n_cell_ph=7;}

    if (photon_coox>-1.425 && photon_coox<1.425 && photon_cooy<4.275 && photon_cooy>1.425) {n_cell_ph=8;}

    if (photon_coox>1.425 && photon_coox<4.275 && photon_cooy<4.275 && photon_cooy>1.425) {n_cell_ph=9;}

    if (photon_coox>4.275 && photon_coox<7.125 && photon_cooy<4.275 && photon_cooy>1.425) {n_cell_ph=10;}
    

    
    if (photon_coox>-7.125 && photon_coox<-4.275 && photon_cooy<1.425 && photon_cooy>-1.425) {n_cell_ph=11;}
    
    if (photon_coox>-4.275 && photon_coox<-1.425 && photon_cooy<1.425 && photon_cooy>-1.425) {n_cell_ph=12;}

    if (photon_coox>-1.425 && photon_coox<1.425 && photon_cooy<1.425 && photon_cooy>-1.425) {n_cell_ph=13;}

    if (photon_coox>1.425 && photon_coox<4.275 && photon_cooy<1.425 && photon_cooy>-1.425) {n_cell_ph=14;}

    if (photon_coox>4.275 && photon_coox<7.125 && photon_cooy<1.425 && photon_cooy>-1.425) {n_cell_ph=15;}
    
    
    
    
    if (photon_coox>-7.125 && photon_coox<-4.275 && photon_cooy<-1.425 && photon_cooy>-4.275) {n_cell_ph=16;}
    
    if (photon_coox>-4.275 && photon_coox<-1.425 && photon_cooy<-1.425 && photon_cooy>-4.275) {n_cell_ph=17;}

    if (photon_coox>-1.425 && photon_coox<1.425 && photon_cooy<-1.425 && photon_cooy>-4.275) {n_cell_ph=18;}

    if (photon_coox>1.425 && photon_coox<4.275 && photon_cooy<-1.425 && photon_cooy>-4.275) {n_cell_ph=19;}

    if (photon_coox>4.275 && photon_coox<7.125 && photon_cooy<-1.425 && photon_cooy>-4.275) {n_cell_ph=20;}
    
    
    
    if (photon_coox>-7.125 && photon_coox<-4.275 && photon_cooy<-4.275 && photon_cooy>-7.125) {n_cell_ph=21;}
    
    if (photon_coox>-4.275 && photon_coox<-1.425 && photon_cooy<-4.275 && photon_cooy>-7.125) {n_cell_ph=22;}

    if (photon_coox>-1.425 && photon_coox<1.425 && photon_cooy<-4.275 && photon_cooy>-7.125) {n_cell_ph=23;}

    if (photon_coox>1.425 && photon_coox<4.275 && photon_cooy<-4.275 && photon_cooy>-7.125) {n_cell_ph=24;}

    if (photon_coox>4.275 && photon_coox<7.125 && photon_cooy<-4.275 && photon_cooy>-7.125) {n_cell_ph=25;} 

        
return n_cell_ph;
}


void FastSim::LoadKineVars(const PxPyPzEVector & p_mu_in,  const PxPyPzEVector & p_e_in, 
			   const PxPyPzEVector & p_mu_out, const PxPyPzEVector & p_e_out, const TMatrixD & coo, const Double_t & TheINT,
			   MuE::KineVars & kv) {
  
  kv.Ee = p_e_out.E();
  kv.Emu = p_mu_out.E();
  kv.the = 1e3* p_e_out.Theta();
  kv.thmu = 1e3* p_mu_out.Theta();
  kv.phe = p_e_out.Phi();
  kv.phmu = p_mu_out.Phi();
 
    

kv.cooXmu = coo[0][0];
kv.cooYmu = coo[1][0];
kv.cooXe = coo[2][0];
kv.cooYe = coo[3][0];

kv.tar = coo[4][0];
 kv.ThEl_interaction = 1e3*TheINT; // angolo elettrone al momento della produzione, in mrad

  kv.pXmu = p_mu_in.Px();
  kv.pYmu = p_mu_in.Py();
  kv.pZmu = p_mu_in.Pz();
  kv.pXe = p_e_in.Px();
  kv.pYe = p_e_in.Py();
  kv.pZe = p_e_in.Pz();
  kv.pXmu_out = p_mu_out.Px();
  kv.pYmu_out = p_mu_out.Py();
  kv.pZmu_out = p_mu_out.Pz();
  kv.pXe_out = p_e_out.Px();
  kv.pYe_out = p_e_out.Py();
  kv.pZe_out = p_e_out.Pz();
kv.Pmu_out = p_mu_out.P();
kv.Pe_out = p_e_out.P();

TMatrixD def_angle=Def_angle(p_mu_in,p_mu_out,p_e_out);
    
kv.def_angle_mu = def_angle[0][0];
kv.def_angle_e = def_angle[1][0]; 
    
    if (abs(kv.cooXe)<0.07125 && abs(kv.cooYe)<0.07125)
    {kv.n_cell_e = ECALe(kv.cooXe,kv.cooYe);}
    else {kv.n_cell_e=0;}
    
  // Note: here Ebeam is the average beam energy, so tt_e and xt_e are defined under this assumption
  MuE::ElasticState emu_state(Ebeam,mm,me, kv.the);
  kv.tt_e = emu_state.GetT();
  kv.xt_e = emu_state.GetX();
  
  PxPyPzEVector q13 = p_mu_in - p_mu_out;
  PxPyPzEVector q24 = p_e_in - p_e_out;
  
  // instead these are exact
  kv.t13 = q13.M2();
  kv.t24 = q24.M2();
  kv.x13 = emu_state.X(kv.t13); 
  kv.x24 = emu_state.X(kv.t24);
  
  XYZVector nvec_mu_in = p_mu_in.Vect().Unit(); 
  XYZVector nvec_mu_out = p_mu_out.Vect().Unit();
  XYZVector nvec_e_out = p_e_out.Vect().Unit();

  // acoplanarity = deltaPhi as long as the incoming muon is collinear with z-axis
  kv.deltaPhi = kv.phe - kv.phmu;
  if (kv.deltaPhi<0.) kv.deltaPhi = kv.deltaPhi + 2.*TMath::Pi();
  kv.deltaPhi = kv.deltaPhi - TMath::Pi();
  
  double dotProduct = nvec_mu_out.Dot(nvec_e_out);
  kv.openingAngle = std::abs(dotProduct)<1. ? 1000.*acos(dotProduct) : 0.;

  XYZVector crossProduct = nvec_mu_out.Cross(nvec_e_out);
  kv.tripleProduct = nvec_mu_in.Dot(crossProduct);
}


      
void FastSim::LoadPhoton(const MuE::Event & event, MuE::Photon & photon,const PxPyPzEVector & p_mu_in,const Double_t & tar,const Double_t & xin,const Double_t & yin) {
  // by now at most one photon
  auto n_photons = event.photons.size();
  
  if (n_photons >0) {  
    PxPyPzEVector p_gamma_Lab = {
                 event.photons[0].px, 
				 event.photons[0].py,
				 event.photons[0].pz,
			     event.photons[0].E};

      
PxPyPzEVector p_gamma_Lab_div=RotDiv(p_mu_in,p_gamma_Lab);
PxPyPzEVector p_gamma_CoM = Lorentz_ToCoM(p_gamma_Lab_div);
  
    
    photon.energy    = p_gamma_Lab_div.E();
    photon.theta     = p_gamma_Lab_div.Theta() *1e3;
    photon.phi       = p_gamma_Lab_div.Phi();
    photon.energyCoM = p_gamma_CoM.E(); 

    TMatrixD coo=MCSphoton(tar,p_gamma_Lab,xin,yin);
      
    photon.coox=coo[0][0];
    photon.cooy=coo[0][1];
      
   if (abs(photon.coox)<0.07125 && abs(photon.cooy)<0.07125 && photon.energy>0.2)
   {photon.n_cell_ph = ECALph(photon.coox,photon.cooy);}
    else{photon.n_cell_ph=0;}   
   
    
  }

  else {
    photon.energyCoM = -1;
    photon.energy    = -1;
    photon.theta     = -1;
    photon.phi       =  0;
    photon.coox     = -1;
    photon.cooy     = -1;
    
    photon.n_cell_ph=0;
  }
      
      
}

// synchronize the random number chain to account for events with negligible weight
//  (skipped in the main event loop)
//
void FastSim::RandomNrSync() {
  gRandom->Gaus(0.,1.); // fake smearing as for the outgoing muon theta X
  gRandom->Gaus(0.,1.); // fake smearing as for the outgoing muon theta Y
  gRandom->Gaus(0.,1.); // fake smearing as for the outgoing electron theta X
  gRandom->Gaus(0.,1.); // fake smearing as for the outgoing electron theta Y
}
