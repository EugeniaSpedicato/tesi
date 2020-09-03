#include <cmath>
#include <iostream>
#include "TRandom.h"
#include "TMath.h"
#include "ElasticState.h"
#include "ResolutionModels.h"
#include "FastSim.h"
#include "TMatrixF.h" 
#include "TMatrixD.h"
#include "TMatrixFBase.h"
#include <TMatrixFSym.h>
#include "TString.h"
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
   
   /* LoadKineVars(p_mu_in, p_e_in, p_mu_out, p_e_out, genKin); 
     
    //DETKIN NORMALI
    PxPyPzEVector p_mu_out_smeared, p_e_out_smeared; 
  if (MSopt ==0) {
    p_mu_out_smeared = Smear(p_mu_out);
    p_e_out_smeared = Smear(p_e_out);
  }
  else if (MSopt ==1) {
    p_mu_out_smeared = SmearX(p_mu_out);
    p_e_out_smeared = SmearX(p_e_out);
  }
  else if (MSopt ==2) {
    p_mu_out_smeared = SmearPolar(p_mu_out);
    p_e_out_smeared = SmearPolar(p_e_out);
  }
  else cout<<"\n"<<"*** ERROR : FastSim, undefined detector MS option = "<<MSopt<<endl;
  
  LoadKineVars(p_mu_in, p_e_in, p_mu_out_smeared, p_e_out_smeared, detKin);
  
   */ 
    //DETKIN RUOTATE (const PxPyPzEVector & k, const Double_t tar, const Double_t xx, const Double_t yy, const Double_t thX, const Double_t thY)
TMatrixD muin=RotDivIN(p_mu_in);
PxPyPzEVector p_mu_in_div(muin[0][0],muin[1][0],muin[2][0],muin[3][0]);
//TMatrixD a=MCSin(p_mu_in_div);//effetto MCS, ritorna una matrice con coo, angoli e momento
//PxPyPzEVector p_mu_in_smeared(a[4][0],a[4][1],a[4][2],a[4][3]);

PxPyPzEVector p_e_in_div=(p_e_in);
    
//cioe div della cinematica e smear del MCS
PxPyPzEVector p_mu_out_div=RotDiv(p_mu_in_div,p_mu_out);
PxPyPzEVector p_e_out_div=RotDiv(p_mu_in_div,p_e_out);
  

    
    
//effetto MCS out, ritorna una matrice con coo, angoli e momento per muone ed elettrone
TMatrixD b=MCSout(p_mu_in_div,p_mu_out_div,p_e_out_div,muin[4][0]);
    
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
  
  LoadKineVars(p_mu_in_div, p_e_in_div, p_mu_out_div_smeared, p_e_out_div_smeared, coo, detKinBeamRot);
  
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

// RMS Theta smearing due to Multiple scattering distribution and intrinsic resolution
//
/*
Double_t FastSim::ThetaRMS(const PxPyPzEVector & k) const
{
  Double_t pmom = k.P();
  
  // resolution: gaussian sigma
  Double_t thrms(0); 
  
  if (model == 0) {
    Double_t msc = thickness > 0 ? MuE::Target_thrms(pmom, thickness) : 0; // in mrad
    thrms = sqrt(msc*msc + intrinsic_resolution*intrinsic_resolution);
  }
  else if (model == 1) {
    thrms = MuE::Antonio_thrms(pmom); // in mrad
  }
  else {
    cout << "\n" << "***ERROR: Undefined smearing model = "<<model << endl;
    exit(999);
  }

  return thrms;
}


*/


 TMatrixD FastSim::RotDivIN(const PxPyPzEVector & k) const
 {
     
Double_t const sS    = 3*0.00128; //m spessore silicio
Double_t const x0S = 0.094; // m
Double_t const sB    = 0.015; //m spessore berillio
Double_t const x0B = 0.353; // m

    
Double_t sigSI=(13.6/(k.E()*1000))*sqrt(sS/x0S)*(1+0.038*log(sS/x0S)); //rad
Double_t sigBE=(13.6/(k.E()*1000))*sqrt(sB/x0B)*(1+0.038*log(sB/x0B)); //rad   
Double_t sigBE2in=(13.6/(k.E()*1000))*sqrt(sB/(2*x0B))*(1+0.038*log(sB/(2*x0B))); //rad   


Double_t divthx = gRandom->Gaus(0., 0.00027);
Double_t divthy = gRandom->Gaus(0., 0.00020); 
    
Int_t tar=gRandom->Integer(2);
     
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


TMatrixD pnewdiv(5,1);
pnewdiv[0][0]=px;
pnewdiv[1][0]=py;
pnewdiv[2][0]=pz; 
pnewdiv[3][0]=k.E();
pnewdiv[4][0]=tar;
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
     
TMatrixD p(2,1);
p[0][0]=0;
p[1][0]=0;
return p;
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


TMatrixD FastSim::MCSout(const PxPyPzEVector & kin, const PxPyPzEVector & k, const PxPyPzEVector & ke, const Double_t & tar) const
{
        
    
Double_t const sSin    = 3*0.00128; //m spessore silicio per fascio entrante

Double_t const sS    = 0.00064; //m spessore silicio
Double_t const x0S = 0.094; // m

Double_t const sB    = 0.015; //m spessore berillio
Double_t const x0B = 0.353; // m



Double_t const dSS = 0.01; // m distanza tra i due 2S


Double_t const d = 0.25-0.005;
Double_t const dx = 0.25+0.005; // m la distanza tra coppie di silici è 0.25. Però 0.25 è tra i due della coppia, quindi considero -dSS/2=0.025

// m la distanza tra coppie di silici è 0.25. Però 0.25 è tra i due della coppia, quindi considero -dSS/2=0.005
Double_t const ris = 18e-6; // m considero questa la risoluzione dei silici
Double_t const dCAL = 0.10; // m distanza silicio calorimetro
    

Double_t sigSIinP=(13.6/(kin.E()*1000))*sqrt(sSin/x0S)*(1+0.038*log(sSin/x0S)); //rad
Double_t sigBE2in=(13.6/(kin.E()*1000))*sqrt(sB/(2*x0B))*(1+0.038*log(sB/(2*x0B))); //rad   
Double_t sigBEin=(13.6/(kin.E()*1000))*sqrt(sB/x0B)*(1+0.038*log(sB/x0B)); //rad   
    
    
Double_t sigSImu=(13.6/(k.E()*1000))*sqrt(sS/x0S)*(1+0.038*log(sS/x0S)); //rad
// considero sB/2 per quandp interagisce a metà 
Double_t sigSIe=(13.6/(ke.E()*1000))*sqrt(sS/x0S)*(1+0.038*log(sS/x0S)); //rad
// considero sB/2 per quandp interagisce a metà 
Double_t sigBE2mu=(13.6/(k.E()*1000))*sqrt(sB/(2*x0B))*(1+0.038*log(sB/(2*x0B))); //rad
// considero sB/2 per quandp interagisce a metà 
Double_t sigBE2e=(13.6/(ke.E()*1000))*sqrt(sB/(2*x0B))*(1+0.038*log(sB/(2*x0B))); //rad
// considero sB/2 per quandp interagisce a metà 
Double_t sigBEmu=(13.6/(k.E()*1000))*sqrt(sB/x0B)*(1+0.038*log(sB/x0B)); //rad   
Double_t sigBEe=(13.6/(ke.E()*1000))*sqrt(sB/x0B)*(1+0.038*log(sB/x0B)); //rad   
    

    Double_t THinX[7];
    Double_t THinY[7];
    Double_t inX[7];
    Double_t inY[7];
    
    TMatrixD thetaX(2,7);
    TMatrixD thetaY(2,7);
    TMatrixD x(2,7);
    TMatrixD y(2,7);
    
    TMatrixD thetaXe(2,7);
    TMatrixD thetaYe(2,7);
    TMatrixD xe(2,7);
    TMatrixD ye(2,7);
    
    Double_t anglexin = atan2(kin.Px(), kin.Pz());
    Double_t angleyin = atan2(kin.Py(), kin.Pz()); 
  
    Double_t anglex = atan2(k.Px(), k.Pz());
    Double_t angley = atan2(k.Py(), k.Pz()); 
    Double_t anglexe = atan2(ke.Px(), ke.Pz());
    Double_t angleye = atan2(ke.Py(), ke.Pz()); 
    
    
    TMatrixD coo_in(1,3);
    coo_in[0][0]=0;
    coo_in[0][1]=0;
    coo_in[0][2]=0;
 
    
    if(tar==0)
    {
       // siamo nelle stazioni con il target di berillio. Ora entra nel berillio ad una distanza d=0.25-0.005 m dagli ultimi silici. Qui però non puoi trascurare lo spessore del berillio, cioè dove interagisce? considero a metà (7.5 mm), quindi aggiungo a d il pezzo in cui x non è modificato e poi sommo con il nuovo angolo di scattering
        
                Double_t xin=(1-sSin)*tan(anglexin)+(1/sqrt(3))*sSin*sigSIinP+(1/sqrt(3))*0.0075*sigBE2in;
                Double_t yin=(1-sSin)*tan(angleyin)+(1/sqrt(3))*sSin*sigSIinP+(1/sqrt(3))*0.0075*sigBE2in;
        
                 thetaX[0][0]= gRandom->Gaus(anglex,sigBE2mu);
                 thetaY[0][0]= gRandom->Gaus(angley,sigBE2mu); 
                 thetaXe[0][0]= gRandom->Gaus(anglexe,sigBE2e);
                 thetaYe[0][0]= gRandom->Gaus(angleye,sigBE2e); 

                 x[0][0]=xin+(1/sqrt(3))*0.0075*sigBE2mu;
                 y[0][0]=yin+(1/sqrt(3))*0.0075*sigBE2mu;  
    
                 xe[0][0]=xin+(1/sqrt(3))*0.0075*sigBE2e;
                 ye[0][0]=yin+(1/sqrt(3))*0.0075*sigBE2e; 
            


        
              
//entro nelle 3 coppie di silici
for (Int_t p=1; p<7; p++)  {
                 if(p==1 || p==3 || p==5){
                 thetaX[0][p]= gRandom->Gaus(thetaX[0][p-1],sigSImu);
                 thetaY[0][p]= gRandom->Gaus(thetaY[0][p-1],sigSImu); 
                 thetaXe[0][p]= gRandom->Gaus(thetaXe[0][p-1],sigSIe);
                 thetaYe[0][p]= gRandom->Gaus(thetaYe[0][p-1],sigSIe); 

                x[0][p]=x[0][p-1]+d*tan(thetaX[0][p-1])+(1/sqrt(3))*sS*sigSImu+gRandom->Gaus(0,ris);
                y[0][p]=y[0][p-1]+d*tan(thetaY[0][p-1])+(1/sqrt(3))*sS*sigSImu;
                xe[0][p]=xe[0][p-1]+d*tan(thetaXe[0][p-1])+(1/sqrt(3))*sS*sigSIe+gRandom->Gaus(0,ris);
                ye[0][p]=ye[0][p-1]+d*tan(thetaYe[0][p-1])+(1/sqrt(3))*sS*sigSIe;
                
                thetaX[0][p+1]= gRandom->Gaus(thetaX[0][p],sigSImu);
                thetaY[0][p+1]= gRandom->Gaus(thetaY[0][p],sigSImu);
                thetaXe[0][p+1]= gRandom->Gaus(thetaXe[0][p],sigSIe);
                thetaYe[0][p+1]= gRandom->Gaus(thetaYe[0][p],sigSIe);

                x[0][p+1]=x[0][p]+dSS*tan(thetaX[0][p])+(1/sqrt(3))*sS*sigSImu;
                y[0][p+1]=y[0][p]+dSS*tan(thetaY[0][p])+(1/sqrt(3))*sS*sigSImu+gRandom->Gaus(0,ris);
                xe[0][p+1]=xe[0][p]+dSS*tan(thetaXe[0][p])+(1/sqrt(3))*sS*sigSIe;
                ye[0][p+1]=ye[0][p]+dSS*tan(thetaYe[0][p])+(1/sqrt(3))*sS*sigSIe+gRandom->Gaus(0,ris);
                 }}
        
                 thetaX[1][0]= gRandom->Gaus(thetaX[0][6],sigBEmu);
                 thetaY[1][0]= gRandom->Gaus(thetaY[0][6],sigBEmu); 
                 thetaXe[1][0]= gRandom->Gaus(thetaXe[0][6],sigBEe);
                 thetaYe[1][0]= gRandom->Gaus(thetaYe[0][6],sigBEe); 

                 x[1][0]=x[0][6]+d*tan(thetaX[0][6])+(1/sqrt(3))*sB*sigBEmu;
                 y[1][0]=y[0][6]+d*tan(thetaY[0][6])+(1/sqrt(3))*sB*sigBEmu;  
                 xe[1][0]=xe[0][6]+d*tan(thetaXe[0][6])+(1/sqrt(3))*sB*sigBEe;
                 ye[1][0]=ye[0][6]+d*tan(thetaYe[0][6])+(1/sqrt(3))*sB*sigBEe;  
   
        //entro nelle 3 coppie di silici
for (Int_t p=1; p<7; p++)  {
                 if(p==1 || p==3 || p==5){
                 thetaX[1][p]= gRandom->Gaus(thetaX[1][p-1],sigSImu);
                 thetaY[1][p]= gRandom->Gaus(thetaY[1][p-1],sigSImu); 
                 thetaXe[1][p]= gRandom->Gaus(thetaXe[1][p-1],sigSIe);
                 thetaYe[1][p]= gRandom->Gaus(thetaYe[1][p-1],sigSIe); 

                x[1][p]=x[1][p-1]+d*tan(thetaX[1][p-1])+(1/sqrt(3))*sS*sigSImu+gRandom->Gaus(0,ris);
                y[1][p]=y[1][p-1]+d*tan(thetaY[1][p-1])+(1/sqrt(3))*sS*sigSImu;
                xe[1][p]=xe[1][p-1]+d*tan(thetaXe[1][p-1])+(1/sqrt(3))*sS*sigSIe+gRandom->Gaus(0,ris);
                ye[1][p]=ye[1][p-1]+d*tan(thetaYe[1][p-1])+(1/sqrt(3))*sS*sigSIe;
                
                thetaX[1][p+1]= gRandom->Gaus(thetaX[1][p],sigSImu);
                thetaY[1][p+1]= gRandom->Gaus(thetaY[1][p],sigSImu);
                thetaXe[1][p+1]= gRandom->Gaus(thetaXe[1][p],sigSIe);
                thetaYe[1][p+1]= gRandom->Gaus(thetaYe[1][p],sigSIe);

                x[1][p+1]=x[1][p]+dSS*tan(thetaX[1][p])+(1/sqrt(3))*sS*sigSImu;
                y[1][p+1]=y[1][p]+dSS*tan(thetaY[1][p])+(1/sqrt(3))*sS*sigSImu+gRandom->Gaus(0,ris);
                xe[1][p+1]=xe[1][p]+dSS*tan(thetaXe[1][p])+(1/sqrt(3))*sS*sigSIe;
                ye[1][p+1]=ye[1][p]+dSS*tan(thetaYe[1][p])+(1/sqrt(3))*sS*sigSIe+gRandom->Gaus(0,ris);                 
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
        coo_ang_fin[18][0]=xin; 
        coo_ang_fin[18][1]=yin;  
        
        
    return coo_ang_fin;
    
    }
     
    
    if(tar==1)
    {
        // siamo nelle stazioni con il target di berillio. Ora entra nel berillio ad una distanza d=0.25-0.005 m dagli ultimi silici. Qui però non puoi trascurare lo spessore del berillio, cioè dove interagisce? considero a metà (7.5 mm), quindi aggiungo a d il pezzo in cui x non è modificato e poi sommo con il nuovo angolo di scattering
                
                Double_t xin=(2-2*sSin-sB)*tan(anglexin)+2*(1/sqrt(3))*sSin*sigSIinP+(1/sqrt(3))*sSin*sigBEin+(1/sqrt(3))*0.0075*sigBE2in;
                Double_t yin=(2-2*sSin-sB)*tan(angleyin)+2*(1/sqrt(3))*sSin*sigSIinP+(1/sqrt(3))*sSin*sigBEin+(1/sqrt(3))*0.0075*sigBE2in;
        
               //  THinX[0]= gRandom->Gaus(anglexin,sigBEin);
            //     THinY[0]= gRandom->Gaus(angleyin,sigBEin); 
                 thetaX[0][0]=0;
                 thetaY[0][0]=0;
                 thetaXe[0][0]=0;
                 thetaYe[0][0]=0;

               //  inX[0]=xin+(1/sqrt(3))*sB*sigBEin;
            //     inY[0]=yin+(1/sqrt(3))*sB*sigBEin;
                 x[0][0]=0;
                 y[0][0]=0;  
                 xe[0][0]=0;
                 ye[0][0]=0;
       
//entro nelle 3 coppie di silici
for (Int_t p=1; p<7; p++)  {
                 if(p==1 || p==3 || p==5){
             //   THinX[p]=gRandom->Gaus(THinX[p-1],sigSIin);
            //  THinY[p]=gRandom->Gaus(THinY[p-1],sigSIin);
                 thetaX[0][p]=0;
                 thetaY[0][p]=0;
                 thetaXe[0][p]= 0;
                 thetaYe[0][p]= 0;

              //  inX[p]=inX[p-1]+d*tan(THinX[p-1])+(1/sqrt(3))*sS*sigSIin+gRandom->Gaus(0,ris);
              //  inY[p]=inY[p-1]+d*tan(THinY[p-1])+(1/sqrt(3))*sS*sigSIin;
                x[0][p]=0;
                y[0][p]=0;
                xe[0][p]=0;
                ye[0][p]=0;
                
                
              //  THinX[p+1]= gRandom->Gaus(THinX[p],sigSIin);
              //  THinY[p+1]= gRandom->Gaus(THinY[p],sigSIin);
                thetaX[0][p+1]=0;
                thetaY[0][p+1]=0;
                thetaXe[0][p+1]= 0;
                thetaYe[0][p+1]= 0;

               // inX[p+1]=inX[p]+dSS*tan(THinX[p])+(1/sqrt(3))*sS*sigSIin;
              //  inY[p+1]=inY[p]+dSS*tan(THinY[p])+(1/sqrt(3))*sS*sigSIin+gRandom->Gaus(0,ris);
                x[0][p+1]=0;
                y[0][p+1]=0;
                xe[0][p+1]=0;
                ye[0][p+1]=0; 
                 }}


                 x[0][0]=xin+(1/sqrt(3))*0.0075*sigBE2mu;
                 y[0][0]=yin+(1/sqrt(3))*0.0075*sigBE2mu;  
    
                 xe[0][0]=xin+(1/sqrt(3))*0.0075*sigBE2e;
                 ye[0][0]=yin+(1/sqrt(3))*0.0075*sigBE2e; 
        
        
             //   Double_t thetaX1= gRandom->Gaus(THinX[6],sigBE2in);
              //  Double_t thetaY1= gRandom->Gaus(THinY[6],sigBE2in); inX[6]+d*tan(THinX[6])
                
                thetaX[1][0]=gRandom->Gaus(anglex,sigBE2mu);
                thetaY[1][0]=gRandom->Gaus(angley,sigBE2mu);
                thetaXe[1][0]=gRandom->Gaus(anglexe,sigBE2e);
                thetaYe[1][0]=gRandom->Gaus(angleye,sigBE2e);

                 x[1][0]=xin+(1/sqrt(3))*0.0075*(sigBE2in)+(1/sqrt(3))*0.0075*(sigBE2mu); 
        
                 y[1][0]=yin+(1/sqrt(3))*0.0075*(sigBE2in)+(1/sqrt(3))*0.0075*(sigBE2mu);
        
                 xe[1][0]=xin+(1/sqrt(3))*0.0075*(sigBE2in)+(1/sqrt(3))*0.0075*(sigBE2e); 
        
                 ye[1][0]=yin+(1/sqrt(3))*0.0075*(sigBE2in)+(1/sqrt(3))*0.0075*(sigBE2e);
   
        //entro nelle 3 coppie di silici
for (Int_t p=1; p<7; p++)  {
                 if(p==1 || p==3 || p==5){
                 thetaX[1][p]= gRandom->Gaus(thetaX[1][p-1],sigSImu);
                 thetaY[1][p]= gRandom->Gaus(thetaY[1][p-1],sigSImu); 
                 thetaXe[1][p]= gRandom->Gaus(thetaXe[1][p-1],sigSIe);
                 thetaYe[1][p]= gRandom->Gaus(thetaYe[1][p-1],sigSIe); 

                x[1][p]=x[1][p-1]+d*tan(thetaX[1][p-1])+(1/sqrt(3))*sS*sigSImu+gRandom->Gaus(0,ris);
                y[1][p]=y[1][p-1]+d*tan(thetaY[1][p-1])+(1/sqrt(3))*sS*sigSImu;
                xe[1][p]=xe[1][p-1]+d*tan(thetaXe[1][p-1])+(1/sqrt(3))*sS*sigSIe+gRandom->Gaus(0,ris);
                ye[1][p]=ye[1][p-1]+d*tan(thetaYe[1][p-1])+(1/sqrt(3))*sS*sigSIe;
                
                thetaX[1][p+1]= gRandom->Gaus(thetaX[1][p],sigSImu);
                thetaY[1][p+1]= gRandom->Gaus(thetaY[1][p],sigSImu);
                thetaXe[1][p+1]= gRandom->Gaus(thetaXe[1][p],sigSIe);
                thetaYe[1][p+1]= gRandom->Gaus(thetaYe[1][p],sigSIe);

                x[1][p+1]=x[1][p]+dSS*tan(thetaX[1][p])+(1/sqrt(3))*sS*sigSImu;
                y[1][p+1]=y[1][p]+dSS*tan(thetaY[1][p])+(1/sqrt(3))*sS*sigSImu+gRandom->Gaus(0,ris);
                xe[1][p+1]=xe[1][p]+dSS*tan(thetaXe[1][p])+(1/sqrt(3))*sS*sigSIe;
                ye[1][p+1]=ye[1][p]+dSS*tan(thetaYe[1][p])+(1/sqrt(3))*sS*sigSIe+gRandom->Gaus(0,ris);                 
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
        coo_ang_fin[18][0]=xin; 
        coo_ang_fin[18][1]=yin;  
        
    return coo_ang_fin;
        
    }

else return coo_in; 

}



TMatrixD FastSim::MCSphoton(const Double_t & tar,const Double_t & theta,const Double_t & phi,const Double_t & xin,const Double_t & yin) const

{ Double_t d0 =2.025; //m
  Double_t d1 =1.025;
         TMatrixD coo (1,2);
        coo[0][0]=0;
        coo[0][1]=0;
 
    if(tar==0)
    {
    Double_t d_xy = d0*tan(theta);//vettore nel piano xy
    Double_t xf = xin+d_xy*cos(phi);
    Double_t yf = yin+d_xy*sin(phi);
        TMatrixD coo (1,2);
        coo[0][0]=xf;
        coo[0][1]=yf;

        return coo;
    }

    if(tar==1)
    {
    Double_t d_xy = d1*tan(theta);//vettore nel piano xy
    Double_t xf = xin+d_xy*cos(phi);
    Double_t yf = yin+d_xy*sin(phi);
        TMatrixD coo (1,2);
        coo[0][0]=xf;
        coo[0][1]=yf;

        return coo;
    }
 else return coo;
}


void FastSim::LoadKineVars(const PxPyPzEVector & p_mu_in,  const PxPyPzEVector & p_e_in, 
			   const PxPyPzEVector & p_mu_out, const PxPyPzEVector & p_e_out, const TMatrixD & coo,
			   MuE::KineVars & kv) {
  
  kv.Ee = p_e_out.E();
  kv.Emu = p_mu_out.E();
  kv.the = 1e3* p_e_out.Theta();
  kv.thmu = 1e3* p_mu_out.Theta();
  kv.phe = p_e_out.Phi();
  kv.phmu = p_mu_out.Phi();
    
/*TMatrixD coo_fin=coo(p_mu_out, p_e_out);
kv.cooXe = coo_fin[1][0];
kv.cooXmu = coo_fin[0][0];
kv.cooYe = coo_fin[1][1];
kv.cooYmu = coo_fin[0][1];

kv.tar = coo_fin[0][2];*/
 
kv.cooXmu = coo[0][0];
kv.cooYmu = coo[1][0];
kv.cooXe = coo[2][0];
kv.cooYe = coo[3][0];

kv.tar = coo[4][0];


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
      
    TMatrixD coo=MCSphoton(tar,photon.theta,photon.phi,xin,yin);
      
    photon.coox=coo[0][0];
    photon.cooy=coo[0][1];
      
    
  }

  else {
    photon.energyCoM = -1;
    photon.energy    = -1;
    photon.theta     = -1;
    photon.phi       =  0;
    photon.coox     = -1;
    photon.cooy     = -1;
      
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
