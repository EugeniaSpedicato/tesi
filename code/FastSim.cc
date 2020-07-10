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
   
    LoadKineVars(p_mu_in, p_e_in, p_mu_out, p_e_out, genKin); 
     
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
  
    
    //DETKIN RUOTATE

PxPyPzEVector p_mu_in_div=RotDivIN(p_mu_in);
PxPyPzEVector p_e_in_div=RotDivIN(p_e_in);
    
PxPyPzEVector p_mu_out_div=RotDiv(p_mu_in_div,p_mu_out);
PxPyPzEVector p_e_out_div=RotDiv(p_mu_in_div,p_e_out);
  
PxPyPzEVector p_mu_out_div_smeared, p_e_out_div_smeared; 
  
    if (MSopt ==0) {
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
  else cout<<"\n"<<"*** ERROR : FastSim, undefined detector MS option = "<<MSopt<<endl;
  
  LoadKineVars(p_mu_in_div, p_e_in_div, p_mu_out_div_smeared, p_e_out_div_smeared, detKinBeamRot);
  
  LoadPhoton(event, photon, p_mu_in_div);
    
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





 PxPyPzEVector FastSim::RotDivIN(const PxPyPzEVector & k) const
 {
Double_t divthx = gRandom->Gaus(0., 0.00027);
Double_t divthy = gRandom->Gaus(0., 0.00020); 

Double_t anglex = atan2(k.Px(), k.Pz());
Double_t angley = atan2(k.Py(), k.Pz()); 
     
anglex += divthx;
angley += divthy;
     
    //NB questi Px Py Pz sono del nuovo!!
Double_t pmuin=sqrt(k.Px()*k.Px()+k.Py()*k.Py()+k.Pz()*k.Pz());

Double_t pz=pmuin/(1+tan(anglex)*tan(anglex)+tan(angley)*tan(angley));
Double_t py=pz*tan(angley);
Double_t px=pz*tan(anglex);

Double_t ptz=sqrt(px*px+pz*pz);

PxPyPzEVector pnewdiv(px, py, pz, k.E());     
return pnewdiv;  
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
    R[1][0]=sin(phi)*sin(psi);
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

PxPyPzEVector pnewdiv(pN[0][0], pN[1][0], pN[2][0], k.E());     
return pnewdiv;    
 }   











// apply resolution smearing to particle momentum
//
PxPyPzEVector FastSim::Smear(const PxPyPzEVector & k) const
{
  
  Double_t thrms = ThetaRMS(k);
  //HO AGGIUNTO QUESTO
  /*Double_t divx = gRandom->Gaus(0., 0.3);
  Double_t divy = gRandom->Gaus(0., 0.2);
                                                               
  Double_t smearx = gRandom->Gaus(divx, thrms);
  Double_t smeary = gRandom->Gaus(divy, thrms);*/
    
  Double_t smearx = gRandom->Gaus(0., thrms);
  Double_t smeary = gRandom->Gaus(0., thrms);    
  
  // angles in the xz and yz planes // defined in -pi, +pi, although should be small angles around zero
  Double_t anglex = atan2(k.Px(), k.Pz());
  Double_t angley = atan2(k.Py(), k.Pz()); // small-angle approx ??
  
  // apply smearing
  anglex += 0.001*smearx;
  angley += 0.001*smeary;
  
  Double_t dxdz = tan(anglex); // could approx tan ~ angle
  Double_t dydz = tan(angley);
  
  // assuming z-motion is always forward
  Double_t skz = sqrt(k.P2() / (dxdz*dxdz + dydz*dydz + 1));
  Double_t skx = skz * dxdz;
  Double_t sky = skz * dydz;
  
  PxPyPzEVector psmeared(skx, sky, skz, k.E());
  
  return psmeared;
}

// apply resolution smearing to particle momentum (SMEAR only on the XZ plane)
//
PxPyPzEVector FastSim::SmearX(const PxPyPzEVector & k) const
{
  Double_t thrms = ThetaRMS(k);
  //HO AGGIUNTO QUESTO
 /* Double_t divx = gRandom->Gaus(0., 0.3);
  Double_t divy = gRandom->Gaus(0., 0.2);
    
  Double_t smearx = gRandom->Gaus(divx, thrms);
  */
  Double_t smeary = gRandom->Gaus(0, thrms);
 Double_t smearx = gRandom->Gaus(0, thrms);  
  gRandom->Gaus(0., 1.); // dummy call to preserve the random chain synchronization
  
  // angles in the xz and yz planes // defined in -pi, +pi, although should be small angles around zero
  Double_t anglex = atan2(k.Px(), k.Pz());
  Double_t angley = atan2(k.Py(), k.Pz()); // small-angle approx ??
  
  // apply smearing
  anglex += 0.001*smearx;
  angley += 0.001*smeary;
  
  Double_t dxdz = tan(anglex); // could approx tan ~ angle
  Double_t dydz = tan(angley);
  
  // assuming z-motion is always forward
  Double_t skz = sqrt(k.P2() / (dxdz*dxdz + dydz*dydz + 1));
  Double_t skx = skz * dxdz;
  Double_t sky = skz * dydz;
  
  PxPyPzEVector psmeared(skx, sky, skz, k.E());
  
  return psmeared;
}

// apply resolution smearing to particle momentum
// (theta smearing of rms*sqrt(2) in the plane defined by the vector and the z-axis)
//
PxPyPzEVector FastSim::SmearPolar(const PxPyPzEVector & k) const
{
  Polar3DVector p(k.P(), k.Theta(), k.Phi()); 

  // assumed smearing is sqrt(2) * smearing on a plane
  Double_t thrms = ThetaRMS(k);  
  Double_t smearth = sqrt(2) * gRandom->Gaus(0., thrms);
  gRandom->Gaus(0., 1.); // dummy call to preserve the random chain synchronization
  Double_t thetasm = p.Theta()+0.001*smearth;
  Double_t phism = p.Phi();
  // note that theta is defined positive (0-pi)
  // going negative means changing the azimuth too and phi is defined in (-pi,+pi)
  if (thetasm < 0) {
    thetasm = -thetasm;
    phism = phism > 0 ? phism - TMath::Pi() : phism + TMath::Pi();
  }

  p.SetTheta(thetasm);
  p.SetPhi(phism);

  return PxPyPzEVector(p.X(), p.Y(), p.Z(), k.E());
}


XYZVector FastSim::coo(const Double_t & the, const Double_t & phi) const
{   Double_t theR = the*0.001;//rad
    Double_t phiR = phi;//rad
    Double_t d0=2.10;//m
    Double_t d1=1.10;//m
 
    Double_t x=gRandom->Gaus(0., 0.026);//m
    Double_t y=gRandom->Gaus(0., 0.027);//m
    Double_t z=0.;
    Double_t zf=2.10;
    XYZVector coo_in(x,y,z);
 
    //interazione target 1 o target 2 
    Int_t tar=gRandom->Integer(2);
    
    if (tar==0)
    {
    //generato random
    Double_t d_xy = d0*tan(theR);//vettore nel piano xy
    Double_t xf = x+d_xy*cos(phiR);
    Double_t yf = y+d_xy*sin(phiR);
       
    XYZVector coo_f(xf,yf,zf);
         return coo_f;
    }
    
    if(tar==1)
    {
    Double_t d_xy = d1*tan(theR);//vettore nel piano xy
    Double_t xf = x+d_xy*cos(phiR);
    Double_t yf = y+d_xy*sin(phiR);
        
    XYZVector coo_f(xf,yf,zf);
     return coo_f;}
       // cout << "Second target ";
        
else return coo_in;
}


void FastSim::LoadKineVars(const PxPyPzEVector & p_mu_in,  const PxPyPzEVector & p_e_in, 
			   const PxPyPzEVector & p_mu_out, const PxPyPzEVector & p_e_out,
			   MuE::KineVars & kv) {
  
  kv.Ee = p_e_out.E();
  kv.Emu = p_mu_out.E();
  kv.the = 1e3* p_e_out.Theta();
  kv.thmu = 1e3* p_mu_out.Theta();
  kv.phe = p_e_out.Phi();
  kv.phmu = p_mu_out.Phi();
/* kv.pXmu = p_mu_in.Px();
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
*/


XYZVector coo_fin_mu=coo(kv.thmu,kv.phmu);
XYZVector coo_fin_e=coo(kv.the,kv.phe);
kv.cooXe = coo_fin_e.X();
kv.cooXmu = coo_fin_mu.X();
kv.cooYe = coo_fin_e.Y();
kv.cooYmu = coo_fin_mu.Y();

    
    
    
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

void FastSim::LoadPhoton(const MuE::Event & event, MuE::Photon & photon,const PxPyPzEVector & p_mu_in_div) {
  // by now at most one photon
  auto n_photons = event.photons.size();
  
  if (n_photons >0) {  /*
    PxPyPzEVector p_gamma_Lab = {event.photons[0].px, 
				 event.photons[0].py,
				 event.photons[0].pz,
			         event.photons[0].E};
    PxPyPzEVector p_gamma_CoM = Lorentz_ToCoM(p_gamma_Lab);
    
    photon.energy    = p_gamma_Lab.E();
    photon.theta     = p_gamma_Lab.Theta() *1e3;
    photon.phi       = p_gamma_Lab.Phi();
    photon.energyCoM = p_gamma_CoM.E();  
  }

  else {
    photon.energyCoM = -1;
    photon.energy    = -1;
    photon.theta     = -1;
    photon.phi       =  0;
  }
     */
      
      PxPyPzEVector p_ph(event.photons[0].px,event.photons[0].py,event.photons[0].pz,event.photons[0].E);
      PxPyPzEVector p_ph_div=RotDiv(p_mu_in_div,p_ph);   
      
     
    PxPyPzEVector p_gamma_Lab = {p_ph_div.Px(), 
                                 p_ph_div.Py(),
                                 p_ph_div.Pz(),
			                     event.photons[0].E};
      
    PxPyPzEVector p_gamma_CoM = Lorentz_ToCoM(p_gamma_Lab);
    
    photon.energy    = p_gamma_Lab.E();
    photon.theta     = p_gamma_Lab.Theta() *1e3;
    photon.phi       = p_gamma_Lab.Phi();
    photon.energyCoM = p_gamma_CoM.E(); 
      
    XYZVector coo_fin_ph=coo(photon.theta,photon.phi);

      photon.cooXph = coo_fin_ph.X();
      photon.cooYph = coo_fin_ph.Y();
    
  }

  else {
    photon.energyCoM = -1;
    photon.energy    = -1;
    photon.theta     = -1;
    photon.phi       =  0;
    photon.cooXph    = 0;
    photon.cooYph    = 0;
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
