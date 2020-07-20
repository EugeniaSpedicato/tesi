#include <cmath>
#include <iostream>
#include "TRandom.h"
#include "TMath.h"
#include "ElasticState.h"
#include "ResolutionModels.h"
#include "FastSim.h"
#include "Math/Vector3D.h"

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
  
  LoadPhoton(event, photon);
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

// apply resolution smearing to particle momentum
//
PxPyPzEVector FastSim::Smear(const PxPyPzEVector & k) const
{
  
  Double_t thrms = ThetaRMS(k);
  
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
  
  Double_t smearx = gRandom->Gaus(0., thrms);
  gRandom->Gaus(0., 1.); // dummy call to preserve the random chain synchronization
  Double_t smeary = 0;
  
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

void FastSim::LoadKineVars(const PxPyPzEVector & p_mu_in,  const PxPyPzEVector & p_e_in, 
			   const PxPyPzEVector & p_mu_out, const PxPyPzEVector & p_e_out,
			   MuE::KineVars & kv) {
  
  kv.Ee = p_e_out.E();
  kv.Emu = p_mu_out.E();
  kv.the = 1e3* p_e_out.Theta();
  kv.thmu = 1e3* p_mu_out.Theta();
  kv.phe = p_e_out.Phi();
  kv.phmu = p_mu_out.Phi();
    
    
 kv.pXmu = p_mu_in.Px();
  kv.pYmu = p_mu_in.Py();
  kv.pZmu = p_mu_in.Pz();
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
  
  XYZVector crossProduct = nvec_mu_out.Cross(nvec_e_out);
  kv.openingAngle = 1000.*asin(crossProduct.R());
  kv.tripleProduct = nvec_mu_in.Dot(crossProduct);
}

void FastSim::LoadPhoton(const MuE::Event & event, MuE::Photon & photon) {
  // by now at most one photon
  auto n_photons = event.photons.size();
  
  if (n_photons >0) { 
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
