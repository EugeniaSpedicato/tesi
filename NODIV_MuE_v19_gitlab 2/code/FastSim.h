#ifndef FastSim_H
#define FastSim_H

///////////////////////////////////////////////
// Fast Simulation of MuE scattering
//
// G.Abbiendi  4/Sep/2018 
///////////////////////////////////////////////
#include "Math/Vector4D.h"
#include "MuEtree.h"
#include "MuEana.h"
#include "Inputs.h"

namespace MuE {

  class FastSim {

  public:
    FastSim(const MCpara & pargen, const FS_Input & fsi, bool _debug_=false);
    virtual ~FastSim(){};

    void Process(const Event & event);

    const KineVars & GetGenKin() const {return genKin;}
    const KineVars & GetDetKin() const {return detKin;}
    const Photon & GetPhoton() const {return photon;}

    void RandomNrSync();

  private:
    typedef ROOT::Math::PxPyPzEVector PxPyPzEVector;

    Double_t P_2bodies_CoM(Double_t Mass, Double_t mass_mu, Double_t mass_e) const;

    PxPyPzEVector Lorentz_ToCoM(const PxPyPzEVector & plab) const;
    PxPyPzEVector Lorentz_ToLab(const PxPyPzEVector & pcm) const;

    Double_t ThetaRMS(const PxPyPzEVector & p) const; 
    PxPyPzEVector Smear(const PxPyPzEVector & p) const; 
    PxPyPzEVector SmearX(const PxPyPzEVector & p) const; 
    PxPyPzEVector SmearPolar(const PxPyPzEVector & p) const; 

    void LoadKineVars(const PxPyPzEVector & p_mu_in,  const PxPyPzEVector & p_e_in, 
		      const PxPyPzEVector & p_mu_out, const PxPyPzEVector & p_e_out, 
		      KineVars & kv);
    void LoadPhoton(const Event & event, Photon & photon);

    static const Double_t mm_PDG; // PDG muon mass 
    static const Double_t me_PDG; // PDG electron mass

    const Double_t & mm; // muon mass
    const Double_t & me; // electron mass
    const Double_t & Ebeam; // average beam energy
    const Double_t & EbeamRMS; // beam energy spread

    const Int_t & model; // model for detector smearing (gaussian resolution)
    const Int_t & MSopt; // options for multiple scattering (default/only Xplane/only polar)
    const Double_t & thickness; // material thickness (in X0) for model_=0
    const Double_t & intrinsic_resolution; // intrinsic resolution (in mrad) for model_=0
    bool debug;

    PxPyPzEVector p_system; // mu-e centre-of-mass system fourmomentum
    Double_t Minv; // event invariant mass = sqrt(s)
    KineVars genKin; // kinematic variables at Gen-level for e and mu track
    KineVars detKin; // kinematic variables at Detector-level for e and mu track
    Photon photon; // photon variables at Gen-level

  };
}

#endif
