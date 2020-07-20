#ifndef MuEana_H
#define MuEana_H

///////////////////////////////////////////////
// Classes defining MuE analysis variables
//
// G.Abbiendi  4/Dec/2018 
///////////////////////////////////////////////

namespace MuE {

  // Analysis variables in the output tree
  
  class KineVars {

  public:
    Double_t t13; // Mandelstam t (muon leg)
    Double_t t24; // Mandelstam t (electron leg)
    Double_t x13; // Feynman x (muon leg)
    Double_t x24; // Feynman x (electron leg)
    Double_t tt_e; // t computed from electron angle with LO formulas
    Double_t xt_e; // x computed from electron angle with LO formulas

    Double_t Ee; // electron energy
    Double_t Emu; // muon energy
    Double_t the; // electron theta (in mrad)
    Double_t thmu; // muon theta (in mrad)
    Double_t phe; // electron phi (from -pi to +pi)
    Double_t phmu; // muon phi (from -pi to +pi) 

    Double_t deltaPhi; // acoplanarity (deltaPhi)
    Double_t openingAngle; // opening angle mu-e out in the Lab
    Double_t tripleProduct; // triple product btw normalized vectors i . mu x e
      
      
      
      
Double_t pXmu;
Double_t pYmu;
Double_t pZmu;
Double_t pXmu_out;
Double_t pYmu_out;
Double_t pZmu_out;
Double_t pXe_out;
Double_t pYe_out;
Double_t pZe_out;
Double_t Pmu_out;
Double_t Pe_out;

    KineVars():
    t13(0),t24(0),x13(0),x24(0),tt_e(0),xt_e(0),Ee(0),Emu(0),the(0),thmu(0),phe(0),phmu(0),deltaPhi(0),openingAngle(0),tripleProduct(0),pXmu(0),pYmu(0),pZmu(0),pXmu_out(0),pYmu_out(0),pZmu_out(0),pXe_out(0),pYe_out(0),pZe_out(0),Pmu_out(0),Pe_out(0)
    {};

    virtual ~KineVars(){};

    ClassDef(KineVars,1)
  };

  class Photon {

  public:
    Double_t energy;    // photon energy in the Lab frame
    Double_t theta;     //   "    theta in the Lab frame (in mrad)
    Double_t phi;       //   "    phi in the Lab frame (in rad)
    Double_t energyCoM; // photon energy in the Centre-of-Mass frame
    
  Photon():
    energy(-1),theta(-1),phi(0),energyCoM(-1)
      {};
    
    virtual ~Photon(){};

    ClassDef(Photon,1)
  };

  class MuEana {

   public:
    UInt_t RunNr;
    Long64_t EventNr;
    Double_t wgt_full, wgt_norun, wgt_lep, wgt_LO;   // event weights 
    Double_t E_mu_in;  // incoming muon energy
    KineVars genKin;   // kinematic variables at Generator-level for e and mu tracks
    KineVars detKin;   // kinematic variables at Detector-level for e and mu tracks
    Photon photon;     // photon kinematic variables at Gen-level
    
    MuEana():
     RunNr(0),EventNr(0),wgt_full(0),wgt_norun(0),wgt_lep(0),wgt_LO(0),E_mu_in(0)
       {};
    
     virtual ~MuEana(){};   

     ClassDef(MuEana,1)
  };
}

#endif
