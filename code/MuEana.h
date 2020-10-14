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
     Double_t cooXe;
      Double_t cooXmu;
      Double_t cooYe;
      Double_t cooYmu;
   Double_t pXmu;
Double_t pYmu;
Double_t pZmu;
Double_t pXe;
Double_t pYe;
Double_t pZe;
Double_t pXmu_out;
Double_t pYmu_out;
Double_t pZmu_out;
Double_t pXe_out;
Double_t pYe_out;
Double_t pZe_out;
Double_t Pmu_out;
Double_t Pe_out;
Double_t tar;
Double_t ThEl_interaction;
Double_t def_angle_mu;
Double_t def_angle_e; 
Double_t n_cell_e;
Double_t E1_e;
Double_t E2_e;
Double_t E3_e;
Double_t E4_e;
Double_t E5_e;
Double_t E6_e;
Double_t E7_e;
Double_t E8_e;
Double_t E9_e;
Double_t E10_e;
Double_t E11_e;
Double_t E12_e;
Double_t E13_e;
Double_t E14_e;
Double_t E15_e;
Double_t E16_e;
Double_t E17_e;
Double_t E18_e;
Double_t E19_e;
Double_t E20_e;
Double_t E21_e;
Double_t E22_e;
Double_t E23_e;
Double_t E24_e;
Double_t E25_e;
      

    KineVars():
    t13(0),t24(0),x13(0),x24(0),tt_e(0),xt_e(0),Ee(0),Emu(0),the(0),thmu(0),phe(0),phmu(0),deltaPhi(0),openingAngle(0),tripleProduct(0),cooXe(0),cooXmu(0),cooYe(0),cooYmu(0),pXmu(0),pYmu(0),pZmu(0),pXe(0),pYe(0),pZe(0),pXmu_out(0),pYmu_out(0),pZmu_out(0),pXe_out(0),pYe_out(0),pZe_out(0),Pmu_out(0),Pe_out(0),tar(0),ThEl_interaction(0),def_angle_mu(0),def_angle_e(0),n_cell_e(0),E1_e(0),E2_e(0),E3_e(0),E4_e(0),E5_e(0),E6_e(0),E7_e(0),E8_e(0),E9_e(0),E10_e(0),E11_e(0),E12_e(0),E13_e(0),E14_e(0),E15_e(0),E16_e(0),E17_e(0),E18_e(0),E19_e(0),E20_e(0),E21_e(0),E22_e(0),E23_e(0),E24_e(0),E25_e(0)
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
   Double_t coox;
   Double_t cooy;
Double_t n_cell_ph;
Double_t E1_ph;
Double_t E2_ph;
Double_t E3_ph;
Double_t E4_ph;
Double_t E5_ph;
Double_t E6_ph;
Double_t E7_ph;
Double_t E8_ph;
Double_t E9_ph;
Double_t E10_ph;
Double_t E11_ph;
Double_t E12_ph;
Double_t E13_ph;
Double_t E14_ph;
Double_t E15_ph;
Double_t E16_ph;
Double_t E17_ph;
Double_t E18_ph;
Double_t E19_ph;
Double_t E20_ph;
Double_t E21_ph;
Double_t E22_ph;
Double_t E23_ph;
Double_t E24_ph;
Double_t E25_ph;
    
  Photon():
    energy(-1),theta(-1),phi(0),energyCoM(-1),coox(-1),cooy(-1),n_cell_ph(0),E1_ph(0),E2_ph(0),E3_ph(0),E4_ph(0),E5_ph(0),E6_ph(0),E7_ph(0),E8_ph(0),E9_ph(0),E10_ph(0),E11_ph(0),E12_ph(0),E13_ph(0),E14_ph(0),E15_ph(0),E16_ph(0),E17_ph(0),E18_ph(0),E19_ph(0),E20_ph(0),E21_ph(0),E22_ph(0),E23_ph(0),E24_ph(0),E25_ph(0)
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
   // KineVars genKin;   // kinematic variables at Generator-level for e and mu tracks
   // KineVars detKin; // kinematic variables at Detector-level for e and mu tracks
    KineVars detKinBeamRot; 
     Photon photon;     // photon kinematic variables at Gen-level
    
    MuEana():
     RunNr(0),EventNr(0),wgt_full(0),wgt_norun(0),wgt_lep(0),wgt_LO(0),E_mu_in(0)
       {};
    
     virtual ~MuEana(){};   

     ClassDef(MuEana,1)
  };
}
  

#endif
