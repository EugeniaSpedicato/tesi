//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul  6 16:08:12 2020 by ROOT version 6.20/00
// from TTree atree/Analysis output tree
// found on file: outtree.root
//////////////////////////////////////////////////////////

#ifndef atree_h
#define atree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class atree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
 //MuE::MuEana     *event;
   UInt_t          RunNr;
   Long64_t        EventNr;
   Double_t        wgt_full;
   Double_t        wgt_norun;
   Double_t        wgt_lep;
   Double_t        wgt_LO;
   Double_t        E_mu_in;
   Double_t        genKin_t13;
   Double_t        genKin_t24;
   Double_t        genKin_x13;
   Double_t        genKin_x24;
   Double_t        genKin_tt_e;
   Double_t        genKin_xt_e;
   Double_t        genKin_Ee;
   Double_t        genKin_Emu;
   Double_t        genKin_the;
   Double_t        genKin_thmu;
   Double_t        genKin_phe;
   Double_t        genKin_phmu;
   Double_t        genKin_deltaPhi;
   Double_t        genKin_openingAngle;
   Double_t        genKin_tripleProduct;
   Double_t        genKin_cooXe;
   Double_t        genKin_cooXmu;
   Double_t        genKin_cooYe;
   Double_t        genKin_cooYmu;
   Double_t        detKin_t13;
   Double_t        detKin_t24;
   Double_t        detKin_x13;
   Double_t        detKin_x24;
   Double_t        detKin_tt_e;
   Double_t        detKin_xt_e;
   Double_t        detKin_Ee;
   Double_t        detKin_Emu;
   Double_t        detKin_the;
   Double_t        detKin_thmu;
   Double_t        detKin_phe;
   Double_t        detKin_phmu;
   Double_t        detKin_deltaPhi;
   Double_t        detKin_openingAngle;
   Double_t        detKin_tripleProduct;
   Double_t        detKin_cooXe;
   Double_t        detKin_cooXmu;
   Double_t        detKin_cooYe;
   Double_t        detKin_cooYmu;
   Double_t        detKinBeamRot_t13;
   Double_t        detKinBeamRot_t24;
   Double_t        detKinBeamRot_x13;
   Double_t        detKinBeamRot_x24;
   Double_t        detKinBeamRot_tt_e;
   Double_t        detKinBeamRot_xt_e;
   Double_t        detKinBeamRot_Ee;
   Double_t        detKinBeamRot_Emu;
   Double_t        detKinBeamRot_the;
   Double_t        detKinBeamRot_thmu;
   Double_t        detKinBeamRot_phe;
   Double_t        detKinBeamRot_phmu;
   Double_t        detKinBeamRot_deltaPhi;
   Double_t        detKinBeamRot_openingAngle;
   Double_t        detKinBeamRot_tripleProduct;
   Double_t        detKinBeamRot_cooXe;
   Double_t        detKinBeamRot_cooXmu;
   Double_t        detKinBeamRot_cooYe;
   Double_t        detKinBeamRot_cooYmu;
   Double_t        photon_energy;
   Double_t        photon_theta;
   Double_t        photon_phi;
   Double_t        photon_energyCoM;
   Double_t        photon_cooXph;
   Double_t        photon_cooYph;

   // List of branches
   TBranch        *b_event_RunNr;   //!
   TBranch        *b_event_EventNr;   //!
   TBranch        *b_event_wgt_full;   //!
   TBranch        *b_event_wgt_norun;   //!
   TBranch        *b_event_wgt_lep;   //!
   TBranch        *b_event_wgt_LO;   //!
   TBranch        *b_event_E_mu_in;   //!
   TBranch        *b_event_genKin_t13;   //!
   TBranch        *b_event_genKin_t24;   //!
   TBranch        *b_event_genKin_x13;   //!
   TBranch        *b_event_genKin_x24;   //!
   TBranch        *b_event_genKin_tt_e;   //!
   TBranch        *b_event_genKin_xt_e;   //!
   TBranch        *b_event_genKin_Ee;   //!
   TBranch        *b_event_genKin_Emu;   //!
   TBranch        *b_event_genKin_the;   //!
   TBranch        *b_event_genKin_thmu;   //!
   TBranch        *b_event_genKin_phe;   //!
   TBranch        *b_event_genKin_phmu;   //!
   TBranch        *b_event_genKin_deltaPhi;   //!
   TBranch        *b_event_genKin_openingAngle;   //!
   TBranch        *b_event_genKin_tripleProduct;   //!
   TBranch        *b_event_genKin_cooXe;   //!
   TBranch        *b_event_genKin_cooXmu;   //!
   TBranch        *b_event_genKin_cooYe;   //!
   TBranch        *b_event_genKin_cooYmu;   //!
   TBranch        *b_event_detKin_t13;   //!
   TBranch        *b_event_detKin_t24;   //!
   TBranch        *b_event_detKin_x13;   //!
   TBranch        *b_event_detKin_x24;   //!
   TBranch        *b_event_detKin_tt_e;   //!
   TBranch        *b_event_detKin_xt_e;   //!
   TBranch        *b_event_detKin_Ee;   //!
   TBranch        *b_event_detKin_Emu;   //!
   TBranch        *b_event_detKin_the;   //!
   TBranch        *b_event_detKin_thmu;   //!
   TBranch        *b_event_detKin_phe;   //!
   TBranch        *b_event_detKin_phmu;   //!
   TBranch        *b_event_detKin_deltaPhi;   //!
   TBranch        *b_event_detKin_openingAngle;   //!
   TBranch        *b_event_detKin_tripleProduct;   //!
   TBranch        *b_event_detKin_cooXe;   //!
   TBranch        *b_event_detKin_cooXmu;   //!
   TBranch        *b_event_detKin_cooYe;   //!
   TBranch        *b_event_detKin_cooYmu;   //!
   TBranch        *b_event_detKinBeamRot_t13;   //!
   TBranch        *b_event_detKinBeamRot_t24;   //!
   TBranch        *b_event_detKinBeamRot_x13;   //!
   TBranch        *b_event_detKinBeamRot_x24;   //!
   TBranch        *b_event_detKinBeamRot_tt_e;   //!
   TBranch        *b_event_detKinBeamRot_xt_e;   //!
   TBranch        *b_event_detKinBeamRot_Ee;   //!
   TBranch        *b_event_detKinBeamRot_Emu;   //!
   TBranch        *b_event_detKinBeamRot_the;   //!
   TBranch        *b_event_detKinBeamRot_thmu;   //!
   TBranch        *b_event_detKinBeamRot_phe;   //!
   TBranch        *b_event_detKinBeamRot_phmu;   //!
   TBranch        *b_event_detKinBeamRot_deltaPhi;   //!
   TBranch        *b_event_detKinBeamRot_openingAngle;   //!
   TBranch        *b_event_detKinBeamRot_tripleProduct;   //!
   TBranch        *b_event_detKinBeamRot_cooXe;   //!
   TBranch        *b_event_detKinBeamRot_cooXmu;   //!
   TBranch        *b_event_detKinBeamRot_cooYe;   //!
   TBranch        *b_event_detKinBeamRot_cooYmu;   //!
   TBranch        *b_event_photon_energy;   //!
   TBranch        *b_event_photon_theta;   //!
   TBranch        *b_event_photon_phi;   //!
   TBranch        *b_event_photon_energyCoM;   //!
   TBranch        *b_event_photon_cooXph;   //!
   TBranch        *b_event_photon_cooYph;   //!

   atree(TTree *tree=0);
   virtual ~atree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef atree_cxx
atree::atree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("outtree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("outtree.root");
      }
      f->GetObject("atree",tree);

   }
   Init(tree);
}

atree::~atree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t atree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t atree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void atree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNr", &RunNr, &b_event_RunNr);
   fChain->SetBranchAddress("EventNr", &EventNr, &b_event_EventNr);
   fChain->SetBranchAddress("wgt_full", &wgt_full, &b_event_wgt_full);
   fChain->SetBranchAddress("wgt_norun", &wgt_norun, &b_event_wgt_norun);
   fChain->SetBranchAddress("wgt_lep", &wgt_lep, &b_event_wgt_lep);
   fChain->SetBranchAddress("wgt_LO", &wgt_LO, &b_event_wgt_LO);
   fChain->SetBranchAddress("E_mu_in", &E_mu_in, &b_event_E_mu_in);
   fChain->SetBranchAddress("genKin.t13", &genKin_t13, &b_event_genKin_t13);
   fChain->SetBranchAddress("genKin.t24", &genKin_t24, &b_event_genKin_t24);
   fChain->SetBranchAddress("genKin.x13", &genKin_x13, &b_event_genKin_x13);
   fChain->SetBranchAddress("genKin.x24", &genKin_x24, &b_event_genKin_x24);
   fChain->SetBranchAddress("genKin.tt_e", &genKin_tt_e, &b_event_genKin_tt_e);
   fChain->SetBranchAddress("genKin.xt_e", &genKin_xt_e, &b_event_genKin_xt_e);
   fChain->SetBranchAddress("genKin.Ee", &genKin_Ee, &b_event_genKin_Ee);
   fChain->SetBranchAddress("genKin.Emu", &genKin_Emu, &b_event_genKin_Emu);
   fChain->SetBranchAddress("genKin.the", &genKin_the, &b_event_genKin_the);
   fChain->SetBranchAddress("genKin.thmu", &genKin_thmu, &b_event_genKin_thmu);
   fChain->SetBranchAddress("genKin.phe", &genKin_phe, &b_event_genKin_phe);
   fChain->SetBranchAddress("genKin.phmu", &genKin_phmu, &b_event_genKin_phmu);
   fChain->SetBranchAddress("genKin.deltaPhi", &genKin_deltaPhi, &b_event_genKin_deltaPhi);
   fChain->SetBranchAddress("genKin.openingAngle", &genKin_openingAngle, &b_event_genKin_openingAngle);
   fChain->SetBranchAddress("genKin.tripleProduct", &genKin_tripleProduct, &b_event_genKin_tripleProduct);
   fChain->SetBranchAddress("genKin.cooXe", &genKin_cooXe, &b_event_genKin_cooXe);
   fChain->SetBranchAddress("genKin.cooXmu", &genKin_cooXmu, &b_event_genKin_cooXmu);
   fChain->SetBranchAddress("genKin.cooYe", &genKin_cooYe, &b_event_genKin_cooYe);
   fChain->SetBranchAddress("genKin.cooYmu", &genKin_cooYmu, &b_event_genKin_cooYmu);
   fChain->SetBranchAddress("detKin.t13", &detKin_t13, &b_event_detKin_t13);
   fChain->SetBranchAddress("detKin.t24", &detKin_t24, &b_event_detKin_t24);
   fChain->SetBranchAddress("detKin.x13", &detKin_x13, &b_event_detKin_x13);
   fChain->SetBranchAddress("detKin.x24", &detKin_x24, &b_event_detKin_x24);
   fChain->SetBranchAddress("detKin.tt_e", &detKin_tt_e, &b_event_detKin_tt_e);
   fChain->SetBranchAddress("detKin.xt_e", &detKin_xt_e, &b_event_detKin_xt_e);
   fChain->SetBranchAddress("detKin.Ee", &detKin_Ee, &b_event_detKin_Ee);
   fChain->SetBranchAddress("detKin.Emu", &detKin_Emu, &b_event_detKin_Emu);
   fChain->SetBranchAddress("detKin.the", &detKin_the, &b_event_detKin_the);
   fChain->SetBranchAddress("detKin.thmu", &detKin_thmu, &b_event_detKin_thmu);
   fChain->SetBranchAddress("detKin.phe", &detKin_phe, &b_event_detKin_phe);
   fChain->SetBranchAddress("detKin.phmu", &detKin_phmu, &b_event_detKin_phmu);
   fChain->SetBranchAddress("detKin.deltaPhi", &detKin_deltaPhi, &b_event_detKin_deltaPhi);
   fChain->SetBranchAddress("detKin.openingAngle", &detKin_openingAngle, &b_event_detKin_openingAngle);
   fChain->SetBranchAddress("detKin.tripleProduct", &detKin_tripleProduct, &b_event_detKin_tripleProduct);
   fChain->SetBranchAddress("detKin.cooXe", &detKin_cooXe, &b_event_detKin_cooXe);
   fChain->SetBranchAddress("detKin.cooXmu", &detKin_cooXmu, &b_event_detKin_cooXmu);
   fChain->SetBranchAddress("detKin.cooYe", &detKin_cooYe, &b_event_detKin_cooYe);
   fChain->SetBranchAddress("detKin.cooYmu", &detKin_cooYmu, &b_event_detKin_cooYmu);
   fChain->SetBranchAddress("detKinBeamRot.t13", &detKinBeamRot_t13, &b_event_detKinBeamRot_t13);
   fChain->SetBranchAddress("detKinBeamRot.t24", &detKinBeamRot_t24, &b_event_detKinBeamRot_t24);
   fChain->SetBranchAddress("detKinBeamRot.x13", &detKinBeamRot_x13, &b_event_detKinBeamRot_x13);
   fChain->SetBranchAddress("detKinBeamRot.x24", &detKinBeamRot_x24, &b_event_detKinBeamRot_x24);
   fChain->SetBranchAddress("detKinBeamRot.tt_e", &detKinBeamRot_tt_e, &b_event_detKinBeamRot_tt_e);
   fChain->SetBranchAddress("detKinBeamRot.xt_e", &detKinBeamRot_xt_e, &b_event_detKinBeamRot_xt_e);
   fChain->SetBranchAddress("detKinBeamRot.Ee", &detKinBeamRot_Ee, &b_event_detKinBeamRot_Ee);
   fChain->SetBranchAddress("detKinBeamRot.Emu", &detKinBeamRot_Emu, &b_event_detKinBeamRot_Emu);
   fChain->SetBranchAddress("detKinBeamRot.the", &detKinBeamRot_the, &b_event_detKinBeamRot_the);
   fChain->SetBranchAddress("detKinBeamRot.thmu", &detKinBeamRot_thmu, &b_event_detKinBeamRot_thmu);
   fChain->SetBranchAddress("detKinBeamRot.phe", &detKinBeamRot_phe, &b_event_detKinBeamRot_phe);
   fChain->SetBranchAddress("detKinBeamRot.phmu", &detKinBeamRot_phmu, &b_event_detKinBeamRot_phmu);
   fChain->SetBranchAddress("detKinBeamRot.deltaPhi", &detKinBeamRot_deltaPhi, &b_event_detKinBeamRot_deltaPhi);
   fChain->SetBranchAddress("detKinBeamRot.openingAngle", &detKinBeamRot_openingAngle, &b_event_detKinBeamRot_openingAngle);
   fChain->SetBranchAddress("detKinBeamRot.tripleProduct", &detKinBeamRot_tripleProduct, &b_event_detKinBeamRot_tripleProduct);
   fChain->SetBranchAddress("detKinBeamRot.cooXe", &detKinBeamRot_cooXe, &b_event_detKinBeamRot_cooXe);
   fChain->SetBranchAddress("detKinBeamRot.cooXmu", &detKinBeamRot_cooXmu, &b_event_detKinBeamRot_cooXmu);
   fChain->SetBranchAddress("detKinBeamRot.cooYe", &detKinBeamRot_cooYe, &b_event_detKinBeamRot_cooYe);
   fChain->SetBranchAddress("detKinBeamRot.cooYmu", &detKinBeamRot_cooYmu, &b_event_detKinBeamRot_cooYmu);
   fChain->SetBranchAddress("photon.energy", &photon_energy, &b_event_photon_energy);
   fChain->SetBranchAddress("photon.theta", &photon_theta, &b_event_photon_theta);
   fChain->SetBranchAddress("photon.phi", &photon_phi, &b_event_photon_phi);
   fChain->SetBranchAddress("photon.energyCoM", &photon_energyCoM, &b_event_photon_energyCoM);
   fChain->SetBranchAddress("photon.cooXph", &photon_cooXph, &b_event_photon_cooXph);
   fChain->SetBranchAddress("photon.cooYph", &photon_cooYph, &b_event_photon_cooYph);
   Notify();
}

Bool_t atree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void atree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t atree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef atree_cxx
