//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep  3 12:15:09 2020 by ROOT version 6.20/02
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
   Double_t        detKinBeamRot_pXmu;
   Double_t        detKinBeamRot_pYmu;
   Double_t        detKinBeamRot_pZmu;
   Double_t        detKinBeamRot_pXe;
   Double_t        detKinBeamRot_pYe;
   Double_t        detKinBeamRot_pZe;
   Double_t        detKinBeamRot_pXmu_out;
   Double_t        detKinBeamRot_pYmu_out;
   Double_t        detKinBeamRot_pZmu_out;
   Double_t        detKinBeamRot_pXe_out;
   Double_t        detKinBeamRot_pYe_out;
   Double_t        detKinBeamRot_pZe_out;
   Double_t        detKinBeamRot_Pmu_out;
   Double_t        detKinBeamRot_Pe_out;
   Double_t        detKinBeamRot_tar;
   Double_t        photon_energy;
   Double_t        photon_theta;
   Double_t        photon_phi;
   Double_t        photon_energyCoM;
   Double_t        photon_coox;
   Double_t        photon_cooy;

   // List of branches
   TBranch        *b_event_RunNr;   //!
   TBranch        *b_event_EventNr;   //!
   TBranch        *b_event_wgt_full;   //!
   TBranch        *b_event_wgt_norun;   //!
   TBranch        *b_event_wgt_lep;   //!
   TBranch        *b_event_wgt_LO;   //!
   TBranch        *b_event_E_mu_in;   //!
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
   TBranch        *b_event_detKinBeamRot_pXmu;   //!
   TBranch        *b_event_detKinBeamRot_pYmu;   //!
   TBranch        *b_event_detKinBeamRot_pZmu;   //!
   TBranch        *b_event_detKinBeamRot_pXe;   //!
   TBranch        *b_event_detKinBeamRot_pYe;   //!
   TBranch        *b_event_detKinBeamRot_pZe;   //!
   TBranch        *b_event_detKinBeamRot_pXmu_out;   //!
   TBranch        *b_event_detKinBeamRot_pYmu_out;   //!
   TBranch        *b_event_detKinBeamRot_pZmu_out;   //!
   TBranch        *b_event_detKinBeamRot_pXe_out;   //!
   TBranch        *b_event_detKinBeamRot_pYe_out;   //!
   TBranch        *b_event_detKinBeamRot_pZe_out;   //!
   TBranch        *b_event_detKinBeamRot_Pmu_out;   //!
   TBranch        *b_event_detKinBeamRot_Pe_out;   //!
   TBranch        *b_event_detKinBeamRot_tar;   //!
   TBranch        *b_event_photon_energy;   //!
   TBranch        *b_event_photon_theta;   //!
   TBranch        *b_event_photon_phi;   //!
   TBranch        *b_event_photon_energyCoM;   //!
   TBranch        *b_event_photon_coox;   //!
   TBranch        *b_event_photon_cooy;   //!

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
   fChain->SetBranchAddress("detKinBeamRot.pXmu", &detKinBeamRot_pXmu, &b_event_detKinBeamRot_pXmu);
   fChain->SetBranchAddress("detKinBeamRot.pYmu", &detKinBeamRot_pYmu, &b_event_detKinBeamRot_pYmu);
   fChain->SetBranchAddress("detKinBeamRot.pZmu", &detKinBeamRot_pZmu, &b_event_detKinBeamRot_pZmu);
   fChain->SetBranchAddress("detKinBeamRot.pXe", &detKinBeamRot_pXe, &b_event_detKinBeamRot_pXe);
   fChain->SetBranchAddress("detKinBeamRot.pYe", &detKinBeamRot_pYe, &b_event_detKinBeamRot_pYe);
   fChain->SetBranchAddress("detKinBeamRot.pZe", &detKinBeamRot_pZe, &b_event_detKinBeamRot_pZe);
   fChain->SetBranchAddress("detKinBeamRot.pXmu_out", &detKinBeamRot_pXmu_out, &b_event_detKinBeamRot_pXmu_out);
   fChain->SetBranchAddress("detKinBeamRot.pYmu_out", &detKinBeamRot_pYmu_out, &b_event_detKinBeamRot_pYmu_out);
   fChain->SetBranchAddress("detKinBeamRot.pZmu_out", &detKinBeamRot_pZmu_out, &b_event_detKinBeamRot_pZmu_out);
   fChain->SetBranchAddress("detKinBeamRot.pXe_out", &detKinBeamRot_pXe_out, &b_event_detKinBeamRot_pXe_out);
   fChain->SetBranchAddress("detKinBeamRot.pYe_out", &detKinBeamRot_pYe_out, &b_event_detKinBeamRot_pYe_out);
   fChain->SetBranchAddress("detKinBeamRot.pZe_out", &detKinBeamRot_pZe_out, &b_event_detKinBeamRot_pZe_out);
   fChain->SetBranchAddress("detKinBeamRot.Pmu_out", &detKinBeamRot_Pmu_out, &b_event_detKinBeamRot_Pmu_out);
   fChain->SetBranchAddress("detKinBeamRot.Pe_out", &detKinBeamRot_Pe_out, &b_event_detKinBeamRot_Pe_out);
   fChain->SetBranchAddress("detKinBeamRot.tar", &detKinBeamRot_tar, &b_event_detKinBeamRot_tar);
   fChain->SetBranchAddress("photon.energy", &photon_energy, &b_event_photon_energy);
   fChain->SetBranchAddress("photon.theta", &photon_theta, &b_event_photon_theta);
   fChain->SetBranchAddress("photon.phi", &photon_phi, &b_event_photon_phi);
   fChain->SetBranchAddress("photon.energyCoM", &photon_energyCoM, &b_event_photon_energyCoM);
   fChain->SetBranchAddress("photon.coox", &photon_coox, &b_event_photon_coox);
   fChain->SetBranchAddress("photon.cooy", &photon_cooy, &b_event_photon_cooy);
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
