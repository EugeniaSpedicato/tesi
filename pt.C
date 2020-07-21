#define atree_cxx
#include "LO.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
void atree::Loop()
{
    TH1::SetDefaultSumw2();
    
    
TH1F* ptmu=new TH1F("a", "PT", 150,0,160);
TH1F* ptBRmu=new TH1F("a", "PT BR", 150,0,160);
    
TH1F* pte=new TH1F("ae", "PT", 150,0,160);
TH1F* ptBRe=new TH1F("ae", "PT BR", 150,0,160);


 if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
     Double_t PTmu=sqrt(detKin_pXmu_out*detKin_pXmu_out+detKin_pYmu_out*detKin_pYmu_out);
     Double_t PTe=sqrt(detKin_pXe_out*detKin_pXe_out+detKin_pYe_out*detKin_pYe_out);
       
     Double_t PTBRmu=sqrt(detKinBeamRot_pXmu_out*detKinBeamRot_pXmu_out+detKinBeamRot_pYmu_out*detKinBeamRot_pYmu_out);
     Double_t PTBRe=sqrt(detKinBeamRot_pXe_out*detKinBeamRot_pXe_out+detKinBeamRot_pYe_out*detKinBeamRot_pYe_out);
    
  
       ptmu->Fill(PTmu,wgt_full);
       ptBRmu->Fill(PTBRmu,wgt_full);
       
       pte->Fill(PTe,wgt_full);
       ptBRe->Fill(PTBRe,wgt_full);
       }
           
    TCanvas * ptq= new TCanvas("ptq","pqt",400,10,1500,1000);
    ptq->Divide(2,1);
    ptq->cd(1);
    ptmu->SetLineColor(40);
    ptmu->Draw("HIST");
    ptBRmu->SetLineColor(46);
    ptBRmu->Draw("HIST same");
    
    cd->(2);
    pte->SetLineColor(40);
    pte->Draw("HIST");
    ptBRe->SetLineColor(46);
    ptBRe->Draw("HIST same");
    
    ptq->SaveAs("pT.png");
    }