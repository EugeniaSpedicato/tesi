#define atree_cxx
#include "LO.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
void atree::Loop()
{
    TH1::SetDefaultSumw2();
    
    
TH1F* ptmu=new TH1F("a", "PT", 150,0,0.25);
TH1F* ptBRmu=new TH1F("a", "PT BR", 150,0,0.25);
    
TH1F* pte=new TH1F("ae", "PT", 150,0,0.25);
TH1F* ptBRe=new TH1F("ae", "PT BR", 150,0,0.25);
    
TH1F* thmu=new TH1F("h3bNj", "theta", 150,0,0.002);
TH1F* thBRmu=new TH1F("h3bNn", "theta BR", 150,0,0.002);
    
TH1F* the=new TH1F("h3bNj", "theta", 150,0,0.05);
TH1F* thBRe=new TH1F("h3bNn", "theta BR", 150,0,0.05);


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
       
        thmu->Fill(detKin_thmu*0.001,wgt_full);
       thBRmu->Fill(detKinBeamRot_thmu*0.001,wgt_full);
       
        the->Fill(detKin_the*0.001,wgt_full);
       thBRe->Fill(detKinBeamRot_the*0.001,wgt_full);
       }
           
    TCanvas * ptq= new TCanvas("ptq","pqt",400,10,1500,1000);
    ptq->Divide(2,2);
    ptq->cd(1);
    ptmu->SetLineColor(40);
    ptmu->Draw("HIST");
    
    ptq->cd(2);
    ptBRmu->SetLineColor(40);
    ptBRmu->Draw("HIST");
    
    ptq->cd(3);
    ptBRe->SetLineColor(46);
    ptBRe->Draw("HIST");
    
    ptq->cd(4);
    pte->SetLineColor(46);
    pte->Draw("HIST");
    
    ptq->SaveAs("pT.png");
    
    
    TCanvas * t= new TCanvas("t","t",400,10,1500,1000);
    t->Divide(2,2);
    t->cd(1);
    thmu->SetLineColor(40);
    thmu->Draw("HIST");
    
    t->cd(2);
    thBRmu->SetLineColor(40);
    thBRmu->Draw("HIST");
    
    t->cd(3);
    the->SetLineColor(46);
    the->Draw("HIST");
    t->cd(4);
    thBRe->SetLineColor(46);
    thBRe->Draw("HIST");
    
    t->SaveAs("TH.png");
    
    }