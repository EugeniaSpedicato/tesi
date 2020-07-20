#define atree_cxx
#include "LO.h"
#include <TH2.h>
#include <TH1.h>

#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{
    TH1::SetDefaultSumw2();
    
Double_t d0=2.10; //meters
Double_t d1=1.10; //meters 
    
           
TH1F* Npx_mu_out=new TH1F("h1a", "pX_out muon BR", 150,-0.3,0.3);
TH1F* Npy_mu_out=new TH1F("h2a", "pY_out muon BR", 150,-0.3,0.3);
TH1F* Npz_mu_out=new TH1F("h3a", "pZ_out muon BR", 150,0,180);
TH1F* Npx_mu_outNO=new TH1F("h1aN", "pX_out muon BR", 150,-0.3,0.3);
TH1F* Npy_mu_outNO=new TH1F("h2aN", "pY_out muon  BR", 150,-0.3,0.3);
TH1F* Npz_mu_outNO=new TH1F("h3aN", "pZ_out muon N BR", 150,0,180);
    
TH1F* Npx_e_out=new TH1F("h1b", "pX_out electron BR", 150,-0.3,0.3);
TH1F* Npy_e_out=new TH1F("h2b", "pY_out electron BR", 150,-0.3,0.3);
TH1F* Npz_e_out=new TH1F("h3b", "pZ_out electron BR", 150,0,5);
TH1F* Npx_e_outNO=new TH1F("h1bN", "pX_out electron  BR", 150,-0.3,0.3);
TH1F* Npy_e_outNO=new TH1F("h2bN", "pY_out electron  BR", 150,-0.3,0.3);
TH1F* Npz_e_outNO=new TH1F("h3bN", "pZ_out electron  BR", 150,0,5);
    
TH1F* px_mu_out=new TH1F("h1a", "pX_out muon NO DIV", 150,-0.3,0.3);
TH1F* py_mu_out=new TH1F("h2a", "pY_out muon NO DIV", 150,-0.3,0.3);
TH1F* pz_mu_out=new TH1F("h3a", "pZ_out muon NO DIV", 150,0,180);
TH1F* px_mu_outNO=new TH1F("h1aN", "pX_out muon NO DIV", 150,-0.3,0.3);
TH1F* py_mu_outNO=new TH1F("h2aN", "pY_out muon NO DIV", 150,-0.3,0.3);
TH1F* pz_mu_outNO=new TH1F("h3aN", "pZ_out muon NO DIV", 150,0,180);
    
TH1F* px_e_out=new TH1F("h1b", "pX_out electron NO DIV", 150,-0.3,0.3);
TH1F* py_e_out=new TH1F("h2b", "pY_out electron NO DIV", 150,-0.3,0.3);
TH1F* pz_e_out=new TH1F("h3b", "pZ_out electron NO DIV", 150,0,5);
TH1F* px_e_outNO=new TH1F("h1bN", "pX_out electron NO DIV", 150,-0.3,0.3);
TH1F* py_e_outNO=new TH1F("h2bN", "pY_out electron NO DIV", 150,-0.3,0.3);
TH1F* pz_e_outNO=new TH1F("h3bN", "pZ_out electron NO DIV", 150,0,5);
    
TH1F* th=new TH1F("h3bN", "theta", 150,-0.1,0.1);
TH1F* thBR=new TH1F("h3bN", "theta BR", 150,-0.1,0.1);
    
TH1F* thXZe_one=new TH1F("h1b", "theta electron 1", 150,-0.3,0.3);
TH1F* thYZe_one=new TH1F("h2b", "theta electron 1", 150,-0.3,0.3);
    
TH1F* thXZmu_one=new TH1F("h1b", "theta muon 1", 150,-0.3,0.3);
TH1F* thYZmu_one=new TH1F("h2b", "theta muon 1", 150,-0.3,0.3);
    
TH1F* thXZe_two=new TH1F("h1b", "theta electron 2", 150,-0.3,0.3);
TH1F* thYZe_two=new TH1F("h2b", "theta electron 2", 150,-0.3,0.3);
    
TH1F* thXZmu_two=new TH1F("h1b", "theta muon 2", 150,-0.3,0.3);
TH1F* thYZmu_two=new TH1F("h2b", "theta muon 2", 150,-0.3,0.3);
    
    
      
    if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
       
       
       Npx_mu_out->Fill(detKinBeamRot_pXmu_out,wgt_LO);
       Npy_mu_out->Fill(detKinBeamRot_pYmu_out,wgt_LO);
       Npz_mu_out->Fill(detKinBeamRot_pZmu_out,wgt_LO);
       Npx_mu_outNO->Fill(detKinBeamRot_pXmu_out,wgt_full);
       Npy_mu_outNO->Fill(detKinBeamRot_pYmu_out,wgt_full);
       Npz_mu_outNO->Fill(detKinBeamRot_pZmu_out,wgt_full);
       
       Npx_e_out->Fill(detKinBeamRot_pXe_out,wgt_LO);
       Npy_e_out->Fill(detKinBeamRot_pYe_out,wgt_LO);
       Npz_e_out->Fill(detKinBeamRot_pZe_out,wgt_LO);
       Npx_e_outNO->Fill(detKinBeamRot_pXe_out,wgt_full);
       Npy_e_outNO->Fill(detKinBeamRot_pYe_out,wgt_full);
       Npz_e_outNO->Fill(detKinBeamRot_pZe_out,wgt_full);
       
        px_mu_out->Fill(detKin_pXmu_out,wgt_LO);
       py_mu_out->Fill(detKin_pYmu_out,wgt_LO);
       pz_mu_out->Fill(detKin_pZmu_out,wgt_LO);
       px_mu_outNO->Fill(detKin_pXmu_out,wgt_full);
       py_mu_outNO->Fill(detKin_pYmu_out,wgt_full);
       pz_mu_outNO->Fill(detKin_pZmu_out,wgt_full);
       
       px_e_out->Fill(detKin_pXe_out,wgt_LO);
       py_e_out->Fill(detKin_pYe_out,wgt_LO);
       pz_e_out->Fill(detKin_pZe_out,wgt_LO);
       px_e_outNO->Fill(detKin_pXe_out,wgt_full);
       py_e_outNO->Fill(detKin_pYe_out,wgt_full);
       pz_e_outNO->Fill(detKin_pZe_out,wgt_full);
       
       th->Fill(detKin_thmu,wgt_full);
        thBR->Fill(detKinBeamRot_thmu,wgt_full);
       
            if (detKinBeamRot_tar==0)
       {
        Double_t th_xze=atan2(detKinBeamRot_cooXe,d0);
        Double_t th_xzmu=atan2(detKinBeamRot_cooXmu,d0);
        Double_t th_yze=atan2(detKinBeamRot_cooYe,d0);
        Double_t th_yzmu=atan2(detKinBeamRot_cooYmu,d0);
                
             thXZe_one->Fill(th_xze,wgt_full);
             thXZmu_one->Fill(th_xzmu,wgt_full);
             thYZe_one->Fill(th_yze,wgt_full);
             thYZmu_one->Fill(th_yzmu,wgt_full);
        
       }
       
       if (detKinBeamRot_tar==1)
       {

        Double_t th_xze=atan2(detKinBeamRot_cooXe,d1);
        Double_t th_xzmu=atan2(detKinBeamRot_cooXmu,d1);
        Double_t th_yze=atan2(detKinBeamRot_cooYe,d1);
        Double_t th_yzmu=atan2(detKinBeamRot_cooYmu,d1);  
           
             thXZe_two->Fill(th_xze,wgt_full);
             thXZmu_two->Fill(th_xzmu,wgt_full);
             thYZe_two->Fill(th_yze,wgt_full);
             thYZmu_two->Fill(th_yzmu,wgt_full);
       }
       
       
   }
    
    TCanvas * diffP= new TCanvas("diff","diff",400,10,1500,1000);
    diffP->Divide(3,2);
    diffP->cd(1);
    Npx_e_outNO->SetLineColor(kRed);
    Npx_e_outNO->Draw("HIST");
    Npx_e_out->SetLineColor(kBlack);
    Npx_e_out->Draw("HIST same");

    
    
    
    diffP->cd(2);
    Npy_e_outNO->SetLineColor(kRed);
    Npy_e_outNO->Draw("HIST");
    Npy_e_out->SetLineColor(kBlack);
    Npy_e_out->Draw("HIST same");

    
    
    
    
    diffP->cd(3);
    Npz_e_outNO->SetLineColor(kRed);
    Npz_e_outNO->Draw("HIST");
    Npz_e_out->SetLineColor(kBlack);
    Npz_e_out->Draw("HIST same");

    
    
    
    diffP->cd(4);
    Npx_mu_outNO->SetLineColor(40);
    Npx_mu_outNO->Draw("HIST");
    Npx_mu_out->SetLineColor(46);
    Npx_mu_out->Draw("HIST same");

    
    
    diffP->cd(5);
    Npy_mu_outNO->SetLineColor(40);
    Npy_mu_outNO->Draw("HIST");
    Npy_mu_out->SetLineColor(46);
    Npy_mu_out->Draw("HIST same");

    
    
    diffP->cd(6);
    Npz_mu_outNO->SetLineColor(40);
    Npz_mu_outNO->Draw("HIST");
    Npz_mu_out->SetLineColor(46);
    Npz_mu_out->Draw("HIST same");

    
    
    diffP->SaveAs("LO-NLO.png");
    
    
        TCanvas * ndiffP= new TCanvas("diff","diff",400,10,1500,1000);
   ndiffP->Divide(3,2);
    ndiffP->cd(1);
    px_e_outNO->SetLineColor(kRed);
    px_e_outNO->Draw("HIST");
    px_e_out->SetLineColor(kBlack);
    px_e_out->Draw("HIST same");

    
    
    
    ndiffP->cd(2);
    py_e_outNO->SetLineColor(kRed);
    py_e_outNO->Draw("HIST");
    py_e_out->SetLineColor(kBlack);
    py_e_out->Draw("HIST same");

    
    
    
    
    ndiffP->cd(3);
    pz_e_outNO->SetLineColor(kRed);
    pz_e_outNO->Draw("HIST");
    pz_e_out->SetLineColor(kBlack);
    pz_e_out->Draw("HIST same");

    
    
    
    ndiffP->cd(4);
    px_mu_outNO->SetLineColor(40);
    px_mu_outNO->Draw("HIST");
    px_mu_out->SetLineColor(46);
    px_mu_out->Draw("HIST same");

    
    
    ndiffP->cd(5);
    py_mu_outNO->SetLineColor(40);
    py_mu_outNO->Draw("HIST");
    py_mu_out->SetLineColor(46);
    py_mu_out->Draw("HIST same");

    
    
    ndiffP->cd(6);
    pz_mu_outNO->SetLineColor(40);
    pz_mu_outNO->Draw("HIST");
    pz_mu_out->SetLineColor(46);
    pz_mu_out->Draw("HIST same");

    
    
    ndiffP->SaveAs("LO-NLOnodiv.png");
    
    TCanvas * t= new TCanvas("t","t",400,10,1500,1000);
    th->SetLineColor(40);
    th->Draw("HIST");
    thBR->SetLineColor(46);
    thBR->Draw("HIST same");
    
    t->SaveAs("theta.png");
    
    
        TCanvas * tar= new TCanvas("tar","tar",400,10,1500,1000);
    tar->Divide(2,2);
    
    tar->cd(1);
  thXZmu_one->SetLineColor(30);
    thXZmu_one->Draw("HIST");
      thXZe_one->SetLineColor(38);
    thXZe_one->Draw("HIST same");
    
    tar->cd(2);
      thYZmu_one->SetLineColor(30);
    thYZmu_one->Draw("HIST");
      thYZe_one->SetLineColor(38);
    thYZe_one->Draw("HIST same");
    
    tar->cd(3);
      thXZmu_two->SetLineColor(30);
    thXZmu_two->Draw("HIST");
      thXZe_two->SetLineColor(38);
    thXZe_two->Draw("HIST same");
    
    tar->cd(4);
  thYZmu_two->SetLineColor(30);
    thYZmu_two->Draw("HIST");
      thYZe_two->SetLineColor(38);
    thYZe_two->Draw("HIST same");
    

    
  tar->SaveAs("thXZYZ.png");
    
    
        
    
}