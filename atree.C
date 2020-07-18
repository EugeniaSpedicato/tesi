#define atree_cxx
#include "atree.h"
#include <TH2.h>
#include <TH1.h>

#include <TStyle.h>
#include <TCanvas.h>

void atree::Loop()
{
    TH1::SetDefaultSumw2();

TH1F* E_mu=new TH1F("h1", "Energy muon", 150,0,200);
TH1F* E_e=new TH1F("h2", "Energy electron", 150,0,200);

TH1F* px_mu=new TH1F("h1", "pX_in muon", 190,-0.3,0.3);
TH1F* py_mu=new TH1F("h2", "pY_in muon", 190,-0.3,0.3);
TH1F* pz_mu=new TH1F("h3", "pZ_in muon", 190,50,180);
    
TH1F* px_mu_out=new TH1F("h1a", "pX_out muon", 150,-0.3,0.3);
TH1F* py_mu_out=new TH1F("h2a", "pY_out muon", 150,-0.3,0.3);
TH1F* pz_mu_out=new TH1F("h3a", "pZ_out muon", 150,0,180);
TH1F* px_mu_outNO=new TH1F("h1aN", "pX_out muon NO DIV", 150,-0.3,0.3);
TH1F* py_mu_outNO=new TH1F("h2aN", "pY_out muon NO DIV", 150,-0.3,0.3);
TH1F* pz_mu_outNO=new TH1F("h3aN", "pZ_out muon NO DIV", 150,0,180);
    
TH1F* px_e_out=new TH1F("h1b", "pX_out electron", 150,-0.3,0.3);
TH1F* py_e_out=new TH1F("h2b", "pY_out electron", 150,-0.3,0.3);
TH1F* pz_e_out=new TH1F("h3b", "pZ_out electron", 150,0,5);
TH1F* px_e_outNO=new TH1F("h1bN", "pX_out electron NO DIV", 150,-0.3,0.3);
TH1F* py_e_outNO=new TH1F("h2bN", "pY_out electron NO DIV", 150,-0.3,0.3);
TH1F* pz_e_outNO=new TH1F("h3bN", "pZ_out electron NO DIV", 150,0,5);
    
TH1F* diffePX=new TH1F("h1N", "pX_out diff electron div-NO DIV", 150,-0.3,0.3);
TH1F* diffePY=new TH1F("h2N", "pY_out diff electron div-NO DIV", 150,-0.3,0.3);
TH1F* diffePZ=new TH1F("h3N", "pZ_out diff electron div-NO DIV", 150,0,5);
        
    
TH1F* coox_mu=new TH1F("h1", "Coo X mu", 140,-0.1,0.1);
TH1F* cooy_mu=new TH1F("h2", "Coo Y mu", 140,-0.1,0.1);
TH1F* coox_e=new TH1F("h1e", "Coo X e", 140,-0.1,0.1);
TH1F* cooy_e=new TH1F("h2e", "Coo Y e", 140,-0.1,0.1);
 
TH1F* tarONEXmu=new TH1F("h1a", "Coo X mu tar1  ", 140,-0.1,0.1);
TH1F* tarONEYmu=new TH1F("h2a", "Coo Y mu tar1 ", 140,-0.1,0.1);
TH1F* tarONEXe=new TH1F("h1ea", "Coo X e tar1 ", 140,-0.1,0.1);
TH1F* tarONEYe=new TH1F("h2ea", "Coo Y e tar1 ", 140,-0.1,0.1);
    
TH1F* tarTWOXmu=new TH1F("h1a", "Coo X mu tar2  ", 140,-0.1,0.1);
TH1F* tarTWOYmu=new TH1F("h2a", "Coo Y mu tar2 ", 140,-0.1,0.1);
TH1F* tarTWOXe=new TH1F("h1ea", "Coo X e tar2 ", 140,-0.1,0.1);
TH1F* tarTWOYe=new TH1F("h2ea", "Coo Y e tar2 ", 140,-0.1,0.1);
 
TH1F* diffX_mue=new TH1F("h", "DiffCoo X mu and e-", 140,-0.3,0.3);
TH1F* diffY_mue=new TH1F("h", "DiffCoo Y mu and e-", 140,-0.3,0.3);
 
TH2F  *X_Y_mu  = new TH2F("h2d" , " X  Vs. y of the muon",140,-0.3,-0.3,100,0,40);
TH2F  *X_Y_e  = new TH2F("h2da" , " X  Vs. y of the electron",140,-0.3,-0.3,100,0,40);
 
 

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       Double_t pmuOUT=sqrt(detKinBeamRot_pXmu_out*detKinBeamRot_pXmu_out+detKinBeamRot_pYmu_out*detKinBeamRot_pYmu_out+detKinBeamRot_pZmu_out*detKinBeamRot_pZmu_out);
       Double_t peOUT=sqrt(detKinBeamRot_pXe_out*detKinBeamRot_pXe_out+detKinBeamRot_pYe_out*detKinBeamRot_pYe_out+detKinBeamRot_pZe_out*detKinBeamRot_pZe_out);
       Double_t DE_mu=sqrt(pmuOUT*pmuOUT+(105.6583745 *0.001)*(105.6583745 *0.001));
       Double_t DE_e=sqrt(peOUT*peOUT+(0.5109989461 *0.001)*(0.5109989461 *0.001));
       
       E_mu->Fill(DE_mu,wgt_full);
       E_e->Fill(DE_e,wgt_full);
    
       
       px_mu->Fill(detKinBeamRot_pXmu,wgt_full);
       py_mu->Fill(detKinBeamRot_pYmu,wgt_full);
       pz_mu->Fill(detKinBeamRot_pZmu,wgt_full);
       
       px_mu_out->Fill(detKinBeamRot_pXmu_out,wgt_full);
       py_mu_out->Fill(detKinBeamRot_pYmu_out,wgt_full);
       pz_mu_out->Fill(detKinBeamRot_pZmu_out,wgt_full);
       px_mu_outNO->Fill(detKin_pXmu_out,wgt_full);
       py_mu_outNO->Fill(detKin_pYmu_out,wgt_full);
       pz_mu_outNO->Fill(detKin_pZmu_out,wgt_full);
       
       px_e_out->Fill(detKinBeamRot_pXe_out,wgt_full);
       py_e_out->Fill(detKinBeamRot_pYe_out,wgt_full);
       pz_e_out->Fill(detKinBeamRot_pZe_out,wgt_full);
       px_e_outNO->Fill(detKin_pXe_out,wgt_full);
       py_e_outNO->Fill(detKin_pYe_out,wgt_full);
       pz_e_outNO->Fill(detKin_pZe_out,wgt_full);
       
        coox_mu->Fill(detKinBeamRot_cooXmu,wgt_full);
       cooy_mu->Fill(detKinBeamRot_cooYmu,wgt_full);
       
       coox_e->Fill(detKinBeamRot_cooXe,wgt_full);
       cooy_e->Fill(detKinBeamRot_cooYe,wgt_full);
       
       Double_t diffX=detKinBeamRot_cooXe-detKinBeamRot_cooXmu;
       diffX_mue->Fill(diffX,wgt_full);
       Double_t diffY=detKinBeamRot_cooYe-detKinBeamRot_cooYmu;
       diffY_mue->Fill(diffY,wgt_full);
       
       Double_t diffPX=detKinBeamRot_pXe_out-detKin_pXe_out;
       diffePX->Fill(diffPX,wgt_full);
       Double_t diffPY=detKinBeamRot_pYe_out-detKin_pYe_out;
       diffePY->Fill(diffPY,wgt_full);
       Double_t diffPZ=detKinBeamRot_pZe_out-detKin_pZe_out;
       diffePZ->Fill(diffPZ,wgt_full);
       
     X_Y_mu ->Fill(detKinBeamRot_cooXmu, detKinBeamRot_cooYmu,wgt_full);
     X_Y_e ->Fill(detKinBeamRot_cooXe, detKinBeamRot_cooYe,wgt_full);
       
       if (detKinBeamRot_tar==0)
       {
         tarONEXmu->Fill(detKinBeamRot_cooXmu,wgt_full);
         tarONEYmu->Fill(detKinBeamRot_cooYmu,wgt_full);
         tarONEXe->Fill(detKinBeamRot_cooXe,wgt_full);
         tarONEYe->Fill(detKinBeamRot_cooYe,wgt_full);
       }
       
       if (detKinBeamRot_tar==1)
       {
         tarTWOXmu->Fill(detKinBeamRot_cooXmu,wgt_full);
         tarTWOYmu->Fill(detKinBeamRot_cooYmu,wgt_full);
         tarTWOXe->Fill(detKinBeamRot_cooXe,wgt_full);
         tarTWOYe->Fill(detKinBeamRot_cooYe,wgt_full);
           
       }
       
   }

    px_mu->Fit("gaus");
    py_mu->Fit("gaus");
    pz_mu->Fit("gaus");
    
    TCanvas * diffP= new TCanvas("diff","diff",400,10,1500,1000);
    diffP->Divide(3,2);
    diffP->cd(1);
    px_e_outNO->SetLineColor(kOrange);
    px_e_outNO->Draw("HIST");
    px_e_out->SetLineColor(kBlack);
    px_e_out->Draw("HIST same");

    
    
    
    diffP->cd(2);
    py_e_outNO->SetLineColor(kOrange);
    py_e_outNO->Draw("HIST");
    py_e_out->SetLineColor(kBlack);
    py_e_out->Draw("HIST same");

    
    
    
    
    diffP->cd(3);
    pz_e_outNO->SetLineColor(kOrange);
    pz_e_outNO->Draw("HIST");
    pz_e_out->SetLineColor(kBlack);
    pz_e_out->Draw("HIST same");

    
    
    
    diffP->cd(4);
    px_mu_outNO->SetLineColor(40);
    px_mu_outNO->Draw("HIST");
    px_mu_out->SetLineColor(46);
    px_mu_out->Draw("HIST same");

    
    
    diffP->cd(5);
    px_mu_outNO->SetLineColor(40);
    px_mu_outNO->Draw("HIST");
    px_mu_out->SetLineColor(46);
    px_mu_out->Draw("HIST same");

    
    
    diffP->cd(6);
    px_mu_outNO->SetLineColor(40);
    px_mu_outNO->Draw("HIST");
    px_mu_out->SetLineColor(46);
    px_mu_out->Draw("HIST same");

    
    
    diffP->SaveAs("divNodiv.jpg");
    
        
        
        
     TCanvas * d= new TCanvas("d","d",400,10,1500,1000);
    d->Divide(3,1);
    d->cd(1);
    diffePX->Draw("HIST");
    d->cd(2);
    diffePY->Draw("HIST");
    d->cd(3);
    diffePZ->Draw("HIST");
    d->SaveAs("diffEp.jpg");
    
    
    
    
    TCanvas * E= new TCanvas("E","E",400,10,1500,1000);
    E->Divide(2,1);
    E->cd(1);
    E_mu->Draw("HIST");
    E->cd(2);
    E_e->SetLineColor(kRed);
    E_e->Draw("HIST");
    
    TCanvas * Pin= new TCanvas("Pin","Pin",400,10,1500,1000);
    Pin->Divide(1,3);
    Pin->cd(1);
    px_mu->Draw("HIST");
    Pin->cd(2);
    py_mu->Draw("HIST");
    Pin->cd(3);
    pz_mu->Draw("HIST");
        
    TCanvas * diff= new TCanvas("diff","diff",400,10,1500,1000);
    diff->Divide(2,2);
    diff->cd(1);
    px_e_out->SetLineColor(kRed);
    px_e_out->Draw("HIST");
    px_mu_out->Draw("HIST same");

    
    
    
    diff->cd(2);
    py_e_out->SetLineColor(kRed);
    py_e_out->Draw("HIST");
    py_mu_out->Draw("HIST same");

    
    
    
    diff->cd(3);
    pz_mu_out->Draw("HIST");
    
    diff->cd(4);
    
    pz_e_out->Draw("HIST");
    diff->SaveAs("diff.jpg");


     
    TCanvas * cooX= new TCanvas("cooX","cooX",400,10,600,400);
    cooX->Divide(2,2);
    cooX->cd(1);
    coox_mu->Draw("HIST");
    coox_e->SetLineColor(kRed);
    coox_e->Draw("HIST same");

    

 
    cooX->cd(2);
    cooy_mu->Draw("HIST");
    cooy_e->SetLineColor(kRed);
    cooy_e->Draw("HIST same");

    
 
    cooX->cd(3);
    diffX_mue->Draw("HIST");
      
    cooX->cd(4);
    diffY_mue->SetLineColor(kRed);
    diffY_mue->Draw("HIST");
  cooX->SaveAs("coo.jpg");
 
    TCanvas * dued= new TCanvas("dued","dued",400,10,600,400);

  X_Y_mu->SetMarkerColor(kBlack);
    X_Y_mu->Draw("HIST");
  X_Y_e->SetMarkerColor(kRed);
    X_Y_e->Draw("HIST same");

    
  dued->SaveAs("duedcoo.jpg");
    
    TCanvas * tar= new TCanvas("tar","tar",400,10,600,400);
    tar->Divide(2,1);
    tar->cd(1);
  tarONEXmu->SetMarkerColor(30);
    tarONEXmu->Draw("HIST");
      tarONEXe->SetMarkerColor(49);
    tarONEXe->Draw("HIST same");
    tar->cd(2);
  tarTWOYmu->SetMarkerColor(50);
    tarTWOYmu->Draw("HIST");
      tarTWOYe->SetMarkerColor(49);
    tarTWOYe->Draw("HIST same");
  tar->SaveAs("tar.jpg");

    
}


