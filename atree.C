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
TH1F* E_muCODE=new TH1F("h1C", "Energy muon", 150,0,200);
TH1F* E_eCODE=new TH1F("h2C", "Energy electron", 150,0,200);
TH1F* p_mu=new TH1F("h1", "momentum out muon", 150,0,200);
TH1F* p_e=new TH1F("h2", "momentum out electron", 150,0,200);
TH1F* p_muCODE=new TH1F("h1C", "momentum k.P() out muon", 150,0,200);
TH1F* p_eCODE=new TH1F("h2C", "momentum k.P() out electron", 150,0,200);
    
TH1F* BRE_mu=new TH1F("h1", "Energy muon", 150,0,200);
TH1F* BRE_e=new TH1F("h2", "Energy electron", 150,0,200);
TH1F* BRE_muCODE=new TH1F("h1C", "Energy muon", 150,0,200);
TH1F* BRE_eCODE=new TH1F("h2C", "Energy electron", 150,0,200);
TH1F* BRp_mu=new TH1F("h1", "momentum out muon", 150,0,200);
TH1F* BRp_e=new TH1F("h2", "momentum out electron", 150,0,200);
TH1F* BRp_muCODE=new TH1F("h1C", "momentum k.P() out muon", 150,0,200);
TH1F* BRp_eCODE=new TH1F("h2C", "momentum k.P() out electron", 150,0,200);    

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
    
TH1F* diffePX=new TH1F("h1N", "pX_out diff electron div-NO DIV", 150,-0.02,0.02);
TH1F* diffePY=new TH1F("h2N", "pY_out diff electron div-NO DIV", 150,-0.02,0.02);
TH1F* diffePZ=new TH1F("h3N", "pZ_out diff electron div-NO DIV", 150,-0.001,0.001);
        
    
TH1F* coox_mu=new TH1F("h1", "Coo X mu", 140,-0.1,0.1);
TH1F* cooy_mu=new TH1F("h2", "Coo Y mu", 140,-0.1,0.1);
TH1F* coox_e=new TH1F("h1e", "Coo X e", 140,-0.1,0.1);
TH1F* cooy_e=new TH1F("h2e", "Coo Y e", 140,-0.1,0.1);
 
TH1F* tarONEXmu=new TH1F("h1a", "Coo X mu tar1  ", 140,-0.15,0.15);
TH1F* tarONEYmu=new TH1F("h2a", "Coo Y mu tar1 ", 140,-0.15,0.15);
TH1F* tarONEXe=new TH1F("h1ea", "Coo X e tar1 ", 140,-0.15,0.15);
TH1F* tarONEYe=new TH1F("h2ea", "Coo Y e tar1 ", 140,-0.15,0.15);
    
TH1F* tarTWOXmu=new TH1F("h1a", "Coo X mu tar2  ", 140,-0.15,0.15);
TH1F* tarTWOYmu=new TH1F("h2a", "Coo Y mu tar2 ", 140,-0.15,0.15);
TH1F* tarTWOXe=new TH1F("h1ea", "Coo X e tar2 ", 140,-0.15,0.15);
TH1F* tarTWOYe=new TH1F("h2ea", "Coo Y e tar2 ", 140,-0.15,0.15);
 
TH1F* diffX_mue=new TH1F("h", "DiffCoo X mu and e-", 140,-0.3,0.3);
TH1F* diffY_mue=new TH1F("h", "DiffCoo Y mu and e-", 140,-0.3,0.3);
 
TH2F  *X_Y_mu  = new TH2F("h2d" , " X  Vs. y of the muon",140,-0.5,-0.5,140,-0.5,0.5);
TH2F  *X_Y_e  = new TH2F("h2da" , " X  Vs. y of the electron",140,-0.5,-0.5,140,-0.5,0.5);
 
 

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
Double_t pmuOUTBR = sqrt(detKinBeamRot_pXmu_out*detKinBeamRot_pXmu_out+detKinBeamRot_pYmu_out*detKinBeamRot_pYmu_out+detKinBeamRot_pZmu_out*detKinBeamRot_pZmu_out);     

Double_t peOUTBR = sqrt(detKinBeamRot_pXe_out*detKinBeamRot_pXe_out+detKinBeamRot_pYe_out*detKinBeamRot_pYe_out+detKinBeamRot_pZe_out*detKinBeamRot_pZe_out);

Double_t pmuOUT = sqrt(detKin_pXmu_out*detKin_pXmu_out+detKin_pYmu_out*detKin_pYmu_out+detKin_pZmu_out*detKin_pZmu_out);
    
Double_t peOUT = sqrt(detKin_pXe_out*detKin_pXe_out+detKin_pYe_out*detKin_pYe_out+detKin_pZe_out*detKin_pZe_out);
 
       

//Double_t DE_eBR =  (2* (0.5109989461 *0.001)*(0.5109989461 *0.001)-  detKinBeamRot_t24)/(2* (0.5109989461 *0.001));
//Double_t DE_e =  (2* (0.5109989461 *0.001)*(0.5109989461 *0.001)-  detKin_t24)/(2* (0.5109989461 *0.001));
       
       
Double_t DE_mu=sqrt(pmuOUT*pmuOUT+(105.6583745 *0.001)*(105.6583745 *0.001));
Double_t DE_e=sqrt(peOUT*peOUT+(0.5109989461 *0.001)*(0.5109989461 *0.001));

Double_t DE_muBR=sqrt(pmuOUTBR*pmuOUTBR+(105.6583745 *0.001)*(105.6583745 *0.001));
Double_t DE_eBR=sqrt(peOUTBR*peOUTBR+(0.5109989461 *0.001)*(0.5109989461 *0.001));
       
         
       E_muCODE->Fill(detKin_Emu,wgt_full);
       E_eCODE->Fill(detKin_Ee,wgt_full);
       E_mu->Fill(DE_mu,wgt_full);
       E_e->Fill(DE_e,wgt_full);
       
       BRE_muCODE->Fill(detKinBeamRot_Emu,wgt_full);
       BRE_eCODE->Fill(detKinBeamRot_Ee,wgt_full);
       BRE_mu->Fill(DE_muBR,wgt_full);
       BRE_e->Fill(DE_eBR,wgt_full);

         
       p_muCODE->Fill(detKin_Pmu_out,wgt_full); //estratto da evento k.P()
       p_eCODE->Fill(detKin_Pe_out,wgt_full); //estratto da evento k.P()
       p_mu->Fill(pmuOUT,wgt_full); //calcolato con px,py,pz
       p_e->Fill(peOUT,wgt_full);   //calcolato con px,py,pz
       
       BRp_muCODE->Fill(detKinBeamRot_Pmu_out,wgt_full); //estratto da evento k.P()
       BRp_eCODE->Fill(detKinBeamRot_Pe_out,wgt_full); //estratto da evento k.P()
       BRp_mu->Fill(pmuOUTBR,wgt_full); //calcolato con px,py,pz
       BRp_e->Fill(peOUTBR,wgt_full);   //calcolato con px,py,pz
       
       
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
    px_e_outNO->SetLineColor(kRed);
    px_e_outNO->Draw("HIST");
    px_e_out->SetLineColor(kBlack);
    px_e_out->Draw("HIST same");

    
    
    
    diffP->cd(2);
    py_e_outNO->SetLineColor(kRed);
    py_e_outNO->Draw("HIST");
    py_e_out->SetLineColor(kBlack);
    py_e_out->Draw("HIST same");

    
    
    
    
    diffP->cd(3);
    pz_e_outNO->SetLineColor(kRed);
    pz_e_outNO->Draw("HIST");
    pz_e_out->SetLineColor(kBlack);
    pz_e_out->Draw("HIST same");

    
    
    
    diffP->cd(4);
    px_mu_outNO->SetLineColor(40);
    px_mu_outNO->Draw("HIST");
    px_mu_out->SetLineColor(46);
    px_mu_out->Draw("HIST same");

    
    
    diffP->cd(5);
    py_mu_outNO->SetLineColor(40);
    py_mu_outNO->Draw("HIST");
    py_mu_out->SetLineColor(46);
    py_mu_out->Draw("HIST same");

    
    
    diffP->cd(6);
    pz_mu_outNO->SetLineColor(40);
    pz_mu_outNO->Draw("HIST");
    pz_mu_out->SetLineColor(46);
    pz_mu_out->Draw("HIST same");

    
    
    diffP->SaveAs("divNodiv.png");
    
        
        
        
     TCanvas * d= new TCanvas("d","d",400,10,1500,1000);
    d->Divide(3,1);
    d->cd(1);
    diffePX->Draw("HIST");
    d->cd(2);
    diffePY->Draw("HIST");
    d->cd(3);
    diffePZ->Draw("HIST");
    d->SaveAs("diffEp.png");
    
    
    
    
    TCanvas * E= new TCanvas("E","E",400,10,1500,1000);
    E->Divide(2,1);
    E->cd(1);
    E_mu->SetLineColor(kBlue);
    E_mu->Draw("HIST");
    E_muCODE->SetLineColor(kBlack);
    E_muCODE->Draw("HIST same");
    BRE_mu->SetLineColor(kRed);
    BRE_mu->Draw("HIST same");
    BRE_muCODE->SetLineColor(kOrange);
    BRE_muCODE->Draw("HIST same");
    E->cd(2);
    E_e->SetLineColor(kBlue);
    E_e->Draw("HIST");
    E_eCODE->SetLineColor(kBlack);
    E_eCODE->Draw("HIST same");
    BRE_e->SetLineColor(kRed);
    BRE_e->Draw("HIST same");
    BRE_eCODE->SetLineColor(kOrange);
    BRE_eCODE->Draw("HIST same");    
    
    E->SaveAs("energy.png");
    
    TCanvas * P= new TCanvas("P","P",400,10,1500,1000);
    P->Divide(2,1);
    P->cd(1);
    p_mu->SetLineColor(kBlue);
    p_mu->Draw("HIST");
    p_muCODE->SetLineColor(kBlack);
    p_muCODE->Draw("HIST same");
    BRp_mu->SetLineColor(kRed);
    BRp_mu->Draw("HIST same");
    BRp_muCODE->SetLineColor(kOrange);
    BRp_muCODE->Draw("HIST same");
    P->cd(2);
    p_e->SetLineColor(kBlue);
    p_e->Draw("HIST");
    p_eCODE->SetLineColor(kBlack);
    p_eCODE->Draw("HIST same");
    BRp_e->SetLineColor(kRed);
    BRp_e->Draw("HIST same");
    BRp_eCODE->SetLineColor(kOrange);
    BRp_eCODE->Draw("HIST same");    
    
    P->SaveAs("momentumTOT.png");
    
    

    
    TCanvas * Pin= new TCanvas("Pin","Pin",400,10,1500,1000);
    Pin->Divide(1,3);
    Pin->cd(1);
    px_mu->Draw("HIST");
    Pin->cd(2);
    py_mu->Draw("HIST");
    Pin->cd(3);
    pz_mu->Draw("HIST");
    
   Pin->SaveAs("p_in.png");
        
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
    diff->SaveAs("diff.png");


     
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
  cooX->SaveAs("coo.png");
 
    TCanvas * dued= new TCanvas("dued","dued",1000,100,1500,2000);

  X_Y_mu->SetMarkerColor(kBlack);
    X_Y_mu->Draw("HIST");
  X_Y_e->SetMarkerColor(kRed);
    X_Y_e->Draw("HIST same");

    
  dued->SaveAs("duedcoo.png");
    
    TCanvas * tar= new TCanvas("tar","tar",400,10,1500,1000);
    tar->Divide(2,2);
    tar->cd(1);
  tarONEXmu->SetLineColor(42);
    tarONEXmu->Draw("HIST");
      tarONEXe->SetLineColor(46);
    tarONEXe->Draw("HIST same");
    tar->cd(2);
      tarONEYmu->SetLineColor(42);
    tarONEXmu->Draw("HIST");
      tarONEYe->SetLineColor(46);
    tarONEYe->Draw("HIST same");
    tar->cd(3);
      tarTWOXmu->SetLineColor(42);
    tarTWOXmu->Draw("HIST");
      tarTWOXe->SetLineColor(46);
    tarTWOXe->Draw("HIST same");
    tar->cd(4);
  tarTWOYmu->SetLineColor(42);
    tarTWOYmu->Draw("HIST");
      tarTWOYe->SetLineColor(46);
    tarTWOYe->Draw("HIST same");
    
  tar->SaveAs("tar.png");

    
}


