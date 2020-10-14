#include <TF2.h>
#include <TMath.h>
#include <cmath>
using namespace std;

int main(){

Double_t Rm = 1.959 ; //raggio di Moliere in centimetri    

    
TF2 *fxy_e = new TF2("fxy_e","(1/(2*pi*[0]*[1]))*(exp(-((x-[2])*(x-[2])+(y-[3])*(y-[3]))/(2*[0]*[1])))",-7.125,7.125,-7.125,7.125 );
    
TF2 *fxy_ph = new TF2("fxy_ph","(1/(2*pi*[0]*[1]))*(exp(-((x-[2])*(x-[2])+(y-[3])*(y-[3]))/(2*[0]*[1])))",-7.125,7.125,-7.125,7.125 );

    
TMatrixD E(3,25);

    
    if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();


    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       

        if (photon_coox!=-1 && photon_cooy!=-1)
    {  
     E_CAL=detKinBeamRot_Ee+photon_energy;}
    else 
    {  E_CAL=detKinBeamRot_Ee;}
       
       
    Double_t d_e_ph=sqrt( (detKinBeamRot_cooXe-photon_coox)*(detKinBeamRot_cooXe-photon_coox)+(detKinBeamRot_cooYe-photon_cooy)*(detKinBeamRot_cooYe-photon_cooy) ); 
       
    Double_t d_e_mu=sqrt( (detKinBeamRot_cooXe-detKinBeamRot_cooXmu)*(detKinBeamRot_cooXe-detKinBeamRot_cooXmu)+(detKinBeamRot_cooYe-detKinBeamRot_cooYmu)*(detKinBeamRot_cooYe-detKinBeamRot_cooYmu) );
       
   
       
       
//if (E_CAL>1){
if (abs(detKinBeamRot_cooXe)<7.125 && abs(detKinBeamRot_cooYe)<7.125)
    
{ 
    fxy_e->SetParameter(0,Rm);
    fxy_e->SetParameter(1,Rm);
    fxy_e->SetParameter(2,detKinBeamRot_cooXe);
    fxy_e->SetParameter(3,detKinBeamRot_cooYe);
    
    
    //cout << "cella elettrone:" << n_cell << endl; 
       
E[0][0]= detKinBeamRot_Ee*fxy_e->Integral(-7.125,-4.275,4.275,7.125); //cella 1
E[0][1]= detKinBeamRot_Ee*fxy_e->Integral(-4.275,-1.425,4.275,7.125); //cella 2
E[0][2]= detKinBeamRot_Ee*fxy_e->Integral(-1.425,1.425,4.275,7.125); //cella 3
E[0][3]= detKinBeamRot_Ee*fxy_e->Integral(1.425,4.275,4.275,7.125); //cella 4
E[0][4]= detKinBeamRot_Ee*fxy_e->Integral(4.275,7.125,4.275,7.125); //cella 5
E[0][5]= detKinBeamRot_Ee*fxy_e->Integral(-7.125,-4.275,1.425,4.275); //cella 6
E[0][6]= detKinBeamRot_Ee*fxy_e->Integral(-4.275,-1.425,1.425,4.275); //cella 7
E[0][7]= detKinBeamRot_Ee*fxy_e->Integral(-1.425,1.425,1.425,4.275); //cella 8
E[0][8]= detKinBeamRot_Ee*fxy_e->Integral(1.425,4.275,1.425,4.275); //cella 9
E[0][9]= detKinBeamRot_Ee*fxy_e->Integral(4.275,7.125,1.425,4.275); //cella 10
E[0][10]= detKinBeamRot_Ee*fxy_e->Integral(-7.125,-4.275,-1.425,1.425); //cella 11
E[0][11]= detKinBeamRot_Ee*fxy_e->Integral(-4.275,-1.425,-1.425,1.425); //cella 12
E[0][12]= detKinBeamRot_Ee*fxy_e->Integral(-1.425,1.425,-1.425,1.425); //cella 13
E[0][13]= detKinBeamRot_Ee*fxy_e->Integral(1.425,4.275,-1.425,1.425); //cella 14
E[0][14]= detKinBeamRot_Ee*fxy_e->Integral(4.275,7.125,-1.425,1.425); //cella 15
E[0][15]= detKinBeamRot_Ee*fxy_e->Integral(-7.125,-4.275,-4.275,-1.425); //cella 16
E[0][16]= detKinBeamRot_Ee*fxy_e->Integral(-4.275,-1.425,-4.275,-1.425); //cella 17
E[0][17]= detKinBeamRot_Ee*fxy_e->Integral(-1.425,1.425,-4.275,-1.425); //cella 18
E[0][18]= detKinBeamRot_Ee*fxy_e->Integral(1.425,4.275,-4.275,-1.425); //cella 19
E[0][19]= detKinBeamRot_Ee*fxy_e->Integral(4.275,7.125,-4.275,-1.425); //cella 20
E[0][20]= detKinBeamRot_Ee*fxy_e->Integral(-7.125,-4.275,-7.125,-4.275); //cella 21
E[0][21]= detKinBeamRot_Ee*fxy_e->Integral(-4.275,-1.425,-7.125,-4.275); //cella 22
E[0][22]= detKinBeamRot_Ee*fxy_e->Integral(-1.425,1.425,-7.125,-4.275); //cella 23
E[0][23]= detKinBeamRot_Ee*fxy_e->Integral(1.425,4.275,-7.125,-4.275); //cella 24
E[0][24]= detKinBeamRot_Ee*fxy_e->Integral(4.275,7.125,-7.125,-4.275); //cella 25
       
if (abs(photon_coox)<7.125 && abs(photon_cooy)<7.125 && photon_energy>0.2)
    
{   fxy_ph->SetParameter(0,Rm);
    fxy_ph->SetParameter(1,Rm);
    fxy_ph->SetParameter(2,photon_coox);
    fxy_ph->SetParameter(3,photon_cooy);
    
E[1][0]= photon_energy*fxy_ph->Integral(-7.125,-4.275,4.275,7.125); //cella 1
E[1][1]= photon_energy*fxy_ph->Integral(-4.275,-1.425,4.275,7.125); //cella 2
E[1][2]= photon_energy*fxy_ph->Integral(-1.425,1.425,4.275,7.125); //cella 3
E[1][3]= photon_energy*fxy_ph->Integral(1.425,4.275,4.275,7.125); //cella 4
E[1][4]= photon_energy*fxy_ph->Integral(4.275,7.125,4.275,7.125); //cella 5
E[1][5]= photon_energy*fxy_ph->Integral(-7.125,-4.275,1.425,4.275); //cella 6
E[1][6]= photon_energy*fxy_ph->Integral(-4.275,-1.425,1.425,4.275); //cella 7
E[1][7]= photon_energy*fxy_ph->Integral(-1.425,1.425,1.425,4.275); //cella 8
E[1][8]= photon_energy*fxy_ph->Integral(1.425,4.275,1.425,4.275); //cella 9
E[1][9]= photon_energy*fxy_ph->Integral(4.275,7.125,1.425,4.275); //cella 10
E[1][10]= photon_energy*fxy_ph->Integral(-7.125,-4.275,-1.425,1.425); //cella 11
E[1][11]= photon_energy*fxy_ph->Integral(-4.275,-1.425,-1.425,1.425); //cella 12
E[1][12]= photon_energy*fxy_ph->Integral(-1.425,1.425,-1.425,1.425); //cella 13
E[1][13]= photon_energy*fxy_ph->Integral(1.425,4.275,-1.425,1.425); //cella 14
E[1][14]= photon_energy*fxy_ph->Integral(4.275,7.125,-1.425,1.425); //cella 15
E[1][15]= photon_energy*fxy_ph->Integral(-7.125,-4.275,-4.275,-1.425); //cella 16
E[1][16]= photon_energy*fxy_ph->Integral(-4.275,-1.425,-4.275,-1.425); //cella 17
E[1][17]= photon_energy*fxy_ph->Integral(-1.425,1.425,-4.275,-1.425); //cella 18
E[1][18]= photon_energy*fxy_ph->Integral(1.425,4.275,-4.275,-1.425); //cella 19
E[1][19]= photon_energy*fxy_ph->Integral(4.275,7.125,-4.275,-1.425); //cella 20
E[1][20]= photon_energy*fxy_ph->Integral(-7.125,-4.275,-7.125,-4.275); //cella 21
E[1][21]= photon_energy*fxy_ph->Integral(-4.275,-1.425,-7.125,-4.275); //cella 22
E[1][22]= photon_energy*fxy_ph->Integral(-1.425,1.425,-7.125,-4.275); //cella 23
E[1][23]= photon_energy*fxy_ph->Integral(1.425,4.275,-7.125,-4.275); //cella 24
E[1][24]= photon_energy*fxy_ph->Integral(4.275,7.125,-7.125,-4.275); //cella 25 
  
// ENERGIA TOTALE NELLE CELLE QUANDO CI SONO ELETTRONE E FOTONE
E[2][0]= E[0][0]+E[1][0]; //cella 1
E[2][1]= E[0][1]+E[1][1]; //cella 2
E[2][2]= E[0][2]+E[1][2]; //cella 3
E[2][3]= E[0][3]+E[1][3]; //cella 4
E[2][4]= E[0][4]+E[1][4]; //cella 5
E[2][5]= E[0][5]+E[1][5]; //cella 6
E[2][6]= E[0][6]+E[1][6]; //cella 7
E[2][7]= E[0][7]+E[1][7]; //cella 8
E[2][8]= E[0][8]+E[1][8]; //cella 9
E[2][9]= E[0][9]+E[1][9]; //cella 10
E[2][10]= E[0][10]+E[1][10]; //cella 11
E[2][11]= E[0][11]+E[1][11]; //cella 12
E[2][12]= E[0][12]+E[1][12]; //cella 13
E[2][13]= E[0][13]+E[1][13]; //cella 14
E[2][14]= E[0][14]+E[1][14]; //cella 15
E[2][15]= E[0][15]+E[1][15]; //cella 16
E[2][16]= E[0][16]+E[1][16]; //cella 17
E[2][17]= E[0][17]+E[1][17]; //cella 18
E[2][18]= E[0][18]+E[1][18]; //cella 19
E[2][19]= E[0][19]+E[1][19]; //cella 20
E[2][20]= E[0][20]+E[1][20]; //cella 21
E[2][21]= E[0][21]+E[1][21]; //cella 22
E[2][22]= E[0][22]+E[1][22]; //cella 23
E[2][23]= E[0][23]+E[1][23]; //cella 24
E[2][24]= E[0][24]+E[1][24]; //cella 25 
            

            } 
    else { 
// ENERGIA TOTALE NELLE CELLE QUANDO C'E' SOLO L'ELETTRONE
E[2][0]= E[0][0]; //cella 1
E[2][1]= E[0][1]; //cella 2
E[2][2]= E[0][2]; //cella 3
E[2][3]= E[0][3]; //cella 4
E[2][4]= E[0][4]; //cella 5
E[2][5]= E[0][5]; //cella 6
E[2][6]= E[0][6]; //cella 7
E[2][7]= E[0][7]; //cella 8
E[2][8]= E[0][8]; //cella 9
E[2][9]= E[0][9]; //cella 10
E[2][10]= E[0][10]; //cella 11
E[2][11]= E[0][11]; //cella 12
E[2][12]= E[0][12]; //cella 13
E[2][13]= E[0][13]; //cella 14
E[2][14]= E[0][14]; //cella 15
E[2][15]= E[0][15]; //cella 16
E[2][16]= E[0][16]; //cella 17
E[2][17]= E[0][17]; //cella 18
E[2][18]= E[0][18]; //cella 19
E[2][19]= E[0][19]; //cella 20
E[2][20]= E[0][20]; //cella 21
E[2][21]= E[0][21]; //cella 22
E[2][22]= E[0][22]; //cella 23
E[2][23]= E[0][23]; //cella 24
E[2][24]= E[0][24]; //cella 25 
        }

    
   }
       

       
       
}
    
}