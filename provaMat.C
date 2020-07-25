#include "TVectorD.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/RotationZYX.h"
#include "TRandom.h"
#include "TMatrixF.h" 
#include "TMatrixFBase.h"
#include <TMatrixFSym.h>
#include "TString.h"
#include <TApplication.h>

#include <cmath>
#include <iostream>

using namespace std;
using namespace ROOT::Math;


  
 PxPyPzEVector RotDiv(const PxPyPzEVector & k)
 {
Double_t divthx = gRandom->Gaus(0., 0.00027);
Double_t divthy = gRandom->Gaus(0., 0.00020); 

Double_t anglex = atan2(k.Px(), k.Pz());
Double_t angley = atan2(k.Py(), k.Pz()); 
     
anglex += divthx;
angley += divthy;
     
    //NB questi Px Py Pz sono del nuovo!!
Double_t pmuin=sqrt(k.Px()*k.Px()+k.Py()*k.Py()+k.Pz()*k.Pz());
Double_t pz=pmuin/(1+tan(anglex)*tan(anglex)+tan(angley)*tan(angley));
Double_t py=pz*tan(angley);
Double_t px=pz*tan(anglex);
Double_t pt=sqrt(px*px+py*py);
Double_t ptz=sqrt(pz*pz+px*px);
     


XYZVector p_div(px,py,pz);

Double_t pOx=k.Px();
Double_t pOy=k.Py();
Double_t pOz=k.Pz();

    
XYZVector p_old(pOx,pOy,pOz);
/*  
Double_t psi=atan2(px,pz); 
Double_t phi=-atan2(py,ptz);  
     
//ruotare il sistema di riferimento con r2 e poi r1 equivale a ruotare il vettore p vecchio con r1 e poi r2. L'angolo poi deve essere opposto
RotationY r1(psi);
RotationX r2(phi); 

XYZVector p_1=r1*p_old;
XYZVector p_new=r2*p_1;
Double_t pNx=p_new.X();
Double_t pNy=p_new.Y();
Double_t pNz=p_new.Z();
*/
cout<<" pxN= "<< px <<" pyN= " <<py<<" pzN= " << pz << endl; 

 
 
 
Double_t psi=atan2(px,pz); 
Double_t phi=atan2(py,ptz);  
     
TMatrixD R(3,3);
R[0][0]=cos(psi);
R[0][1]=0;
R[0][2]=sin(psi);
    R[1][0]=-sin(phi)*sin(psi);
    R[1][1]=cos(phi);
    R[1][2]=sin(phi)*cos(psi);
        R[2][0]=-cos(phi)*sin(psi);
        R[2][1]=-sin(phi);
        R[2][2]=cos(phi)*cos(psi);
     

TMatrixD pO(3,1);
    pO[0][0]=k.Px();
    pO[1][0]=k.Py();
    pO[2][0]=k.Pz();

Double_t det;

TMatrixD pN(R, TMatrixD::kMult,pO);

PxPyPzEVector pnewdiv(pN[0][0], pN[1][0], pN[2][0], k.E()); 
     
cout<<" pxN= "<< pN[0][0] <<" pyN= " <<pN[1][0]<<" pzN= " << pN[2][0] << endl; 
     
cout<<" Tx= "<< anglex <<" Ty= "<< angley << " Dx= "<< divthx <<endl;
return pnewdiv;    

     
     
     
     
     
     
     

    

//PxPyPzEVector pnewdiv(pNx, pNy, pNz, k.E());

cout<< divthx <<" " <<divthy<<endl;




 }   
/*  

//  devi prendere il theta Mu e theta E, metti sotto loadkine vars, sarà forse genKin.the o .thmu e .phe phmu
TMatrixD coo(const PxPyPzEVector & k,const Double_t & the, const Double_t & phi) 
{   Double_t theR = the;//rad
    Double_t phiR = phi;//rad
 
Double_t anglex = atan2(k.Px(), k.Pz());
Double_t angley = atan2(k.Py(), k.Pz()); 
    
    Double_t d0=2.10;//m
    Double_t d1=1.10;//m
 
    Double_t x=gRandom->Gaus(0., 0.026);//m
    Double_t y=gRandom->Gaus(0., 0.027);//m
    Double_t z=0.;
    Double_t zf=2.10;
    TMatrixD coo_in(1,3);
 coo_in[0][0]=0;
 coo_in[0][1]=0;
 coo_in[0][2]=0;
 
    //interazione target 1 o target 2 
    Int_t tar=gRandom->Integer(2);
    
    if (tar==0)
    {

        //primo modo
    Double_t d_xy = d0*tan(theR);//vettore nel piano xy
    Double_t xf = x+d_xy*cos(phiR);
    Double_t yf = y+d_xy*sin(phiR);
        //secondo modo
    Double_t xfe = x+d0*tan(anglex);
    Double_t yfe = y+d0*tan(angley);
    
    TMatrixD coo_f(2,3);
    coo_f [0][0]=xf;  
    coo_f [0][1]=yf;  
    coo_f [0][2]=tar;  
    coo_f [1][0]=xfe;  
    coo_f [1][1]=yfe;  
    coo_f [1][2]=tar;  
        

         return coo_f;
    }
    
    if(tar==1)
    {
        //primo modo
    Double_t d_xy = d1*tan(theR);//vettore nel piano xy
    Double_t xf = x+d_xy*cos(phiR);
    Double_t yf = y+d_xy*sin(phiR);
        //secondo modo
    Double_t xfe = x+d1*tan(anglex);
    Double_t yfe = y+d1*tan(angley);

    TMatrixD coo_f (2,3);
    coo_f [0][0]=xf;  
    coo_f [0][1]=yf;  
    coo_f [0][2]=tar;  
    coo_f [1][0]=xfe;  
    coo_f [1][1]=yfe;  
    coo_f [1][2]=tar;  
        

    return coo_f;

    }
      
        
else return coo_in;
}



TMatrixD MCSin(const PxPyPzEVector & k)
{   
Double_t const energy   = 150000; //MeV

Double_t const sS    = 0.00032; //m spessore silicio
Double_t const x0S = 0.094; // m

Double_t const sB    = 0.015; //m spessore berillio
Double_t const x0B = 0.353; // m



Double_t const dSS = 0.01; // m distanza tra i due 2S


Double_t const d = 0.25-0.005;
Double_t const dx = 0.25+0.005; // m la distanza tra coppie di silici è 0.25. Però 0.25 è tra i due della coppia, quindi considero -dSS/2=0.005

// m la distanza tra coppie di silici è 0.25. Però 0.25 è tra i due della coppia, quindi considero -dSS/2=0.005
Double_t const ris = 18e-6; // m considero questa la risoluzione dei silici


Double_t const sigX=0.026; //m
Double_t const sigY=0.027; //m
        
    
Double_t sigSI=(13.6/energy)*sqrt(sS/x0S)*(1+0.038*log(sS/x0S)); //rad
    
    
    Double_t thetaX[7];
    Double_t thetaY[7];
    Double_t x[7];
    Double_t y[7];
    
    Double_t anglex = atan2(k.Px(), k.Pz());
    Double_t angley = atan2(k.Py(), k.Pz()); 
    
 
        
    //questi parametri sono dati dallo spread intrinseco del fascio
           
            thetaX[0]=anglex;
            x[0]=gRandom->Gaus(0., sigX);
            
            thetaY[0]=angley;
            y[0]=gRandom->Gaus(0., sigY);
   

// ora entra nel primo silicio, considero una distanza dall' origine d=0.25-0.005 m. Mi calcolo ThetaX e ThetaY all'uscita, x e y sono all'uscita e considero anche la risoluzione ris data dal silicio che ora misura x, quindi solo sulla x
                thetaX[1]= gRandom->Gaus(thetaX[0],sigSI);
                thetaY[1]= gRandom->Gaus(thetaY[0],sigSI); 

                x[1]=x[0]+(1/sqrt(3))*sS*thetaX[1]+gRandom->Gaus(0,ris);
                y[1]=y[0]+(1/sqrt(3))*sS*thetaY[1];
                
// ora secondo silicio della coppia che misura la y quindi ris va aggiunta alla y
                
                thetaX[2]= gRandom->Gaus(thetaX[1],sigSI);
                thetaY[2]= gRandom->Gaus(thetaY[1],sigSI);

                x[2]=x[1]+dSS*tan(thetaX[1])+(1/sqrt(3))*sS*thetaX[2];
                y[2]=y[1]+dSS*tan(thetaY[1])+(1/sqrt(3))*sS*thetaY[2]+gRandom->Gaus(0,ris);
                
// ora entra nella seconda coppia distante d
                
                thetaX[3]= gRandom->Gaus(thetaX[2],sigSI);
                thetaY[3]= gRandom->Gaus(thetaY[2],sigSI);

                x[3]=x[2]+d*tan(thetaX[2])+(1/sqrt(3))*sS*thetaX[3]+gRandom->Gaus(0,ris);
                y[3]=y[2]+d*tan(thetaY[2])+(1/sqrt(3))*sS*thetaY[3];
                
                thetaX[4]= gRandom->Gaus(thetaX[3],sigSI);
                thetaY[4]= gRandom->Gaus(thetaY[3],sigSI);

                x[4]=x[3]+dSS*tan(thetaX[3])+(1/sqrt(3))*sS*thetaX[4];
                y[4]=y[3]+dSS*tan(thetaY[3])+(1/sqrt(3))*sS*thetaY[4]+gRandom->Gaus(0,ris);
                
// terza coppia
                
                thetaX[5]= gRandom->Gaus(thetaX[4],sigSI);
                thetaY[5]= gRandom->Gaus(thetaY[4],sigSI);

                x[5]=x[4]+d*tan(thetaX[4])+(1/sqrt(3))*sS*thetaX[5]+gRandom->Gaus(0,ris);
                y[5]=y[4]+d*tan(thetaY[4])+(1/sqrt(3))*sS*thetaY[5];
                
                thetaX[6]= gRandom->Gaus(thetaX[5],sigSI);
                thetaY[6]= gRandom->Gaus(thetaY[5],sigSI);

                x[6]=x[5]+dSS*tan(thetaX[5])+(1/sqrt(3))*sS*thetaX[6];
                y[6]=y[5]+dSS*tan(thetaY[5])+(1/sqrt(3))*sS*thetaY[6]+gRandom->Gaus(0,ris);
                  
    TMatrixD coo_ang_in(4,7);
    
    for (Int_t i=0; i<7; i++)
    {
        coo_ang_in[0][i]=x[i]; //prima riga ha le coo x
        coo_ang_in[1][i]=y[i]; //seconda riga ha le coo y
        coo_ang_in[2][i]=thetaX[i];
        coo_ang_in[3][i]=thetaY[i];
    }

        return coo_ang_in;
}



TMatrixD MCSout(const PxPyPzEVector & k, const Double_t tar, const Double_t xx, const Double_t yy, const Double_t thX, const Double_t thY)
{   
Double_t const energy   = 150000; //MeV

Double_t const sS    = 0.00032; //m spessore silicio
Double_t const x0S = 0.094; // m

Double_t const sB    = 0.015; //m spessore berillio
Double_t const x0B = 0.353; // m



Double_t const dSS = 0.01; // m distanza tra i due 2S


Double_t const d = 0.25-0.005;
Double_t const dx = 0.25+0.005; // m la distanza tra coppie di silici è 0.25. Però 0.25 è tra i due della coppia, quindi considero -dSS/2=0.005

// m la distanza tra coppie di silici è 0.25. Però 0.25 è tra i due della coppia, quindi considero -dSS/2=0.005
Double_t const ris = 18e-6; // m considero questa la risoluzione dei silici
Double_t const dCAL = 0.10; // m distanza silicio calorimetro
        
Double_t sigSI=(13.6/energy)*sqrt(sS/x0S)*(1+0.038*log(sS/x0S)); //rad
// considero sB/2 per quandp interagisce a metà 
Double_t sigBE2=(13.6/energy)*sqrt(sB/2*x0B)*(1+0.038*log(sB/2*x0B)); //rad
// considero sB/2 per quandp interagisce a metà 
Double_t sigBE=(13.6/energy)*sqrt(sB/2*x0B)*(1+0.038*log(sB/2*x0B)); //rad    

    
    
    TMatrixD thetaX(2,7);
    TMatrixD thetaY(2,7);
    TMatrixD x(2,7);
    TMatrixD y(2,7);
    
    Double_t anglex = atan2(k.Px(), k.Pz());
    Double_t angley = atan2(k.Py(), k.Pz()); 
    
    if(tar==0)
    {
       // siamo nelle stazioni con il target di berillio. Ora entra nel berillio ad una distanza d=0.25-0.005 m dagli ultimi silici. Qui però non puoi trascurare lo spessore del berillio, cioè dove interagisce? considero a metà (7.5 mm), quindi aggiungo a d il pezzo in cui x non è modificato e poi sommo con il nuovo angolo di scattering
                
                Double_t thetaX1= gRandom->Gaus(thX,sigBE2);
                Double_t thetaY1= gRandom->Gaus(thY,sigBE2); 
                anglex += thetaX1;
                angley += thetaY1;
        
                 thetaX[0][0]= gRandom->Gaus(anglex,sigBE2);
                 thetaY[0][0]= gRandom->Gaus(angley,sigBE2); 

                 x[0][0]=xx+d*tan(thX)+(1/sqrt(3))*0.0075*(thetaX1)+(1/sqrt(3))*0.0075*(thetaX[0][0]);
                 y[0][0]=yy+d*tan(thY)+(1/sqrt(3))*0.0075*(thetaY1)+(1/sqrt(3))*0.0075*(thetaY[0][0]);  
              
//entro nelle 3 coppie di silici
for (Int_t p=1; p<7; p++)  {
                 if(p==1 || p==3 || p==5){
                 thetaX[0][p]= gRandom->Gaus(thetaX[0][p-1],sigSI);
                 thetaY[0][p]= gRandom->Gaus(thetaY[0][p-1],sigSI); 

                x[0][p]=x[0][p-1]+d*tan(thetaX[0][p-1])+(1/sqrt(3))*sS*thetaX[0][p]+gRandom->Gaus(0,ris);
                y[0][p]=y[0][p-1]+d*tan(thetaY[0][p-1])+(1/sqrt(3))*sS*thetaY[0][p];
                
                thetaX[0][p+1]= gRandom->Gaus(thetaX[0][p],sigSI);
                thetaY[0][p+1]= gRandom->Gaus(thetaY[0][p],sigSI);

                x[0][p+1]=x[0][p]+dSS*tan(thetaX[0][p])+(1/sqrt(3))*sS*thetaX[0][p+1];
                y[0][p+1]=y[0][p]+dSS*tan(thetaY[0][p])+(1/sqrt(3))*sS*thetaY[0][p+1]+gRandom->Gaus(0,ris);                 
                 }}
        
                 thetaX[1][0]= gRandom->Gaus(thetaX[0][6],sigBE);
                 thetaY[1][0]= gRandom->Gaus(thetaY[0][6],sigBE); 

                 x[1][0]=x[0][6]+d*tan(thetaX[0][6])+(1/sqrt(3))*sB*(thetaX[1][0]);
                 y[1][0]=y[0][6]+d*tan(thetaY[0][6])+(1/sqrt(3))*sB*(thetaY[1][0]);  
   
        //entro nelle 3 coppie di silici
for (Int_t p=1; p<7; p++)  {
                 if(p==1 || p==3 || p==5){
                 thetaX[1][p]= gRandom->Gaus(thetaX[1][p-1],sigSI);
                 thetaY[1][p]= gRandom->Gaus(thetaY[1][p-1],sigSI); 

                x[1][p]=x[1][p-1]+d*tan(thetaX[1][p-1])+(1/sqrt(3))*sS*thetaX[1][p]+gRandom->Gaus(0,ris);
                y[1][p]=y[1][p-1]+d*tan(thetaY[1][p-1])+(1/sqrt(3))*sS*thetaY[1][p];
                
                thetaX[1][p+1]= gRandom->Gaus(thetaX[1][p],sigSI);
                thetaY[1][p+1]= gRandom->Gaus(thetaY[1][p],sigSI);

                x[1][p+1]=x[1][p]+dSS*tan(thetaX[1][p])+(1/sqrt(3))*sS*thetaX[1][p+1];
                y[1][p+1]=y[1][p]+dSS*tan(thetaY[1][p])+(1/sqrt(3))*sS*thetaY[1][p+1]+gRandom->Gaus(0,ris);                 
                 }}
    
    }
     
    
    if(tar==1)
    {
        // siamo nelle stazioni con il target di berillio. Ora entra nel berillio ad una distanza d=0.25-0.005 m dagli ultimi silici. Qui però non puoi trascurare lo spessore del berillio, cioè dove interagisce? considero a metà (7.5 mm), quindi aggiungo a d il pezzo in cui x non è modificato e poi sommo con il nuovo angolo di scattering
                

        
                 thetaX[0][0]= gRandom->Gaus(thX,sigBE);
                 thetaY[0][0]= gRandom->Gaus(thY,sigBE); 

                 x[0][0]=xx+d*tan(thX)+(1/sqrt(3))*sB*(thetaX[0][0]);
                 y[0][0]=yy+d*tan(thY)+(1/sqrt(3))*sB*(thetaY[0][0]);
       
//entro nelle 3 coppie di silici
for (Int_t p=1; p<7; p++)  {
                 if(p==1 || p==3 || p==5){
                 thetaX[0][p]= gRandom->Gaus(thetaX[0][p-1],sigSI);
                 thetaY[0][p]= gRandom->Gaus(thetaY[0][p-1],sigSI); 

                x[0][p]=x[0][p-1]+d*tan(thetaX[0][p-1])+(1/sqrt(3))*sS*thetaX[0][p]+gRandom->Gaus(0,ris);
                y[0][p]=y[0][p-1]+d*tan(thetaY[0][p-1])+(1/sqrt(3))*sS*thetaY[0][p];
                
                thetaX[0][p+1]= gRandom->Gaus(thetaX[0][p],sigSI);
                thetaY[0][p+1]= gRandom->Gaus(thetaY[0][p],sigSI);

                x[0][p+1]=x[0][p]+dSS*tan(thetaX[0][p])+(1/sqrt(3))*sS*thetaX[0][p+1];
                y[0][p+1]=y[0][p]+dSS*tan(thetaY[0][p])+(1/sqrt(3))*sS*thetaY[0][p+1]+gRandom->Gaus(0,ris);                 
                 }}
                                 
                Double_t thetaX1= gRandom->Gaus(thetaX[0][6],sigBE2);
                Double_t thetaY1= gRandom->Gaus(thetaY[0][6],sigBE2); 
        
                anglex += thetaX1;
                angley += thetaY1;
                
                thetaX[1][0]=gRandom->Gaus(anglex,sigBE2);
                thetaY[1][0]=gRandom->Gaus(angley,sigBE2);

                 x[1][0]=x[0][6]+d*tan(thetaX[0][6])+(1/sqrt(3))*0.0075*(thetaX1)+(1/sqrt(3))*0.0075*(thetaX[1][0]); 
        
                 y[1][0]=y[0][6]+d*tan(thetaY[0][6])+(1/sqrt(3))*0.0075*(thetaY1)+(1/sqrt(3))*0.0075*(thetaY[1][0]);
   
        //entro nelle 3 coppie di silici
for (Int_t p=1; p<7; p++)  {
                 if(p==1 || p==3 || p==5){
                 thetaX[1][p]= gRandom->Gaus(thetaX[1][p-1],sigSI);
                 thetaY[1][p]= gRandom->Gaus(thetaY[1][p-1],sigSI); 

                x[1][p]=x[1][p-1]+d*tan(thetaX[1][p-1])+(1/sqrt(3))*sS*thetaX[1][p]+gRandom->Gaus(0,ris);
                y[1][p]=y[1][p-1]+d*tan(thetaY[1][p-1])+(1/sqrt(3))*sS*thetaY[1][p];
                
                thetaX[1][p+1]= gRandom->Gaus(thetaX[1][p],sigSI);
                thetaY[1][p+1]= gRandom->Gaus(thetaY[1][p],sigSI);

                x[1][p+1]=x[1][p]+dSS*tan(thetaX[1][p])+(1/sqrt(3))*sS*thetaX[1][p+1];
                y[1][p+1]=y[1][p]+dSS*tan(thetaY[1][p])+(1/sqrt(3))*sS*thetaY[1][p+1]+gRandom->Gaus(0,ris);                 
                 }}
        // sul calorimetro
        
    Double_t xf = x[1][6]+dCAL*tan(thetaX[1][6]);
    Double_t yf = y[1][6]+dCAL*tan(thetaY[1][6]);
    
    TMatrixD coo_ang_fin(9,7);
        
                for (Int_t i=0; i<7; i++)
    {
        coo_ang_fin[0][i]=x[0][i]; 
        coo_ang_fin[1][i]=y[0][i]; 
        coo_ang_fin[2][i]=thetaX[0][i]; 
        coo_ang_fin[3][i]=thetaY[0][i]; 
                    
        coo_ang_fin[4][i]=x[1][i]; 
        coo_ang_fin[5][i]=y[1][i]; 
        coo_ang_fin[6][i]=thetaX[1][i]; 
        coo_ang_fin[7][i]=thetaY[1][i];         
    }
        
        coo_ang_fin[8][0]=xf;   
        coo_ang_fin[8][1]=yf;  
        coo_ang_fin[8][2]=thetaX[1][6];      
        coo_ang_fin[8][3]=thetaY[1][6]; 
        coo_ang_fin[8][4]=tar;      
        coo_ang_fin[8][5]=0;      
        coo_ang_fin[8][6]=0;      
    
return coo_ang_fin;
        
    }

}
*/

void provaMat()
{

Double_t ppx=0;
Double_t ppy=0;
Double_t ppz=150;
Double_t Ein=200;
    
Double_t px=10;
Double_t py=30;
Double_t pz=110;

    
PxPyPzEVector p_in(ppx,ppy,ppz,Ein);
PxPyPzEVector p_out(px,py,pz,Ein);
    
    
PxPyPzEVector p_in_div=RotDiv(p_in);
    
 
Double_t theta=p_in_div.Theta();
Double_t phi=p_in_div.Phi();

/*
Int_t tar=gRandom->Integer(2);  
TMatrixD coo_in=MCSin(p_in_div);
TMatrixD coo_f=MCSout(p_out,tar,coo_in[0][6],coo_in[1][6],coo_in[2][6],coo_in[3][6]);
    
 
    for (Int_t i=0; i<7; i++)
    {
         cout  << " x1= " <<
        coo_f[0][i] << " y1= " <<
        coo_f[1][i] << " thX1= " <<
        coo_f[2][i] << " thX1= " <<
        coo_f[3][i] << " x2= " <<
                    
        coo_f[4][i] <<  " y2= " <<
        coo_f[5][i] << " thX2= " <<
        coo_f[6][i] <<  " thY2 " <<
        coo_f[7][i] <<endl;       
    }
 cout << " xf MC= " << coo_f[8][0]<< " yf MC= " << coo_f[8][1] << " ThXf MC= " << coo_f[8][2] << " ThYf MC= " << coo_f[8][3] << " tarf= " << coo_f[8][4] <<endl;  
}
    */
    
    
    
}
    
    
    
    
    
