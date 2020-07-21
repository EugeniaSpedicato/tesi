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




void provaMat()
{

Double_t ppx=0;
Double_t ppy=0;
Double_t ppz=12342342342;
Double_t Ein=200;
    
PxPyPzEVector p_in(ppx,ppy,ppz,Ein);
    
PxPyPzEVector p_in_div=RotDiv(p_in);
    

 
Double_t theta=p_in_div.Theta();
Double_t phi=p_in_div.Phi();


TMatrixD coo_f=coo(p_in_div,theta,phi);
    
 cout << " x1= " << coo_f[0][0] << " y1= " << coo_f[0][1] << " x2= " << coo_f[1][0] << " y2= " <<coo_f[1][1] <<endl;   
}

    
    
    
    
    
