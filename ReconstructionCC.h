//fonction pour reconstruire les donnees de la CC
#include <TMatrixD.h>
#include <TMatrix.h>

Double_t* ReconstructionCC(Double_t x1, Double_t y1, Double_t z1,Double_t x2, Double_t y2, Double_t z2, Double_t edep1, Double_t edep2, Double_t E0_in, Double_t test)
{

///-------------------oooooOOOOO00000OOOOOooooo---------------------#
//                                                                 #
//      DEFINITION DU TYPE DE DONNEES 				               #
//       	 			   			                               #
//                                                                 #
//-------------------oooooOOOOO00000OOOOOooooo---------------------#

 Double_t costheta;
Double_t E0=0;
 Double_t E0_sum;
 Double_t mec2=0.511;
Double_t E0_ini =E0_in;
Double_t test_reconstruction=2;// 2 pour une valeur au hasard.
	Double_t egality=0;

	 E0_sum=edep1+edep2;
     test_reconstruction= test;
	if(E0_ini == E0_sum) egality =1;


	if(test_reconstruction ==0) E0= E0_sum;
	if(test_reconstruction ==1) E0 = E0_ini;

	//cout<< "E0 "<<E0<<" test " <<test<< " E0_sum "<<E0_sum<< " E0_ini "<<E0_ini<<endl;

 if(E0==0 ||edep2==0)
 {
 	cout<<"Erreur pendant la reconstruction, E0 "<<E0<<" edep1 "<<edep1<<" edep2 "<<edep2<<" pos1 "<<x1<<" "<<y1<<" "<<z1<<" pos2 "<<x2<<" "<<y2<<" "<<z2<<endl;
	return NULL;
 }
 else costheta=1-mec2*(1/edep2-1/E0);

//cout<<"Cos theta calculation ok"<<endl;


 //hodoscope
 //
 TMatrixD ph(3,1);
 ph(0,0)=0;
 ph(1,0)=1;
 ph(2,0)=0;
 TMatrixD dh(3,1);
 dh(0,0)=0;
 dh(1,0)=1;
 dh(2,0)=0;
 TMatrixD dht(1,3);
 dht.Transpose(dh);

 //vertex et axe du cone
 //
 TMatrixD pv(3,1);
 pv(0,0)=x1;
 pv(1,0)=y1;
 pv(2,0)=z1;
// cout<<pv(0,0)<<" "<<pv(1,0)<<" "<<pv(2,0)<<endl;

 TMatrixD diff(3,1); diff=ph-pv;
 //cout<<"diff "<<diff(0,0)<<" "<<diff(1,0)<<" "<<diff(2,0)<<endl;

 TMatrixD difft(1,3);
 difft.Transpose(diff);
 //cout<<"diff T"<<difft(0,0)<<" "<<difft(0,1)<<" "<<difft(0,2)<<endl;

 Double_t c1,c2,c3;
 c1=(x1-x2)/sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
 c2=(y1-y2)/sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
 c3=(z1-z2)/sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
 TMatrixD c(3,1);
 c(0,0)=c1;
 c(1,0)=c2;
 c(2,0)=c3;

 //Matrice
 //
 TMatrixD M(3,3);
 M(0,0)=pow(c1,2)-pow(costheta,2);
 M(0,1)=c1*c2;
 M(0,2)=c1*c3;
 M(1,0)=c2*c1;
 M(1,1)=pow(c2,2)-pow(costheta,2);
 M(1,2)=c2*c3;
 M(2,0)=c3*c1;
 M(2,1)=c3*c2;
 M(2,2)=pow(c3,2)-pow(costheta,2);

 /*M(0,0)=1;
 M(0,1)=0;
 M(0,2)=0;
 M(1,0)=0;
 M(1,1)=1;
 M(1,2)=0;
 M(2,0)=0;
 M(2,1)=0;
 M(2,2)=1;*/

 //coeffs de l'equation du 2nd degre
 TMatrix A0(difft,TMatrix::kMult,M);
 TMatrix A01(A0,TMatrix::kMult,diff);
 Double_t a0=A01(0,0);

 TMatrix A1(dht,TMatrix::kMult,M);
 TMatrix A11(A1,TMatrix::kMult,diff);
 Double_t a1=A11(0,0);

 TMatrix A2(dht,TMatrix::kMult,M);
 TMatrix A21(A2,TMatrix::kMult,dh);
 Double_t a2=A21(0,0);

 //cout<<a0<<" "<<a1<<" "<<a2<<endl;

//-------------------oooooOOOOO00000OOOOOooooo---------------------#
//                                                                 #
//   		   DEBUT DU PROGRAMME D'ANALYSE			   #
//                                                                 #
//-------------------oooooOOOOO00000OOOOOooooo---------------------#

 Double_t Delta=pow(a1,2)-a0*a2;


 Double_t tplus,tmoins;
 Double_t* y = NULL;


 if(a2!=0 && Delta>0)
 {
 	tplus=(-a1+sqrt(Delta))/a2;
	tmoins=(-a1-sqrt(Delta))/a2;
	//cout<<"2 solutions "<<" "<<tplus<<" "<<tmoins<<endl;



		y=new Double_t[2];

	 //if(test_reconstruction==1 && egality == 1)
	// {
		y[0]=tplus;
		y[1]=tmoins;
	//}
	// else
	// {
		// y[0]=y[1]=300;
	// }


 }
 //else cout<<"Reconstruction impossible"<<endl;




 if(y==NULL)
 {
 	//cout <<"costheta "<<costheta<<" a0 "<<a0<<" a1 "<<a1<<" a2 "<<a2<<" tplus "<<tplus<<" tmoins "<<tmoins<<" Delta "<<Delta<<endl;
	//cout<<" x1 "<<x1<<" x2 "<<x2<<" y1 "<<y1<<" y2 "<<y2<<" z1 "<<z1<<" z2 "<<z2<<"edep1 "<<edep1<<" edep2 "<<edep2<<endl;
}


 return y;

 }
