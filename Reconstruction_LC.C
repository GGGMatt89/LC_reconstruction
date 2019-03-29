#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TFrame.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "TNtuple.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "Riostream.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TApplication.h"
#include <TGraph.h>

#include "ReconstructionCC.h"


void Reconstruction_LC()
{

  Int_t nfile = 0;
  Double_t Nion = 100*pow(10.,6);
  Double_t edep1,edep2,x1,y1,z1,x2,y2,z2,time1,time2;
  //Donn�es geometrie de la camera
  Double_t d0,d1,d2,l1,e1,l2,e2;	//en mm
  d0=200;//100;
  d1=10;
  d2=400;
  l1=90;//80;
  e1=2;
  l2=380;//300;
  e2=30;//25;
  Double_t c=3*pow(10.,8)*pow(10.,3)/pow(10.,9);	// en mm.ns-1
  Double_t* y = NULL;
  y=new Double_t[2];
  Int_t nbin=340;//68;
  Int_t limitym_histo_inf=-180;
  Int_t limitym_histo_sup=160;

  TH1* ym_histo = new TH1F("ym","",nbin,limitym_histo_inf,limitym_histo_sup);
  double pce=3.65/1000;           // pair creation energy in kev
  double F=0.115;                 // Fano factor
  double ENC=600;                 // Equivalent Noise Charge
  std::ifstream file;

  std::string line, line1;
  double eg,ee,deg,dee;
  int x =0;
  Bool_t new_event = kFALSE;
  std:string name;
  std::ostringstream ff;
  std::string fr;
  std::string path = "/Users/fontana/PhD/Work/Simulations/Geant4/ComptonCamera_Simulation/CC_Hadrontherapy_Monitoring/reconstruction/LC_reconstruction/data/7.5_10_8_events/30realizations/";
  TFile *f;

  for(Int_t fl = 0; fl<30; fl ++){
  ym_histo->Reset();
  ff.str("");
  ff.clear();
  ff << fl;
  fr = ff.str();
  name = path + "reconstruction" + fr + "_MLEM.txt";

  file.open(name);
  while (!file.eof()){
    file>>line>>line1;
    //std::cout<<line<<" "<<line1<<std::endl;
    //new_event = (line.compare(0,19,"---Nouvel evenement")==0);
    file >> line ; file >> line ; file >> x1 >> y1 >> z1;
    file >> line ; file >> line ; file >> edep1;    // initial values in MeV
    //std::cout<<x1<<"\t"<<y1<<"\t"<<z1<<"\t"<<edep1<<std::endl;
    file >> line ; file >> line ; file >> x2 >> y2 >> z2;
    file >> line ; file >> line ; file >> edep2;
    //std::cout<<x2<<"\t"<<y2<<"\t"<<z2<<"\t"<<edep2<<std::endl;

    dee=pce*sqrt(ENC*ENC+F*edep1/pce);     // These MarieHelene: Eg_FWHM=2.355*formuleFano
    deg=(8./100.)/2.35*sqrt(edep2);

    y=ReconstructionCC(x1,y1,z1,x2,y2,z2,edep1,edep2,0,0);
    //std::cout<<"Reconstruction ok"<<endl;

    if(y!=NULL)
    {
      //std::cout<<"Reconstructed positions : "<<y[0]<<"\t"<<y[1]<<std::endl;
      if(y[0] >-300 && y[0] <300 || y[1] >-300 && y[1] <300)
      {
        ym_histo->Fill(y[0],1);
        ym_histo->Fill(y[1],1);
      }
    }else{continue;}

    x1=0.;
    y1=0.;
    z1=0.;
    edep1=0.;
    x2=0.;
    y2=0.;
    z2=0.;
    edep2=0.;
    y = NULL;
  }

  ym_histo->Scale(1/Nion);
  ym_histo->SetLineColor(kBlack);
  ym_histo->SetLineWidth(2);
  ym_histo->SetLineStyle(2);
  f = new TFile(Form("/Users/fontana/PhD/Work/Simulations/Geant4/ComptonCamera_Simulation/CC_Hadrontherapy_Monitoring/reconstruction/LC_reconstruction/results/7.5_10_8_events/reconstructionLC_%d.root", fl+70),"RECREATE");
  ym_histo->Write("reconstructionLC");
  f->Close();
  file.close();
}
  return ;
}
