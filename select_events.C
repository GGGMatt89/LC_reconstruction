#include "Riostream.h"
#include "string.h"
#include "sstream"
#include "sys/stat.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMatrixD.h"
#include "TMatrix.h"
#include "TMath.h"
#include "TApplication.h"
#include "TRandom.h"

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

void select_events(string filename, int seed){
  TString out_file;
  int count = 0;
  string line;
  ifstream filein(filename);
  ofstream output;
  while(!filein.eof()){
    getline(filein, line);
    count++;
    //if(count%100==0){cout<<count<<" lines read "<<endl;}

  }
  cout<<"Total number of events for the MLEM reconstruction for 10^10 incident protons "<<count/5<<endl;
  filein.clear();
  filein.seekg(0, filein.beg);
  cout<<count<<endl;
  long proton_N = 100000000;
  cout<<proton_N<<endl;
  Double_t normal = 0.;
  normal = (double)(5e10/proton_N);
  TRandom *generator = new TRandom(seed);
  Int_t number_events = 0;
  int extraction = 0;
  int ev_counter = 0;
  int x = 0;

  for(Int_t ll = 0; ll<70; ll++){

    ev_counter = 0;
    number_events = 0;
    x=0;

    out_file = Form("./results/10_8_events/reconstruction%d_MLEM.txt", ll);
    output.open(out_file);
    number_events = generator->Gaus((double)(count/normal), 10.);

    //cout<<number_events<<endl;



  while(ev_counter<number_events){
    cout<<"Entered in first while "<<endl;
    filein.clear();
    filein.seekg(0, filein.beg);
    while((ev_counter<number_events)&&(!filein.eof())){
      extraction = generator->Integer(11);
      //cout<<extraction<<endl;
      for(x=0; x<5; x++){
        getline(filein, line);
        if(extraction>=5){output<<line<<endl;}
      }
      if(extraction>=5)ev_counter++;
    }
  }

  cout<<"Number of events in new file "<<ev_counter<<endl;
  output.close();

  }

  cout<<"You can start the MLEM reconstruction!"<<endl;
  filein.close();
  return;
}
