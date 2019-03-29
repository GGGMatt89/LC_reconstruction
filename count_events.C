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

int count_events(string filename){

  int count = 0;
  string line;
  ifstream filein(filename);
  while(!filein.eof()){
    getline(filein, line);
    count++;
  }
  cout<<"Total number of events for the MLEM reconstruction"<<endl;
  return count/5;
  }
