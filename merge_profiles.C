#include <TCanvas.h>

#include <TROOT.h>
#include "TH1.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include <math.h>
#include <Riostream.h>


void merge_profiles(){

TFile *out_root;
TFile *in_root_1;

out_root = new TFile("./results/profile_complete_5_10_9_events.root", "RECREATE");
TH1F profiles[100];

  for(Int_t x = 0; x<100; x++){
    in_root_1 = new TFile(Form("./results/5_10_9_events/reconstructionLC_%d.root", x), "READ");
    if(in_root_1->IsOpen()){
      std::cout<<"File number "<<x<<" correctly opened! "<<std::endl;
      profiles[x] = (TH1F*)in_root_1->Get("reconstructionLC");
      out_root->cd();
      profiles[x]->Write(Form("TrueGamma_%d", x));
      std::cout<<"Profile "<<x<<" correctly copied to the merged file! "<<std::endl;
      in_root_1->Close();
    }
  }
out_root->Close();
}
