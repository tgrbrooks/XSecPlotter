#include <string>
#include <vector>
#include <TChain.h>
#include <TTree.h>

#include <map>
#include <cmath>
#include <iostream>
#include <ctime>
#include <algorithm>

#include <TFile.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLine.h>
#include <TVector3.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TH3D.h>
#include <THStack.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCanvas.h>
#include <sstream>
#include <fstream>
#include <math.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

bool InFiducial(double vtx_x, double vtx_y, double vtx_z){
  if(vtx_x < -200 || vtx_x > 200 ||
         vtx_y < -200 || vtx_x > 200 ||
         vtx_z < 0 || vtx_x > 500){ 
        return false;
  }
  return true;
}

// Main
void TotalRate(){

  // Read in fake data
  TFile data_file("../Trees/xsectree_pi.root");
  TTreeReader tree_reader("XSecTree/interaction", &data_file);
  TTreeReaderValue<double> vtx_x(tree_reader, "vtx_x");
  TTreeReaderValue<double> vtx_y(tree_reader, "vtx_y");
  TTreeReaderValue<double> vtx_z(tree_reader, "vtx_z");
  TTreeReaderValue<int>    nu_pdg(tree_reader, "reco_nu_pdg");
  double rate = 0;
  std::vector<bool> used;
  while(tree_reader.Next()){
    if(*nu_pdg != 14){ used.push_back(false); continue; }
    if(!InFiducial(*vtx_x, *vtx_y, *vtx_z)){ used.push_back(false); continue; }
    rate += 1;
    used.push_back(true);
  }

  std::vector<double> rate_rwt;
  for(size_t i = 0; i < 100; i++) rate_rwt.push_back(0);

  // Get the reweighting from file
  TTreeReader weight_reader("XSecTree/weight", &data_file);
  TTreeReaderArray<double> fw(weight_reader, "flux_weights");
  int index = 0;
  while(weight_reader.Next()){
    if(!used[index]){ index++; continue;}
    for(size_t j = 0; j < 100; j++){
      if(fw[j] > 0 && fw[j] < 100){
        rate_rwt[j] += fw[j];
      }
    }
    index++;
  }

  double mean = 0;
  for(size_t ns = 0; ns < rate_rwt.size(); ns++){
    mean += rate_rwt[ns];
  }
  mean /= rate_rwt.size();
  double std_dev = 0;
  for(size_t ns = 0; ns < rate_rwt.size(); ns++){
    std_dev += std::pow(rate_rwt[ns] - mean, 2.);
  }
  std_dev = std::sqrt(std_dev/(rate_rwt.size()-1));

  std::cout<<"Rate = "<<rate<<"\n"
           <<"Flux reweight: mean = "<<mean<<" std dev = "<<std_dev<<"\n"
           <<"% of mean = "<<100*std_dev/mean<<" of rate = "<<100*std_dev/rate<<"\n";
}
