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
  if(vtx_x < -190 || vtx_x > 190 ||
         vtx_y < -190 || vtx_x > 180 ||
         vtx_z < 15 || vtx_x > 450 ||
         (vtx_x > -5 && vtx_x < 5) ||
         (vtx_z > 247.5 && vtx_z < 252.5)){ 
        return false;
  }
  return true;
}

std::vector<double> ChangeBinning(TH1D* hist, double max, double err){
  std::vector<double> bin_edges;
  for(size_t i = 0; i < hist->GetNbinsX(); i++){
    bin_edges.push_back(hist->GetBinLowEdge(i+1));
  }
  bin_edges.push_back(max);

  for(size_t i = 0; i < hist->GetNbinsX(); i++){
    if(hist->GetBinError(i+1)/hist->GetBinContent(i+1) > err || hist->GetBinContent(i+1) == 0){

      if(i != hist->GetNbinsX()-1){
        bin_edges.erase(bin_edges.begin()+i+1);
        double edges_array[bin_edges.size()];
        std::copy(bin_edges.begin(), bin_edges.end(), edges_array);
        TH1D* new_hist = (TH1D*)hist->Rebin(bin_edges.size()-1, "new", edges_array);
        return ChangeBinning(new_hist, max, err);
      }

      else{
        bin_edges.erase(bin_edges.begin()+(bin_edges.size()-2));
        return bin_edges;
      }

    }
  }
  return bin_edges;
}

std::vector<double> GetBinning(std::vector<double> data, double err, std::string var, double scale, bool twoD=false){
  std::map<std::string, double> mins = {{"mom", 0}, {"theta", -1}, {"phi", -180}, {"vise", 0}, {"ntracks", 0}};
  std::map<std::string, double> maxs = {{"mom", 2}, {"theta", 1}, {"phi", 180}, {"vise", 3}, {"ntracks", 7}};
  std::map<std::string, double> nbins = {{"mom", 50}, {"theta", 50}, {"phi", 50}, {"vise", 50}, {"ntracks", 7}};
  double nb = nbins[var];
  if(twoD && var!="ntracks") nb = nb * 2./5.;
  TH1D* hist = new TH1D("temp_hist", "", nb, mins[var], maxs[var]);
  for(auto const& d : data){
    hist->Fill(d);
  }
  for(size_t i = 1; i <= hist->GetNbinsX(); i++){
    hist->SetBinContent(i, hist->GetBinContent(i)*scale);
  }
  std::vector<double> binning = ChangeBinning(hist, maxs[var], err);
  delete hist;
  return binning;
}

std::vector<std::vector<double>> ChangeBinning2D(TH2D* hist, std::vector<double> xbin_edges, double err){
    
  std::vector<std::vector<double>> all_edges;
  // Loop over slices in one dimension
  for(size_t i = 1; i <= hist->GetNbinsY(); i++){
    TH1D* temp_hist = (TH1D*) hist->ProjectionX("temp_slice", i, i);
    // Rebin slice so errors below maximum
    std::vector<double> bin_edges_new = ChangeBinning(temp_hist, xbin_edges.back(), err);
    // Push new bin edges to vector
    all_edges.push_back(bin_edges_new);
    delete temp_hist;
  }
  return all_edges;
}

std::vector<std::vector<double>> GetBinning2D(std::vector<double> data_i, std::vector<double> data_j, std::vector<double> xbins, std::vector<double> ybins, double err, double scale){

  // Perform two dimensional rebinning
  double xedges_array[xbins.size()];
  std::copy(xbins.begin(), xbins.end(), xedges_array);
  double yedges_array[ybins.size()];
  std::copy(ybins.begin(), ybins.end(), yedges_array);
  // Create 2D histogram with the calculated binning
  TH2D *temp_hist = new TH2D("temp_hist", "", xbins.size()-1, xedges_array, ybins.size()-1, yedges_array);
  // Fill with the data
  for(size_t i = 0; i < data_i.size(); i++){
    temp_hist->Fill(data_i[i], data_j[i]);
  }
  for(size_t i = 1; i <= temp_hist->GetNbinsX(); i++){
    for(size_t j = 1; j <= temp_hist->GetNbinsY(); j++){
      temp_hist->SetBinContent(i, j, temp_hist->GetBinContent(i, j)*scale);
    }
  }
  std::vector<std::vector<double>> all_bin_edges = ChangeBinning2D(temp_hist, xbins, err);
  delete temp_hist;
  
  return all_bin_edges;
}

TH2D* Covariance(std::vector<TH1D*> universes){

  int nbins = universes[0]->GetNbinsX();
  TString name = TString(universes[0]->GetName())+"cov";
  TH2D* covariance = new TH2D(name, "", nbins, 1, nbins+1, nbins, 1, nbins+1);

  std::vector<double> means;
  for(size_t n = 1; n <= nbins; n++){
    double mean = 0;
    for(size_t ns = 0; ns < universes.size(); ns++){
      mean += universes[ns]->GetBinContent(n);
    }
    mean /= universes.size();
    means.push_back(mean);
    double std_dev = 0;
    for(size_t ns = 0; ns < universes.size(); ns++){
      std_dev += std::pow(universes[ns]->GetBinContent(n) - mean, 2.);
    }
    std_dev = std::sqrt(std_dev/(universes.size()-1));
  }

  // Calculate the covariance and correlation over all universes
  for(size_t i = 1; i <= nbins; i++){
    double cv_i = means[i-1];
    for(size_t j = 1; j <= nbins; j++){
      double cv_j = means[j-1];
      double E_ij = 0;
      for(size_t ns = 0; ns < universes.size(); ns++){
        E_ij += (universes[ns]->GetBinContent(i)-cv_i)*(universes[ns]->GetBinContent(j)-cv_j);
      }
      E_ij /= universes.size();
      covariance->SetBinContent(i, j, E_ij);
    }
  }
  return covariance;
}

TH2D* Covariance2D(std::vector<TH2Poly*> universes){

  int nbins = universes[0]->GetNumberOfBins();
  TString name = TString(universes[0]->GetName())+"cov";
  TH2D* covariance = new TH2D(name, "", nbins, 1, nbins+1, nbins, 1, nbins+1);

  std::vector<double> means;
  for(size_t n = 1; n <= nbins; n++){
    double mean = 0;
    for(size_t ns = 0; ns < universes.size(); ns++){
      mean += universes[ns]->GetBinContent(n);
    }
    mean /= universes.size();
    means.push_back(mean);
    double std_dev = 0;
    for(size_t ns = 0; ns < universes.size(); ns++){
      std_dev += std::pow(universes[ns]->GetBinContent(n) - mean, 2.);
    }
    std_dev = std::sqrt(std_dev/(universes.size()-1));
  }

  // Calculate the covariance and correlation over all universes
  for(size_t i = 1; i <= nbins; i++){
    double cv_i = means[i-1];
    for(size_t j = 1; j <= nbins; j++){
      double cv_j = means[j-1];
      double E_ij = 0;
      for(size_t ns = 0; ns < universes.size(); ns++){
        E_ij += (universes[ns]->GetBinContent(i)-cv_i)*(universes[ns]->GetBinContent(j)-cv_j);
      }
      E_ij /= universes.size();
      covariance->SetBinContent(i, j, E_ij);
    }
  }
  return covariance;
}

double ChiSquare(TH1D* data, TH1D* mc, TH2D* covariance){
  int nbins = data->GetNbinsX();
  TMatrix cov;
  cov.Clear();
  cov.ResizeTo(nbins, nbins);

  for (int i = 1; i <= nbins; i ++) {
    for (int j = 1; j <= nbins; j ++) {
      cov[i-1][j-1] = covariance->GetBinContent(i, j);
    }
  }

  TMatrix cov_inv = cov.Invert();

  double chi = 0;
  for(size_t i = 1; i <= nbins; i++){
    double data_i = data->GetBinContent(i);
    double mc_i = mc->GetBinContent(i);
    for(size_t j = 1; j <= nbins; j++){
      double data_j = data->GetBinContent(j);
      double mc_j = mc->GetBinContent(j);
      chi += (data_i-mc_i)*cov_inv[i-1][j-1]*(data_j-mc_j);
    }
  }
  return chi;
}

double ChiSquare2D(TH2Poly* data, TH2Poly* mc, TH2D* covariance){
  int nbins = data->GetNumberOfBins();
  TMatrix cov;
  cov.Clear();
  cov.ResizeTo(nbins, nbins);

  for (int i = 1; i <= nbins; i ++) {
    for (int j = 1; j <= nbins; j ++) {
      cov[i-1][j-1] = covariance->GetBinContent(i, j);
    }
  }

  TMatrix cov_inv = cov.Invert();

  double chi = 0;
  for(size_t i = 1; i <= nbins; i++){
    double data_i = data->GetBinContent(i);
    double mc_i = mc->GetBinContent(i);
    for(size_t j = 1; j <= nbins; j++){
      double data_j = data->GetBinContent(j);
      double mc_j = mc->GetBinContent(j);
      chi += (data_i-mc_i)*cov_inv[i-1][j-1]*(data_j-mc_j);
    }
  }
  return chi;
}

// Main
void Variables(){

  gStyle->SetPalette(kBlueGreenYellow);

  // Read in fake data
  std::vector<double> data_mom;
  std::vector<double> data_theta;
  std::vector<double> data_phi;
  std::vector<double> data_vise;
  std::vector<double> data_ntracks;
  TFile data_file("../Trees/xsectree_v3rwt.root");
  TTreeReader tree_reader("XSecTree/interaction", &data_file);
  TTreeReaderValue<double> vtx_x(tree_reader, "vtx_x");
  TTreeReaderValue<double> vtx_y(tree_reader, "vtx_y");
  TTreeReaderValue<double> vtx_z(tree_reader, "vtx_z");
  TTreeReaderValue<bool>   particles_contained(tree_reader, "reco_particles_contained");
  TTreeReaderValue<int>    nu_pdg(tree_reader, "reco_nu_pdg");
  TTreeReaderValue<double> lep_mom(tree_reader, "reco_lep_mom");
  TTreeReaderValue<double> lep_theta(tree_reader, "reco_lep_theta");
  TTreeReaderValue<double> lep_phi(tree_reader, "reco_lep_phi");
  TTreeReaderValue<double> vise(tree_reader, "reco_nu_energy");
  TTreeReaderValue<unsigned int> npr(tree_reader, "reco_n_pr");
  TTreeReaderValue<unsigned int> npi(tree_reader, "reco_n_pipm");
  while(tree_reader.Next()){
    if(*nu_pdg != 14) continue;
    if(!(*particles_contained)) continue;
    if(!InFiducial(*vtx_x, *vtx_y, *vtx_z)) continue;
    std::vector<double> data_v;
    data_mom.push_back(*lep_mom);
    data_theta.push_back(cos(*lep_theta));
    data_phi.push_back(*lep_phi*180./TMath::Pi());
    data_vise.push_back(std::sqrt(std::pow(*lep_mom, 2) + std::pow(0.10566, 2))+*vise);
    data_ntracks.push_back((double)(*npr + *npi));
  }
  std::vector<std::vector<double>> data = {data_mom, data_theta, data_phi, data_vise, data_ntracks};

  // Get fake data POT from tree
  TTreeReader pot_reader("XSecTree/metadata", &data_file);
  TTreeReaderValue<double> pot_val(pot_reader, "pot");
  double data_pot = 0;
  while(pot_reader.Next()){
    data_pot += *pot_val;
  }

  double data_scale = 6.6e20/data_pot;

  // Bin fake data in 1D according to percentage statistical error per bin
  // Create an array of 1D histograms
  std::vector<std::string> vars = {"mom", "theta", "phi", "vise", "ntracks"};
  std::vector<TString> va = {"P [GeV]", "cos #theta", "#phi [rad]", "E_{visible} [GeV]", "N tracks"};
  std::vector<TString> vs = {"P", "c#theta", "#phi", "E", "NT"};
  std::vector<TH1D*> hists;
  for(size_t i = 0; i < vars.size(); i++){
    std::vector<double> bin_edges = GetBinning(data[i], 0.01, vars[i], data_scale);
    double edges_array[bin_edges.size()];
    std::copy(bin_edges.begin(), bin_edges.end(), edges_array);
    TH1D* hist = new TH1D(vars[i].c_str()+TString("_hist"), "", bin_edges.size()-1, edges_array);
    hists.push_back(hist);
  }

  std::vector<TH2Poly*> hists_2D;
  std::vector<std::pair<size_t, size_t>> ind_2D;
  for(size_t i = 0; i < vars.size(); i++){
    std::vector<double> xbins = GetBinning(data[i], 0.01/3, vars[i], data_scale, true);
    for(size_t j = 0; j < vars.size(); j++){
      if(i==j) continue;
      std::vector<double> ybins = GetBinning(data[j], 0.01/3, vars[j], data_scale, true);;
      std::vector<std::vector<double>> xedges = GetBinning2D(data[i], data[j], xbins, ybins, 0.01, data_scale);
      TH2Poly* hist = new TH2Poly();
      hist->SetName((vars[i]+vars[j]).c_str());
      hist->SetTitle("");
      for(size_t k = 0; k < ybins.size()-1; k++){
        std::vector<double> xbin_edges = xedges[k];
        for(size_t l = 0; l < xbin_edges.size()-1; l++){
          hist->AddBin(xbin_edges[l], ybins[k], xbin_edges[l+1], ybins[k+1]);
        }
      }
      for(size_t n = 0; n < data[i].size(); n++){
        hist->Fill(data[i][n], data[j][n]);
      }
      for(size_t x = 1; x <= hist->GetNumberOfBins(); x++){
        hist->SetBinContent(x, hist->GetBinContent(x)*data_scale);
      }
      hists_2D.push_back(hist);
      ind_2D.push_back(std::make_pair(i, j));
    }
  }
  
  // Loop over fake data and fill the histograms
  for(size_t i = 0; i < hists.size(); i++){
    for(auto const& d : data[i]){
      hists[i]->Fill(d);
    }
    for(size_t x = 1; x <= hists[i]->GetNbinsX(); x++){
      hists[i]->SetBinContent(x, hists[i]->GetBinContent(x)*data_scale);
    }
  }

  // Read in MC and fill equivalent histograms
  std::vector<TH1D*> mc_hists;
  for(size_t i = 0; i < hists.size(); i++){
    TH1D* mc_hist = (TH1D*) hists[i]->Clone();
    mc_hist->Reset();
    mc_hists.push_back(mc_hist);
  }

  std::vector<TH2Poly*> mc_hists_2D;
  for(size_t i = 0; i < hists_2D.size(); i++){
    TH2Poly* mc_hist = (TH2Poly*) hists_2D[i]->Clone();
    mc_hist->ClearBinContents();
    mc_hists_2D.push_back(mc_hist);
  }
  
  // Read MC data from file
  std::vector<bool> mc_used;
  std::vector<double> mc_mom;
  std::vector<double> mc_theta;
  std::vector<double> mc_phi;
  std::vector<double> mc_vise;
  std::vector<double> mc_ntracks;
  TFile mc_file("../Trees/xsectree_v2rwt.root");
  TTreeReader mc_reader("XSecTree/interaction", &mc_file);
  TTreeReaderValue<double> mc_vtx_x(mc_reader, "vtx_x");
  TTreeReaderValue<double> mc_vtx_y(mc_reader, "vtx_y");
  TTreeReaderValue<double> mc_vtx_z(mc_reader, "vtx_z");
  TTreeReaderValue<bool>   mc_particles_contained(mc_reader, "reco_particles_contained");
  TTreeReaderValue<int>    mc_nu_pdg(mc_reader, "reco_nu_pdg");
  TTreeReaderValue<double> mc_lep_mom(mc_reader, "reco_lep_mom");
  TTreeReaderValue<double> mc_lep_theta(mc_reader, "reco_lep_theta");
  TTreeReaderValue<double> mc_lep_phi(mc_reader, "reco_lep_phi");
  TTreeReaderValue<double> mc_rvise(mc_reader, "reco_nu_energy");
  TTreeReaderValue<unsigned int> mc_npr(mc_reader, "reco_n_pr");
  TTreeReaderValue<unsigned int> mc_npi(mc_reader, "reco_n_pipm");
  while(mc_reader.Next()){
    if(*mc_nu_pdg != 14) {mc_used.push_back(false); continue;}
    if(!(*mc_particles_contained)) {mc_used.push_back(false); continue;}
    if(!InFiducial(*mc_vtx_x, *mc_vtx_y, *mc_vtx_z)) {mc_used.push_back(false); continue;}
    mc_used.push_back(true);
    mc_mom.push_back(*mc_lep_mom);
    mc_theta.push_back(cos(*mc_lep_theta));
    mc_phi.push_back(*mc_lep_phi*180./TMath::Pi());
    mc_vise.push_back(std::sqrt(std::pow(*mc_lep_mom, 2) + std::pow(0.10566, 2)) + *mc_rvise);
    mc_ntracks.push_back((double)(*mc_npr + *mc_npi));
  }
  std::vector<std::vector<double>> mc = {mc_mom, mc_theta, mc_phi, mc_vise, mc_ntracks};

  for(size_t n = 0; n < mc[0].size(); n++){
    for(size_t i = 0; i < hists.size(); i++){
      mc_hists[i]->Fill(mc[i][n]);
    }
    for(size_t i = 0; i < hists_2D.size(); i++){
      mc_hists_2D[i]->Fill(mc[ind_2D[i].first][n], mc[ind_2D[i].second][n]);
    }
  }

  // Get MC POT from file
  TTreeReader mc_pot_reader("XSecTree/metadata", &mc_file);
  TTreeReaderValue<double> mc_pot_val(mc_pot_reader, "pot");
  double mc_pot = 0;
  while(mc_pot_reader.Next()){
    mc_pot += *mc_pot_val;
  }
  std::cout<<"data pot = "<<data_pot<<" mc = "<<mc_pot<<" Ratio = "<<data_pot/mc_pot<<"\n";
  double mc_scale = 6.6e20/mc_pot;

  // Scale MC to data POT 
  for(size_t i = 0; i < hists.size(); i++){
    for(size_t j = 0; j <= hists[i]->GetNbinsX(); j++){
      mc_hists[i]->SetBinContent(j, mc_hists[i]->GetBinContent(j)*mc_scale);
    }
  }
  for(size_t i = 0; i < mc_hists_2D.size(); i++){
    for(size_t j = 0; j <= mc_hists_2D[i]->GetNumberOfBins(); j++){
      mc_hists_2D[i]->SetBinContent(j, mc_hists_2D[i]->GetBinContent(j)*mc_scale);
    }
  }
  
  // Calculate systematic errors for the MC
  std::vector<TH2D*> covariances;
  std::vector<std::vector<TH1D*>> genie;
  std::vector<std::vector<TH1D*>> flux;
  std::vector<std::vector<TH1D*>> det;
  for(size_t i = 0; i < hists.size(); i++){
    int nbins = hists[i]->GetNbinsX();
    TH2D* cov = new TH2D(Form("cov%i", (int)i), "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    // Add up all the constant errors
    for(size_t j = 1; j <= nbins; j++){
      // Percent errors
      double perror = pow(mc_hists[i]->GetBinContent(j)*(0.02+0.01+0.01),2.);
      if(vars[i]!="mom" && vars[i]!="theta"){
        perror = pow(mc_hists[i]->GetBinContent(j)*(0.02+0.01+0.01+0.01),2.);
      }
      // Stat errors
      double serror = pow(hists[i]->GetBinError(j), 2.);
      //cov->SetBinContent(j, j, perror+serror);
      cov->SetBinContent(j, j, serror);
    }
    // Genie and flux reweighting systematics
    std::vector<TH1D*> genie_v;
    std::vector<TH1D*> flux_v;
    for(size_t j = 0; j < 100; j++){
      TH1D* genie_tmp = (TH1D*)hists[i]->Clone("genie");
      genie_tmp->Reset();
      genie_v.push_back(genie_tmp);
      TH1D* flux_tmp = (TH1D*)hists[i]->Clone("flux");
      flux_tmp->Reset();
      flux_v.push_back(flux_tmp);
    }
    // Detector variation systematics
    std::vector<TH1D*> det_v;
    for(size_t j = 0; j < 50; j++){
      TH1D* det_tmp = (TH1D*)hists[i]->Clone("det");
      det_tmp->Reset();
      det_v.push_back(det_tmp);
    }
    genie.push_back(genie_v);
    flux.push_back(flux_v);
    det.push_back(det_v);
    covariances.push_back(cov);
  }
  // Calculate systematic errors for the MC
  std::vector<TH2D*> covariances_2D;
  std::vector<std::vector<TH2Poly*>> genie_2D;
  std::vector<std::vector<TH2Poly*>> flux_2D;
  std::vector<std::vector<TH2Poly*>> det_2D;
  for(size_t i = 0; i < hists_2D.size(); i++){
    int nbins = hists_2D[i]->GetNumberOfBins();
    TH2D* cov = new TH2D(Form("cov%i", (int)(i+1)*10), "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    // Add up all the constant errors
    for(size_t j = 1; j <= nbins; j++){
      // Percent errors
      double perror = pow(mc_hists_2D[i]->GetBinContent(j)*(0.02+0.01+0.01),2.);
      if(vars[i]!="mom" && vars[i]!="theta"){
        perror = pow(mc_hists_2D[i]->GetBinContent(j)*(0.02+0.01+0.01+0.01),2.);
      }
      // Stat errors
      double serror = pow(hists_2D[i]->GetBinError(j), 2.);
      //cov->SetBinContent(j, j, perror+serror);
      cov->SetBinContent(j, j, serror);
    }
    // Genie and flux reweighting systematics
    std::vector<TH2Poly*> genie_v;
    std::vector<TH2Poly*> flux_v;
    for(size_t j = 0; j < 100; j++){
      TH2Poly* genie_tmp = (TH2Poly*)hists_2D[i]->Clone("genie2D");
      genie_tmp->ClearBinContents();
      genie_v.push_back(genie_tmp);
      TH2Poly* flux_tmp = (TH2Poly*)hists_2D[i]->Clone("flux2D");
      flux_tmp->ClearBinContents();
      flux_v.push_back(flux_tmp);
    }
    // Detector variation systematics
    std::vector<TH2Poly*> det_v;
    for(size_t j = 0; j < 50; j++){
      TH2Poly* det_tmp = (TH2Poly*)hists_2D[i]->Clone("det2D");
      det_tmp->ClearBinContents();
      det_v.push_back(det_tmp);
    }
    genie_2D.push_back(genie_v);
    flux_2D.push_back(flux_v);
    det_2D.push_back(det_v);
    covariances_2D.push_back(cov);
  }
/*
  // Get the reweighting from file
  TTreeReader weight_reader("XSecTree/weight", &mc_file);
  TTreeReaderArray<double> gw(weight_reader, "genie_weights");
  TTreeReaderArray<double> fw(weight_reader, "flux_weights");
  int index = 0;
  int data_i = 0;
  while(weight_reader.Next()){
    if(!mc_used[index]){ index++; continue;}
    for(size_t j = 0; j < 100; j++){
      for(size_t i = 0; i < hists.size(); i++){
        if(gw[j] > 0 && gw[j] < 100){
          genie[i][j]->Fill(mc[i][data_i], gw[j]);
        }
        if(fw[j] > 0 && fw[j] < 100){
          flux[i][j]->Fill(mc[i][data_i], fw[j]);
        }
      }
      for(size_t i = 0; i < hists_2D.size(); i++){
        if(gw[j] > 0 && gw[j] < 100){
          genie_2D[i][j]->Fill(mc[ind_2D[i].first][data_i], mc[ind_2D[i].second][data_i], gw[j]);
        }
        if(fw[j] > 0 && fw[j] < 100){
          flux_2D[i][j]->Fill(mc[ind_2D[i].first][data_i], mc[ind_2D[i].second][data_i], fw[j]);
        }
      }
    }
    index++;
    data_i++;
  }
  
  // Get the detector variations from file
  TTreeReader det_reader("XSecTree/detsys", &mc_file);
  TTreeReaderValue<double> ds_vtx_x(det_reader, "ds_vtx_x");
  TTreeReaderValue<double> ds_vtx_y(det_reader, "ds_vtx_y");
  TTreeReaderValue<double> ds_vtx_z(det_reader, "ds_vtx_z");
  TTreeReaderArray<bool> ds_particles_contained(det_reader, "ds_particles_contained");
  TTreeReaderArray<int> ds_nu_pdg(det_reader, "ds_nu_pdg");
  TTreeReaderArray<double> ds_lep_mom(det_reader, "ds_lep_mom");
  TTreeReaderArray<double> ds_lep_theta(det_reader, "ds_lep_theta");
  TTreeReaderArray<double> ds_lep_phi(det_reader, "ds_lep_phi");
  TTreeReaderArray<double> ds_vise(det_reader, "ds_vise");
  TTreeReaderArray<unsigned int> ds_ntracks(det_reader, "ds_ntracks");
  while(det_reader.Next()){
    if(!InFiducial(*ds_vtx_x, *ds_vtx_y, *ds_vtx_z)) continue;
    for(size_t j = 0; j < 50; j++){
      if(ds_nu_pdg[j] != 14) continue;
      if(!(ds_particles_contained[j])) continue;
      std::vector<double> ds_data = {ds_lep_mom[j], cos(ds_lep_theta[j]), std::sqrt(std::pow(ds_lep_mom[j], 2) + std::pow(0.10566, 2)) + ds_lep_phi[j], ds_vise[j], (double)ds_ntracks[j]};
      for(size_t i = 0; i < hists.size(); i++){
        det[i][j]->Fill(ds_data[i]);
      }
      for(size_t i = 0; i < hists_2D.size(); i++){
        det_2D[i][j]->Fill(ds_data[ind_2D[i].first], ds_data[ind_2D[i].second]);
      }
    }
  }

  // Scale each varied histogram to the data POT
  for(size_t i = 0; i < mc_hists.size(); i++){
    // Scale hists first
    for(size_t j = 0; j < 100; j++){
      for(size_t k = 1; k <= genie[i][j]->GetNbinsX(); k++){
        genie[i][j]->SetBinContent(k, genie[i][j]->GetBinContent(k)*mc_scale);
        flux[i][j]->SetBinContent(k, flux[i][j]->GetBinContent(k)*mc_scale);
      }
    }
    for(size_t j = 0; j < 50; j++){
      for(size_t k = 1; k <= det[i][j]->GetNbinsX(); k++){
        det[i][j]->SetBinContent(k, det[i][j]->GetBinContent(k)*mc_scale);
      }
    }
    // Calculate covariance matrices and add up
    TH2D* genie_cov = Covariance(genie[i]);
    covariances[i]->Add(genie_cov);
    delete genie_cov;
    TH2D* flux_cov = Covariance(flux[i]);
    covariances[i]->Add(flux_cov);
    delete flux_cov;
    TH2D* det_cov = Covariance(det[i]);
    covariances[i]->Add(det_cov);
    delete det_cov;
  }
  // Scale each varied histogram to the data POT
  for(size_t i = 0; i < hists_2D.size(); i++){
    // Scale hists first
    for(size_t j = 0; j < 100; j++){
      for(size_t k = 1; k <= genie_2D[i][j]->GetNumberOfBins(); k++){
        genie_2D[i][j]->SetBinContent(k, genie_2D[i][j]->GetBinContent(k)*mc_scale);
        flux_2D[i][j]->SetBinContent(k, flux_2D[i][j]->GetBinContent(k)*mc_scale);
      }
    }
    for(size_t j = 0; j < 50; j++){
      for(size_t k = 1; k <= det_2D[i][j]->GetNumberOfBins(); k++){
        det_2D[i][j]->SetBinContent(k, det_2D[i][j]->GetBinContent(k)*mc_scale);
      }
    }
    // Calculate covariance matrices and add up
    TH2D* genie_cov = Covariance2D(genie_2D[i]);
    covariances_2D[i]->Add(genie_cov);
    delete genie_cov;
    TH2D* flux_cov = Covariance2D(flux_2D[i]);
    covariances_2D[i]->Add(flux_cov);
    delete flux_cov;
    TH2D* det_cov = Covariance2D(det_2D[i]);
    covariances_2D[i]->Add(det_cov);
    delete det_cov;
  }
*/
  // Calculate chi2 between data and MC for each histogram
  std::vector<double> chis;
  for(size_t i = 0; i < hists.size(); i++){
    double chi = ChiSquare(hists[i], mc_hists[i], covariances[i]);
    chis.push_back(chi);
    std::cout<<"var = "<<vars[i]<<" chi = "<<chi<<"\n";

    // Plot all the histogram variations
    TCanvas *c1 = new TCanvas(Form("canvas%i", (int)i), "", 900, 600);
    TH1D* error_hist = (TH1D*)mc_hists[i]->Clone();
    for(size_t j = 1; j <= mc_hists[i]->GetNbinsX(); j++){
      error_hist->SetBinError(j, sqrt(covariances[i]->GetBinContent(j, j)-pow(hists[i]->GetBinError(j),2)));
    }
    mc_hists[i]->SetLineColor(42);
    mc_hists[i]->SetMarkerSize(0);
    mc_hists[i]->Scale(1., "width");
    mc_hists[i]->GetXaxis()->SetTitle(va[i]);
    mc_hists[i]->GetYaxis()->SetTitle("Events (/bin width)");
    mc_hists[i]->Draw("HIST");
    error_hist->SetLineWidth(0);
    error_hist->SetMarkerStyle(0);
    error_hist->SetFillColor(15);
    error_hist->SetFillStyle(3001);
    error_hist->Scale(1., "width");
    error_hist->Draw("E2 SAME");
    hists[i]->SetLineColor(46);
    hists[i]->SetMarkerSize(0);
    hists[i]->Scale(1., "width");
    hists[i]->Draw("HIST E1 SAME");
    int maxbin = mc_hists[i]->GetMaximumBin();
    double ymax = mc_hists[i]->GetBinContent(maxbin);
    mc_hists[i]->GetYaxis()->SetRangeUser(0, 1.1*ymax);
    TLegend* legend = new TLegend(0.7, 0.73, 0.92, 0.89);
    legend->SetFillStyle(0);
    legend->AddEntry(hists[i], "Fake data", "l");
    legend->AddEntry(mc_hists[i], "Simulation", "l");
    legend->Draw();
    c1->SaveAs(TString("Plots/")+vars[i].c_str()+TString("_hist.png"));
  }

  std::vector<double> chis_2D;
  for(size_t i = 0; i < hists_2D.size(); i++){
    double chi = ChiSquare2D(hists_2D[i], mc_hists_2D[i], covariances_2D[i]);
    chis_2D.push_back(chi);
    std::cout<<"var1 = "<<vars[ind_2D[i].first]<<" var2 = "<<vars[ind_2D[i].second]<<" chi = "<<chi<<"\n";

    // Plot all the histogram variations
    TCanvas *c1 = new TCanvas(Form("canvas%i", (int)(i+1)*10), "", 900, 600);
    c1->SetRightMargin(0.16);
    //mc_hists_2D[i]->Scale(1., "width");
    for(auto const& obj : *mc_hists_2D[i]->GetBins()){
      TH2PolyBin *bin = (TH2PolyBin*)obj;
      double wy = abs(bin->GetYMax() - bin->GetYMin());
      double wx = abs(bin->GetXMax() - bin->GetXMin());
      double width = wy*wx;
      int j = bin->GetBinNumber();
      mc_hists_2D[i]->SetBinContent(j, mc_hists_2D[i]->GetBinContent(j)/width);
    }
    mc_hists_2D[i]->GetXaxis()->SetTitle(va[ind_2D[i].first]);
    mc_hists_2D[i]->GetYaxis()->SetTitle(va[ind_2D[i].second]);
    mc_hists_2D[i]->Draw("COLZ");
    c1->SaveAs(TString("Plots/")+vars[ind_2D[i].first].c_str()+TString(vars[ind_2D[i].second].c_str())+TString("_hist.png"));
  }

  // Plot chi2 vs percentage statistical error per bin
  TCanvas *canvas = new TCanvas("canvas", "", 900, 600);
  canvas->SetRightMargin(0.12);
  std::vector<double> var_i = {0, 1, 2, 3, 4};
  TGraph *graph = new TGraph(chis.size(), &var_i[0], &chis[0]);
  graph->GetXaxis()->SetTitle("Maximum % error on bin");
  graph->GetYaxis()->SetTitle("#chi^{2}");
  graph->SetMarkerColor(46);
  graph->Draw("AP");
  canvas->SaveAs("Plots/variables.png");

  TH2D *all = new TH2D("all", "", vars.size(), 0, vars.size(), vars.size(), 0, vars.size());
  int pind = 0;
  for(size_t i = 0; i < vars.size(); i++){
    for(size_t j = 0; j < vars.size(); j++){
      if(i==j) all->SetBinContent(i+1, j+1, chis[i]);
      else{
        all->SetBinContent(i+1, j+1, chis_2D[pind]);
        pind++;
      }
    }
  }
  TCanvas *c2 = new TCanvas("c2", "", 900, 900);
  c2->SetRightMargin(0.16);
  c2->SetLeftMargin(0.16);
  c2->SetTopMargin(0.16);
  c2->SetBottomMargin(0.16);
  gStyle->SetPaintTextFormat("4.1f");
  all->SetMarkerSize(1.5);
  all->SetMarkerColor(10);
  all->GetXaxis()->CenterLabels();
  all->GetYaxis()->CenterLabels();
  for(size_t i = 1; i <= all->GetNbinsX(); i++){
    all->GetXaxis()->SetBinLabel(i, vs[i-1]);
    all->GetYaxis()->SetBinLabel(i, vs[i-1]);
  }
  all->Draw("COLZ TEXT");
  c2->SaveAs("Plots/all_variables.png");

}
