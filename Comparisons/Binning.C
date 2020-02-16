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

std::vector<double> GetBinning(std::vector<double> data, double err){
  TH1D* hist = new TH1D("temp_hist", "", 100, 0, 2);
  for(auto const& d : data){
    hist->Fill(d);
  }
  std::vector<double> binning = ChangeBinning(hist, 2, err);
  delete hist;
  return binning;
}

TH2D* Covariance(std::vector<TH1D*> universes){

  TH1D *mean_syst = (TH1D*) universes[0]->Clone();
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
    mean_syst->SetBinContent(n, mean);
    mean_syst->SetBinError(n, std_dev);
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

// Main
void Binning(){

  // Read in fake data
  std::vector<double> data;
  TFile data_file("../Trees/xsectree_v3rwt.root");
  TTreeReader tree_reader("XSecTree/interaction", &data_file);
  TTreeReaderValue<double> vtx_x(tree_reader, "vtx_x");
  TTreeReaderValue<double> vtx_y(tree_reader, "vtx_y");
  TTreeReaderValue<double> vtx_z(tree_reader, "vtx_z");
  TTreeReaderValue<bool>   particles_contained(tree_reader, "reco_particles_contained");
  TTreeReaderValue<int>    nu_pdg(tree_reader, "reco_nu_pdg");
  TTreeReaderValue<double> lep_mom(tree_reader, "reco_lep_mom");
  while(tree_reader.Next()){
    if(*nu_pdg != 14) continue;
    if(!(*particles_contained)) continue;
    if(!InFiducial(*vtx_x, *vtx_y, *vtx_z)) continue;
    data.push_back(*lep_mom);
  }

  // Get fake data POT from tree
  TTreeReader pot_reader("XSecTree/metadata", &data_file);
  TTreeReaderValue<double> pot_val(pot_reader, "pot");
  double data_pot = 0;
  while(pot_reader.Next()){
    data_pot += *pot_val;
  }

  // Bin fake data in 1D according to percentage statistical error per bin
  std::vector<double> errs = {0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02};
  // Create an array of 1D histograms
  std::vector<TH1D*> hists;
  std::vector<std::vector<double>> edges;
  for(size_t i = 0; i < errs.size(); i++){
    std::vector<double> bin_edges = GetBinning(data, errs[i]);
    double edges_array[bin_edges.size()];
    std::copy(bin_edges.begin(), bin_edges.end(), edges_array);
    TH1D* hist = new TH1D(Form("hist%i", (int)i), "", bin_edges.size()-1, edges_array);
    hists.push_back(hist);
    edges.push_back(bin_edges);
  }
  
  // Loop over fake data and fill the histograms
  for(auto const& d : data){
    for(size_t i = 0; i < hists.size(); i++){
      hists[i]->Fill(d);
    }
  }

  // Read in MC and fill equivalent histograms
  std::vector<bool> mc_used;
  std::vector<double> mc;
  std::vector<TH1D*> mc_hists;
  for(size_t i = 0; i < hists.size(); i++){
    double edges_array[edges[i].size()];
    std::copy(edges[i].begin(), edges[i].end(), edges_array);
    TH1D* mc_hist = new TH1D(Form("mc%i", (int)i), "", edges[i].size()-1, edges_array);
    mc_hists.push_back(mc_hist);
  }
  
  // Read MC data from file
  TFile mc_file("../Trees/xsectree_v2rwt.root");
  TTreeReader mc_reader("XSecTree/interaction", &mc_file);
  TTreeReaderValue<double> mc_vtx_x(mc_reader, "vtx_x");
  TTreeReaderValue<double> mc_vtx_y(mc_reader, "vtx_y");
  TTreeReaderValue<double> mc_vtx_z(mc_reader, "vtx_z");
  TTreeReaderValue<bool>   mc_particles_contained(mc_reader, "reco_particles_contained");
  TTreeReaderValue<int>    mc_nu_pdg(mc_reader, "reco_nu_pdg");
  TTreeReaderValue<double> mc_lep_mom(mc_reader, "reco_lep_mom");
  while(mc_reader.Next()){
    if(*mc_nu_pdg != 14) {mc_used.push_back(false); continue;}
    if(!(*mc_particles_contained)) {mc_used.push_back(false); continue;}
    if(!InFiducial(*mc_vtx_x, *mc_vtx_y, *mc_vtx_z)) {mc_used.push_back(false); continue;}
    mc_used.push_back(true);
    mc.push_back(*mc_lep_mom);
    for(size_t i = 0; i < hists.size(); i++){
      mc_hists[i]->Fill(*mc_lep_mom);
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

  // Scale MC to data POT 
  for(size_t i = 0; i < hists.size(); i++){
    for(size_t j = 0; j < hists[i]->GetNbinsX(); j++){
      mc_hists[i]->SetBinContent(j, mc_hists[i]->GetBinContent(j)*data_pot/mc_pot);
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
      // Stat errors
      double serror = pow(hists[i]->GetBinError(j), 2.);
      cov->SetBinContent(j, j, perror+serror);
    }
    // Genie and flux reweighting systematics
    std::vector<TH1D*> genie_v;
    std::vector<TH1D*> flux_v;
    for(size_t j = 0; j < 100; j++){
      TH1D* genie_tmp = (TH1D*)hists[i]->Clone("genie");
      genie_v.push_back(genie_tmp);
      TH1D* flux_tmp = (TH1D*)hists[i]->Clone("flux");
      flux_v.push_back(flux_tmp);
    }
    // Detector variation systematics
    std::vector<TH1D*> det_v;
    for(size_t j = 0; j < 50; j++){
      TH1D* det_tmp = (TH1D*)hists[i]->Clone("det");
      det_v.push_back(det_tmp);
    }
    genie.push_back(genie_v);
    flux.push_back(flux_v);
    det.push_back(det_v);
    covariances.push_back(cov);
  }

  // Get the reweighting from file
  TTreeReader weight_reader("XSecTree/weight", &mc_file);
  TTreeReaderArray<double> gw(weight_reader, "genie_weights");
  TTreeReaderArray<double> fw(weight_reader, "flux_weights");
  int index = 0;
  int data_i = 0;
  while(weight_reader.Next()){
    if(!mc_used[index]){ index++; continue;}
    for(size_t i = 0; i < hists.size(); i++){
      for(size_t j = 0; j < 100; j++){
        if(gw[j] > 0 && gw[j] < 100){
          genie[i][j]->Fill(mc[data_i], gw[j]);
        }
        if(fw[j] > 0 && fw[j] < 100){
          flux[i][j]->Fill(mc[data_i], fw[j]);
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
  while(det_reader.Next()){
    if(!InFiducial(*ds_vtx_x, *ds_vtx_y, *ds_vtx_z)) continue;
    for(size_t i = 0; i < hists.size(); i++){
      for(size_t j = 0; j < 50; j++){
        if(ds_nu_pdg[j] != 14) continue;
        if(!(ds_particles_contained[j])) continue;
        det[i][j]->Fill(ds_lep_mom[j]);
      }
    }
  }

  // Scale each varied histogram to the data POT
  for(size_t i = 0; i < hists.size(); i++){
    // Scale hists first
    for(size_t j = 0; j < 100; j++){
      for(size_t k = 1; k <= genie[i][j]->GetNbinsX(); k++){
        genie[i][j]->SetBinContent(k, genie[i][j]->GetBinContent(k)*data_pot/mc_pot);
        flux[i][j]->SetBinContent(k, flux[i][j]->GetBinContent(k)*data_pot/mc_pot);
      }
    }
    for(size_t j = 0; j < 50; j++){
      for(size_t k = 1; k <= det[i][j]->GetNbinsX(); k++){
        det[i][j]->SetBinContent(k, det[i][j]->GetBinContent(k)*data_pot/mc_pot);
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

  // Calculate chi2 between data and MC for each histogram
  std::vector<double> chis;
  for(size_t i = 0; i < hists.size(); i++){
    double chi = ChiSquare(hists[i], mc_hists[i], covariances[i]);
    chis.push_back(chi);
    std::cout<<"error = "<<errs[i]<<" chi = "<<chi<<"\n";

    // Plot all the histogram variations
    TCanvas *c1 = new TCanvas(Form("canvas%i", (int)i), "", 900, 600);
    TH1D* error_hist = (TH1D*)hists[i]->Clone();
    error_hist->Reset();
    for(size_t j = 1; j <= hists[i]->GetNbinsX(); j++){
      error_hist->SetBinContent(j, mc_hists[i]->GetBinContent(j));
      error_hist->SetBinError(j, sqrt(covariances[i]->GetBinContent(j, j)-pow(hists[i]->GetBinError(j),2)));
    }
    error_hist->SetLineWidth(0);
    error_hist->SetMarkerStyle(0);
    error_hist->SetFillColor(15);
    error_hist->SetFillStyle(3001);
    error_hist->Scale(1., "width");
    error_hist->GetXaxis()->SetTitle("P [GeV]");
    error_hist->GetYaxis()->SetTitle("Events (/bin width)");
    error_hist->Draw("E2");
    mc_hists[i]->SetLineColor(42);
    mc_hists[i]->SetMarkerSize(0);
    mc_hists[i]->Scale(1., "width");
    mc_hists[i]->Draw("HIST SAME");
    hists[i]->SetLineColor(46);
    hists[i]->SetMarkerSize(0);
    hists[i]->Scale(1., "width");
    hists[i]->Draw("HIST E1 SAME");
    TLegend* legend = new TLegend(0.7, 0.73, 0.92, 0.89);
    legend->SetFillStyle(0);
    legend->AddEntry(hists[i], "Fake data", "l");
    legend->AddEntry(mc_hists[i], "Simulation", "l");
    legend->Draw();
    c1->SaveAs(Form("Plots/hist_%i.png", (int)i));
  }

  // Plot chi2 vs percentage statistical error per bin
  TCanvas *canvas = new TCanvas("canvas", "", 900, 600);
  canvas->SetRightMargin(0.14);
  TGraph *graph = new TGraph(chis.size(), &errs[0], &chis[0]);
  graph->GetXaxis()->SetTitle("Maximum % error on bin");
  graph->GetYaxis()->SetTitle("#chi^{2}");
  graph->SetMarkerColor(46);
  graph->Draw("AP");
  canvas->SaveAs("Plots/binning.png");

}
