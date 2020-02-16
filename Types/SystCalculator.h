#ifndef SYSTCALCULATOR_H
#define SYSTCALCULATOR_H

#include "Configuration.h"
#include "HistManager.h"
#include "DataManager.h"
#include "Systematics.h"
#include "Selection.h"

// Structure for holding interaction information
class SystCalculator
{
  public:

  Configuration *config;
  HistManager *histman;
  DataManager *dataman;

  size_t file_i;

  // Constructor
  SystCalculator(Configuration *c, HistManager *h, DataManager *d, size_t f)
  {
    config = c;
    histman = h;
    dataman = d;
    file_i = f;

    // Work out which systematics to calculate
    bool calc_genie = false;
    bool calc_flux = false;
    bool calc_detector = false;
    bool calc_background = false;
    bool calc_constant = false;
    for(auto const& syst : config->systematics){
      if(syst == "genie") calc_genie = true;
      if(syst == "flux") calc_flux = true;
      if(syst == "detector" && config->stage=="reco") calc_detector = true;
      if(syst == "background") calc_background = true;
      if(syst == "constant") calc_constant = true;
    }

    // Calculate the systematics
    if(calc_genie || calc_flux) GetReweightSysts(calc_genie, calc_flux);
    if(calc_detector) GetDetectorSysts();
    if(calc_background) GetBackgroundSysts();
    if(calc_constant) GetConstantSysts();

    // Calculate the total systematics for each histogram
    // Total
    histman->total->systematics->GetTotal();
    histman->total->PrintSummary();
    // 2D histograms
    for(auto& kv : histman->histos_1D){
      kv.second->systematics->GetTotal();
    }
    // 2D histograms
    for(auto& kv : histman->histos_2D){
      kv.second->systematics->GetTotal();
    }

  }

  // Apply constant systematic errors to each bin
  void GetConstantSysts(){

    // Total systematics
    double err = config->constant_syst;
    for(size_t n = 1; n <= histman->total->total_hist->GetNbinsX(); n++){
      histman->total->systematics->constant->mean_syst->SetBinError(n, err*histman->total->total_hist->GetBinContent(n));
    }

    // 1D systematics
    for(auto& kv1D : histman->histos_1D){
      for(size_t n = 1; n <= kv1D.second->total_hist->GetNbinsX(); n++){
        kv1D.second->systematics->constant->mean_syst->SetBinError(n, err*kv1D.second->total_hist->GetBinContent(n));
      }
    }

    // 2D systematics
    for(auto const& kv2D : histman->histos_2D){
      for(size_t i = 1; i <= kv2D.second->total_hist->GetNbinsX(); i++){
        for(size_t j = 1; j <= kv2D.second->total_hist->GetNbinsY(); j++){
          histman->histos_2D[kv2D.first]->systematics->constant->mean_syst->SetBinError(i, j, err*kv2D.second->total_hist->GetBinContent(i,j));
        }
      }
    }

    histman->CalculateSyst("constant");
    
  }

  // Work out the error on the amount of background subtracted from a bin
  double BkgSubtractionError(double mid, double width, TH1D* bkg){

    int bin = bkg->GetXaxis()->FindBin(mid);
    double bin_width = bkg->GetBinWidth(bin);
    double scale = (width/bin_width) * (config->pot[file_i]*config->pot_scale_fac[file_i]/6.6e20);
    double percent_error = bkg->GetBinError(bin)/bkg->GetBinContent(bin);
    double subtracted = bkg->GetBinContent(bin) * scale;
    double subtraction_error = subtracted * percent_error;
    if(std::isnan(subtraction_error)) subtraction_error = 0.;

    return subtraction_error;
  }

  // Work out total error on amount of background subtracted
  double TotalBkgError(TH1D* bkg){

    double scale = config->pot[file_i]*config->pot_scale_fac[file_i]/6.6e20;
    double total_error = 0;
    double total_content = bkg->IntegralAndError(0, bkg->GetNbinsX()+1, total_error);
    for(size_t i = 0; i <= bkg->GetNbinsX()+1; i++){
      double err = std::pow(bkg->GetBinError(i),2);
      if(std::isnan(err)) err = 0;
      total_error += err;
    }
    double percent_error = std::sqrt(total_error)/total_content;
    double subtracted = total_content * scale;
    double subtraction_error = subtracted * percent_error;
    if(std::isnan(subtraction_error)) subtraction_error = 0.;

    return subtraction_error;
  }

  // Calculate external background systematic uncertainties
  void GetBackgroundSysts(){

    // Background templates should all be scales to 6.6e20
    TFile *bkg_file = new TFile("BackgroundTemplates.root");
    TH1D* hMomCos = (TH1D*)bkg_file->Get("hMomCosErr");
    TH1D* hCosThetaCos = (TH1D*)bkg_file->Get("hCosThetaCosErr");
    TH2D* hMomCosThetaCos = (TH2D*)bkg_file->Get("hMomCosThetaCosErr");
    TH1D* hMomDirt = (TH1D*)bkg_file->Get("hMomDirtErr");
    TH1D* hCosThetaDirt = (TH1D*)bkg_file->Get("hCosThetaDirtErr");
    TH2D* hMomCosThetaDirt = (TH2D*)bkg_file->Get("hMomCosThetaDirtErr");

    // Total systematics
    double tot_dirt_esq = std::pow(TotalBkgError(hMomDirt), 2);
    double tot_cos_esq = std::pow(TotalBkgError(hMomCos), 2);
    double total_err = std::sqrt(tot_dirt_esq + tot_cos_esq);
    for(size_t n = 1; n <= histman->total->total_hist->GetNbinsX(); n++){
      histman->total->systematics->background->mean_syst->SetBinError(n, total_err);
    }

    // 1D systematics
    for(auto& kv1D : histman->histos_1D){
      for(size_t n = 1; n <= kv1D.second->total_hist->GetNbinsX(); n++){
        double mid = kv1D.second->total_hist->GetBinCenter(n);
        double width = kv1D.second->total_hist->GetBinWidth(n);
        // Determine the plotting variable
        if(kv1D.first=="lep_mom"){
          // Determine the bin of the background template
          double cos_sub_err = BkgSubtractionError(mid, width, hMomCos);
          double dirt_sub_err = BkgSubtractionError(mid, width, hMomDirt);
          double tot_err = std::sqrt(std::pow(cos_sub_err, 2.)+std::pow(dirt_sub_err, 2.));
          kv1D.second->systematics->background->mean_syst->SetBinError(n, tot_err);
        }
        else if(kv1D.first=="cos_lep_theta"){
          double cos_sub_err = BkgSubtractionError(mid, width, hCosThetaCos);
          double dirt_sub_err = BkgSubtractionError(mid, width, hCosThetaDirt);
          double tot_err = std::sqrt(std::pow(cos_sub_err, 2.)+std::pow(dirt_sub_err, 2.));
          kv1D.second->systematics->background->mean_syst->SetBinError(n, tot_err);
        }
        else{
          kv1D.second->systematics->background->mean_syst->SetBinError(n, 0.01*kv1D.second->total_hist->GetBinContent(n));
        }
      }
    }

    // 2D systematics TODO Add 2D systematics properly
    for(auto& kv2D : histman->histos_2D){
      for(size_t i = 1; i <= kv2D.second->total_hist->GetNbinsX(); i++){
        for(size_t j = 1; j <= kv2D.second->total_hist->GetNbinsY(); j++){
          kv2D.second->systematics->background->mean_syst->SetBinError(i, j, 0.01*kv2D.second->total_hist->GetBinContent(i,j));
        }
      }
    }

    histman->CalculateSyst("background");

  }

  // Calculate the detector systematic uncertainties
  void GetDetectorSysts(){

    // Number of universes
    int nsims = 50;

    // Create empty histograms for universes
    histman->CreateUniverses("detector", nsims);

    Selection sel(config);

    // Read in variation data from file
    TFile data_file(config->input_file[file_i]);

    TTreeReader tree_reader("XSecTree/detsyst", &data_file);
    TTreeReaderValue<double> vtx_x(tree_reader, "ds_vtx_x");
    TTreeReaderValue<double> vtx_y(tree_reader, "ds_vtx_y");
    TTreeReaderValue<double> vtx_z(tree_reader, "ds_vtx_z");
    TTreeReaderArray<bool> particles_contained(tree_reader, "ds_particles_contained");
    TTreeReaderArray<bool> lep_contained(tree_reader, "ds_lep_contained");
    TTreeReaderArray<int> cc(tree_reader, "ds_cc");
    TTreeReaderArray<int> nu_pdg(tree_reader, "ds_nu_pdg");
    TTreeReaderArray<double> lep_mom(tree_reader, "ds_lep_mom");
    TTreeReaderArray<double> lep_theta(tree_reader, "ds_lep_theta");

    // Loop over the tree entries
    while (tree_reader.Next()) {

      // Apply fiducial volume cut
      if(!sel.InFiducial(*vtx_x, *vtx_y, *vtx_z)) continue;

      // Loop over the universes
      for(size_t ns = 0; ns < nsims; ns++){
        // Apply selection
        if(!sel.IsSelected(nu_pdg[ns], true, lep_contained[ns], particles_contained[ns], -1, -1, -1, -1)) continue;
        // Total
        histman->total->systematics->detector->universes[ns]->Fill(1.);
        // 1D histograms
        for(auto& kv : histman->histos_1D){
          if(kv.first == "lep_mom") kv.second->systematics->detector->universes[ns]->Fill(lep_mom[ns]);
          if(kv.first == "lep_theta") kv.second->systematics->detector->universes[ns]->Fill(lep_theta[ns]);
          if(kv.first == "cos_lep_theta") kv.second->systematics->detector->universes[ns]->Fill(cos(lep_theta[ns]));
        }
        // 2D histograms
        for(auto& kv : histman->histos_2D){
          if(kv.first.first == "lep_mom"){
            if(kv.first.second == "lep_theta") kv.second->systematics->detector->universes[ns]->Fill(lep_mom[ns], lep_theta[ns]);
            if(kv.first.second == "cos_lep_theta") kv.second->systematics->detector->universes[ns]->Fill(lep_mom[ns], cos(lep_theta[ns]));
          }
          if(kv.first.first == "lep_theta"){
            if(kv.first.second == "lep_mom") kv.second->systematics->detector->universes[ns]->Fill(lep_theta[ns], lep_mom[ns]);
            if(kv.first.second == "cos_lep_theta") kv.second->systematics->detector->universes[ns]->Fill(lep_theta[ns], cos(lep_theta[ns]));
          }
          if(kv.first.first == "cos_lep_theta"){
            if(kv.first.second == "lep_theta") kv.second->systematics->detector->universes[ns]->Fill(cos(lep_theta[ns]), lep_theta[ns]);
            if(kv.first.second == "lep_mom") kv.second->systematics->detector->universes[ns]->Fill(cos(lep_mom[ns]), lep_mom[ns]);
          }
        }
      }

    }

    // Apply scale factors
    histman->ScaleUniverses("detector");
    // Calculate means and covariances
    histman->CalculateSyst("detector");

  }

  // Calculate reweighting systematics (genie and flux)
  void GetReweightSysts(bool calc_genie, bool calc_flux){
    int nsims = 100;

    // Create empty histograms for universes
    histman->CreateUniverses("genie", nsims);
    histman->CreateUniverses("flux", nsims);

    // Read in data from file
    TFile data_file(config->input_file[file_i]);

    //Read in TTree
    TTreeReader tree_reader("XSecTree/weight", &data_file);
    TTreeReaderArray<double> genie_weight(tree_reader, "genie_weights");
    TTreeReaderArray<double> flux_weight(tree_reader, "flux_weights");
    
    // Loop over tree
    int index = 0;
    int data_i = 0;
    std::cout<<"This takes...\n";
    while (tree_reader.Next()) {
      // Determine if interaction was selected
      if(!dataman->data_used[index]){ index++; continue; }
      // Loop over the number of universes
      for(size_t ns = 0; ns < nsims; ns++){
        // Total
        if(genie_weight[ns] > 0 && genie_weight[ns] < 100 && calc_genie){
          histman->total->systematics->genie->universes[ns]->Fill(1., genie_weight[ns]);
        }
        if(flux_weight[ns] > 0 && flux_weight[ns] < 100 && calc_flux){
          histman->total->systematics->flux->universes[ns]->Fill(1., flux_weight[ns]);
        }
        // 1D histograms
        for(size_t i = 0; i < config->plot_variables.size(); i++){
          TString key = config->plot_variables[i];
          if(genie_weight[ns] > 0 && genie_weight[ns] < 100 && calc_genie){
            histman->histos_1D[key]->systematics->genie->universes[ns]->Fill(dataman->total_data[data_i][i], genie_weight[ns]);
          }
          if(flux_weight[ns] > 0 && flux_weight[ns] < 100 && calc_flux){
            histman->histos_1D[key]->systematics->flux->universes[ns]->Fill(dataman->total_data[data_i][i], flux_weight[ns]);
          }
          // 2D histograms
          for(size_t j = 0; j < config->plot_variables.size(); j++){
            if(i==j) continue;
            std::pair<TString, TString> key2D = std::make_pair(key, config->plot_variables[j]);
            if(genie_weight[ns] > 0 && genie_weight[ns] < 100 && calc_genie){
              histman->histos_2D[key2D]->systematics->genie->universes[ns]->Fill(dataman->total_data[data_i][i], dataman->total_data[data_i][j], genie_weight[ns]);
            }
            if(flux_weight[ns] > 0 && flux_weight[ns] < 100 && calc_flux){
              histman->histos_2D[key2D]->systematics->flux->universes[ns]->Fill(dataman->total_data[data_i][i], dataman->total_data[data_i][j], flux_weight[ns]);
            }
          }
        }
      }
      index++;
      data_i++;
    }
    std::cout<<"... Aaages\n";

    // Scale universes
    if(calc_genie){
      histman->ScaleUniverses("genie");
      histman->CalculateSyst("genie");
    }
    // Calculate means and covariances
    if(calc_flux){
      histman->ScaleUniverses("flux");
      histman->CalculateSyst("flux");
    }

  }


};

#endif
