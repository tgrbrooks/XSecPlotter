#ifndef SYSTCALCULATOR_H
#define SYSTCALCULATOR_H

#include "Configuration.h"
#include "HistManager.h"
#include "HistSystematics.h"
#include "DataManager.h"
#include "Systematics.h"
#include "SystHolder.h"
#include "Selection.h"

// Structure for holding interaction information
class SystCalculator
{
  public:

  Configuration *config;
  HistManager *histman;
  DataManager *dataman;

  size_t file_i;

  SystCalculator(Configuration *c, HistManager *h, DataManager *d, size_t f)
  {
    config = c;
    histman = h;
    dataman = d;
    file_i = f;

    /*HistSystematics* genie = new HistSystematics(histman, "genie");
    HistSystematics* flux = new HistSystematics(histman, "flux");
    HistSystematics* detector = new HistSystematics(histman, "detector");
    HistSystematics* background = new HistSystematics(histman, "background");
    HistSystematics* constant = new HistSystematics(histman, "constant");*/

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

    if(calc_genie || calc_flux) GetReweightSysts(calc_genie, calc_flux);
    if(calc_detector) GetDetectorSysts();
    if(calc_background) GetBackgroundSysts();
    if(calc_constant) GetConstantSysts();

    // Calculate the systematics for each histogram
    // Total
    histman->total->systematics->GetTotal();
    histman->total->PrintSummary();
    // 2D histograms
    for(auto const& kv : histman->histos_1D){
      histman->histos_1D[kv.first]->systematics->GetTotal();
    }
    // 2D histograms
    for(auto const& kv : histman->histos_2D){
      histman->histos_2D[kv.first]->systematics->GetTotal();
    }
/*
    // Calculate the systematics for each histogram
    // Total
    std::cout<<"Here\n"; 
    std::cout<<genie->total->mean_syst->GetName()<<"\n";
    SystSummary* syst = new SystSummary(genie->total, flux->total, detector->total, background->total, constant->total);
    std::cout<<"1\n"; 
    histman->SetTotalSystematics(syst);
    histman->total->systematics->GetTotal();
    // 2D histograms
    for(auto const& kv : histman->histos_1D){
      SystSummary* syst1D = new SystSummary(genie->h1D[kv.first], flux->h1D[kv.first], detector->h1D[kv.first], background->h1D[kv.first], constant->h1D[kv.first]);
      histman->Set1DSystematics(kv.first, syst1D);
      histman->histos_1D[kv.first]->systematics->GetTotal();
    }
    std::cout<<"2\n"; 
    // 2D histograms
    for(auto const& kv : histman->histos_2D){
      SystSummary2D* syst2D = new SystSummary2D(genie->h2D[kv.first], flux->h2D[kv.first], detector->h2D[kv.first], background->h2D[kv.first], constant->h2D[kv.first]);
      histman->Set2DSystematics(kv.first.first, kv.first.second, syst2D);
      histman->histos_2D[kv.first]->systematics->GetTotal();
    }
    std::cout<<"3\n"; 
*/
    /*
    // 1D histograms
    for(size_t i = 0; i < config->plot_variables.size(); i++){
      TString key = config->plot_variables[i];
      // Get the reweighting systematics
      TH1D* empty = (TH1D*)histman->Get1DHist(key)->Clone();
      empty->Reset();
      Systematics* no_syst = new Systematics(empty);
      Systematics* genie = no_syst;
      Systematics* flux = no_syst;
      Systematics* detector = no_syst;
      Systematics* background = no_syst;
      Systematics* constant = no_syst;

      bool calc_genie = false;
      bool calc_flux = false;
      bool calc_detector = false;
      bool calc_background = false;
      bool calc_constant = false;
      for(auto const& syst : config->systematics){
        if(syst == "genie") calc_genie = true;
        if(syst == "flux") calc_flux = true;
        if(syst == "detector") calc_detector = true;
        if(syst == "background") calc_background = true;
        if(syst == "constant") calc_constant = true;
      }

      if(calc_genie || calc_flux){
        std::pair<std::vector<TH1D*>, std::vector<TH1D*>> rw_systs = GetReweightSysts(i);
        if(calc_genie) genie = new Systematics(rw_systs.first);
        if(calc_flux) flux = new Systematics(rw_systs.second);
      }
      if(calc_detector){
        std::vector<TH1D*> det_systs = GetDetSysts(i);
        detector = new Systematics(det_systs);
      }
      if(calc_background){
        TH1D* bkg_syst = GetBkgSyst(i);
        background = new Systematics(bkg_syst);
      }
      if(calc_constant){
        TH1D* const_syst = GetConstSyst(i, config->constant_syst);
        constant = new Systematics(const_syst);
      }
      SystSummary* syst = new SystSummary(genie, flux, detector, background, constant, config->systematics);
      histman->Set1DSystematics(key, syst);
    }

    // 2D histograms
    */

  }

  void GetConstantSysts(){

    // Total systematics
    double err = config->constant_syst;
    for(size_t n = 1; n <= histman->total->total_hist->GetNbinsX(); n++){
      histman->total->systematics->constant->mean_syst->SetBinError(n, err*histman->total->total_hist->GetBinContent(n));
    }

    // 1D systematics
    for(auto const& kv1D : histman->histos_1D){
      for(size_t n = 1; n <= kv1D.second->total_hist->GetNbinsX(); n++){
        histman->histos_1D[kv1D.first]->systematics->constant->mean_syst->SetBinError(n, err*kv1D.second->total_hist->GetBinContent(n));
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
    
  }
/*
  TH1D* GetConstSyst(size_t var_i, double err){

    TH1D *syst_hist = (TH1D*) histman->Get1DHist(config->plot_variables[var_i])->Clone();
    // Loop over the actual binning
    for(size_t n = 1; n <= syst_hist->GetNbinsX(); n++){
      syst_hist->SetBinError(n, err*syst_hist->GetBinContent(n));
    }

    return syst_hist;

  }
*/
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

  void GetBackgroundSysts(){

    // Background templates should all be scales to 6.6e20
    TFile *bkg_file = new TFile("BackgroundTemplates.root");
    TH1D* hMomCos = (TH1D*)bkg_file->Get("hMomCosErr");
    TH1D* hCosThetaCos = (TH1D*)bkg_file->Get("hCosThetaCosErr");
    TH1D* hMomDirt = (TH1D*)bkg_file->Get("hMomDirtErr");
    TH1D* hCosThetaDirt = (TH1D*)bkg_file->Get("hCosThetaDirtErr");

    // Total systematics
    double tot_dirt_esq = std::pow(TotalBkgError(hMomDirt), 2);
    double tot_cos_esq = std::pow(TotalBkgError(hMomCos), 2);
    double total_err = std::sqrt(tot_dirt_esq + tot_cos_esq);
    for(size_t n = 1; n <= histman->total->total_hist->GetNbinsX(); n++){
      histman->total->systematics->background->mean_syst->SetBinError(n, total_err);
    }

    // 1D systematics
    for(auto const& kv1D : histman->histos_1D){
      for(size_t n = 1; n <= kv1D.second->total_hist->GetNbinsX(); n++){
        double mid = kv1D.second->total_hist->GetBinCenter(n);
        double width = kv1D.second->total_hist->GetBinWidth(n);
        // Determine the plotting variable
        if(kv1D.first=="lep_mom"){
          // Determine the bin of the background template
          double cos_sub_err = BkgSubtractionError(mid, width, hMomCos);
          double dirt_sub_err = BkgSubtractionError(mid, width, hMomDirt);
          double tot_err = std::sqrt(std::pow(cos_sub_err, 2.)+std::pow(dirt_sub_err, 2.));
          histman->histos_1D[kv1D.first]->systematics->background->mean_syst->SetBinError(n, tot_err);
        }
        else if(kv1D.first=="cos_lep_theta"){
          double cos_sub_err = BkgSubtractionError(mid, width, hCosThetaCos);
          double dirt_sub_err = BkgSubtractionError(mid, width, hCosThetaDirt);
          double tot_err = std::sqrt(std::pow(cos_sub_err, 2.)+std::pow(dirt_sub_err, 2.));
          histman->histos_1D[kv1D.first]->systematics->background->mean_syst->SetBinError(n, tot_err);
        }
        else{
          histman->histos_1D[kv1D.first]->systematics->background->mean_syst->SetBinError(n, 0.01*kv1D.second->total_hist->GetBinContent(n));
        }
      }
    }

    // 2D systematics
    for(auto const& kv2D : histman->histos_2D){
      for(size_t i = 1; i <= kv2D.second->total_hist->GetNbinsX(); i++){
        for(size_t j = 1; j <= kv2D.second->total_hist->GetNbinsY(); j++){
          histman->histos_2D[kv2D.first]->systematics->background->mean_syst->SetBinError(i, j, 0.01*kv2D.second->total_hist->GetBinContent(i,j));
        }
      }
    }

  }

/*
  TH1D* GetBkgSyst(size_t var_i){

    std::vector<double> bin_edges = histman->Get1DBinning(config->plot_variables[var_i]); 
    double edges_array[bin_edges.size()];
    std::copy(bin_edges.begin(), bin_edges.end(), edges_array);
    TH1D *syst_hist = new TH1D(config->plot_variables[var_i]+"_bkgsyst", "", bin_edges.size()-1, edges_array);

    // Background templates should all be scales to 6.6e20
    TFile *bkg_file = new TFile("BackgroundTemplates.root");
    TH1D* hMomCos = (TH1D*)bkg_file->Get("hMomCosErr");
    TH1D* hCosThetaCos = (TH1D*)bkg_file->Get("hCosThetaCosErr");
    TH1D* hMomDirt = (TH1D*)bkg_file->Get("hMomDirtErr");
    TH1D* hCosThetaDirt = (TH1D*)bkg_file->Get("hCosThetaDirtErr");

    double total;

    // Loop over the actual binning
    for(size_t n = 1; n <= syst_hist->GetNbinsX(); n++){
      double mid = syst_hist->GetBinCenter(n);
      double width = syst_hist->GetBinWidth(n);
      // Determine the plotting variable
      if(config->plot_variables[var_i]=="lep_mom"){
        // Determine the bin of the background template
        double cos_sub_err = BkgSubtractionError(mid, width, hMomCos);
        double dirt_sub_err = BkgSubtractionError(mid, width, hMomDirt);
        double tot_err = std::sqrt(std::pow(cos_sub_err, 2.)+std::pow(dirt_sub_err, 2.));
        syst_hist->SetBinError(n, tot_err);
        total += tot_err;
      }
      else if(config->plot_variables[var_i]=="cos_lep_theta"){
        double cos_sub_err = BkgSubtractionError(mid, width, hCosThetaCos);
        double dirt_sub_err = BkgSubtractionError(mid, width, hCosThetaDirt);
        double tot_err = std::sqrt(std::pow(cos_sub_err, 2.)+std::pow(dirt_sub_err, 2.));
        syst_hist->SetBinError(n, tot_err);
        total += tot_err;
      }
    }

    std::cout<<"Total background error = "<<total<<"\n";

    return syst_hist;

  }
*/

  void GetDetectorSysts(){

    int nsims = 50;

    //SystHolder *holder = new SystHolder(histman, nsims);
    histman->CreateUniverses("detector", nsims);

    Selection sel(config);

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

    while (tree_reader.Next()) {

      if(!sel.InFiducial(*vtx_x, *vtx_y, *vtx_z)) continue;

      for(size_t ns = 0; ns < nsims; ns++){
        if(!sel.IsSelected(nu_pdg[ns], cc[ns], lep_contained[ns], particles_contained[ns], 0, 0, 0, 0)) continue;
        //holder->Fill(1., ns);
        histman->total->systematics->detector->universes[ns]->Fill(1.);
        for(auto& kv : histman->histos_1D){
          //if(kv.first == "lep_mom") holder->Fill1D(kv.first, lep_mom[ns], 1., ns);
          if(kv.first == "lep_mom") kv.second->systematics->detector->universes[ns]->Fill(lep_mom[ns]);
          if(kv.first == "lep_theta") kv.second->systematics->detector->universes[ns]->Fill(lep_theta[ns]);
          if(kv.first == "cos_lep_theta") kv.second->systematics->detector->universes[ns]->Fill(cos(lep_theta[ns]));
          //if(kv.first == "lep_theta") holder->Fill1D(kv.first, lep_theta[ns], 1., ns);
          //if(kv.first == "cos_lep_theta") holder->Fill1D(kv.first, cos(lep_theta[ns]), 1., ns);
        }
        for(auto& kv : histman->histos_2D){
          if(kv.first.first == "lep_mom"){
            if(kv.first.second == "lep_theta") kv.second->systematics->detector->universes[ns]->Fill(lep_mom[ns], lep_theta[ns]);
            if(kv.first.second == "cos_lep_theta") kv.second->systematics->detector->universes[ns]->Fill(lep_mom[ns], cos(lep_theta[ns]));
            //if(kv.first.second == "lep_theta") holder->Fill2D(kv.first.first, lep_mom[ns], kv.first.second, lep_theta[ns], 1., ns);
            //if(kv.first.second == "cos_lep_theta") holder->Fill2D(kv.first.first, lep_mom[ns], kv.first.second, cos(lep_theta[ns]), 1., ns);
          }
          if(kv.first.first == "lep_theta"){
            if(kv.first.second == "lep_mom") kv.second->systematics->detector->universes[ns]->Fill(lep_theta[ns], lep_mom[ns]);
            if(kv.first.second == "cos_lep_theta") kv.second->systematics->detector->universes[ns]->Fill(lep_theta[ns], cos(lep_theta[ns]));
            //if(kv.first.second == "lep_mom") holder->Fill2D(kv.first.first, lep_theta[ns], kv.first.second, lep_mom[ns], 1., ns);
            //if(kv.first.second == "cos_lep_theta") holder->Fill2D(kv.first.first, lep_theta[ns], kv.first.second, cos(lep_theta[ns]), 1., ns);
          }
          if(kv.first.first == "cos_lep_theta"){
            if(kv.first.second == "lep_theta") kv.second->systematics->detector->universes[ns]->Fill(cos(lep_theta[ns]), lep_theta[ns]);
            if(kv.first.second == "lep_mom") kv.second->systematics->detector->universes[ns]->Fill(cos(lep_mom[ns]), lep_mom[ns]);
            //if(kv.first.second == "lep_theta") holder->Fill2D(kv.first.first, cos(lep_theta[ns]), kv.first.second, lep_theta[ns], 1., ns);
            //if(kv.first.second == "lep_mom") holder->Fill2D(kv.first.first, cos(lep_theta[ns]), kv.first.second, lep_mom[ns], 1., ns);
          }
        }
      }

    }

    //holder->Scale(config, file_i);
    histman->ScaleUniverses("detector");
    histman->CalculateSyst("detector");

/*
    // Total systematics
    systs->total->Recalculate(holder->total);
    // 1D systematics
    for(auto const& kv1D : systs->h1D){
      systs->h1D[kv1D.first]->Recalculate(holder->hists_1D[kv1D.first]);
    }
    // 2D systematics
    for(auto const& kv2D : systs->h2D){
      systs->h2D[kv2D.first]->Recalculate(holder->hists_2D[kv2D.first]);
    }
    */

  }

/*
  std::vector<TH1D*> GetDetSysts(size_t var_i){

    int nsims = 50;

    std::vector<double> bin_edges = histman->Get1DBinning(config->plot_variables[var_i]); 
    double edges_array[bin_edges.size()];
    std::copy(bin_edges.begin(), bin_edges.end(), edges_array);

    std::vector<double> total;
    std::vector<TH1D*> detsyst_hists;
    for(size_t ns = 0; ns < nsims; ns++){
      total.push_back(0);
      TH1D *detsyst_hist = new TH1D(Form(config->plot_variables[var_i]+"_detsyst%i",ns), "", bin_edges.size()-1, edges_array);
      detsyst_hists.push_back(detsyst_hist);
    }

    Selection sel(config);

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

    while (tree_reader.Next()) {

      if(!sel.InFiducial(*vtx_x, *vtx_y, *vtx_z)) continue;

      for(size_t ns = 0; ns < nsims; ns++){
        if(!sel.IsSelected(nu_pdg[ns], cc[ns], lep_contained[ns], particles_contained[ns], 0, 0, 0, 0)) continue;
        total[ns]++;
        if(config->plot_variables[var_i] == "lep_mom") detsyst_hists[ns]->Fill(lep_mom[ns]);
        if(config->plot_variables[var_i] == "lep_theta") detsyst_hists[ns]->Fill(lep_theta[ns]);
        if(config->plot_variables[var_i] == "cos_lep_theta") detsyst_hists[ns]->Fill(cos(lep_theta[ns]));
      }

    }

    GetTotalRateSyst(total);
    ScaleUniverses(detsyst_hists, var_i);

    return detsyst_hists;

  }
*/

  // Get the total (unstacked) histogram
  void GetReweightSysts(bool calc_genie, bool calc_flux){
    int nsims = 100;

    //SystHolder genie_holder(histman, nsims);
    //SystHolder flux_holder(histman, nsims);
    histman->CreateUniverses("genie", nsims);
    histman->CreateUniverses("flux", nsims);

    TFile data_file(config->input_file[file_i]);

    //Read in TTree
    TTreeReader tree_reader(config->weight_path, &data_file);
    TTreeReaderArray<double> genie_weight(tree_reader, "genie_weights");
    TTreeReaderArray<double> flux_weight(tree_reader, "flux_weights");
    
    int index = 0;
    int data_i = 0;
    while (tree_reader.Next()) {
      if(!dataman->data_used[index]){ index++; continue; }
      for(size_t ns = 0; ns < nsims; ns++){
        // Total
        if(genie_weight[ns] > 0 && genie_weight[ns] < 100 && calc_genie){
          //genie_holder.Fill(genie_weight[ns], ns);
          histman->total->systematics->genie->universes[ns]->Fill(1., genie_weight[ns]);
        }
        if(flux_weight[ns] > 0 && flux_weight[ns] < 100 && calc_flux){
          //flux_holder.Fill(flux_weight[ns], ns);
          histman->total->systematics->flux->universes[ns]->Fill(1., flux_weight[ns]);
        }
        // 1D
        for(size_t i = 0; i < config->plot_variables.size(); i++){
          TString key = config->plot_variables[i];
          if(genie_weight[ns] > 0 && genie_weight[ns] < 100 && calc_genie){
            //genie_holder.Fill1D(key, dataman->total_data[data_i][i], genie_weight[ns], ns);
            histman->histos_1D[key]->systematics->genie->universes[ns]->Fill(dataman->total_data[data_i][i], genie_weight[ns]);
          }
          if(flux_weight[ns] > 0 && flux_weight[ns] < 100 && calc_flux){
            //flux_holder.Fill1D(key, dataman->total_data[data_i][i], flux_weight[ns], ns);
            histman->histos_1D[key]->systematics->flux->universes[ns]->Fill(dataman->total_data[data_i][i], flux_weight[ns]);
          }
          // 2D
          for(size_t j = 0; j < config->plot_variables.size(); j++){
            if(i==j) continue;
            std::pair<TString, TString> key2D = std::make_pair(key, config->plot_variables[j]);
            if(genie_weight[ns] > 0 && genie_weight[ns] < 100 && calc_genie){
              //genie_holder.Fill2D(key2D.first, dataman->total_data[data_i][i], key2D.second, dataman->total_data[data_i][j], genie_weight[ns], ns);
              histman->histos_2D[key2D]->systematics->genie->universes[ns]->Fill(dataman->total_data[data_i][i], dataman->total_data[data_i][j], genie_weight[ns]);
            }
            if(flux_weight[ns] > 0 && flux_weight[ns] < 100 && calc_flux){
              //flux_holder.Fill2D(key2D.first, dataman->total_data[data_i][i], key2D.second, dataman->total_data[data_i][j], flux_weight[ns], ns);
              histman->histos_2D[key2D]->systematics->flux->universes[ns]->Fill(dataman->total_data[data_i][i], dataman->total_data[data_i][j], flux_weight[ns]);
            }
          }
        }
      }
      index++;
      data_i++;
    }

    if(calc_genie){
      histman->ScaleUniverses("genie");
      histman->CalculateSyst("genie");
    }
    if(calc_flux){
      histman->ScaleUniverses("flux");
      histman->CalculateSyst("flux");
    }
    /*
    genie_holder.Scale(config, file_i);
    flux_holder.Scale(config, file_i);

    // Total systematics
    if(calc_genie) geniesysts->total->Recalculate(genie_holder.total);
    if(calc_flux) fluxsysts->total->Recalculate(flux_holder.total);
    // 1D systematics
    for(auto const& kv1D : geniesysts->h1D){
      if(calc_genie) geniesysts->h1D[kv1D.first]->Recalculate(genie_holder.hists_1D[kv1D.first]);
    }
    for(auto const& kv1D : fluxsysts->h1D){
      if(calc_flux) fluxsysts->h1D[kv1D.first]->Recalculate(flux_holder.hists_1D[kv1D.first]);
    }
    // 2D systematics
    for(auto const& kv2D : geniesysts->h2D){
      if(calc_genie) geniesysts->h2D[kv2D.first]->Recalculate(genie_holder.hists_2D[kv2D.first]);
    }
    for(auto const& kv2D : fluxsysts->h2D){
      if(calc_flux) fluxsysts->h2D[kv2D.first]->Recalculate(flux_holder.hists_2D[kv2D.first]);
    }
  */

  }

/*
  // Get the total (unstacked) histogram
  std::pair<std::vector<TH1D*>, std::vector<TH1D*>> GetReweightSysts(size_t var_i){
      
    int nsims = 100;

    std::vector<double> bin_edges = histman->Get1DBinning(config->plot_variables[var_i]); 
    double edges_array[bin_edges.size()];
    std::copy(bin_edges.begin(), bin_edges.end(), edges_array);

    std::vector<double> genie_total;
    std::vector<double> flux_total;
    std::vector<TH1D*> geniesyst_hists;
    std::vector<TH1D*> fluxsyst_hists;
    for(size_t ns = 0; ns < nsims; ns++){
      genie_total.push_back(0);
      flux_total.push_back(0);
      TH1D *geniesyst_hist = new TH1D(Form(config->plot_variables[var_i]+"_geniesyst%i",ns), "", bin_edges.size()-1, edges_array);
      TH1D *fluxsyst_hist = new TH1D(Form(config->plot_variables[var_i]+"_fluxsyst%i",ns), "", bin_edges.size()-1, edges_array);
      geniesyst_hists.push_back(geniesyst_hist);
      fluxsyst_hists.push_back(fluxsyst_hist);
    }

    TFile data_file(config->input_file[file_i]);

    //Read in TTree
    TTreeReader tree_reader(config->weight_path, &data_file);
    TTreeReaderArray<double> genie_weight(tree_reader, "genie_weights");
    TTreeReaderArray<double> flux_weight(tree_reader, "flux_weights");
    
    int index = 0;
    int data_i = 0;
    while (tree_reader.Next()) {
      if(!dataman->data_used[index]){ index++; continue; }
      for(size_t ns = 0; ns < nsims; ns++){
        if(genie_weight[ns] > 0 && genie_weight[ns] < 100){
          geniesyst_hists[ns]->Fill(dataman->total_data[data_i][var_i], genie_weight[ns]);
          genie_total[ns] += genie_weight[ns];
        }
        if(flux_weight[ns] > 0 && flux_weight[ns] < 100){
          fluxsyst_hists[ns]->Fill(dataman->total_data[data_i][var_i], flux_weight[ns]);
          flux_total[ns] += flux_weight[ns];
        }
      }
      index++;
      data_i++;
    }

    ScaleUniverses(geniesyst_hists, var_i);
    ScaleUniverses(fluxsyst_hists, var_i);

    //GetTotalRateSyst(genie_total);
    //GetTotalRateSyst(flux_total);

    return std::make_pair(geniesyst_hists, fluxsyst_hists);

  }
*/
  void ScaleUniverses(std::vector<TH1D*> hists, size_t var_i){
    // Include scale factor
    for(size_t ns = 0; ns < hists.size(); ns++){
      hists[ns]->Scale(config->pot_scale_fac[file_i]);
    }

    // If plotting cross section convert from rate
    if(config->plot_xsec){
      double xsec_scale = 1e38/(config->flux[file_i] * config->targets);
      for(size_t ns = 0; ns < hists.size(); ns++){
        hists[ns]->Scale(xsec_scale, "width");
      }
    }
    // Else if max error used divide each bin by width
    else if (config->max_error > 0 || config->bin_edges[var_i].size()>1){
      for(size_t ns = 0; ns < hists.size(); ns++){
        hists[ns]->Scale(1, "width");
      }
    }
  }

  // Get total systematic errors from variations
  void  GetTotalRateSyst(std::vector<double> total){
    double mean = 0;
    for(size_t ns = 0; ns < total.size(); ns++){
      mean += total[ns];
    }
    mean /= total.size();
    double std_dev = 0;
    for(size_t ns = 0; ns < total.size(); ns++){
      std_dev += std::pow(total[ns] - mean, 2.);
    }
    std_dev = std::sqrt(std_dev/(total.size()-1));
    std::cout<<"Total rate = "<<mean<<" +/- "<<std_dev<<" "<<100*std_dev/mean<<"%\n";
  }


};

#endif
