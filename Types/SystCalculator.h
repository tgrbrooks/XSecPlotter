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

  SystCalculator(Configuration *c, HistManager *h, DataManager *d, size_t f)
  {
    config = c;
    histman = h;
    dataman = d;
    file_i = f;

    // Calculate the systematics for each histogram
    // 1D histograms
    for(size_t i = 0; i < config->plot_variables.size(); i++){
      TString key = config->plot_variables[i];
      // Get the reweighting systematics
      std::pair<std::vector<TH1D*>, std::vector<TH1D*>> rw_systs = GetReweightSysts(i);
      Systematics* genie = new Systematics(rw_systs.first);
      Systematics* flux = new Systematics(rw_systs.second);
      std::vector<TH1D*> det_systs = GetDetSysts(i);
      Systematics* detector = new Systematics(det_systs);
      TH1D* bkg_syst = GetBkgSyst(i);
      Systematics* background = new Systematics(bkg_syst);
      TH1D* const_syst = GetConstSyst(i, config->constant_syst);
      Systematics* constant = new Systematics(const_syst);
      SystSummary* syst = new SystSummary(genie, flux, detector, background, constant);
      histman->Set1DSystematics(key, syst);
    }

    // 2D histograms

  }

  TH1D* GetConstSyst(size_t var_i, double err){

    TH1D *syst_hist = (TH1D*) histman->Get1DHist(config->plot_variables[var_i])->Clone();
    // Loop over the actual binning
    for(size_t n = 1; n <= syst_hist->GetNbinsX(); n++){
      syst_hist->SetBinError(n, err*syst_hist->GetBinContent(n));
    }

    return syst_hist;

  }

  double BkgSubtractionError(double mid, double width, TH1D* bkg){

    int bin = bkg->GetXaxis()->FindBin(mid);
    double bin_width = bkg->GetBinWidth(bin);
    double scale = (width/bin_width) * (config->pot_scale/6.6e20);
    double percent_error = bkg->GetBinError(bin)/bkg->GetBinContent(bin);
    double subtracted = bkg->GetBinContent(bin) * scale;
    double subtraction_error = subtracted * percent_error;
    if(std::isnan(subtraction_error)) subtraction_error = 0.;

    return subtraction_error;
  }

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

    // Loop over the actual binning
    for(size_t n = 1; n <= syst_hist->GetNbinsX(); n++){
      double mid = syst_hist->GetBinCenter(n);
      double width = syst_hist->GetBinWidth(n);
      // Determine the plotting variable
      if(config->plot_variables[var_i]=="lep_mom"){
        // Determine the bin of the background template
        double cos_sub_err = BkgSubtractionError(mid, width, hMomCos);
        double dirt_sub_err = BkgSubtractionError(mid, width, hMomDirt);
        syst_hist->SetBinError(n, std::sqrt(std::pow(cos_sub_err, 2.)+std::pow(dirt_sub_err, 2.)));
      }
      else if(config->plot_variables[var_i]=="cos_lep_theta"){
        double cos_sub_err = BkgSubtractionError(mid, width, hCosThetaCos);
        double dirt_sub_err = BkgSubtractionError(mid, width, hCosThetaDirt);
        syst_hist->SetBinError(n, std::sqrt(std::pow(cos_sub_err, 2.)+std::pow(dirt_sub_err, 2.)));
      }
    }

    return syst_hist;

  }

  std::vector<TH1D*> GetDetSysts(size_t var_i){

    int nsims = 50;

    std::vector<double> bin_edges = histman->Get1DBinning(config->plot_variables[var_i]); 
    double edges_array[bin_edges.size()];
    std::copy(bin_edges.begin(), bin_edges.end(), edges_array);

    std::vector<TH1D*> detsyst_hists;
    for(size_t ns = 0; ns < nsims; ns++){
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
        if(config->plot_variables[var_i] == "lep_mom") detsyst_hists[ns]->Fill(lep_mom[ns]);
        if(config->plot_variables[var_i] == "lep_theta") detsyst_hists[ns]->Fill(lep_theta[ns]);
        if(config->plot_variables[var_i] == "cos_lep_theta") detsyst_hists[ns]->Fill(cos(lep_theta[ns]));
      }

    }

    ScaleUniverses(detsyst_hists, var_i);

    return detsyst_hists;

  }

  // Get the total (unstacked) histogram
  std::pair<std::vector<TH1D*>, std::vector<TH1D*>> GetReweightSysts(size_t var_i){
      
    int nsims = 100;

    std::vector<double> bin_edges = histman->Get1DBinning(config->plot_variables[var_i]); 
    double edges_array[bin_edges.size()];
    std::copy(bin_edges.begin(), bin_edges.end(), edges_array);

    std::vector<TH1D*> geniesyst_hists;
    std::vector<TH1D*> fluxsyst_hists;
    for(size_t ns = 0; ns < nsims; ns++){
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
        }
        if(flux_weight[ns] > 0 && flux_weight[ns] < 100){
          fluxsyst_hists[ns]->Fill(dataman->total_data[data_i][var_i], flux_weight[ns]);
        }
      }
      index++;
      data_i++;
    }

    ScaleUniverses(geniesyst_hists, var_i);
    ScaleUniverses(fluxsyst_hists, var_i);

    return std::make_pair(geniesyst_hists, fluxsyst_hists);

  }

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


};

#endif
