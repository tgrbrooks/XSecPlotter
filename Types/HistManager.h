#ifndef HISTMANAGER_H
#define HISTMANAGER_H

#include "Configuration.h"
#include "Titles.h"
#include "DataManager.h"
#include "BinManager.h"
#include "XSecCalculator.h"
#include "Histo1D.h"
#include "Histo2D.h"

// Structure for holding all of the histograms
class HistManager
{
  public:

  Configuration *config; // Global configuration
  Titles *titles;        // Histogram titles
  DataManager *dataman;  // Stored data
  BinManager *binman;    // Histogram binnings
  XSecCalculator *xsec;  // Cross section calculator functions

  size_t file_i; // File index

  Histo1D* total;                                            // Total histogram
  std::map<TString, Histo1D*> histos_1D;                     // All 1D histograms
  std::map<std::pair<TString, TString>, Histo2D*> histos_2D; // All 2D histograms

  // Constructor
  HistManager(Configuration *c, Titles *t, DataManager *d, BinManager *b, XSecCalculator *x, size_t f)
  {
    config = c;
    titles = t;
    dataman = d;
    binman = b;
    xsec = x;
    file_i = f;

    // Total event rate
    TH1D* hist0D = GetTotal();
    total = new Histo1D(hist0D);
    // Calculate the cross section if plotting
    if(config->plot_xsec){
      SetBackground();
      SetEfficiency();
      SetXSec();
    }

    // 1D histograms
    for(size_t i = 0; i < config->plot_variables.size(); i++){
      TString key1D = config->plot_variables[i];

      // Make the histograms
      TH1D* hist1D = GetTotalHist1D(i);
      std::pair<THStack*, TLegend*> stack1D = StackHist1D(i);

      // Create the histo object
      Histo1D* histo1D = new Histo1D(hist1D, stack1D);
      histos_1D[key1D] = histo1D;

      // Calculate efficiency and purity if plotting
      if(config->plot_xsec) SetBackground(i);
      if(config->plot_eff_pur || config->plot_xsec) SetEfficiency(i);
      if(config->plot_response || config->plot_xsec) SetResponse(i);
      if(config->plot_xsec) SetXSec(i);

      // 2D histograms
      for(size_t j = 0; j < config->plot_variables.size(); j++){
        if(i==j) continue;
        std::pair<TString, TString> key2D = std::make_pair(key1D, config->plot_variables[j]);

        // Make the histograms
        TH2Poly* hist2D = GetTotalHist2D(i, j);
        std::pair<std::vector<THStack*>, TLegend*> stack2D = StackHist2D(i, j);

        // Get the binning
        std::vector<double> ybins = binman->Get1DBinning(config->plot_variables[j]);
        std::vector<std::vector<double>> xbins = binman->Get2DBinning(config->plot_variables[i], config->plot_variables[j]);

        // Create the histo object
        Histo2D* histo2D = new Histo2D(hist2D, stack2D, ybins, xbins);
        histos_2D[key2D] = histo2D;

        // Calculate efficiency and response if plotting
        if(config->plot_xsec) SetBackground(i, j);
        if(config->plot_eff_pur || config->plot_xsec) SetEfficiency(i, j);
        if(config->plot_response || config->plot_xsec) SetResponse(i, j, hist2D);
        if(config->plot_xsec) SetXSec(i, j);
      }
    }
    
  }

  // Accessors for 1D histograms
  TH1D* Get1DHist(TString plot_var){
    return histos_1D[plot_var]->total_hist;
  }
  TH1D* Get1DHist(size_t var_i){
    return histos_1D[config->plot_variables[var_i]]->total_hist;
  }
  Histo1D* GetHisto1D(TString plot_var){
    return histos_1D[plot_var];
  }
  Histo1D* GetHisto1D(size_t var_i){
    return histos_1D[config->plot_variables[var_i]];
  }

  // Accessors for 1D binning
  std::vector<double> Get1DBinning(TString plot_var){
    return binman->Get1DBinning(plot_var);
  }

  // Accessor for 2D histograms
  Histo2D* GetHisto2D(TString plot_var1, TString plot_var2){
    std::pair<TString, TString> key = std::make_pair(plot_var1, plot_var2);
    return histos_2D[key];
  }
  Histo2D* GetHisto2D(size_t var_i, size_t var_j){
    std::pair<TString, TString> key = std::make_pair(config->plot_variables[var_i], config->plot_variables[var_j]);
    return histos_2D[key];
  }

  // Accessors for 2D binning
  std::vector<std::vector<double>> Get2DBinning(TString plot_var1, TString plot_var2){
    return binman->Get2DBinning(plot_var1, plot_var2);
  }

  // Return the number of histograms that would be plotted
  int GetNHists(){
    int n_1D = histos_1D.size();
    int n_2D = histos_2D.size();
    if(n_2D == 0) return n_1D;
    n_2D = 1;
    int n_1D_slice = 0;
    bool first = true;
    for(auto const& kv : histos_2D){
      if(!first) continue;
      n_1D_slice += kv.second->ybins.size();
      first = false;
    }
    return n_1D + n_2D + n_1D_slice;
  }

  // ------------------------------------------------------------------------------------------------------------------
  //                                        HANDLING SYSTEMATIC UNIVERSES
  // ------------------------------------------------------------------------------------------------------------------

  // Create empty histograms for systematic variations
  void CreateUniverses(TString syst, int n){
    total->systematics->GetSyst(syst)->CreateUniverses(n);
    for(auto& kv : histos_1D){
      kv.second->systematics->GetSyst(syst)->CreateUniverses(n);
    }
    for(auto& kv : histos_2D){
      kv.second->systematics->GetSyst(syst)->CreateUniverses(n);
    }
  }

  // Scale each universe by pot, xsec, bin width where appropriate
  void ScaleUniverses(TString syst){
    total->systematics->GetSyst(syst)->ScaleUniverses(config, file_i);
    for(auto& kv : histos_1D){
      kv.second->systematics->GetSyst(syst)->ScaleUniverses(config, file_i);
    }
    for(auto& kv : histos_2D){
      kv.second->systematics->GetSyst(syst)->ScaleUniverses(config, file_i);
    }
  }

  // Calculate mean, covariance and correlation for universe variations
  void CalculateSyst(TString syst){
    if(config->plot_xsec){
      total->systematics->GetSyst(syst)->Calculate(total->xsec_hist);
    }
    else{
      total->systematics->GetSyst(syst)->Calculate(total->total_hist);
    }
    for(auto& kv : histos_1D){
      if(config->plot_xsec){
        kv.second->systematics->GetSyst(syst)->Calculate(kv.second->xsec_hist);
      }
      else{
        kv.second->systematics->GetSyst(syst)->Calculate(kv.second->total_hist);
      }
    }
    for(auto& kv : histos_2D){
      if(config->plot_xsec){
        kv.second->systematics->GetSyst(syst)->Calculate(kv.second->xsec_hist);
      }
      else{
        kv.second->systematics->GetSyst(syst)->Calculate(kv.second->total_hist);
      }
    }
  }

  // ------------------------------------------------------------------------------------------------------------------
  //                                      HANDLING XSEC SYSTEMATIC UNIVERSES
  // ------------------------------------------------------------------------------------------------------------------
  
  // Create cross section universes
  void CreateXSecUni(TString syst, int n){
    total->systematics->GetSyst(syst)->CreateXSecUni(n);
    for(auto& kv : histos_1D){
      kv.second->systematics->GetSyst(syst)->CreateXSecUni(n);
    }
    for(auto& kv : histos_2D){
      kv.second->systematics->GetSyst(syst)->CreateXSecUni(n);
    }
  }

  // Fill cross section universes with reweighting systematic variations for all histograms
  void FillXSecUni(TString syst, int index, int uni, double weight){

    // Is the true interaction part of selection
    bool is_true = dataman->interactions[index].true_selected;
    // Is the reconstructed interaction selected
    bool sel = dataman->interactions[index].selected;

    // Fill histograms needed for folded cross section calculation
    // Total
    if(is_true){
      // Fill true generated histograms
      total->systematics->GetSyst(syst)->xsecuni[uni]->generated->Fill(1., weight);
    }
    if(sel){
      // Fill true selected histograms
      if(is_true){
        total->systematics->GetSyst(syst)->xsecuni[uni]->selected->Fill(1., weight);
      }
      // Fill background histograms
      else{
        total->systematics->GetSyst(syst)->xsecuni[uni]->background->Fill(1., weight);
      }
    }

    // 1D
    for(size_t i = 0; i < config->plot_variables.size(); i++){
      TString k1D = config->plot_variables[i];
      double data_i = dataman->interactions[index].variables[i];
      double true_data_i = dataman->interactions[index].true_variables[i];
      if(is_true){
        // Fill true generated histograms
        histos_1D[k1D]->systematics->GetSyst(syst)->xsecuni[uni]->generated->Fill(true_data_i, weight);
        // Fill migration matrices
        if(data_i != -99999){
          histos_1D[k1D]->systematics->GetSyst(syst)->xsecuni[uni]->migration->Fill(true_data_i, data_i);
        }
      }
      if(sel){
        // Fill true selected histograms
        if(is_true){
          histos_1D[k1D]->systematics->GetSyst(syst)->xsecuni[uni]->selected->Fill(true_data_i, weight);
        }
        // Fill background histograms
        else{
          histos_1D[k1D]->systematics->GetSyst(syst)->xsecuni[uni]->background->Fill(data_i, weight);
        }
      }

      for(size_t j = 0; j < config->plot_variables.size(); j++){
        if(i==j) continue;
        std::pair<TString, TString> k2D = std::make_pair(k1D, config->plot_variables[j]);
        double data_j = dataman->interactions[index].variables[j];
        double true_data_j = dataman->interactions[index].true_variables[j];
        if(is_true){
          // Fill true generated histograms
          histos_2D[k2D]->systematics->GetSyst(syst)->xsecuni[uni]->generated->Fill(true_data_i, true_data_j, weight);
          // Fill migration matrices
          if(data_i != -99999 && data_j != -99999){
            int true_bin = histos_2D[k2D]->total_hist->FindBin(true_data_i, true_data_j);
            int reco_bin = histos_2D[k2D]->total_hist->FindBin(data_i, data_j);
            histos_2D[k2D]->systematics->GetSyst(syst)->xsecuni[uni]->migration->Fill(true_bin-0.5, reco_bin-0.5);
          }
        }
        if(sel){
          // Fill true selected histograms
          if(is_true){
            histos_2D[k2D]->systematics->GetSyst(syst)->xsecuni[uni]->selected->Fill(true_data_i, true_data_j, weight);
          }
          // Fill background histograms
          else{
            histos_2D[k2D]->systematics->GetSyst(syst)->xsecuni[uni]->background->Fill(data_i, data_j, weight);
          }
        }
      }
    }

  }

  // Calculate the cross section from universe variations for all histograms
  void CalculateXSec(TString syst){
    int nuni = total->systematics->GetSyst(syst)->xsecuni.size();
    for(int i = 0; i < nuni; i++){
      TH1D *xsec_total_temp = xsec->ToXSec(total->systematics->GetSyst(syst)->xsecuni[i], total->total_hist, syst, i);
      total->systematics->GetSyst(syst)->universes[i]->Add(xsec_total_temp);
      for(auto& kv : histos_1D){
        TH1D *xsec_1D_temp = xsec->ToXSec(kv.second->systematics->GetSyst(syst)->xsecuni[i], kv.second->total_hist, syst, i);
        kv.second->systematics->GetSyst(syst)->universes[i]->Add(xsec_1D_temp); 
      }
      for(auto& kv : histos_2D){
        TH2Poly *xsec_2D_temp = xsec->ToXSec(kv.second->systematics->GetSyst(syst)->xsecuni[i], kv.second->total_hist, syst, i);
        for(int j = 0; j <= xsec_2D_temp->GetNumberOfBins()+1; j++){
          kv.second->systematics->GetSyst(syst)->universes[i]->SetBinContent(j, xsec_2D_temp->GetBinContent(j)); 
        }
      }
    }
  }

  // ------------------------------------------------------------------------------------------------------------------
  //                                        BACKGROUND HISTOGRAM CREATION
  // ------------------------------------------------------------------------------------------------------------------

  // Calculate the background for 1D histograms
  void SetBackground(){

    for (auto const& in : dataman->interactions){
      // Denominator of efficiency plot is all interactions that are selected in truth
      if(in.selected && !in.true_selected){ 
        total->bkg_hist->Fill(1);
      }
    }
    total->bkg_hist->SetBinContent(1, total->bkg_hist->GetBinContent(1)*config->pot_scale_fac[file_i]);

  }

  // Calculate the background for 1D histograms
  void SetBackground(size_t var_i){

    TString plot_var = config->plot_variables[var_i];

    for (auto const& in : dataman->interactions){
      // Denominator of efficiency plot is all interactions that are selected in truth
      if(in.selected && !in.true_selected){ 
        histos_1D[plot_var]->bkg_hist->Fill(in.true_variables[var_i]);
      }
    }
    for(int n = 0; n <= histos_1D[plot_var]->bkg_hist->GetNbinsX(); n++){
      histos_1D[plot_var]->bkg_hist->SetBinContent(n, histos_1D[plot_var]->bkg_hist->GetBinContent(n)*config->pot_scale_fac[file_i]);
    }
    // If plotting cross section convert from rate
    histos_1D[plot_var]->bkg_hist->Scale(1, "width");

  }

  // Calculate the background for 2D histograms
  void SetBackground(size_t var_i, size_t var_j){

    std::pair<TString, TString> key = std::make_pair(config->plot_variables[var_i], config->plot_variables[var_j]);

    for (auto const& in : dataman->interactions){
      // Denominator of efficiency plot is all interactions that are selected in truth
      if(in.selected && !in.true_selected){ 
        histos_2D[key]->bkg_hist->Fill(in.true_variables[var_i], in.true_variables[var_j]);
      }
    }
    for(int n = 0; n <= histos_2D[key]->bkg_hist->GetNumberOfBins(); n++){
      histos_2D[key]->bkg_hist->SetBinContent(n, histos_2D[key]->bkg_hist->GetBinContent(n)*config->pot_scale_fac[file_i]);
    }
    // If plotting cross section convert from rate
    histos_2D[key]->bkg_hist->Scale(1, "width");

  }

  // ------------------------------------------------------------------------------------------------------------------
  //                                      EFFICIENCY AND PURITY CALCULATION
  // ------------------------------------------------------------------------------------------------------------------

  // Calculate the efficiency and purity for 1D histograms
  void SetEfficiency(){

    for (auto const& in : dataman->interactions){
      // Denominator of efficiency plot is all interactions that are selected in truth
      if(in.true_selected){ 
         total->efficiency.second->Fill(1);
      }
      // Denominator of purity plot is all interaction that are selected after reconstruction
      if(in.selected){ 
        total->purity.second->Fill(1);
      }
      // Numerator of efficiency and purity plots is all interactions that are selected in both truth and reco
      if(in.selected && in.true_selected){ 
        total->efficiency.first->Fill(1);
        total->purity.first->Fill(1);
      }
    }

  }

  // Calculate the efficiency and purity for 1D histograms
  void SetEfficiency(size_t var_i){

    TString plot_var = config->plot_variables[var_i];

    for (auto const& in : dataman->interactions){
      // Denominator of efficiency plot is all interactions that are selected in truth
      if(in.true_selected){ 
        histos_1D[plot_var]->efficiency.second->Fill(in.true_variables[var_i]);
      }
      // Denominator of purity plot is all interaction that are selected after reconstruction
      if(in.selected){ 
        histos_1D[plot_var]->purity.second->Fill(in.variables[var_i]);
      }
      // Numerator of efficiency and purity plots is all interactions that are selected in both truth and reco
      if(in.selected && in.true_selected){ 
        histos_1D[plot_var]->efficiency.first->Fill(in.true_variables[var_i]);
        histos_1D[plot_var]->purity.first->Fill(in.variables[var_i]);
      }
    }

  }

  // Calculate the efficiency and purity for 2D histograms
  void SetEfficiency(size_t var_i, size_t var_j){

    std::pair<TString, TString> key = std::make_pair(config->plot_variables[var_i], config->plot_variables[var_j]);

    for (auto const& in : dataman->interactions){
      // Denominator of efficiency plot is all interactions that are selected in truth
      if(in.true_selected){ 
        histos_2D[key]->efficiency.second->Fill(in.true_variables[var_i], in.true_variables[var_j]);
      }
      // Denominator of purity plot is all interaction that are selected after reconstruction
      if(in.selected){ 
        histos_2D[key]->purity.second->Fill(in.variables[var_i], in.variables[var_j]);
      }
      // Numerator of efficiency and purity plots is all interactions that are selected in both truth and reco
      if(in.selected && in.true_selected){ 
        histos_2D[key]->efficiency.first->Fill(in.true_variables[var_i], in.true_variables[var_j]);
        histos_2D[key]->purity.first->Fill(in.variables[var_i], in.variables[var_j]);
      }
    }

  }

  // ------------------------------------------------------------------------------------------------------------------
  //                                        RESPONSE MATRIX CREATION
  // ------------------------------------------------------------------------------------------------------------------

  // Calculate the response matrix for 1D histograms
  void SetResponse(size_t var_i){

    TString key = config->plot_variables[var_i];
    // Create a histogram and fill it with truth and reco pairs
    // Truth in X and reco in Y
    std::vector<double> bin_edges = Get1DBinning(key);
    double edges_array[bin_edges.size()];
    std::copy(bin_edges.begin(), bin_edges.end(), edges_array);
    TH2D *temp_response = new TH2D("temp_response", "", bin_edges.size()-1, edges_array, bin_edges.size()-1, edges_array);
    int index = 0;
    for(auto const& in : dataman->interactions){
      // Response matrix for all reconstructed events that are selected in truth
      if(in.true_selected && in.variables[var_i] != -99999){
        temp_response->Fill(in.true_variables[var_i], in.variables[var_i]);
      }
    }

    // Create the response matrix 
    // Loop over the true bins from the first to last
    for(size_t bin_j = 1; bin_j <= bin_edges.size(); bin_j++){
      // Calculate the total number of events generated in bin j
      double total = 0;
      // Loop over the reco bins including underflow and overflow
      for(size_t bin_i = 0; bin_i <= bin_edges.size(); bin_i++){
        total += (double)temp_response->GetBinContent(bin_j, bin_i);
      }   
      // Loop over the reco bins again and fill with number of events in reco bin from true bin / total in true bin
      for(size_t bin_i = 1; bin_i <= bin_edges.size(); bin_i++){
        histos_1D[key]->response->SetBinContent(bin_j, bin_i, (double)temp_response->GetBinContent(bin_j, bin_i)/total);
        if(temp_response->GetBinContent(bin_j, bin_i)==0){
          histos_1D[key]->response->SetBinContent(bin_j, bin_i, 0.00001);
        } 
      }   
    } 
    delete temp_response;

  }

  // Calculate the response matrix for 2D histograms
  void SetResponse(size_t var_i, size_t var_j, TH2Poly* hist){

    std::pair<TString, TString> key = std::make_pair(config->plot_variables[var_i], config->plot_variables[var_j]);
    int nbins = hist->GetNumberOfBins();

    TH2D *response = new TH2D("response", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    for(auto const& in : dataman->interactions){
      // Response matrix defined for all reconstructed events match to true event in FV
      int true_bin = hist->FindBin(in.true_variables[var_i], in.true_variables[var_j]);
      int reco_bin = hist->FindBin(in.variables[var_i], in.variables[var_j]);
      if(in.true_selected && in.variables[var_i] != -99999 && in.variables[var_j] != -99999){
        response->Fill(true_bin-0.5, reco_bin-0.5);
      }
    }

    // Create the response matrix 
    // Loop over the true bins from the first to last
    for(int bin_j = 1; bin_j <= nbins; bin_j++){
      // Calculate the total number of events generated in bin j
      double total = 0;
      // Loop over the reco bins including underflow and overflow
      for(int bin_i = 0; bin_i <= nbins; bin_i++){
        total += (double)response->GetBinContent(bin_j, bin_i);
      }   
      // Loop over the reco bins again and fill with number of events in reco bin from true bin / total in true bin
      for(int bin_i = 1; bin_i <= nbins; bin_i++){
        histos_2D[key]->response->SetBinContent(bin_j, bin_i, (double)response->GetBinContent(bin_j, bin_i)/total);
        if(response->GetBinContent(bin_j, bin_i)==0){
          histos_2D[key]->response->SetBinContent(bin_j, bin_i, 0.000001);
        } 
      }   
    } 
    delete response;

  }

  // ------------------------------------------------------------------------------------------------------------------
  //                                        CROSS SECTION CONVERSION
  // ------------------------------------------------------------------------------------------------------------------

  // Calculate the total cross section
  void SetXSec(){
    total->xsec_hist = xsec->ToXSec(total);
    // Set the total systematic histogram
    for(int i = 0; i <= total->xsec_hist->GetNbinsX()+1; i++){
      total->systematics->GetSyst("total")->mean_syst->SetBinContent(i, total->xsec_hist->GetBinContent(i));
    }
  }

  // Calculate the 1D cross section
  void SetXSec(int var_i){
    TString key = config->plot_variables[var_i];
    histos_1D[key]->xsec_hist = xsec->ToXSec(histos_1D[key]);

    // Do a little formatting
    histos_1D[key]->xsec_hist->SetFillColor(config->cols[0]);
    histos_1D[key]->xsec_hist->SetLineColor(config->cols[0]);
    if(file_i != 0){
      histos_1D[key]->xsec_hist->SetFillStyle(config->fsty[file_i]);
      histos_1D[key]->xsec_hist->SetLineStyle(config->lsty[file_i]);
    }
    if(!config->plot_filled){
      histos_1D[key]->xsec_hist->SetFillColor(0);
      histos_1D[key]->xsec_hist->SetLineWidth(3);
    } 

    // Set the total systematic histogram
    for(int i = 0; i <= histos_1D[key]->xsec_hist->GetNbinsX()+1; i++){
      histos_1D[key]->systematics->GetSyst("total")->mean_syst->SetBinContent(i, histos_1D[key]->xsec_hist->GetBinContent(i));
    }
  }

  // Calculate the 2D cross section
  void SetXSec(int var_i, int var_j){
    std::pair<TString, TString> key = std::make_pair(config->plot_variables[var_i], config->plot_variables[var_j]);
    std::pair<TH2Poly*, TH2Poly*> xsec_conv = xsec->ToXSec(histos_2D[key]);
    histos_2D[key]->xsec_hist = xsec_conv.first;
    histos_2D[key]->xsec_err = xsec_conv.second;

    // Do a little formatting
    histos_2D[key]->xsec_hist->SetFillColor(config->cols[0]);
    histos_2D[key]->xsec_hist->SetLineColor(config->cols[0]);
    if(file_i != 0){
      histos_2D[key]->xsec_hist->SetFillStyle(config->fsty[file_i]);
      histos_2D[key]->xsec_hist->SetLineStyle(config->lsty[file_i]);
    }
    if(!config->plot_filled){
      histos_2D[key]->xsec_hist->SetFillColor(0);
      histos_2D[key]->xsec_hist->SetLineWidth(3);
    }

    // Set the total systematic histogram
    for(int i = 0; i <= histos_2D[key]->xsec_hist->GetNumberOfBins()+1; i++){
      histos_2D[key]->systematics->GetSyst("total")->mean_syst->SetBinContent(i, histos_2D[key]->xsec_hist->GetBinContent(i));
    }
  }
  
  // ------------------------------------------------------------------------------------------------------------------
  //                                        TOTAL RATE HISTOGRAM CREATION
  // ------------------------------------------------------------------------------------------------------------------
  
  // Get the total (unstacked) histogram
  TH1D* GetTotal(){

    TString tune = "";
    if(config->tune_name.size() > 1) tune = config->tune_name[file_i];
    double total;
    TH1D *total_hist = new TH1D(tune+"total", "", 1, 0, 2);

    for (size_t n = 0; n < dataman->total_data.size(); n++){
      total_hist->Fill(1);
      total++;
    }
    total_hist->SetBinContent(1, total_hist->GetBinContent(1)*config->pot_scale_fac[file_i]);

    return total_hist;
    
  }


  // ------------------------------------------------------------------------------------------------------------------
  //                                        1D RATE HISTOGRAM CREATION
  // ------------------------------------------------------------------------------------------------------------------
  
  // Create stacked histogram and legend from data
  std::pair<THStack*, TLegend*> StackHist1D(size_t var_i){

    TString tune = "";
    if(config->tune_name.size() > 1) tune = config->tune_name[file_i];
    std::vector<double> bin_edges = binman->Get1DBinning(config->plot_variables[var_i]);
    double edges_array[bin_edges.size()];
    std::copy(bin_edges.begin(), bin_edges.end(), edges_array);
    THStack *hstack = new THStack(tune+config->plot_variables[var_i]+"_stack", "");
    double x1 = 0.14; double x2 = 0.94; double y1 = 0.; double y2 = 0.06;
    if(!config->show_error_band){
      x1 = 0.24; x2 = 0.96; y1 = 0.85; y2 = 0.91;
    }
    TLegend *legend = new TLegend(x1, y1, x2, y2);

    int index = 0;
    for(auto const& dat: dataman->stack_data){
      TH1D* hist = new TH1D(tune+config->plot_variables[var_i]+dat.first.c_str(), "", bin_edges.size()-1, edges_array);
      for(size_t n = 0; n < dat.second.size(); n++){
        hist->Fill(dat.second[n][var_i]);
      }
      hist->Scale(config->pot_scale_fac[file_i]);
      hist->Scale(1, "width");
      hist->SetFillColor(config->cols[index]);
      hist->SetLineColor(config->cols[index]);
      if(config->input_file.size() > 1){
        hist->SetFillStyle(config->fsty[file_i]);
        hist->SetLineStyle(config->lsty[file_i]);
      }
      if(!config->plot_filled){
        hist->SetFillColor(0);
        hist->SetLineWidth(3);
      }
      hstack->Add(hist);
      legend->AddEntry(hist, dat.first.c_str(), "lf");
      index++;
    }

    return std::make_pair(hstack, legend);
  }

  // Get the total (unstacked) histogram
  TH1D* GetTotalHist1D(size_t var_i){
      
    TString tune = "";
    if(config->tune_name.size() > 1) tune = config->tune_name[file_i];
    std::vector<double> bin_edges = binman->Get1DBinning(config->plot_variables[var_i]);
    double edges_array[bin_edges.size()];
    std::copy(bin_edges.begin(), bin_edges.end(), edges_array);
    TH1D *total_hist = new TH1D(tune+config->plot_variables[var_i], "", bin_edges.size()-1, edges_array);

    for (size_t n = 0; n < dataman->total_data.size(); n++){
      total_hist->Fill(dataman->total_data[n][var_i]);
    }
    // Include scale factor bin by bin as ->Scale() won't change errors
    for(int n = 0; n <= total_hist->GetNbinsX(); n++){
      total_hist->SetBinContent(n, total_hist->GetBinContent(n)*config->pot_scale_fac[file_i]);
    }
    // If plotting cross section convert from rate
    total_hist->Scale(1, "width");
    total_hist->SetLineColor(config->cols[0]);
    if(config->input_file.size() > 1) total_hist->SetLineStyle(config->lsty[file_i]);
    return total_hist;
  }


  // ------------------------------------------------------------------------------------------------------------------
  //                                        2D RATE HISTOGRAM CREATION
  // ------------------------------------------------------------------------------------------------------------------

  // Create stacked histogram and legend from data
  std::pair<std::vector<THStack*>, TLegend*> StackHist2D(size_t var_i, size_t var_j){

    // Include name of tune if more than one file used
    TString tune = "";
    if(config->tune_name.size() > 1) tune = config->tune_name[file_i];
    TString name = tune+config->plot_variables[var_i]+"_"+config->plot_variables[var_j];

    // Get the bin edges
    std::vector<double> ybin_edges = binman->Get1DBinning(config->plot_variables[var_j]);
    std::vector<std::vector<double>> bin_edges = binman->Get2DBinning(config->plot_variables[var_i], config->plot_variables[var_j]);

    // Create a stacked hist for each Y bin
    std::vector<THStack*> hstacks;
    for(size_t bin_j = 0; bin_j < ybin_edges.size()-1; bin_j++){
      TString name = tune+config->plot_variables[var_i]+"_"+config->plot_variables[var_j]+Form("_%.1f_%.1f", ybin_edges[bin_j], ybin_edges[bin_j+1]);
      TString title = titles->hist_titles[var_j] + Form(": [%.2f, %.2f] ", ybin_edges[bin_j], ybin_edges[bin_j+1]) + titles->units[var_j];
      THStack *hstack = new THStack(name, title);
      hstacks.push_back(hstack);
    }
    double x1 = 0.14; double x2 = 0.94; double y1 = 0.; double y2 = 0.06;
    if(!config->show_error_band){
      x1 = 0.24; x2 = 0.96; y1 = 0.85; y2 = 0.91;
    }
    TLegend *legend = new TLegend(x1, y1, x2, y2);

    for(size_t bin_j = 0; bin_j < ybin_edges.size()-1; bin_j++){
      int index = 0;
      std::vector<double> xbin_edges = bin_edges[bin_j];
      double xedges_array[xbin_edges.size()];
      std::copy(xbin_edges.begin(), xbin_edges.end(), xedges_array);

      for(auto const& dat: dataman->stack_data){
        TString name_1D = Form(tune+config->plot_variables[var_i]+config->plot_variables[var_j]+dat.first.c_str()+"%i", bin_j);
        TH1D* hist = new TH1D(name_1D, "", xbin_edges.size()-1, xedges_array);
        for(size_t n = 0; n < dat.second.size(); n++){
          if(dat.second[n][var_j] >= ybin_edges[bin_j] && dat.second[n][var_j] < ybin_edges[bin_j+1]){
            hist->Fill(dat.second[n][var_i]);
          }
        }
        hist->Scale(config->pot_scale_fac[file_i]);
        // If plotting cross section convert from rate
        double width = (ybin_edges[bin_j+1] - ybin_edges[bin_j]);
        hist->Scale(1/width, "width");
        hist->SetFillColor(config->cols[index]);
        hist->SetLineColor(config->cols[index]);
        if(config->input_file.size() > 1){
          hist->SetFillStyle(config->fsty[file_i]);
          hist->SetLineStyle(config->lsty[file_i]);
        }
        if(!config->plot_filled){
          hist->SetFillColor(0);
          hist->SetLineWidth(3);
        }
        hstacks[bin_j]->Add(hist);
        if(bin_j==0) legend->AddEntry(hist, dat.first.c_str(), "lf");
        index++;
      }
    }
  
    return std::make_pair(hstacks, legend);
  }

  // Get the total (unstacked) histogram
  TH2Poly* GetTotalHist2D(size_t var_i, size_t var_j){
      
    // Include name of tune if more than one file used
    TString tune = "";
    if(config->tune_name.size() > 1) tune = config->tune_name[file_i];
    TString name = tune+config->plot_variables[var_i]+"_"+config->plot_variables[var_j];

    // Get the bin edges
    std::vector<double> ybin_edges = binman->Get1DBinning(config->plot_variables[var_j]);
    std::vector<std::vector<double>> bin_edges = binman->Get2DBinning(config->plot_variables[var_i], config->plot_variables[var_j]);

    // Create a 2d poly hist
    TH2Poly* total_hist = new TH2Poly();
    total_hist->SetName(name);
    total_hist->SetTitle("");

    // Add all the bins
    for(size_t i = 0; i < ybin_edges.size()-1; i++){
      std::vector<double> xbin_edges = bin_edges[i];
      for(size_t j = 0; j < xbin_edges.size()-1; j++){
        total_hist->AddBin(xbin_edges[j], ybin_edges[i], xbin_edges[j+1], ybin_edges[i+1]);
      }
    }

    for (size_t n = 0; n < dataman->total_data.size(); n++){
      total_hist->Fill(dataman->total_data[n][var_i], dataman->total_data[n][var_j]);
    }
    // Include scale factor bin by bin as ->Scale() won't change errors
    for(int i = 0; i <= total_hist->GetNumberOfBins(); i++){
      total_hist->SetBinContent(i, total_hist->GetBinContent(i)*config->pot_scale_fac[file_i]);
    }
    // If plotting cross section convert from rate
    total_hist->Scale(1, "width");
    total_hist->SetLineColor(config->cols[0]);
    if(config->plot_variables.size() > 1) total_hist->SetLineStyle(config->lsty[file_i]);
    return total_hist;
  }

};

#endif
