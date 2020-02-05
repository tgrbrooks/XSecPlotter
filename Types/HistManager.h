#ifndef HISTMANAGER_H
#define HISTMANAGER_H

#include "Configuration.h"
#include "Titles.h"
#include "DataManager.h"
#include "BinManager.h"
#include "Histo1D.h"
#include "Histo2D.h"

// Structure for holding interaction information
class HistManager
{
  public:

  Configuration *config;
  Titles *titles;
  DataManager *dataman;
  BinManager *binman;
  size_t file_i;

  Histo1D* total;
  std::map<TString, Histo1D*> histos_1D;
  std::map<std::pair<TString, TString>, Histo2D*> histos_2D;

  HistManager(Configuration *c, Titles *t, DataManager *d, BinManager *b, size_t f)
  {
    config = c;
    titles = t;
    dataman = d;
    binman = b;
    file_i = f;

    // Total event rate
    TH1D* hist0D = GetTotal();
    total = new Histo1D(hist0D);
    // 1D histograms
    for(size_t i = 0; i < config->plot_variables.size(); i++){
      TString key1D = config->plot_variables[i];
      TH1D* hist1D = GetTotalHist1D(i);
      std::pair<THStack*, TLegend*> stack1D = StackHist1D(i);
      Histo1D* histo1D = new Histo1D(hist1D, stack1D);
      histos_1D[key1D] = histo1D;
      // Calculate efficiency and purity
      if(config->plot_eff_pur) SetEfficiency(i);
      // 2D histograms
      for(size_t j = 0; j < config->plot_variables.size(); j++){
        if(i==j) continue;
        std::pair<TString, TString> key2D = std::make_pair(key1D, config->plot_variables[j]);
        TH2D* hist2D = GetTotalHist2D(i, j);
        std::pair<std::vector<THStack*>, TLegend*> stack2D = StackHist2D(i, j);
        Histo2D* histo2D = new Histo2D(hist2D, stack2D);
        histos_2D[key2D] = histo2D;
        if(config->plot_eff_pur) SetEfficiency(i, j);
      }
    }
    
  }

  TH1D* Get1DHist(TString plot_var){
    return histos_1D[plot_var]->total_hist;
  }
  TH1D* Get1DErrorBand(TString plot_var){
    return GetErrorBand(histos_1D[plot_var]->total_hist);
  }
  THStack* Get1DStack(TString plot_var){
    return histos_1D[plot_var]->stacked_hist;
  }
  TLegend* Get1DLegend(TString plot_var){
    return histos_1D[plot_var]->legend;
  }

  Histo1D* GetHisto1D(TString plot_var){
    return histos_1D[plot_var];
  }

  std::vector<double> Get1DBinning(TString plot_var){
    return binman->Get1DBinning(plot_var);
  }

  void SetTotalSystematics(SystSummary* syst){
    total->systematics = syst;
  }

  void Set1DSystematics(TString plot_var, SystSummary* syst){
    histos_1D[plot_var]->systematics = syst;
  }

  void Set2DSystematics(TString plot_var1, TString plot_var2, SystSummary2D* syst){
    std::pair<TString, TString> key = std::make_pair(plot_var1, plot_var2);
    histos_2D[key]->systematics = syst;
  }

  TH1D* Get1DHistSlice(TString plot_var, TString slice_var, size_t slice_bin){
    std::pair<TString, TString> key = std::make_pair(plot_var, slice_var);
    return histos_2D[key]->total_hist->ProjectionX(Form(key.first+key.second+"total%i",slice_bin), slice_bin, slice_bin);
  }
  TH1D* Get1DSliceErrorBand(TString plot_var, TString slice_var, size_t slice_bin){
    std::pair<TString, TString> key = std::make_pair(plot_var, slice_var);
    return GetErrorBand(histos_2D[key]->total_hist->ProjectionX(Form(key.first+key.second+"total%i",slice_bin), slice_bin, slice_bin));
  }
  THStack* Get1DStackSlice(TString plot_var, TString slice_var, size_t slice_bin){
    std::pair<TString, TString> key = std::make_pair(plot_var, slice_var);
    return histos_2D[key]->stacked_hist[slice_bin];
  }
  TH2D* Get2DHist(TString plot_var, TString slice_var){
    std::pair<TString, TString> key = std::make_pair(plot_var, slice_var);
    return histos_2D[key]->total_hist;
  }
  TLegend* Get2DLegend(TString plot_var, TString slice_var){
    std::pair<TString, TString> key = std::make_pair(plot_var, slice_var);
    return histos_2D[key]->legend;
  }

  Histo2D* GetHisto2D(TString plot_var, TString slice_var){
    std::pair<TString, TString> key = std::make_pair(plot_var, slice_var);
    return histos_2D[key];
  }

  int GetNHists(){
    int n_1D = histos_1D.size();
    int n_2D = histos_2D.size();
    if(n_2D == 0) return n_1D;
    int n_1D_slice = 0;
    bool first = true;
    for(auto const& kv : histos_2D){
      if(!first) continue;
      n_1D_slice += kv.second->total_hist->GetNbinsX();
      //n_1D_slice += kv.second->total_hist->GetNbinsY();
      first = false;
    }
    return n_1D + n_2D + n_1D_slice;
  }

  void CreateUniverses(TString syst, int n){
    if(syst=="detector"){
      total->systematics->detector->CreateUniverses(n);
      for(auto const& kv : histos_1D){
        histos_1D[kv.first]->systematics->detector->CreateUniverses(n);
      }
      for(auto const& kv : histos_2D){
        histos_2D[kv.first]->systematics->detector->CreateUniverses(n);
      }
    }
    else if(syst=="genie"){
      total->systematics->genie->CreateUniverses(n);
      for(auto const& kv : histos_1D){
        histos_1D[kv.first]->systematics->genie->CreateUniverses(n);
      }
      for(auto const& kv : histos_2D){
        histos_2D[kv.first]->systematics->genie->CreateUniverses(n);
      }
    }
    else if(syst=="flux"){
      total->systematics->flux->CreateUniverses(n);
      for(auto const& kv : histos_1D){
        histos_1D[kv.first]->systematics->flux->CreateUniverses(n);
      }
      for(auto const& kv : histos_2D){
        histos_2D[kv.first]->systematics->flux->CreateUniverses(n);
      }
    }
  }

  void ScaleUniverses(TString syst){
    if(syst=="detector"){
      total->systematics->detector->ScaleUniverses(config, file_i);
      for(auto& kv : histos_1D){
        kv.second->systematics->detector->ScaleUniverses(config, file_i);
      }
      for(auto& kv : histos_2D){
        kv.second->systematics->detector->ScaleUniverses(config, file_i);
      }
    }
    else if(syst=="genie"){
      total->systematics->genie->ScaleUniverses(config, file_i);
      for(auto& kv : histos_1D){
        kv.second->systematics->genie->ScaleUniverses(config, file_i);
      }
      for(auto& kv : histos_2D){
        kv.second->systematics->genie->ScaleUniverses(config, file_i);
      }
    }
    else if(syst=="flux"){
      total->systematics->flux->ScaleUniverses(config, file_i);
      for(auto& kv : histos_1D){
        kv.second->systematics->flux->ScaleUniverses(config, file_i);
      }
      for(auto& kv : histos_2D){
        kv.second->systematics->flux->ScaleUniverses(config, file_i);
      }
    }
  }

  void CalculateSyst(TString syst){
    if(syst=="detector"){
      total->systematics->detector->Calculate();
      for(auto& kv : histos_1D){
        kv.second->systematics->detector->Calculate();
      }
      for(auto& kv : histos_2D){
        kv.second->systematics->detector->Calculate();
      }
    }
    else if(syst=="genie"){
      total->systematics->genie->Calculate();
      for(auto& kv : histos_1D){
        kv.second->systematics->genie->Calculate();
      }
      for(auto& kv : histos_2D){
        kv.second->systematics->genie->Calculate();
      }
    }
    else if(syst=="flux"){
      total->systematics->flux->Calculate();
      for(auto& kv : histos_1D){
        kv.second->systematics->flux->Calculate();
      }
      for(auto& kv : histos_2D){
        kv.second->systematics->flux->Calculate();
      }
    }
  }

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

  void SetEfficiency(size_t var_i, size_t var_j){

    std::pair<TString, TString> key = std::make_pair(config->plot_variables[var_i], config->plot_variables[var_j]);
    /*TH2D *eff_numerator = (TH2D*) histos_2D[key]->total_hist->Clone();
    eff_numerator->Reset();
    TH2D *eff_denom = (TH2D*) histos_2D[key]->total_hist->Clone();
    eff_denom->Reset();
    TH2D *pur_numerator = (TH2D*) histos_2D[key]->total_hist->Clone();
    pur_numerator->Reset();
    TH2D *pur_denom = (TH2D*) histos_2D[key]->total_hist->Clone();
    pur_denom->Reset();*/

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

    //histos_2D[key]->efficiency = std::make_pair(eff_numerator, eff_denom);
    //histos_2D[key]->purity = std::make_pair(pur_numerator, pur_denom);
    //histos_2D[key]->eff_calculated = true;
  }

/*
  // Create stacked histogram and legend from data
  std::pair<THStack*, TLegend*> StackHist1D(const std::map<std::string, std::vector<std::vector<double>>> &data, TString name, TString title, std::vector<std::vector<double>> bin_edges, int i, int j = -1, int bin_j = -1){

    double edges_array[bin_edges[i].size()];
    std::copy(bin_edges[i].begin(), bin_edges[i].end(), edges_array);
    THStack *hstack = new THStack(name, title);
    TLegend *legend = new TLegend(0.14, 0., 0.94, 0.06);

    int index = 0;
    for(auto const& dat: data){
      TH1D* hist = new TH1D(name+dat.first.c_str(), title, bin_edges[i].size()-1, edges_array);
      for(size_t n = 0; n < dat.second.size(); n++){
        if(j == -1){
          hist->Fill(dat.second[n][i]);
        }
        else if(dat.second[n][j] >= bin_edges[j][bin_j] && dat.second[n][j] < bin_edges[j][bin_j+1]){
          hist->Fill(dat.second[n][i]);
        }
      }
      hist->Scale(config->pot_scale_fac[0]);
      // If plotting cross section convert from rate
      if(config->plot_xsec){
        double width = 1;
        if(j != -1) width = width * (bin_edges[j][bin_j+1] - bin_edges[j][bin_j]);
        double xsec_scale = 1e38/(width * config->flux[0] * config->targets);
        hist->Scale(xsec_scale, "width");
      }
      // Else if max error used divide each bin by width
      else if (config->max_error > 0 || config->bin_edges[i].size()>1){
        hist->Scale(1, "width");
      }
      hist->SetFillColor(config->cols[index]);
      hist->SetLineColor(config->cols[index]);
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
  TH1D* GetTotalHist(const std::vector<std::vector<double>> &data, TString name, std::vector<std::vector<double>> bin_edges, int tune_i, int i, int j = -1, int bin_j = -1){
      
    double edges_array[bin_edges[i].size()];
    std::copy(bin_edges[i].begin(), bin_edges[i].end(), edges_array);
    TH1D *total_hist = new TH1D("total"+name+config->tune_name[tune_i], "hist", bin_edges[i].size()-1, edges_array);

    for (int n = 0; n < data.size(); n++){
      if(j == -1){
        total_hist->Fill(data[n][i]);
      }
      else if(data[n][j] >= bin_edges[j][bin_j] && data[n][j] < bin_edges[j][bin_j+1]){
        total_hist->Fill(data[n][i]);
      }
    }
    // Include scale factor bin by bin as ->Scale() won't change errors
    for(size_t n = 0; n <= total_hist->GetNbinsX(); n++){
      total_hist->SetBinContent(n, total_hist->GetBinContent(n)*config->pot_scale_fac[tune_i]);
    }
    // If plotting cross section convert from rate
    if(config->plot_xsec){
      double width = 1;
      if(j != -1) width = width * (bin_edges[j][bin_j+1] - bin_edges[j][bin_j]);
      double xsec_scale = 1e38/(width * config->flux[tune_i] * config->targets);
      total_hist->Scale(xsec_scale, "width");
    }
    // Else if max error used divide each bin by width FIXME is this ok for stat errors?
    else if (config->max_error > 0 || config->bin_edges[i].size()>1){
      total_hist->Scale(1, "width");
    }
    return total_hist;
  }
*/
  // Create stacked histogram and legend from data
  std::pair<THStack*, TLegend*> StackHist1D(size_t var_i){

    std::vector<double> bin_edges = binman->Get1DBinning(config->plot_variables[var_i]);
    double edges_array[bin_edges.size()];
    std::copy(bin_edges.begin(), bin_edges.end(), edges_array);
    THStack *hstack = new THStack(config->plot_variables[var_i]+"_stack", "");
    TLegend *legend = new TLegend(0.14, 0., 0.94, 0.06);

    int index = 0;
    for(auto const& dat: dataman->stack_data){
      TH1D* hist = new TH1D(config->plot_variables[var_i]+dat.first.c_str(), "", bin_edges.size()-1, edges_array);
      for(size_t n = 0; n < dat.second.size(); n++){
        hist->Fill(dat.second[n][var_i]);
      }
      hist->Scale(config->pot_scale_fac[0]);
      // If plotting cross section convert from rate
      if(config->plot_xsec){
        double xsec_scale = 1e38/(config->flux[0] * config->targets);
        hist->Scale(xsec_scale, "width");
      }
      // Else if max error used divide each bin by width
      else if (config->max_error > 0 || config->bin_edges[var_i].size()>1){
        hist->Scale(1, "width");
      }
      hist->SetFillColor(config->cols[index]);
      hist->SetLineColor(config->cols[index]);
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
  TH1D* GetTotal(){

    double total;
    TH1D *total_hist = new TH1D("total", "", 1, 0, 2);

    for (int n = 0; n < dataman->total_data.size(); n++){
      total_hist->Fill(1);
      total++;
    }
    total_hist->SetBinContent(1, total_hist->GetBinContent(1)*config->pot_scale_fac[file_i]);
    // If plotting cross section convert from rate
    if(config->plot_xsec){
      double xsec_scale = 1e38/(config->flux[file_i] * config->targets);
      total_hist->Scale(xsec_scale);
    }
    return total_hist;
    
  }

  // Get the total (unstacked) histogram
  TH1D* GetTotalHist1D(size_t var_i){
      
    std::vector<double> bin_edges = binman->Get1DBinning(config->plot_variables[var_i]);
    double edges_array[bin_edges.size()];
    std::copy(bin_edges.begin(), bin_edges.end(), edges_array);
    TH1D *total_hist = new TH1D(config->plot_variables[var_i], "", bin_edges.size()-1, edges_array);

    for (int n = 0; n < dataman->total_data.size(); n++){
      total_hist->Fill(dataman->total_data[n][var_i]);
    }
    // Include scale factor bin by bin as ->Scale() won't change errors
    for(size_t n = 0; n <= total_hist->GetNbinsX(); n++){
      total_hist->SetBinContent(n, total_hist->GetBinContent(n)*config->pot_scale_fac[file_i]);
    }
    // If plotting cross section convert from rate
    if(config->plot_xsec){
      double xsec_scale = 1e38/(config->flux[file_i] * config->targets);
      total_hist->Scale(xsec_scale, "width");
    }
    // Else if max error used divide each bin by width FIXME is this ok for stat errors?
    else if (config->max_error > 0 || config->bin_edges[var_i].size()>1){
      total_hist->Scale(1, "width");
    }
    return total_hist;
  }

  // Create stacked histogram and legend from data
  std::pair<std::vector<THStack*>, TLegend*> StackHist2D(size_t var_i, size_t var_j){

    std::pair<std::vector<double>, std::vector<double>> bin_edges = binman->Get2DBinning(config->plot_variables[var_i], config->plot_variables[var_j]);
    double xedges_array[bin_edges.first.size()];
    std::copy(bin_edges.first.begin(), bin_edges.first.end(), xedges_array);
    double yedges_array[bin_edges.second.size()];
    std::copy(bin_edges.second.begin(), bin_edges.second.end(), yedges_array);

    std::vector<THStack*> hstacks;
    for(size_t bin_j = 0; bin_j < bin_edges.second.size()-1; bin_j++){
      TString title = titles->hist_titles[var_j] + Form(": [%.2f, %.2f] ", yedges_array[bin_j], yedges_array[bin_j+1]) + titles->units[var_j];
      THStack *hstack = new THStack(Form(config->plot_variables[var_i]+config->plot_variables[var_j]+"_stack%i", bin_j), title);
      hstacks.push_back(hstack);
    }
    TLegend *legend = new TLegend(0.14, 0., 0.94, 0.06);

    for(size_t bin_j = 0; bin_j < bin_edges.second.size()-1; bin_j++){
      int index = 0;
      for(auto const& dat: dataman->stack_data){
        TH1D* hist = new TH1D(Form(config->plot_variables[var_i]+config->plot_variables[var_j]+dat.first.c_str()+"%i", bin_j), "", bin_edges.first.size()-1, xedges_array);
        for(size_t n = 0; n < dat.second.size(); n++){
          if(dat.second[n][var_j] >= bin_edges.second[bin_j] && dat.second[n][var_j] < bin_edges.second[bin_j+1]){
            hist->Fill(dat.second[n][var_i]);
          }
        }
        hist->Scale(config->pot_scale_fac[0]);
        // If plotting cross section convert from rate
        if(config->plot_xsec){
          double width = (bin_edges.second[bin_j+1] - bin_edges.second[bin_j]);
          double xsec_scale = 1e38/(width * config->flux[0] * config->targets);
          hist->Scale(xsec_scale, "width");
        }
        // Else if max error used divide each bin by width
        else if (config->max_error > 0 || config->bin_edges[var_i].size()>1){
          double width = (bin_edges.second[bin_j+1] - bin_edges.second[bin_j]);
          hist->Scale(1/width, "width");
        }
        hist->SetFillColor(config->cols[index]);
        hist->SetLineColor(config->cols[index]);
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
  TH2D* GetTotalHist2D(size_t var_i, size_t var_j){
      
    std::pair<std::vector<double>, std::vector<double>> bin_edges = binman->Get2DBinning(config->plot_variables[var_i], config->plot_variables[var_j]);
    double xedges_array[bin_edges.first.size()];
    std::copy(bin_edges.first.begin(), bin_edges.first.end(), xedges_array);
    double yedges_array[bin_edges.second.size()];
    std::copy(bin_edges.second.begin(), bin_edges.second.end(), yedges_array);
    TH2D *total_hist = new TH2D(config->plot_variables[var_i]+"_"+config->plot_variables[var_j], "", bin_edges.first.size()-1, xedges_array, bin_edges.second.size()-1, yedges_array);

    for (int n = 0; n < dataman->total_data.size(); n++){
      total_hist->Fill(dataman->total_data[n][var_i], dataman->total_data[n][var_j]);
    }
    // Include scale factor bin by bin as ->Scale() won't change errors
    for(size_t i = 0; i <= total_hist->GetNbinsX(); i++){
      for(size_t j = 0; j <= total_hist->GetNbinsX(); j++){
        total_hist->SetBinContent(i, j, total_hist->GetBinContent(i, j)*config->pot_scale_fac[file_i]);
      }
    }
    // If plotting cross section convert from rate
    if(config->plot_xsec){
      double xsec_scale = 1e38/(config->flux[file_i] * config->targets);
      total_hist->Scale(xsec_scale, "width");
    }
    // Else if max error used divide each bin by width FIXME is this ok for stat errors?
    else if (config->max_error > 0 || config->bin_edges[var_i].size()>1){
      total_hist->Scale(1, "width");
    }
    return total_hist;
  }

  // Get the percentage statistical error per bin
  TH1D* GetErrorBand(TH1D* total_hist){

    TH1D *error_band = (TH1D*)total_hist->Clone();

    // Set the bin errors on seperate plot
    for (int n = 1; n <= total_hist->GetNbinsX(); n++){
     error_band->SetBinContent(n, 0);
     if (total_hist->GetBinContent(n) > 0)
       error_band->SetBinContent(n, 100*total_hist->GetBinError(n)/total_hist->GetBinContent(n));
    }
    return error_band;

  }

/*
  TH2D* Get2DHist(const std::vector<std::vector<double>>& data, std::vector<std::vector<double>> bin_edges, TString name, int d_i, int d_j){
    
    double edges_array[bin_edges[d_i].size()];
    std::copy(bin_edges[d_i].begin(), bin_edges[d_i].end(), edges_array);

    double edges_array2[bin_edges[d_j].size()];
    std::copy(bin_edges[d_j].begin(), bin_edges[d_j].end(), edges_array2);

    TH2D* hist_2D = new TH2D(name+config->plot_variables[d_j], "", bin_edges[d_i].size()-1, edges_array, bin_edges[d_j].size()-1, edges_array2);

    for (int n = 0; n < data.size(); n++){
      hist_2D->Fill(data[n][d_i], data[n][d_j]);
    }

    return hist_2D;
  }
*/
  
};

#endif
