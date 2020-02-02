#ifndef HISTMANAGER_H
#define HISTMANAGER_H

#include "Configuration.h"
#include "DataManager.h"
#include "BinManager.h"
#include "Histo1D.h"

// Structure for holding interaction information
class HistManager
{
  public:

  Configuration *config;
  DataManager *dataman;
  BinManager *binman;

  std::map<TString, TH1D*> hists_1D;
  std::map<std::pair<TString, TString>, TH2D*> hists_2D;
  std::map<std::pair<TString, TString>, std::vector<TH1D*>> hists_1D_slice_2D;
  std::map<std::pair<TString, TString>, std::vector<std::pair<THStack*, TLegend*>>> stacks_1D_slice_2D;
  std::map<std::vector<TString>, std::vector<std::vector<TH1D*>>> hists_1D_slice_3D;
  std::map<std::vector<TString>, std::vector<std::vector<std::pair<THStack*, TLegend*>>>> stacks_1D_slice_3D;
  std::map<std::vector<TString>, std::vector<TH2D*>> hists_2D_slice_3D;

  HistManager(Configuration *c, DataManager *d, BinManager *b)
  {
    config = c;
    dataman = d;
    binman = b;

    /*
    for(size_t i = 0; i < config->plot_variables.size(); i++){
      TString name_1D = config->plot_variables[d_i];
      TString key1D = config->plot_variables[i];
      hists_1D[key1D] = GetTotalHist(dataman->total_data, name_1D, bin_edges, file_i, d_i);
      stacks_1D[key1D] = StackHist1D(dataman->stack_data, name_1D, title_1D, bin_edges, d_i);
    }
    */
  }

  TH1D* Get1DHist(TString plot_var){
    return hists_1D[plot_var];
  }
   
  std::pair<THStack*, TLegend*> Get1DStack(TString plot_var){
    return stacks_1D[plot_var];
  }

  TH1D* Get1DHistFrom2DSlice(TString plot_var, TString slice_var, size_t slice_bin){
    std::pair<TString, TString> key = std::make_pair(plot_var, slice_var);
    return hists_1D_slice_2D[key][slice_bin];
  }

  std::pair<THStack*, TLegend*> Get1DStackFrom2DSlice(TString plot_var, TString slice_var, size_t slice_bin){
    std::pair<TString, TString> key = std::make_pair(plot_var, slice_var);
    return stacks_1D_slice_2D[key][slice_bin];
  }

  TH1D* Get1DHistFrom3DSlice(TString plot_var, TString slice_var1, size_t slice_bin1, TString slice_var2, size_t slice_bin2){
    std::vector<TString> key = {plot_var, slice_var1, slice_var2};
    return hists_1D_slice_3D[key][slice_bin1][slice_bin2];
  }

  std::pair<THStack*, TLegend*> Get1DStackFrom3DSlice(TString plot_var, TString slice_var1, size_t slice_bin1, TString slice_var2, size_t slice_bin2){
    std::vector<TString> key = {plot_var, slice_var1, slice_var2};
    return stacks_1D_slice_3D[key][slice_bin1][slice_bin2];
  }

  TH2D* Get2DHist(TString plot_var1, TString plot_var2){
    std::pair<TString, TString> key = std::make_pair(plot_var1, plot_var2);
    return hists_2D[key];
  }

  TH2D* Get2DHistFrom3DSlice(TString plot_var1, TString plot_var2, TString slice_var, size_t slice_bin){
    std::vector<TString> key = {plot_var1, plot_var2, slice_var};
    return hists_2D_slice_3D[key][slice_bin];
  }

  // Create stacked histogram and legend from data
  std::pair<THStack*, TLegend*> StackHist1D(const std::map<std::string, std::vector<std::vector<double>>> &data, TString name, TString title, std::vector<std::vector<double>> bin_edges, int i, int j = -1, int bin_j = -1, int k = -1, int bin_k = -1){

    double edges_array[bin_edges[i].size()];
    std::copy(bin_edges[i].begin(), bin_edges[i].end(), edges_array);
    THStack *hstack = new THStack(name, title);
    TLegend *legend = new TLegend(0.14, 0., 0.94, 0.06);

    int index = 0;
    for(auto const& dat: data){
      TH1D* hist = new TH1D(name+dat.first.c_str(), title, bin_edges[i].size()-1, edges_array);
      for(size_t n = 0; n < dat.second.size(); n++){
        if(j==-1 && k == -1){
          hist->Fill(dat.second[n][i]);
        }
        else if(dat.second[n][j] >= bin_edges[j][bin_j] && dat.second[n][j] < bin_edges[j][bin_j+1]){
          if(k == -1){
            hist->Fill(dat.second[n][i]);
          }
          else{
            if(dat.second[n][k] >= bin_edges[k][bin_k] && dat.second[n][k] < bin_edges[k][bin_k+1]){
              hist->Fill(dat.second[n][i]);
            }
          }
        }
      }
      hist->Scale(config->pot_scale_fac[0]);
      // If plotting cross section convert from rate
      if(config->plot_xsec){
        double width = 1;
        if(j != -1) width = width * (bin_edges[j][bin_j+1] - bin_edges[j][bin_j]);
        if(k != -1) width = width * (bin_edges[k][bin_k+1] - bin_edges[k][bin_k]);
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
  TH1D* GetTotalHist(const std::vector<std::vector<double>> &data, TString name, std::vector<std::vector<double>> bin_edges, int tune_i, int i, int j = -1, int bin_j = -1, int k = -1, int bin_k = -1){
      
    double edges_array[bin_edges[i].size()];
    std::copy(bin_edges[i].begin(), bin_edges[i].end(), edges_array);
    TH1D *total_hist = new TH1D("total"+name+config->tune_name[tune_i], "hist", bin_edges[i].size()-1, edges_array);

    for (int n = 0; n < data.size(); n++){
      if(j == -1 && k == -1){
        total_hist->Fill(data[n][i]);
      }
      else if(data[n][j] >= bin_edges[j][bin_j] && data[n][j] < bin_edges[j][bin_j+1]){
        if(k == -1)
          total_hist->Fill(data[n][i]);
        else if(data[n][k] >= bin_edges[k][bin_k] && data[n][k] < bin_edges[k][bin_k+1])
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
      if(k != -1) width = width * (bin_edges[k][bin_k+1] - bin_edges[k][bin_k]);
      double xsec_scale = 1e38/(width * config->flux[tune_i] * config->targets);
      total_hist->Scale(xsec_scale, "width");
    }
    // Else if max error used divide each bin by width FIXME is this ok for stat errors?
    else if (config->max_error > 0 || config->bin_edges[i].size()>1){
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
  
};

#endif
