#ifndef SYSTHOLDER_H
#define SYSTHOLDER_H

#include "HistManager.h"

// Structure for holding interaction information
class SystHolder
{
  public:

  std::vector<TH1D*> total;
  std::map<TString, std::vector<TH1D*>> hists_1D;
  std::map<std::pair<TString, TString>, std::vector<TH2D*>> hists_2D;
  size_t universes;

  SystHolder(HistManager* histman, size_t u)
  {
    universes = u;

    for(size_t u = 0; u < universes; u++){
      TH1D* total_temp = (TH1D*)histman->total->total_hist->Clone();
      total_temp->Reset();
      total.push_back(total_temp);
      for(auto const& kv1D : histman->histos_1D){
        TH1D* temp_1D = (TH1D*)kv1D.second->total_hist->Clone();
        temp_1D->Reset();
        hists_1D[kv1D.first].push_back(temp_1D);
      }
      for(auto const& kv2D : histman->histos_2D){
        TH2D* temp_2D = (TH2D*)kv2D.second->total_hist->Clone();
        temp_2D->Reset();
        hists_2D[kv2D.first].push_back(temp_2D);
      }
    }
  }

  void Fill(double weight, size_t u){
    total[u]->Fill(1., weight);
  }

  void Fill1D(TString plot_var, double val, double weight, size_t u){
    if(hists_1D.find(plot_var) != hists_1D.end())
      hists_1D[plot_var][u]->Fill(val, weight);
  }

  void Fill2D(TString plot_var1, double val1, TString plot_var2, double val2, double weight, size_t u){
    std::pair<TString, TString> key = std::make_pair(plot_var1, plot_var2);
    if(hists_2D.find(key) != hists_2D.end())
      hists_2D[key][u]->Fill(val1, val2, weight);
  }

  void Scale(Configuration *config, size_t file_i){

    double xsec_scale = 1e38/(config->flux[file_i] * config->targets);

    for(size_t u = 0; u < universes; u++){
      total[u]->Scale(config->pot_scale_fac[file_i]);
      if(config->plot_xsec){
        total[u]->Scale(xsec_scale, "width");
      }
      else if(config->max_error > 0 || config->bin_edges[0].size() > 1){
        total[u]->Scale(1, "width");
      }
      for(auto const& kv1D : hists_1D){
        hists_1D[kv1D.first][u]->Scale(config->pot_scale_fac[file_i]);
        if(config->plot_xsec){
          hists_1D[kv1D.first][u]->Scale(xsec_scale, "width");
        }
        else if(config->max_error > 0 || config->bin_edges[0].size() > 1){
          hists_1D[kv1D.first][u]->Scale(1, "width");
        }
      }
      for(auto const& kv2D : hists_2D){
        hists_2D[kv2D.first][u]->Scale(config->pot_scale_fac[file_i]);
        if(config->plot_xsec){
          hists_2D[kv2D.first][u]->Scale(xsec_scale, "width");
        }
        else if(config->max_error > 0 || config->bin_edges[0].size() > 1){
          hists_2D[kv2D.first][u]->Scale(1, "width");
        }
      }
    }
  }
  
};

#endif
