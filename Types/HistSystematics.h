#ifndef HISTSYSTEMATICS_H
#define HISTSYSTEMATICS_H

#include "Systematics.h"
#include "Systematics2D.h"

// Structure for holding interaction information
class HistSystematics
{
  public:

  Systematics* total;
  std::map<TString, Systematics*> h1D;
  std::map<std::pair<TString, TString>, Systematics2D*> h2D;

  HistSystematics(HistManager* histman, TString name)
  {
    TH1D *temp_tot = (TH1D*)histman->total->total_hist->Clone(name);
    total = new Systematics(temp_tot, true);
    for(auto const& kv1D : histman->histos_1D){
      TH1D *temp_1D = (TH1D*)kv1D.second->total_hist->Clone(name+kv1D.first);
      h1D[kv1D.first] = new Systematics(temp_1D, true);
    }
    for(auto const& kv2D : histman->histos_2D){
      TH2D *temp_2D = (TH2D*)kv2D.second->total_hist->Clone(name+kv2D.first.first+kv2D.first.second);
      h2D[kv2D.first] = new Systematics2D(temp_2D, true);
    }
  }
  
};

#endif
