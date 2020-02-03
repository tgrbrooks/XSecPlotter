#ifndef SYSTSUMMARY_H
#define SYSTSUMMARY_H

#include "Systematics.h"

// Structure for holding systematic error information
class SystSummary
{
  public:

  Systematics* genie;
  Systematics* flux;
  Systematics* detector;
  Systematics* background;
  Systematics* constant;
  TH1D* total;

  SystSummary(){}

  SystSummary(Systematics *g, Systematics *f, Systematics *d, Systematics *b, Systematics *c)
  {
    genie = g;
    flux = f;
    detector = d;
    background = b;
    constant = c;

    total = (TH1D*)genie->mean_syst->Clone();
    total->Reset();
    for(auto const& syst : config->systematics){
      if(syst == "genie") AddErrors(total, genie->mean_syst);
      if(syst == "flux") AddErrors(total, flux->mean_syst);
      if(syst == "detector") AddErrors(total, detector->mean_syst);
      if(syst == "background") AddErrors(total, background->mean_syst);
      if(syst == "constant") AddErrors(total, constant->mean_syst);
    }

  }

  // Add histogram errors in quadrature, ignoring bin contents
  void AddErrors(TH1D* syst_hist, TH1D* hist){
    for(size_t i = 1; i <= syst_hist->GetNbinsX(); i++){
      double new_err = std::sqrt(std::pow(syst_hist->GetBinError(i),2)+std::pow(hist->GetBinError(i),2));
      syst_hist->SetBinError(i, new_err);
    }
  }

};

#endif
