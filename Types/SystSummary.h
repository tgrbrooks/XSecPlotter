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

  SystSummary(TH1D* hist){
    genie = new Systematics(hist, "_geniesyst");
    flux = new Systematics(hist, "_fluxsyst");
    detector = new Systematics(hist, "_detsyst");
    background = new Systematics(hist, "_bkgsyst");
    constant = new Systematics(hist, "_constsyst");
    total = (TH1D*)hist->Clone(TString(hist->GetName())+"_totalsyst");
    total->Reset();
  }

  SystSummary(Systematics *g, Systematics *f, Systematics *d, Systematics *b, Systematics *c)
  {
    genie = g;
    flux = f;
    detector = d;
    background = b;
    constant = c;

    total = (TH1D*)genie->mean_syst->Clone();
    total->Reset();
    AddErrors(total, genie->mean_syst);
    AddErrors(total, flux->mean_syst);
    AddErrors(total, detector->mean_syst);
    AddErrors(total, background->mean_syst);
    AddErrors(total, constant->mean_syst);

  }

  void GetTotal(){
    AddErrors(total, genie->mean_syst);
    AddErrors(total, flux->mean_syst);
    AddErrors(total, detector->mean_syst);
    AddErrors(total, background->mean_syst);
    AddErrors(total, constant->mean_syst);
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
