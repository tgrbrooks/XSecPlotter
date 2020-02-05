#ifndef SYSTSUMMARY2D_H
#define SYSTSUMMARY2D_H

#include "Systematics2D.h"

// Structure for holding systematic error information
class SystSummary2D
{
  public:

  Systematics2D* genie;
  Systematics2D* flux;
  Systematics2D* detector;
  Systematics2D* background;
  Systematics2D* constant;
  TH2D* total;

  SystSummary2D(){}

  SystSummary2D(TH2D* hist){
    genie = new Systematics2D(hist, "_geniesyst");
    flux = new Systematics2D(hist, "_fluxsyst");
    detector = new Systematics2D(hist, "_detsyst");
    background = new Systematics2D(hist, "_bkgsyst");
    constant = new Systematics2D(hist, "_constsyst");
    total = (TH2D*)hist->Clone(TString(hist->GetName())+"_totalsyst");
    total->Reset();
  }

  SystSummary2D(Systematics2D *g, Systematics2D *f, Systematics2D *d, Systematics2D *b, Systematics2D *c)
  {
    genie = g;
    flux = f;
    detector = d;
    background = b;
    constant = c;

    total = (TH2D*)genie->mean_syst->Clone();
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
  void AddErrors(TH2D* syst_hist, TH2D* hist){
    for(size_t i = 1; i <= syst_hist->GetNbinsX(); i++){
      for(size_t j = 1; j <= syst_hist->GetNbinsY(); j++){
        double new_err = std::sqrt(std::pow(syst_hist->GetBinError(i, j),2)+std::pow(hist->GetBinError(i, j),2));
        syst_hist->SetBinError(i, j, new_err);
      }
    }
  }

  SystSummary* Slice(size_t i){
    Systematics* genie_s = genie->Slice(i);
    Systematics* flux_s = flux->Slice(i);
    Systematics* detector_s = detector->Slice(i);
    Systematics* background_s = background->Slice(i);
    Systematics* constant_s = constant->Slice(i);
    SystSummary* summary = new SystSummary(genie_s, flux_s, detector_s, background_s, constant_s);
    summary->GetTotal();
    return summary;
  }

};

#endif
