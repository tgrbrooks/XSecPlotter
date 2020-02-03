#ifndef HISTO1D_H
#define HISTO1D_H

#include "SystSummary.h"

// Structure for holding 1D histogram information
class Histo1D
{
  public:

  TH1D* total_hist;
  TH1D* error_band;
  THStack* stacked_hist;
  TLegend* legend;
  std::pair<TH1D*, TH1D*> efficiency;
  std::pair<TH1D*, TH1D*> purity;
  SystSummary* systematics;

  Histo1D(TH1D* hist, std::pair<THStack*, TLegend*> stack1D)
  {
    total_hist = hist;
    stacked_hist = stack1D.first;
    legend = stack1D.second;

    error_band = GetErrorBand(total_hist);
  }

  // Get the percentage statistical error per bin
  TH1D* GetErrorBand(TH1D* total_hist){

    TH1D *error_band = (TH1D*)total_hist->Clone();

    // Set the bin errors on seperate plot
    for (int n = 1; n <= total_hist->GetNbinsX(); n++){
     error_band->SetBinContent(n, 0);
     error_band->SetBinError(n, 0);
     if (total_hist->GetBinContent(n) > 0)
       error_band->SetBinContent(n, 100*total_hist->GetBinError(n)/total_hist->GetBinContent(n));
    }
    return error_band;

  }

};

#endif
