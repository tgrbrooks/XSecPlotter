#ifndef XSECCALCULATOR_H
#define XSECCALCULATOR_H

#include "Configuration.h"
#include "FluxManager.h"
#include "Histo1D.h"
#include "Histo2D.h"

// Structure for holding interaction information
class XSecCalculator
{
  public:

  Configuration *config; // Global configurations
  FluxManager *fluxman;

  // Constructor
  XSecCalculator(Configuration *c, FluxManager *f)
  {
    config = c;
    fluxman = f;
  }

  TH1D* ToXSec(Histo1D* histo){
    TH1D* xsec_hist = (TH1D*)histo->total_hist->Clone(TString(histo->total_hist->GetName())+"_xsec");
    // Subtract background
    xsec_hist->Add(histo->bkg_hist, -1.);
    // Preserve the statistical errors as percentage
    // Correct efficiency with response
    TH1D* eff = GetEfficiency(histo->efficiency, histo->response);
    xsec_hist->Divide(eff);
    double nt = config->targets;
    double flux = fluxman->IntegratedFlux();
    double scale = 1e38/(nt * flux);
    xsec_hist->Scale(scale);
    return xsec_hist;
  }

  // Apply a response matrix to a true histogram to get a reconstructed histogram
  TH1D* ApplyResponse(TH1D* hist, TH2D* response){
    // N_i = R_ij * N_j

    // Create a reconstructed histogram with the same binning as the true one
    TH1D* reco_hist = (TH1D*) hist->Clone();
    reco_hist->Reset();

    // Loop over the bins of the reco histogram
    for(size_t bin_i = 1; bin_i <= reco_hist->GetNbinsX(); bin_i++){
      double content_i = 0;
      // Loop over bins of the true histogram and apply the response matrix
      for(size_t bin_j = 1; bin_j <= hist->GetNbinsX()+1; bin_j++){
        content_i += response->GetBinContent(bin_j, bin_i) * hist->GetBinContent(bin_j);
      }   
      reco_hist->SetBinContent(bin_i, content_i);
    }
    return reco_hist;
  }

  // Apply response matrix to efficiency
  TH1D* GetEfficiency(std::pair<TH1D*, TH1D*> eff, TH2D* response){

    TH1D *selected_resp = ApplyResponse(eff.first, response);
    TH1D *true_resp = ApplyResponse(eff.second, response);

    selected_resp->Divide(true_resp);

    delete true_resp;
    return selected_resp;

  }

};

#endif
