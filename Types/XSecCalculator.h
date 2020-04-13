#ifndef XSECCALCULATOR_H
#define XSECCALCULATOR_H

#include "Configuration.h"
#include "FluxManager.h"
#include "Histo1D.h"
#include "Histo2D.h"
#include "XSecUniverse.h"
#include "XSecUniverse2D.h"

// Calculator for getting cross section from rates
class XSecCalculator
{
  public:

  Configuration *config; // Global configurations
  FluxManager *fluxman;  // Reweighted flux holder

  size_t file_i;  // File index

  // Constructor
  XSecCalculator(Configuration *c, FluxManager *f, size_t fi)
  {
    config = c;
    fluxman = f;
    file_i = fi;
  }

  // Cross section from 1D histogram ("data")
  TH1D* ToXSec(Histo1D* histo){
    
    // Clone the rate
    TH1D* xsec_hist = (TH1D*)histo->total_hist->Clone(TString(histo->total_hist->GetName())+"_xsec");

    // Subtract background
    xsec_hist->Add(histo->bkg_hist, -1.);

    // Correct efficiency with response
    TH1D* eff = GetEfficiency(histo->efficiency, histo->response);
    xsec_hist->Divide(eff);

    // Get and apply the scale factor
    double nt = config->targets;
    double flux = fluxman->IntegratedFlux();
    double scale = 1e38/(nt * flux);
    xsec_hist->Scale(scale);

    if(xsec_hist->GetNbinsX()==1){
      std::cout<<"Total = "<<histo->total_hist->GetBinContent(1)<<"\n"
               <<"Background = "<<histo->bkg_hist->GetBinContent(1)<<"\n"
               <<"Efficiency = "<<eff->GetBinContent(1)<<"\n"
               <<"Targets = "<<nt<<"\n"
               <<"Flux = "<<flux<<"\n"
               <<"XSec = "<<xsec_hist->GetBinContent(1)<<"\n";
    }
    
    return xsec_hist;
  }

  // Cross section from 2D histogram ("data") TODO errors
  std::pair<TH2Poly*, TH2Poly*> ToXSec(Histo2D* histo){

    // Clone the rate
    TH2Poly* xsec_hist = (TH2Poly*)histo->total_hist->Clone(TString(histo->total_hist->GetName())+"_xsec");
    TH2Poly* xsec_err = (TH2Poly*)histo->total_hist->Clone(TString(histo->total_hist->GetName())+"_xsecerr");

    // Store percentage errors
    std::vector<double> perrs;
    for(int i = 1; i <= xsec_hist->GetNumberOfBins(); i++){
      perrs.push_back(xsec_hist->GetBinError(i)/xsec_hist->GetBinContent(i));
    };

    // Subtract background
    for(int i = 0; i <= xsec_hist->GetNumberOfBins()+1; i++){
      xsec_hist->SetBinContent(i, xsec_hist->GetBinContent(i) - histo->bkg_hist->GetBinContent(i));
    }

    // Correct efficiency with response
    TH2Poly* eff = GetEfficiency(histo->efficiency, histo->response);
    for(int i = 0; i <= xsec_hist->GetNumberOfBins()+1; i++){
      xsec_hist->SetBinContent(i, xsec_hist->GetBinContent(i) / eff->GetBinContent(i));
    }

    // Need to manually scale by bin width because TH2Poly is trash
    for(auto const& obj : *xsec_hist->GetBins()){
      TH2PolyBin *bin = (TH2PolyBin*)obj;
      double wy = abs(bin->GetYMax() - bin->GetYMin());
      double wx = abs(bin->GetXMax() - bin->GetXMin());
      double width = wy*wx;
      int j = bin->GetBinNumber();
      xsec_hist->SetBinContent(j, xsec_hist->GetBinContent(j)/width);
    }

    // Get and apply the scale factor
    double nt = config->targets;
    double flux = fluxman->IntegratedFlux();
    double scale = 1e38/(nt * flux);
    xsec_hist->Scale(scale);

    // Reapply percentage errors
    for(int i = 1; i <= xsec_hist->GetNumberOfBins(); i++){
      xsec_err->SetBinContent(i, perrs[i-1]*xsec_hist->GetBinContent(i));
    };
    
    return std::make_pair(xsec_hist, xsec_err);
  }

  // Cross section from 1D universe variation
  TH1D* ToXSec(XSecUniverse* xsecuni, TH1D* rate, TString syst, int uni){

    // Scale xsec universe to the required POT
    xsecuni->Scale(config->pot_scale_fac[file_i]);

    // Clone the total rate from data
    TH1D* xsec_hist = (TH1D*)rate->Clone(Form(TString(rate->GetName())+syst+"_unitemp%i", uni));
    // Subtract background
    xsec_hist->Add(xsecuni->background, -1.);

    // Correct efficiency with response
    std::pair<TH1D*, TH1D*> efficiency = std::make_pair(xsecuni->selected, xsecuni->generated);
    if(xsecuni->selected->GetNbinsX() > 1){
      TH1D* eff = GetEfficiency(efficiency, xsecuni->Response());
      xsec_hist->Divide(eff);
    }
    // Don't apply response for total hist
    else{
      TH1D* eff = (TH1D*)efficiency.first->Clone();
      eff->Divide(efficiency.second);
      xsec_hist->Divide(eff);
    }

    // Get and apply scale
    double nt = config->targets;
    double flux = fluxman->IntegratedFlux();
    if(syst=="flux") flux = fluxman->IntegratedFlux(uni);
    double scale = 1e38/(nt * flux);
    xsec_hist->Scale(scale);

    return xsec_hist;
  }

  // Cross section from 2D universe variation TODO errors
  TH2Poly* ToXSec(XSecUniverse2D* xsecuni, TH2Poly* rate, TString syst, int uni){
    // Scale xsec universe to the required POT
    xsecuni->Scale(config->pot_scale_fac[file_i]);

    // Clone the total rate from "data"
    TH2Poly* xsec_hist = (TH2Poly*)rate->Clone(Form(TString(rate->GetName())+syst+"_unitemp%i", uni));
    // Subtract background
    for(int i = 0; i <= xsec_hist->GetNumberOfBins()+1; i++){
      xsec_hist->SetBinContent(i, xsec_hist->GetBinContent(i) - xsecuni->background->GetBinContent(i));
    }

    // Correct efficiency with response
    std::pair<TH2Poly*, TH2Poly*> efficiency = std::make_pair(xsecuni->selected, xsecuni->generated);
    TH2Poly* eff = GetEfficiency(efficiency, xsecuni->Response());
    for(int i = 0; i <= xsec_hist->GetNumberOfBins()+1; i++){
      xsec_hist->SetBinContent(i, xsec_hist->GetBinContent(i) / eff->GetBinContent(i));
    }

    // Need to manually scale by bin width because TH2Poly is trash
    for(auto const& obj : *xsec_hist->GetBins()){
      TH2PolyBin *bin = (TH2PolyBin*)obj;
      double wy = abs(bin->GetYMax() - bin->GetYMin());
      double wx = abs(bin->GetXMax() - bin->GetXMin());
      double width = wy*wx;
      int j = bin->GetBinNumber();
      xsec_hist->SetBinContent(j, xsec_hist->GetBinContent(j)/width);
    }

    // Get and apply scale
    double nt = config->targets;
    double flux = fluxman->IntegratedFlux();
    if(syst=="flux") flux = fluxman->IntegratedFlux(uni);
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
    for(int bin_i = 1; bin_i <= reco_hist->GetNbinsX(); bin_i++){
      double content_i = 0;
      // Loop over bins of the true histogram and apply the response matrix
      for(int bin_j = 1; bin_j <= hist->GetNbinsX()+1; bin_j++){
        content_i += response->GetBinContent(bin_j, bin_i) * hist->GetBinContent(bin_j);
      }   
      reco_hist->SetBinContent(bin_i, content_i);
    }
    return reco_hist;
  }

  // Apply response matrix to efficiency
  TH1D* GetEfficiency(std::pair<TH1D*, TH1D*> eff, TH2D* response){

    // If looking at total hist don't apply response
    if(eff.first->GetNbinsX() <= 1){
      TH1D* efficiency = (TH1D*)eff.first->Clone();
      efficiency->Divide(eff.second);
      return efficiency;
    }

    TH1D *selected_resp = ApplyResponse(eff.first, response);
    TH1D *true_resp = ApplyResponse(eff.second, response);

    selected_resp->Divide(true_resp);

    delete true_resp;
    return selected_resp;

  }

  // Apply a response matrix to a true histogram to get a reconstructed histogram
  TH2Poly* ApplyResponse(TH2Poly* hist, TH2D* response){
    // N_i = R_ij * N_j

    // Create a reconstructed histogram with the same binning as the true one
    TH2Poly* reco_hist = (TH2Poly*) hist->Clone();
    reco_hist->ClearBinContents();

    // Loop over the bins of the reco histogram
    for(int bin_i = 1; bin_i <= reco_hist->GetNumberOfBins(); bin_i++){
      double content_i = 0;
      // Loop over bins of the true histogram and apply the response matrix
      for(int bin_j = 1; bin_j <= hist->GetNumberOfBins()+1; bin_j++){
        content_i += response->GetBinContent(bin_j, bin_i) * hist->GetBinContent(bin_j);
      }   
      reco_hist->SetBinContent(bin_i, content_i);
    }
    return reco_hist;
  }

  // Apply response matrix to efficiency TODO errors
  TH2Poly* GetEfficiency(std::pair<TH2Poly*, TH2Poly*> eff, TH2D* response){

    if(eff.first->GetNbinsX() <= 1){
      TH2Poly* efficiency = (TH2Poly*)eff.first->Clone();
      for(int i = 0; i <= efficiency->GetNumberOfBins()+1; i++){
        efficiency->SetBinContent(i, efficiency->GetBinContent(i) / eff.first->GetBinContent(i));
      }
      return efficiency;
    }

    TH2Poly *selected_resp = ApplyResponse(eff.first, response);
    TH2Poly *true_resp = ApplyResponse(eff.second, response);

    for(int i = 0; i <= selected_resp->GetNumberOfBins()+1; i++){
      selected_resp->SetBinContent(i, selected_resp->GetBinContent(i) / true_resp->GetBinContent(i));
    }

    delete true_resp;
    return selected_resp;

  }

};

#endif
