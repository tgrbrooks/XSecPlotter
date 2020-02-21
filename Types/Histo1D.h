#ifndef HISTO1D_H
#define HISTO1D_H

#include "SystSummary.h"

// Structure for holding 1D histogram information
class Histo1D
{
  public:

  TString name;
  TH1D* total_hist;
  TH1D* bkg_hist;
  TH1D* empty;
  TH1D* error_band;
  double scale;
  THStack* stacked_hist;
  TLegend* legend;
  std::pair<TH1D*, TH1D*> efficiency;
  std::pair<TH1D*, TH1D*> purity;
  TH2D* response;
  SystSummary* systematics;

  // Constructor (used for total data)
  Histo1D(TH1D* hist)
  {
    total_hist = hist;
    empty = (TH1D*) total_hist->Clone();
    empty->Reset();
    name = TString(total_hist->GetName());
    bkg_hist = (TH1D*) empty->Clone(name+"_bkg");
    GetErrorBand();

    // Initialise efficiency, purity and systematics
    TH1D *eff_numerator = (TH1D*) empty->Clone(name+"_effnum");
    TH1D *eff_denom = (TH1D*) empty->Clone(name+"_effden");
    efficiency = std::make_pair(eff_numerator, eff_denom);

    TH1D *pur_numerator = (TH1D*) empty->Clone(name+"_purnum");
    TH1D *pur_denom = (TH1D*) empty->Clone(name+"_purden");
    purity = std::make_pair(pur_numerator, pur_denom);

    // Initialise empty systematics
    systematics = new SystSummary(total_hist);
  }

  // Constructor (used for normal 1D histograms)
  Histo1D(TH1D* hist, std::pair<THStack*, TLegend*> stack1D)
  {
    total_hist = hist;
    empty = (TH1D*) total_hist->Clone();
    empty->Reset();
    name = TString(total_hist->GetName());
    bkg_hist = (TH1D*) empty->Clone(name+"_bkg");
    stacked_hist = stack1D.first;
    legend = stack1D.second;

    GetErrorBand();

    // Initialise efficiency, purity and systematics
    TH1D *eff_numerator = (TH1D*) empty->Clone(name+"_effnum");
    TH1D *eff_denom = (TH1D*) empty->Clone(name+"_effden");
    efficiency = std::make_pair(eff_numerator, eff_denom);

    TH1D *pur_numerator = (TH1D*) empty->Clone(name+"_purnum");
    TH1D *pur_denom = (TH1D*) empty->Clone(name+"_purden");
    purity = std::make_pair(pur_numerator, pur_denom);

    int nbins = total_hist->GetNbinsX();
    response = new TH2D(name+"_response", "", nbins, 1, nbins+1, nbins, 1, nbins+1);

    // Initialise empty systematics
    systematics = new SystSummary(total_hist);
  }

  // Constructor (used for 1D slices of 2D histograms)
  Histo1D(TH1D* hist, TH1D* bkg, std::pair<THStack*, TLegend*> stack1D, std::pair<TH1D*, TH1D*> eff, std::pair<TH1D*, TH1D*> pur, SystSummary *syst)
  {
    total_hist = hist;
    bkg_hist = bkg;
    name = TString(total_hist->GetName());
    stacked_hist = stack1D.first;
    legend = stack1D.second;
    efficiency = eff;
    purity = pur;
    systematics = syst;

    GetErrorBand();
  }

  // Get the percentage statistical error per bin
  void GetErrorBand(){

    error_band = (TH1D*)total_hist->Clone();

    // Set the bin errors on seperate plot
    for (int n = 1; n <= total_hist->GetNbinsX(); n++){
     error_band->SetBinContent(n, 0);
     error_band->SetBinError(n, 0);
     if (total_hist->GetBinContent(n) > 0)
       error_band->SetBinContent(n, 100*total_hist->GetBinError(n)/total_hist->GetBinContent(n));
    }

  }

  void SystErrorBand(TString systname, bool stat){
    for (int n = 1; n <= total_hist->GetNbinsX(); n++){
      double stat_error = 0;
      if(stat && systname=="total") stat_error = total_hist->GetBinError(n);
      double syst_error = systematics->GetSyst(systname)->mean_syst->GetBinError(n);
      double error = std::sqrt(std::pow(stat_error, 2) + std::pow(syst_error, 2));
      if (total_hist->GetBinContent(n) > 0)
        error_band->SetBinContent(n, 100*error/total_hist->GetBinContent(n));
    }
  }

  // Set the efficiency and purity
  void SetEfficiency(std::pair<TH1D*, TH1D*> eff, std::pair<TH1D*, TH1D*> pur){
    efficiency = eff;
    purity = pur;
  }

  // Print error summary for total data
  void PrintSummary(){
    double total = total_hist->GetBinContent(1);
    double stat_e = total_hist->GetBinError(1);
    double genie_e = systematics->genie->mean_syst->GetBinError(1);
    double flux_e = systematics->flux->mean_syst->GetBinError(1);
    double det_e = systematics->detector->mean_syst->GetBinError(1);
    double bkg_e = systematics->background->mean_syst->GetBinError(1);
    double const_e = systematics->constant->mean_syst->GetBinError(1);
    double syst_e = systematics->total->mean_syst->GetBinError(1);
    double tot_e = std::sqrt(pow(syst_e, 2)+pow(stat_e, 2));
    std::cout<<"Total = "<<total<<"\n"
             <<"-------------------------------------------------\n"
             <<"Statistical error = "<<stat_e<<" ("<<100*stat_e/total<<"%)\n"
             <<"-------------------------------------------------\n"
             <<"Systematic errors:\n"
             <<"Genie      = "<<genie_e<<" ("<<100*genie_e/total<<"%)\n"
             <<"Flux       = "<<flux_e<<" ("<<100*flux_e/total<<"%)\n"
             <<"Detector   = "<<det_e<<" ("<<100*det_e/total<<"%)\n"
             <<"Background = "<<bkg_e<<" ("<<100*bkg_e/total<<"%)\n"
             <<"Constant   = "<<const_e<<" ("<<100*const_e/total<<"%)\n"
             <<"-------------------------------------------------\n"
             <<"Systematic error = "<<syst_e<<" ("<<100*syst_e/total<<"%)\n"
             <<"-------------------------------------------------\n"
             <<"Total error = "<<tot_e<<" ("<<100*tot_e/total<<"%)\n";
  }

};

#endif
