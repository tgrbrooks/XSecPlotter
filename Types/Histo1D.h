#ifndef HISTO1D_H
#define HISTO1D_H

#include "SystSummary.h"

// Structure for holding 1D histogram information
class Histo1D
{
  public:

  TString name;     // Histogram base name

  TH1D* total_hist; // Total rate histogram
  TH1D* bkg_hist;   // Background histogram
  TH1D* xsec_hist;  // Cross section histogram
  TH1D* error_band; // Percentage error histogram

  THStack* stacked_hist;              // Stacked histogram by type
  TLegend* legend;                    // Legend for stacked histogram
  std::pair<TH1D*, TH1D*> efficiency; // Numerator and denominator for efficiency
  std::pair<TH1D*, TH1D*> purity;     // Numerator and denominator for purity
  TH2D* response;                     // Response matrix
  SystSummary* systematics;           // Systematics holder

  // Constructor (used for total data)
  Histo1D(TH1D* hist)
  {
    // Set the rate histogram
    total_hist = hist;
    name = TString(total_hist->GetName());

    // Create an empty histogram for quicker cloning
    TH1D* empty = (TH1D*) total_hist->Clone();
    empty->Reset();

    // Initialise background and cross section histograms as empty
    bkg_hist = (TH1D*) empty->Clone(name+"_bkg");
    xsec_hist = (TH1D*) empty->Clone(name+"_xsec");

    // Calculate the percentag error band
    GetErrorBand();

    // Initialise efficiency, purity as empty histograms
    TH1D *eff_numerator = (TH1D*) empty->Clone(name+"_effnum");
    TH1D *eff_denom = (TH1D*) empty->Clone(name+"_effden");
    efficiency = std::make_pair(eff_numerator, eff_denom);

    TH1D *pur_numerator = (TH1D*) empty->Clone(name+"_purnum");
    TH1D *pur_denom = (TH1D*) empty->Clone(name+"_purden");
    purity = std::make_pair(pur_numerator, pur_denom);

    // Initialise systematics holder
    systematics = new SystSummary(total_hist);

    delete empty;
  }

  // Constructor (used for normal 1D histograms)
  Histo1D(TH1D* hist, std::pair<THStack*, TLegend*> stack1D)
  {
    // Set the rate histogram
    total_hist = hist;
    name = TString(total_hist->GetName());

    // Create and empty histogram for quicker cloning
    TH1D* empty = (TH1D*) total_hist->Clone();
    empty->Reset();

    // Initialise background and cross section histograms as empty
    bkg_hist = (TH1D*) empty->Clone(name+"_bkg");
    xsec_hist = (TH1D*) empty->Clone(name+"_xsec");

    // Set the stacked histogram and legend
    stacked_hist = stack1D.first;
    legend = stack1D.second;

    // Get the percentage error band
    GetErrorBand();

    // Initialise efficiency and purity as empty histograms
    TH1D *eff_numerator = (TH1D*) empty->Clone(name+"_effnum");
    TH1D *eff_denom = (TH1D*) empty->Clone(name+"_effden");
    efficiency = std::make_pair(eff_numerator, eff_denom);

    TH1D *pur_numerator = (TH1D*) empty->Clone(name+"_purnum");
    TH1D *pur_denom = (TH1D*) empty->Clone(name+"_purden");
    purity = std::make_pair(pur_numerator, pur_denom);

    // Create an empty response matrix
    int nbins = total_hist->GetNbinsX();
    response = new TH2D(name+"_response", "", nbins, 1, nbins+1, nbins, 1, nbins+1);

    // Initialise systematics holder
    systematics = new SystSummary(total_hist);

    delete empty;
  }

  // Constructor (used for 1D slices of 2D histograms)
  Histo1D(TH1D* hist, TH1D* bkg, TH1D* xsec, std::pair<THStack*, TLegend*> stack1D, std::pair<TH1D*, TH1D*> eff, std::pair<TH1D*, TH1D*> pur, SystSummary *syst)
  {
    // Set all members from inputs
    total_hist = hist;
    name = TString(total_hist->GetName());

    bkg_hist = bkg;
    xsec_hist = xsec;
    stacked_hist = stack1D.first;
    legend = stack1D.second;
    efficiency = eff;
    purity = pur;
    systematics = syst;

    // Calculate the percentage error band
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

  // Get the percentage statistical and systematic error per bin
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

  // Get the percentage statistical and systematic error per bin for cross sections
  void XSecErrorBand(TString systname, bool stat){
    for (int n = 1; n <= total_hist->GetNbinsX(); n++){
      double stat_error = 0;
      if(stat && (systname=="total" || systname=="")) stat_error = xsec_hist->GetBinError(n);
      double syst_error = 0;
      if(systname != "") syst_error = systematics->GetSyst(systname)->mean_syst->GetBinError(n);
      double error = std::sqrt(std::pow(stat_error, 2) + std::pow(syst_error, 2));
      if (xsec_hist->GetBinContent(n) > 0)
        error_band->SetBinContent(n, 100*error/xsec_hist->GetBinContent(n));
      else error_band->SetBinContent(n, 0);
    }
  }

  // Set the efficiency and purity
  void SetEfficiency(std::pair<TH1D*, TH1D*> eff, std::pair<TH1D*, TH1D*> pur){
    efficiency = eff;
    purity = pur;
  }

  // Print error summary for total data
  void PrintSummary(int bin){
    double total = total_hist->GetBinContent(bin);
    double stat_e = total_hist->GetBinError(bin);
    double genie_e = systematics->genie->mean_syst->GetBinError(bin);
    double flux_e = systematics->flux->mean_syst->GetBinError(bin);
    double det_e = systematics->detector->mean_syst->GetBinError(bin);
    double bkg_e = systematics->background->mean_syst->GetBinError(bin);
    double const_e = systematics->constant->mean_syst->GetBinError(bin);
    double syst_e = systematics->total->mean_syst->GetBinError(bin);
    double tot_e = std::sqrt(pow(syst_e, 2)+pow(stat_e, 2));
    std::cout<<"Total rate = "<<total<<"\n"
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

  // Print error summary for total data
  void PrintXSecSummary(int bin){
    double total = xsec_hist->GetBinContent(bin);
    double stat_e = xsec_hist->GetBinError(bin);
    double genie_e = systematics->genie->mean_syst->GetBinError(bin);
    double flux_e = systematics->flux->mean_syst->GetBinError(bin);
    double det_e = systematics->detector->mean_syst->GetBinError(bin);
    double bkg_e = systematics->background->mean_syst->GetBinError(bin);
    double const_e = systematics->constant->mean_syst->GetBinError(bin);
    double syst_e = systematics->total->mean_syst->GetBinError(bin);
    double tot_e = std::sqrt(pow(syst_e, 2)+pow(stat_e, 2));
    std::cout<<"Total cross section = "<<total<<"\n"
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
