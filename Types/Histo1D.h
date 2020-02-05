#ifndef HISTO1D_H
#define HISTO1D_H

#include "SystSummary.h"

// Structure for holding 1D histogram information
class Histo1D
{
  public:

  TString name;
  TH1D* total_hist;
  TH1D* empty;
  TH1D* error_band;
  THStack* stacked_hist;
  TLegend* legend;
  std::pair<TH1D*, TH1D*> efficiency;
  std::pair<TH1D*, TH1D*> purity;
  SystSummary* systematics;

  Histo1D(TH1D* hist)
  {
    total_hist = hist;
    empty = (TH1D*) total_hist->Clone();
    empty->Reset();
    name = TString(total_hist->GetName());
    error_band = GetErrorBand(total_hist);

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

  Histo1D(TH1D* hist, std::pair<THStack*, TLegend*> stack1D)
  {
    total_hist = hist;
    empty = (TH1D*) total_hist->Clone();
    empty->Reset();
    name = TString(total_hist->GetName());
    stacked_hist = stack1D.first;
    legend = stack1D.second;

    error_band = GetErrorBand(total_hist);

    // Initialise efficiency, purity and systematics
    TH1D *eff_numerator = (TH1D*) empty->Clone(name+"_effnum");
    TH1D *eff_denom = (TH1D*) empty->Clone(name+"_effden");
    efficiency = std::make_pair(eff_numerator, eff_denom);

    TH1D *pur_numerator = (TH1D*) empty->Clone(name+"_purnum");
    TH1D *pur_denom = (TH1D*) empty->Clone(name+"_purden");
    purity = std::make_pair(pur_numerator, pur_denom);

    // Initialise empty systematics
    systematics = new SystSummary(empty);
  }

  Histo1D(TH1D* hist, std::pair<THStack*, TLegend*> stack1D, std::pair<TH1D*, TH1D*> eff, std::pair<TH1D*, TH1D*> pur, SystSummary *syst)
  {
    total_hist = hist;
    name = TString(total_hist->GetName());
    stacked_hist = stack1D.first;
    legend = stack1D.second;
    efficiency = eff;
    purity = pur;
    systematics = syst;

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

  void SetEfficiency(std::pair<TH1D*, TH1D*> eff, std::pair<TH1D*, TH1D*> pur){
    efficiency = eff;
    purity = pur;
  }

  void PrintSummary(){
    double total = total_hist->GetBinContent(1);
    double stat_e = total_hist->GetBinError(1);
    double genie_e = systematics->genie->mean_syst->GetBinError(1);
    double flux_e = systematics->flux->mean_syst->GetBinError(1);
    double det_e = systematics->detector->mean_syst->GetBinError(1);
    double bkg_e = systematics->background->mean_syst->GetBinError(1);
    double const_e = systematics->constant->mean_syst->GetBinError(1);
    double syst_e = systematics->total->GetBinError(1);
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
