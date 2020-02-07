#ifndef HISTO2D_H
#define HISTO2D_H

#include "Histo1D.h"
#include "SystSummary2D.h"

// Structure for holding 2D histogram information
class Histo2D
{
  public:

  TString name;
  TH2D* total_hist;
  TH2D* empty;
  std::vector<THStack*> stacked_hist;
  TLegend* legend;
  std::pair<TH2D*, TH2D*> efficiency;
  std::pair<TH2D*, TH2D*> purity;
  TH2D* response;
  bool eff_calculated = false;
  SystSummary2D* systematics;

  // Constructor
  Histo2D(TH2D* hist, std::pair<std::vector<THStack*>, TLegend*> stack2D)
  {
    total_hist = hist;
    empty = (TH2D*)total_hist->Clone();
    empty->Reset();
    name = TString(total_hist->GetName());
    stacked_hist = stack2D.first;
    legend = stack2D.second;

    TH2D *eff_numerator = (TH2D*) empty->Clone(name+"_effnum");
    TH2D *eff_denom = (TH2D*) empty->Clone(name+"_effden");
    efficiency = std::make_pair(eff_numerator, eff_denom);

    TH2D *pur_numerator = (TH2D*) empty->Clone(name+"_purnum");
    TH2D *pur_denom = (TH2D*) empty->Clone(name+"_purden");
    purity = std::make_pair(pur_numerator, pur_denom);

    int nbins = total_hist->GetNbinsX()*total_hist->GetNbinsY();
    response = new TH2D(name+"_response", "", nbins, 1, nbins+1, nbins, 1, nbins+1);

    systematics = new SystSummary2D(empty);
  }

  // Slice in y axis bins
  Histo1D* Slice(size_t i){
    int bin = i;
    TString slice_name = name + Form("_%.1f_%.1f", total_hist->GetYaxis()->GetBinLowEdge(i), total_hist->GetYaxis()->GetBinLowEdge(i+1));
    TH1D* total = (TH1D*)total_hist->ProjectionX(slice_name, bin, bin)->Clone();
    THStack* stack = stacked_hist[i-1];
    std::pair<TH1D*, TH1D*> eff = std::make_pair(efficiency.first->ProjectionX(slice_name+"_en", bin, bin)
                                               , efficiency.second->ProjectionX(slice_name+"_ed", bin, bin));
    std::pair<TH1D*, TH1D*> pur = std::make_pair(purity.first->ProjectionX(slice_name+"_pn", bin, bin)
                                               , purity.second->ProjectionX(slice_name+"_pd", bin, bin));
    SystSummary* syst1D = systematics->Slice(i);
    Histo1D* sliced_histo = new Histo1D(total, std::make_pair(stack, legend), eff, pur, syst1D);
    return sliced_histo;
  }

};

#endif
