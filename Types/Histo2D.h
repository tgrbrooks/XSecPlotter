#ifndef HISTO2D_H
#define HISTO2D_H

#include "Histo1D.h"
#include "SystSummary2D.h"
#include "../Functions/PolySlicer.h"

// Structure for holding 2D histogram information
class Histo2D
{
  public:

  TString name;
  TH2Poly* total_hist;
  TH2Poly* empty;
  std::vector<THStack*> stacked_hist;
  std::vector<double> ybins;
  std::vector<std::vector<double>> xbins;
  TLegend* legend;
  std::pair<TH2Poly*, TH2Poly*> efficiency;
  std::pair<TH2Poly*, TH2Poly*> purity;
  TH2D* response;
  bool eff_calculated = false;
  SystSummary2D* systematics;

  // Constructor
  Histo2D(TH2Poly* hist, std::pair<std::vector<THStack*>, TLegend*> stack2D, std::vector<double> yb, std::vector<std::vector<double>> xb)
  {
    total_hist = hist;
    empty = (TH2Poly*)total_hist->Clone();
    empty->ClearBinContents();
    name = TString(total_hist->GetName());
    stacked_hist = stack2D.first;
    ybins = yb;
    xbins = xb;
    legend = stack2D.second;

    TH2Poly *eff_numerator = (TH2Poly*) empty->Clone(name+"_effnum");
    TH2Poly *eff_denom = (TH2Poly*) empty->Clone(name+"_effden");
    efficiency = std::make_pair(eff_numerator, eff_denom);

    TH2Poly *pur_numerator = (TH2Poly*) empty->Clone(name+"_purnum");
    TH2Poly *pur_denom = (TH2Poly*) empty->Clone(name+"_purden");
    purity = std::make_pair(pur_numerator, pur_denom);

    int nbins = total_hist->GetNumberOfBins();
    response = new TH2D(name+"_response", "", nbins, 1, nbins+1, nbins, 1, nbins+1);

    systematics = new SystSummary2D(empty, ybins, xbins);
  }


  // Slice in y axis bins
  Histo1D* Slice(size_t i){
    int bin = i;
    TString slice_name = name + Form("_%.1f_%.1f", ybins[i], ybins[i+1]);
    TH1D* total = SlicePoly(total_hist, i, slice_name, ybins, xbins);
    THStack* stack = stacked_hist[i];
    std::pair<TH1D*, TH1D*> eff = std::make_pair(SlicePoly(efficiency.first, i, slice_name+"_en", ybins, xbins)
                                               , SlicePoly(efficiency.second, i, slice_name+"_ed", ybins, xbins));
    std::pair<TH1D*, TH1D*> pur = std::make_pair(SlicePoly(purity.first, i, slice_name+"_pn", ybins, xbins)
                                               , SlicePoly(purity.second, i, slice_name+"_pd", ybins, xbins));
    SystSummary* syst1D = systematics->Slice(i);
    Histo1D* sliced_histo = new Histo1D(total, std::make_pair(stack, legend), eff, pur, syst1D);
    return sliced_histo;
  }

};

#endif
