#ifndef HISTO2D_H
#define HISTO2D_H

#include "Histo1D.h"
#include "SystSummary2D.h"
#include "../Functions/PolySlicer.h"

// Structure for holding 2D histogram information
class Histo2D
{
  public:

  TString name;   // Histogram base name

  TH2Poly* total_hist;  // Total rate histogram
  TH2Poly* bkg_hist;    // Background histogram
  TH2Poly* xsec_hist;   // Cross section histogram

  std::vector<THStack*> stacked_hist;       // Stacked histogram by type
  TLegend* legend;                          // Legend for stacked histogram
  std::vector<double> ybins;                // Binning of slicing variable
  std::vector<std::vector<double>> xbins;   // Binning of other variable in slices

  std::pair<TH2Poly*, TH2Poly*> efficiency; // Numerator and denominator for efficiency
  std::pair<TH2Poly*, TH2Poly*> purity;     // Numerator and denominator for purity
  TH2D* response;                           // Response matrix
  SystSummary2D* systematics;               // Systematics holder

  // Constructor
  Histo2D(TH2Poly* hist, std::pair<std::vector<THStack*>, TLegend*> stack2D, std::vector<double> yb, std::vector<std::vector<double>> xb)
  {
    // Set total rate hist
    total_hist = hist;
    name = TString(total_hist->GetName());

    // Create an empty histogram for quick cloning
    TH2Poly* empty = (TH2Poly*)total_hist->Clone();
    empty->ClearBinContents();

    // Create background and cross section histograms
    bkg_hist = (TH2Poly*) empty->Clone(name+"_bkg");
    xsec_hist = (TH2Poly*) empty->Clone(name+"_xsec");

    // Set stacked histograms and legend
    stacked_hist = stack2D.first;
    legend = stack2D.second;
    // Set binning
    ybins = yb;
    xbins = xb;

    // Create efficiency and purity histograms
    TH2Poly *eff_numerator = (TH2Poly*) empty->Clone(name+"_effnum");
    TH2Poly *eff_denom = (TH2Poly*) empty->Clone(name+"_effden");
    efficiency = std::make_pair(eff_numerator, eff_denom);

    TH2Poly *pur_numerator = (TH2Poly*) empty->Clone(name+"_purnum");
    TH2Poly *pur_denom = (TH2Poly*) empty->Clone(name+"_purden");
    purity = std::make_pair(pur_numerator, pur_denom);

    // Create an empty response matrix
    int nbins = total_hist->GetNumberOfBins();
    response = new TH2D(name+"_response", "", nbins, 1, nbins+1, nbins, 1, nbins+1);

    // Create new systematics holder
    systematics = new SystSummary2D(total_hist, ybins, xbins);
  }


  // Slice in y axis bins
  Histo1D* Slice(size_t i){
    int bin = i;
    // Name the slice
    TString slice_name = name + Form("_%.1f_%.1f", ybins[i], ybins[i+1]);

    // Slice the 2D histograms
    TH1D* total = SlicePoly(total_hist, i, slice_name, ybins, xbins);
    TH1D* bkg = SlicePoly(bkg_hist, i, slice_name+"_bkg", ybins, xbins);
    TH1D* xsec = SlicePoly(xsec_hist, i, slice_name+"_xsec", ybins, xbins);
    xsec->SetLineWidth(3);

    // Get the relevant stacked histogram
    THStack* stack = stacked_hist[i];

    // Slice the efficiency and purity
    std::pair<TH1D*, TH1D*> eff = std::make_pair(SlicePoly(efficiency.first, i, slice_name+"_en", ybins, xbins)
                                               , SlicePoly(efficiency.second, i, slice_name+"_ed", ybins, xbins));
    std::pair<TH1D*, TH1D*> pur = std::make_pair(SlicePoly(purity.first, i, slice_name+"_pn", ybins, xbins)
                                               , SlicePoly(purity.second, i, slice_name+"_pd", ybins, xbins));

    // Slice the systematics
    SystSummary* syst1D = systematics->Slice(i);

    // Create a new 1D histogram
    Histo1D* sliced_histo = new Histo1D(total, bkg, xsec, std::make_pair(stack, legend), eff, pur, syst1D);

    return sliced_histo;
  }

};

#endif
