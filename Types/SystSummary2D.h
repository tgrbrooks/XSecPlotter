#ifndef SYSTSUMMARY2D_H
#define SYSTSUMMARY2D_H

#include "Systematics2D.h"

// Structure for holding all systematic error information associated with 2D histograms
class SystSummary2D
{
  public:

  Systematics2D* genie;      // GENIE generator systematics
  Systematics2D* flux;       // Flux systematics
  Systematics2D* detector;   // Detector performance systematics
  Systematics2D* background; // External background systematics
  Systematics2D* constant;   // Constant systematics
  Systematics2D* total;      // Total systematics

  // Constructor
  SystSummary2D(TH2Poly* hist, std::vector<double> yb, std::vector<std::vector<double>> xb){
    genie = new Systematics2D(hist, "_geniesyst", yb, xb);
    flux = new Systematics2D(hist, "_fluxsyst", yb, xb);
    detector = new Systematics2D(hist, "_detsyst", yb, xb);
    background = new Systematics2D(hist, "_bkgsyst", yb, xb);
    constant = new Systematics2D(hist, "_constsyst", yb, xb);
    total = new Systematics2D(hist, "_totalsyst", yb, xb);
  }

  // Get systematic by name
  Systematics2D* GetSyst(TString name){
    if(name == "genie") return genie;
    if(name == "flux") return flux;
    if(name == "detector") return detector;
    if(name == "background") return background;
    if(name == "constant") return constant;
    return total;
  }

  // Calculate total errors assuming uncorrelated
  void GetTotal(){
    ClearErrors(total);
    AddSyst(total, genie);
    AddSyst(total, flux);
    AddSyst(total, detector);
    AddSyst(total, background);
    AddSyst(total, constant);

    int nbins = total->mean_syst->GetNumberOfBins();

    for(int i = 1; i <= nbins; i++){
      // Central value in bin i
      double cv_i = total->mean_syst->GetBinContent(i);
      // Error on bin i
      double s_ii = std::sqrt(total->covariance->GetBinContent(i, i));

      for(int j = 1; j <= nbins; j++){
        // Covariance of bin i and bin j
        double cov_ij = total->covariance->GetBinContent(i, j);
        // Central value in bin j
        double cv_j = total->mean_syst->GetBinContent(j);
        // Error on bin j
        double s_jj = std::sqrt(total->covariance->GetBinContent(j, j));
        // Fractional covariance and correlation calculation
        total->frac_covariance->SetBinContent(i, j, cov_ij/(cv_i*cv_j));
        total->correlation->SetBinContent(i, j, cov_ij/(s_ii*s_jj));
        if(i==j) total->correlation->SetBinContent(i, j, 1);
      }
    }
  }
  
  // Delete the systematic errors
  void ClearErrors(Systematics2D* s1){
    for(int i = 1; i <= s1->std_syst->GetNumberOfBins(); i++){
      s1->std_syst->SetBinContent(i, 0);
    }
  }

  // Add the errors and covariance matrices
  void AddSyst(Systematics2D* s1, Systematics2D* s2){
    AddErrors(s1->std_syst, s2->std_syst);
    s1->covariance->Add(s2->covariance);
  }

  // Add histogram errors in quadrature, ignoring bin contents
  void AddErrors(TH2Poly* syst_hist, TH2Poly* hist){
    for(int i = 1; i <= syst_hist->GetNumberOfBins(); i++){
      double new_err = std::sqrt(std::pow(syst_hist->GetBinContent(i),2)+std::pow(hist->GetBinContent(i),2));
      syst_hist->SetBinContent(i, new_err);
    }
  }

  // Get 1D systematics for slice in Y bin
  SystSummary* Slice(size_t i, bool xsec){
    Systematics* genie_s = genie->Slice(i, xsec);
    Systematics* flux_s = flux->Slice(i, xsec);
    Systematics* detector_s = detector->Slice(i, xsec);
    Systematics* background_s = background->Slice(i, xsec);
    Systematics* constant_s = constant->Slice(i, xsec);
    Systematics* total_s = total->Slice(i, xsec);
    SystSummary* summary = new SystSummary(genie_s, flux_s, detector_s, background_s, constant_s, total_s);
    return summary;
  }

};

#endif
