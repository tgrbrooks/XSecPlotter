#ifndef SYSTSUMMARY_H
#define SYSTSUMMARY_H

#include "Systematics.h"

// Structure for holding all systematic error information associated with 1D histogram
class SystSummary
{
  public:

  Systematics* genie;      // GENIE generator systematics
  Systematics* flux;       // Flux systematics
  Systematics* detector;   // Detector performance systematics
  Systematics* background; // External background systematics
  Systematics* constant;   // Constant systematics
  Systematics* total;      // Total systematics

  // Constructor
  SystSummary(TH1D* hist){
    genie = new Systematics(hist, "_geniesyst");
    flux = new Systematics(hist, "_fluxsyst");
    detector = new Systematics(hist, "_detsyst");
    background = new Systematics(hist, "_bkgsyst");
    constant = new Systematics(hist, "_constsyst");
    total = new Systematics(hist, "_totalsyst");
  }

  // Constructor (for 1D slices)
  SystSummary(Systematics *g, Systematics *f, Systematics *d, Systematics *b, Systematics *c, Systematics *t)
  {
    genie = g;
    flux = f;
    detector = d;
    background = b;
    constant = c;
    total = t;
  }

  // Return systematic by name
  Systematics* GetSyst(TString name){
    if(name == "genie") return genie;
    if(name == "flux") return flux;
    if(name == "detector") return detector;
    if(name == "background") return background;
    if(name == "constant") return constant;
    return total;
  }
  
  // Calculate the total error assuming uncorrelated
  void GetTotal(){
    ClearErrors(total);
    AddSyst(total, genie);
    AddSyst(total, flux);
    AddSyst(total, detector);
    AddSyst(total, background);
    AddSyst(total, constant);

    // Get the fractional covariance and correlation from covariance
    for(int i = 1; i <= total->mean_syst->GetNbinsX(); i++){
      // Central value in bin i
      double cv_i = total->mean_syst->GetBinContent(i);
      // Error on bin i
      double s_ii = std::sqrt(total->covariance->GetBinContent(i, i));

      for(int j = 1; j <= total->mean_syst->GetNbinsX(); j++){
        // Covariance of bin i and bin j
        double cov_ij = total->covariance->GetBinContent(i, j);
        // Central value in bin i
        double cv_j = total->mean_syst->GetBinContent(j);
        // Error on bin j
        double s_jj = std::sqrt(total->covariance->GetBinContent(j, j));
        // Fractional covariance and correlation calculation
        total->frac_covariance->SetBinContent(i, j, cov_ij/(cv_i*cv_j));
        total->correlation->SetBinContent(i, j, cov_ij/(s_ii*s_jj));
      }
    }
  }

  // Delete the systematic errors
  void ClearErrors(Systematics* s1){
    for(int i = 1; i <= s1->mean_syst->GetNbinsX(); i++){
      s1->mean_syst->SetBinError(i, 0);
    }
  }

  // Add the errors and covariance matrices
  void AddSyst(Systematics* s1, Systematics* s2){
    AddErrors(s1->mean_syst, s2->mean_syst);
    s1->covariance->Add(s2->covariance);
  }

  // Add histogram errors in quadrature, ignoring bin contents
  void AddErrors(TH1D* syst_hist, TH1D* hist){
    for(int i = 1; i <= syst_hist->GetNbinsX(); i++){
      double new_err = std::sqrt(std::pow(syst_hist->GetBinError(i),2)+std::pow(hist->GetBinError(i),2));
      syst_hist->SetBinError(i, new_err);
    }
  }

};

#endif
