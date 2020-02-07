#ifndef SYSTSUMMARY2D_H
#define SYSTSUMMARY2D_H

#include "Systematics2D.h"

// Structure for holding all systematic error information associated with 2D histograms
class SystSummary2D
{
  public:

  Systematics2D* genie;
  Systematics2D* flux;
  Systematics2D* detector;
  Systematics2D* background;
  Systematics2D* constant;
  Systematics2D* total;

  // Constructor
  SystSummary2D(TH2D* hist){
    genie = new Systematics2D(hist, "_geniesyst");
    flux = new Systematics2D(hist, "_fluxsyst");
    detector = new Systematics2D(hist, "_detsyst");
    background = new Systematics2D(hist, "_bkgsyst");
    constant = new Systematics2D(hist, "_constsyst");
    total = new Systematics2D(hist, "_totalsyst");
    //total = (TH2D*)hist->Clone(TString(hist->GetName())+"_totalsyst");
    //total->Reset();
  }

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
    AddSyst(total, genie);
    AddSyst(total, flux);
    AddSyst(total, detector);
    AddSyst(total, background);
    AddSyst(total, constant);

    int nxbins = total->mean_syst->GetNbinsX();
    int nybins = total->mean_syst->GetNbinsY();
    int nbins = nxbins*nybins;

    for(size_t i = 1; i <= nbins; i++){
      int i_x = ceil((double)i/nybins);
      int i_y = i - nybins*(i_x-1);
      double cv_i = total->mean_syst->GetBinContent(i_x, i_y);
      double s_ii = std::sqrt(total->covariance->GetBinContent(i, i));
      for(size_t j = 1; j <= nbins; j++){
        int j_x = ceil((double)j/nybins);
        int j_y = j - nybins*(j_x-1);
        double cov_ij = total->covariance->GetBinContent(i, j);
        double cv_j = total->mean_syst->GetBinContent(j_x, j_y);
        double s_jj = std::sqrt(total->covariance->GetBinContent(j, j));
        total->frac_covariance->SetBinContent(i, j, cov_ij/(cv_i*cv_j));
        total->correlation->SetBinContent(i, j, cov_ij/(s_ii*s_jj));
        if(cov_ij==0||(s_ii*s_jj)==0||isnan(s_ii*s_jj)) total->correlation->SetBinContent(i, j, 0);
        if(i==j) total->correlation->SetBinContent(i, j, 1);
      }
    }
  }

  void AddSyst(Systematics2D* s1, Systematics2D* s2){
    AddErrors(s1->mean_syst, s2->mean_syst);
    s1->covariance->Add(s2->covariance);
  }

  // Add histogram errors in quadrature, ignoring bin contents
  void AddErrors(TH2D* syst_hist, TH2D* hist){
    for(size_t i = 1; i <= syst_hist->GetNbinsX(); i++){
      for(size_t j = 1; j <= syst_hist->GetNbinsY(); j++){
        double new_err = std::sqrt(std::pow(syst_hist->GetBinError(i, j),2)+std::pow(hist->GetBinError(i, j),2));
        syst_hist->SetBinError(i, j, new_err);
      }
    }
  }

  // Get 1D systematics for slice in Y bin
  SystSummary* Slice(size_t i){
    Systematics* genie_s = genie->Slice(i);
    Systematics* flux_s = flux->Slice(i);
    Systematics* detector_s = detector->Slice(i);
    Systematics* background_s = background->Slice(i);
    Systematics* constant_s = constant->Slice(i);
    Systematics* total_s = total->Slice(i);
    SystSummary* summary = new SystSummary(genie_s, flux_s, detector_s, background_s, constant_s, total_s);
    return summary;
  }

};

#endif
