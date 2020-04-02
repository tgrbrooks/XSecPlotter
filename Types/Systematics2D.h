#ifndef SYSTEMATICS2D_H
#define SYSTEMATICS2D_H

#include "../Functions/PolySlicer.h"
#include "XSecUniverse2D.h"

// Structure for holding 2D systematic error information
class Systematics2D
{
  public:

  TString sname;    // Systematics base name

  std::vector<TH2Poly*> universes;      // Universe variations for rate or cross section
  std::vector<XSecUniverse2D*> xsecuni; // Universe variations for cross section components

  TH2Poly* mean_syst;                     // Means of universe variations
  TH2Poly* std_syst;                      // Standard deviations of universe variations
  std::vector<double> ybins;              // Binning of the slicing variables
  std::vector<std::vector<double>> xbins; // Binning of other variables

  TH2D* covariance;      // Covariance matrix  
  TH2D* frac_covariance; // Fractional covariance matrix
  TH2D* correlation;     // Correlation matrix

  // Constructor
  Systematics2D(TH2Poly* hist, TString name, std::vector<double> yb, std::vector<std::vector<double>> xb){

    // Set the mean and standard deviation histograms
    mean_syst = (TH2Poly*)hist->Clone(TString(hist->GetName())+name);
    std_syst = (TH2Poly*)hist->Clone(TString(hist->GetName())+name+"_stddev");

    // Reset if not total
    if(name!="_totalsyst") mean_syst->ClearBinContents();
    std_syst->ClearBinContents();

    // Get the name
    sname = TString(mean_syst->GetName());
    ybins = yb;
    xbins = xb;

    // Create the covariance and fractional covariance matrices
    size_t nbins = mean_syst->GetNumberOfBins();
    covariance = new TH2D(sname+"_covariance", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    frac_covariance = new TH2D(sname+"_frac_covariance", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    correlation = new TH2D(sname+"_correlation", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
  }


  // Create 1D slice in Y bin
  Systematics* Slice(size_t i, bool xsec){
    int bin = i;

    // Universe slices
    std::vector<TH1D*> uni_s;
    // Name the universe
    TString slice_name = sname + Form("_%.1f_%.1f", ybins[i], ybins[i+1]);
    // Get the 1D slice for each universe
    for(size_t u = 0; u < universes.size(); u++){
      uni_s.push_back(SlicePoly(universes[u], i, Form(slice_name+"_uni%i", (int)u), ybins, xbins, xsec));
    }

    // Slice the mean and standard deviation histograms
    TH1D* mean_s = SlicePoly(mean_syst, i, slice_name, ybins, xbins, xsec);
    // Need to get the errors separately
    TH1D* std_s = SlicePoly(std_syst, i, slice_name+"_stddev", ybins, xbins, xsec);
    for(int x = 1; x <= mean_s->GetNbinsX(); x++) mean_s->SetBinError(x, std_s->GetBinContent(x));

    // Create 1D systematics
    Systematics* syst_s = new Systematics(sname, uni_s, mean_s);
    return syst_s;
  }
  

  // Create empty universes for variations
  void CreateUniverses(size_t nuni){
    universes.clear();
    for(size_t i = 0; i < nuni; i++){
      TH2Poly* uni = (TH2Poly*) mean_syst->Clone(Form(sname+"_uni%i",(int)i));
      uni->ClearBinContents();
      universes.push_back(uni);
    }
  }

  // Create empty cross section universes for variations
  void CreateXSecUni(size_t nuni){
    xsecuni.clear();
    for(size_t i = 0; i < nuni; i++){
      XSecUniverse2D* uni = new XSecUniverse2D(mean_syst, Form(sname+"_xsecuni%i",(int)i));
      xsecuni.push_back(uni);
    }
  }

  // Scale by appropriate factor
  void ScaleUniverses(Configuration* config, size_t file_i){

    for(size_t u = 0; u < universes.size(); u++){
      universes[u]->Scale(config->pot_scale_fac[file_i], "width");
    }
  }

  // Calculate mean, covariance and correlation from universe variations
  void Calculate(TH2Poly* cv_hist){

    size_t nbins = mean_syst->GetNumberOfBins();

    // For constant errors covariance is just diagonal variance matrix
    if(universes.size()==0){
      for(size_t i = 1; i <= nbins; i++){
        for(size_t j = 1; j <= nbins; j++){
          if(i==j){ 
            covariance->SetBinContent(i, j, pow(std_syst->GetBinContent(i), 2));
            frac_covariance->SetBinContent(i, j, pow(std_syst->GetBinContent(i),2)/pow(cv_hist->GetBinContent(i),2));
            correlation->SetBinContent(i, j, 1.);
          }
        }
      }
      return;
    }

    // Calculate the mean and standard deviation for each bin over all universes
    std::vector<double> means;
    for(size_t i = 1; i <= nbins; i++){
      double mean = 0;
      for(size_t ns = 0; ns < universes.size(); ns++){
        mean += universes[ns]->GetBinContent(i);
      }
      mean /= universes.size();
      means.push_back(mean);
      double std_dev = 0;
      for(size_t ns = 0; ns < universes.size(); ns++){
        //std_dev += std::pow(universes[ns]->GetBinContent(i) - mean, 2.);
        std_dev += std::pow(universes[ns]->GetBinContent(i) - cv_hist->GetBinContent(i), 2.);
      }
      std_dev = std::sqrt(std_dev/(universes.size()-1));
      mean_syst->SetBinContent(i, mean);
      std_syst->SetBinContent(i, std_dev);
    }


    // Calculate the covariance and correlation over all universes
    for(size_t i = 1; i <= nbins; i++){
      //double cv_i = means[i-1];
      double cv_i = cv_hist->GetBinContent(i);
      for(size_t j = 1; j <= nbins; j++){
        //double cv_j = means[j-1];
        double cv_j = cv_hist->GetBinContent(j);
        double E_ij = 0;
        for(size_t ns = 0; ns < universes.size(); ns++){
          E_ij += (universes[ns]->GetBinContent(i)-cv_i)*(universes[ns]->GetBinContent(j)-cv_j);
        }
        E_ij /= universes.size();
        covariance->SetBinContent(i, j, E_ij);
        frac_covariance->SetBinContent(i, j, E_ij/(cv_i*cv_j));
      }
    }

    for(size_t i = 1; i <= nbins; i++){
      for(size_t j = 1; j <= nbins; j++){
        double corr_bin = covariance->GetBinContent(i, j)/(std::sqrt(covariance->GetBinContent(i,i))*std::sqrt(covariance->GetBinContent(j,j)));
        correlation->SetBinContent(i, j, corr_bin);
      }
    }
  }

};

#endif
