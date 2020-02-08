#ifndef SYSTEMATICS2D_H
#define SYSTEMATICS2D_H

#include "../Functions/PolySlicer.h"

// Structure for holding 2D systematic error information
class Systematics2D
{
  public:

  TString sname;
  std::vector<TH2Poly*> universes;
  TH2Poly* mean_syst;
  std::vector<double> ybins;
  std::vector<std::vector<double>> xbins;
  TH2D* covariance;
  TH2D* frac_covariance;
  TH2D* correlation;

  // Constructor
  Systematics2D(TH2Poly* hist, TString name, std::vector<double> yb, std::vector<std::vector<double>> xb){
    mean_syst = (TH2Poly*)hist->Clone(TString(hist->GetName())+name);
    sname = TString(mean_syst->GetName());
    mean_syst->ClearBinContents();
    ybins = yb;
    xbins = xb;
    size_t nbins = mean_syst->GetNumberOfBins();
    covariance = new TH2D(sname+"_covariance", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    frac_covariance = new TH2D(sname+"_frac_covariance", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    correlation = new TH2D(sname+"_correlation", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
  }


  // Create 1D slice in Y bin
  Systematics* Slice(size_t i){
    int bin = i;
    std::vector<TH1D*> uni_s;
    TString slice_name = sname + Form("_%.1f_%.1f", ybins[i], ybins[i+1]);
    for(size_t u = 0; u < universes.size(); u++){
      uni_s.push_back(SlicePoly(universes[u], i, Form(slice_name+"_uni%i", (int)u), ybins, xbins));
    }
    TH1D* mean_s = SlicePoly(mean_syst, i, slice_name, ybins, xbins);
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

  // Scale by appropriate factor
  void ScaleUniverses(Configuration* config, size_t file_i){
    double xsec_scale = 1e38/(config->flux[file_i] * config->targets);

    for(size_t u = 0; u < universes.size(); u++){
      universes[u]->Scale(config->pot_scale_fac[file_i]);
      if(config->plot_xsec){
        universes[u]->Scale(xsec_scale, "width");
      }
      else if(config->max_error > 0 || config->bin_edges[0].size() > 1){
        universes[u]->Scale(1, "width");
      }
    }
  }

  // Calculate mean, covariance and correlation from universe variations
  void Calculate(){

    size_t nbins = mean_syst->GetNumberOfBins();

    // For constant errors covariance is just diagonal variance matrix
    if(universes.size()==0){
      for(size_t i = 1; i <= nbins; i++){
        for(size_t j = 1; j <= nbins; j++){
          if(i==j){ 
            covariance->SetBinContent(i, j, pow(mean_syst->GetBinError(i),2));
            frac_covariance->SetBinContent(i, j, pow(mean_syst->GetBinError(i),2)/pow(mean_syst->GetBinContent(i),2));
            correlation->SetBinContent(i, j, 1.);
          }
        }
      }
      return;
    }

    mean_syst->ClearBinContents();

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
        std_dev += std::pow(universes[ns]->GetBinContent(i) - mean, 2.);
      }
      std_dev = std::sqrt(std_dev/(universes.size()-1));
      mean_syst->SetBinContent(i, mean);
      mean_syst->SetBinError(i, std_dev);
    }


    // Calculate the covariance and correlation over all universes
    for(size_t i = 1; i <= nbins; i++){
      double cv_i = means[i-1];
      for(size_t j = 1; j <= nbins; j++){
        double cv_j = means[j-1];
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
