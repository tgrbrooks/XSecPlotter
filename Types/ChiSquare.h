#ifndef CHISQUARE_H
#define CHISQUARE_H

#include "TMatrix.h"

#include "Configuration.h"
#include "Systematics.h"
#include "Histo1D.h"
#include "Histo2D.h"

// Class for handling all plotting
class ChiSquare
{
  public:

  // Constructor
  ChiSquare(){}

  // Calculate chi2 between fake data and MC for 1D distributions
  std::pair<double, int> Calculate(TH1D* data, Histo1D* mc){

    int nbins = mc->systematics->total->covariance->GetNbinsX();

    TMatrix cov;
    cov.Clear();
    cov.ResizeTo(nbins, nbins);

    for (int i = 1; i <= nbins; i ++) {
      for (int j = 1; j <= nbins; j ++) {
        cov[i-1][j-1] = mc->systematics->total->covariance->GetBinContent(i, j);
        if(i==j) cov[i-1][j-1] += pow(data->GetBinError(i), 2);
      }
    }

    TMatrix cov_inv = cov.Invert();

    double chi2 = 0.;
    for (int i = 1; i <= nbins; i++) {
      for (int j = 1; j <= nbins; j++) {

        double data_i = data->GetBinContent(i);
        double mc_i = mc->total_hist->GetBinContent(i);

        double data_j = data->GetBinContent(j);
        double mc_j = mc->total_hist->GetBinContent(j);

        chi2 += (data_i - mc_i) * cov_inv[i-1][j-1] * (data_j - mc_j);

      }
    }
    return std::make_pair(chi2, nbins-1);
  }

  // Calculate chi2 between fake data and MC for 2D distributions
  std::pair<double, int> Calculate(TH2Poly* data, Histo2D* mc){

    int nbins = mc->total_hist->GetNumberOfBins();
    double scale = 5000;

    std::vector<int> filled_bins;
    for(int i = 1; i <= nbins; i++){
      // Doesn't work if there is no entry in MC bin
      if(mc->systematics->total->covariance->GetBinContent(i, i)+data->GetBinError(i)< 1e-6
         || mc->total_hist->GetBinContent(i)< 1e-6) continue;
      filled_bins.push_back(i);
    }

    TMatrix cov;
    cov.Clear();
    cov.ResizeTo(filled_bins.size(), filled_bins.size());

    for (size_t i = 0; i < filled_bins.size(); i ++) {
      int bin_i = filled_bins[i];
      for (size_t j = 0; j < filled_bins.size(); j ++) {
        int bin_j = filled_bins[j];
        cov[i][j] = mc->systematics->total->covariance->GetBinContent(bin_i, bin_j);
        if(bin_i == bin_j) cov[i][j] += pow(data->GetBinError(bin_i), 2); 
      }
    }
    //cov.Print();
    std::cout<<"Det = "<<cov.Determinant()<<"\n";

    TMatrix cov_inv = cov.Invert();
    
    double chi2 = 0.;
    for (size_t i = 0; i < filled_bins.size(); i++) {
      int bin_i = filled_bins[i];
      
      for (size_t j = 0; j < filled_bins.size(); j++) {
        int bin_j = filled_bins[j];

        double data_i = data->GetBinContent(bin_i);
        double mc_i = mc->total_hist->GetBinContent(bin_i);

        double data_j = data->GetBinContent(bin_j);
        double mc_j = mc->total_hist->GetBinContent(bin_j);

        chi2 += (data_i - mc_i) * cov_inv[i][j] * (data_j - mc_j);

      }
    }
    return std::make_pair(chi2, filled_bins.size()-1);
  }
  
};

#endif
