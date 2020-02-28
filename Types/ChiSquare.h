#ifndef CHISQUARE_H
#define CHISQUARE_H

#include "TMatrix.h"

#include "Configuration.h"
#include "Systematics.h"
#include "Histo1D.h"
#include "Histo2D.h"

// Class for handling chi2 calculations between models
class ChiSquare
{
  public:

  // Constructor
  ChiSquare(){}

  // Calculate chi2 between fake data and MC for 1D distributions
  std::pair<double, int> Calculate(TH1D* data, Histo1D* mc, bool xsec=false){

    int nbins = mc->systematics->total->covariance->GetNbinsX();

    // Create covariance matrix to be inverted
    TMatrix cov;
    cov.Clear();
    cov.ResizeTo(nbins, nbins);

    // Fill covariance matrix
    for (int i = 1; i <= nbins; i ++) {
      for (int j = 1; j <= nbins; j ++) {
        cov[i-1][j-1] = mc->systematics->total->covariance->GetBinContent(i, j);
        if(i==j) cov[i-1][j-1] += pow(data->GetBinError(i), 2);
      }
    }

    cov *= 1./cov[0][0];
    // Invert the covariance matrix
    TMatrix cov_inv = cov.Invert();
    cov_inv *= cov[0][0];

    // Loop over all combinations of bins
    double chi2 = 0.;
    for (int i = 1; i <= nbins; i++) {
      for (int j = 1; j <= nbins; j++) {

        double data_i = data->GetBinContent(i);
        double mc_i = mc->total_hist->GetBinContent(i);
        if(xsec) mc_i = mc->xsec_hist->GetBinContent(i);

        double data_j = data->GetBinContent(j);
        double mc_j = mc->total_hist->GetBinContent(j);
        if(xsec) mc_j = mc->xsec_hist->GetBinContent(j);

        chi2 += (data_i - mc_i) * cov_inv[i-1][j-1] * (data_j - mc_j);

      }
    }
    return std::make_pair(chi2, nbins-1);
  }

  // Calculate chi2 between fake data and MC for 2D distributions TODO xsec stat errors
  std::pair<double, int> Calculate(TH2Poly* data, Histo2D* mc, bool xsec=false){

    int nbins = mc->xsec_hist->GetNumberOfBins();

    // Find all of the filled bins
    std::vector<int> filled_bins;
    for(int i = 1; i <= nbins; i++){
      // Doesn't work if there is no entry in MC bin
      if(mc->systematics->total->covariance->GetBinContent(i, i)+data->GetBinError(i) < 1e-6
         || mc->total_hist->GetBinContent(i) < 1e-6) continue;
      filled_bins.push_back(i);
    }

    // Create covariance matrix for inversion
    TMatrix cov;
    cov.Clear();
    cov.ResizeTo(filled_bins.size(), filled_bins.size());

    // Fill the covariance matrix
    for (size_t i = 0; i < filled_bins.size(); i ++) {
      int bin_i = filled_bins[i];
      for (size_t j = 0; j < filled_bins.size(); j ++) {
        int bin_j = filled_bins[j];
        cov[i][j] = mc->systematics->total->covariance->GetBinContent(bin_i, bin_j);
        if(bin_i == bin_j) cov[i][j] += pow(data->GetBinError(bin_i), 2); 
      }
    }
    //cov.Print();
    cov *= 1./cov[0][0];
    std::cout<<"Det = "<<cov.Determinant()<<"\n";

    // Invert the covariance matrix
    TMatrix cov_inv = cov.Invert();
    cov_inv += cov[0][0];
    
    // Loop over all combinations of bins
    double chi2 = 0.;
    for (size_t i = 0; i < filled_bins.size(); i++) {
      int bin_i = filled_bins[i];
      
      for (size_t j = 0; j < filled_bins.size(); j++) {
        int bin_j = filled_bins[j];

        double data_i = data->GetBinContent(bin_i);
        double mc_i = mc->total_hist->GetBinContent(bin_i);
        if(xsec) mc_i = mc->xsec_hist->GetBinContent(bin_i);

        double data_j = data->GetBinContent(bin_j);
        double mc_j = mc->total_hist->GetBinContent(bin_j);
        if(xsec) mc_j = mc->xsec_hist->GetBinContent(bin_j);

        chi2 += (data_i - mc_i) * cov_inv[i][j] * (data_j - mc_j);

      }
    }
    return std::make_pair(chi2, filled_bins.size()-1);
  }
  
};

#endif
