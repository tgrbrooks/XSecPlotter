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
  std::pair<double, int> Calculate(TH2D* data, Histo2D* mc){

    int nxbins = mc->total_hist->GetNbinsX();
    int nybins = mc->total_hist->GetNbinsY();
    int nbins = mc->systematics->total->covariance->GetNbinsX();

    std::vector<int> filled_bins;
    for(int i = 1; i <= nbins; i++){
      int i_x = ceil((double)i/nybins);
      int i_y = i - nybins*(i_x-1);
      // Doesn't work if there is no entry in MC bin
      if(mc->systematics->total->covariance->GetBinContent(i, i)+data->GetBinError(i_x, i_y)< 1e-6
         || mc->total_hist->GetBinContent(i_x, i_y)< 1e-6) continue;
      filled_bins.push_back(i);
    }

    TMatrix cov;
    cov.Clear();
    cov.ResizeTo(filled_bins.size(), filled_bins.size());

    for (int i = 0; i < filled_bins.size(); i ++) {
      int bin_i = filled_bins[i];
      int i_x = ceil((double)bin_i/nybins);
      int i_y = bin_i - nybins*(i_x-1);
      for (int j = 0; j < filled_bins.size(); j ++) {
        int bin_j = filled_bins[j];
        cov[i][j] = mc->systematics->total->covariance->GetBinContent(bin_i, bin_j);
        if(bin_i == bin_j) cov[i][j] += pow(data->GetBinError(i_x, i_y), 2); 
      }
    }
    cov.Print();
    std::cout<<"Det = "<<cov.Determinant()<<"\n";

    TMatrix cov_inv = cov.Invert();
    
    double chi2 = 0.;
    for (int i = 0; i < filled_bins.size(); i++) {
      int bin_i = filled_bins[i];
      int i_x = ceil((double)bin_i/nybins);
      int i_y = bin_i - nybins*(i_x-1);
      
      for (int j = 0; j < filled_bins.size(); j++) {
        int bin_j = filled_bins[j];
        int j_x = ceil((double)bin_j/nybins);
        int j_y = bin_j - nybins*(j_x-1);

        double data_i = data->GetBinContent(i_x, i_y);
        double mc_i = mc->total_hist->GetBinContent(i_x, i_y);

        double data_j = data->GetBinContent(j_x, j_y);
        double mc_j = mc->total_hist->GetBinContent(j_x, j_y);

        chi2 += (data_i - mc_i) * cov_inv[i][j] * (data_j - mc_j);

      }
    }
    return std::make_pair(chi2, filled_bins.size()-1);
  }
  
};

#endif
