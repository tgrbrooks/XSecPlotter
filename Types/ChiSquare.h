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

  Configuration *config;
  TRandom3 *randgen;

  // Constructor
  ChiSquare(Configuration *c, TRandom3 *r)
  {
    config = c;
    randgen = r;
  }

  // -------------------------------------------------------------------------------------------------
  //                                    1D HISTOGRAMS
  // -------------------------------------------------------------------------------------------------

  // Calculate chi2 between fake data and MC for 1D distributions
  std::pair<double, int> Calculate(TH1D* data, Histo1D* mc){

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

    double scale = 1./cov[0][0];
    cov *= scale;
    // Invert the covariance matrix
    TMatrix cov_inv = cov.Invert();
    cov_inv *= scale;

    // Loop over all combinations of bins
    double chi2 = 0.;
    for (int i = 1; i <= nbins; i++) {
      for (int j = 1; j <= nbins; j++) {

        double data_i = data->GetBinContent(i);
        double mc_i = mc->total_hist->GetBinContent(i);
        if(config->plot_xsec) mc_i = mc->xsec_hist->GetBinContent(i);

        double data_j = data->GetBinContent(j);
        double mc_j = mc->total_hist->GetBinContent(j);
        if(config->plot_xsec) mc_j = mc->xsec_hist->GetBinContent(j);

        chi2 += (data_i - mc_i) * cov_inv[i-1][j-1] * (data_j - mc_j);

      }
    }
    return std::make_pair(chi2, nbins-1);
  }

  // Calculate chi2 between fake data and MC for 1D distributions with no covariance
  std::pair<double, int> Calculate(TH1D* data, TH1D* mc){

    int nbins = mc->GetNbinsX();

    // Loop over all combinations of bins
    double chi2 = 0.;
    for (int i = 1; i <= nbins; i++) {

      double data_i = data->GetBinContent(i);
      double mc_i = mc->GetBinContent(i);

      chi2 += (data_i - mc_i) * (data_i - mc_i) / data_i;

    }
    return std::make_pair(chi2, nbins-1);
  }

  // -------------------------------------------------------------------------------------------------
  //                                    2D HISTOGRAMS
  // -------------------------------------------------------------------------------------------------

  // Calculate chi2 between fake data and MC for 2D distributions TODO xsec stat errors
  std::pair<double, int> Calculate(TH2Poly* data, Histo2D* mc){

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
    // Scale covariance matrix to around unity to make the inversion not fail
    double scale = 1./cov[0][0];
    cov *= scale;
    //std::cout<<"Det = "<<cov.Determinant()<<"\n";
    // Invert the covariance matrix
    TMatrix cov_inv = cov.Invert();
    // Scale the inverted covariance matrix back
    cov_inv *= scale;
    
    // Loop over all combinations of bins
    double chi2 = 0.;
    for (size_t i = 0; i < filled_bins.size(); i++) {
      int bin_i = filled_bins[i];
      
      for (size_t j = 0; j < filled_bins.size(); j++) {
        int bin_j = filled_bins[j];

        double data_i = data->GetBinContent(bin_i);
        double mc_i = mc->total_hist->GetBinContent(bin_i);
        if(config->plot_xsec) mc_i = mc->xsec_hist->GetBinContent(bin_i);

        double data_j = data->GetBinContent(bin_j);
        double mc_j = mc->total_hist->GetBinContent(bin_j);
        if(config->plot_xsec) mc_j = mc->xsec_hist->GetBinContent(bin_j);

        chi2 += (data_i - mc_i) * cov_inv[i][j] * (data_j - mc_j);

      }
    }
    return std::make_pair(chi2, filled_bins.size()-1);
  }

  // -------------------------------------------------------------------------------------------------
  //                                        1D P VALUES
  // -------------------------------------------------------------------------------------------------

  double PValue(TH1D* data, Histo1D* mc){

    int N = 1000000;
    int n = 0;
    int nbins = mc->total_hist->GetNbinsX();
    // Calculate chi square between data and mc
    std::pair<double, double> data_chi2 = Calculate(data, mc);

    double maxbin = max(3.*nbins, data_chi2.first+10.);
    TH1D* chi_dist = new TH1D("chi_dist", "", 100, 0, maxbin);
    std::vector<double> chi_vec;
    // Loop over N universes
    for(int i = 0; i < N; i++){
      TH1D* uni = (TH1D*) mc->total_hist->Clone("uni");
      uni->Reset();
      // Loop over bins in histogram
      for(int b = 1; b <= nbins; b++){
        // Generate Poisson distributed random number for bin with mean as mc value
        double mean = mc->total_hist->GetBinContent(b);
        double rand = randgen->Poisson(mean);
        // TODO if looking at cross section is random variable a gaussian? what's the standard deviation? stat or syst?
        // Fill histogram for each universe
        uni->SetBinContent(b, rand);
      }
      // Calculate chi square statistic for universe
      std::pair<double, double> chi2 = Calculate(uni, mc->total_hist);
      // Fill histogram of chi squares
      chi_dist->Fill(chi2.first);
      chi_vec.push_back(chi2.first);
      delete uni;
    }

    // Count number of universes with more extreme chi square, n
    for(auto const& chi : chi_vec){
      if(chi >= data_chi2.first) n++;
    }

    // Draw chi2 distribution
    bool draw = true;
    if(draw){
      TCanvas *canvas = new TCanvas("temp", "", 900, 600);
      canvas->SetMargin(0.15, 0.04, 0.15, 0.15);
      chi_dist->GetXaxis()->SetTitle("#chi^{2}");
      chi_dist->GetYaxis()->SetTitle("P(#chi^{2})");
      chi_dist->SetLineWidth(3);
      chi_dist->SetLineColor(46);
      chi_dist->Scale(1./N);
      chi_dist->Draw("HIST C");
      double hmax = chi_dist->GetMaximum();
      TLine *line = new TLine(data_chi2.first, 0, data_chi2.first, hmax);
      line->SetLineStyle(9);
      line->Draw("same");
      TString output_file = config->output_file;
      output_file.ReplaceAll(".","_"+TString(mc->total_hist->GetName())+"_chi2.");
      canvas->SaveAs(output_file);
      delete canvas;
      delete line;
    }

    // Return n/N
    delete chi_dist;
    return (double)n/N;
  }

  
  
};

#endif
