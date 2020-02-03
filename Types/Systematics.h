#ifndef SYSTEMATICS_H
#define SYSTEMATICS_H

// Structure for holding systematic error information
class Systematics
{
  public:

  std::vector<TH1D*> universes;
  TH1D* mean_syst;
  TH2D* covariance;
  TH2D* frac_covariance;
  TH2D* correlation;

  Systematics(){}

  Systematics(TH1D* m){
    mean_syst = m;
  }

  Systematics(std::vector<TH1D*> u){
    universes = u;

    mean_syst = (TH1D*)universes[0]->Clone();
    mean_syst->Reset();

    // Calculate the mean and standard deviation for each bin over all universes
    std::vector<double> means;
    for(size_t n = 1; n <= universes[0]->GetNbinsX(); n++){
      double mean = 0;
      for(size_t ns = 0; ns < u.size(); ns++){
        mean += universes[ns]->GetBinContent(n);
      }
      mean /= universes.size();
      means.push_back(mean);
      double std_dev = 0;
      for(size_t ns = 0; ns < universes.size(); ns++){
        std_dev += std::pow(universes[ns]->GetBinContent(n) - mean, 2.);
      }
      std_dev = std::sqrt(std_dev/(universes.size()-1));
      mean_syst->SetBinContent(n, mean);
      mean_syst->SetBinError(n, std_dev);
    }

    size_t nbins = mean_syst->GetNbinsX();
    TString name = mean_syst->GetName();

    // Calculate the covariance and correlation over all universes
    covariance = new TH2D(name+"_covariance", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    frac_covariance = new TH2D(name+"_frac_covariance", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    correlation = new TH2D(name+"_correlation", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    for(size_t i = 1; i <= universes[0]->GetNbinsX(); i++){
      double cv_i = means[i-1];
      for(size_t j = 1; j <= universes[0]->GetNbinsX(); j++){
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

    for(size_t i = 1; i <= universes[0]->GetNbinsX(); i++){
      for(size_t j = 1; j <= universes[0]->GetNbinsX(); j++){
        double corr_bin = covariance->GetBinContent(i, j)/(std::sqrt(covariance->GetBinContent(i,i))*std::sqrt(covariance->GetBinContent(j,j)));
        correlation->SetBinContent(i, j, corr_bin);
      }
    }
  }

};

#endif
