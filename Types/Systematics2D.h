#ifndef SYSTEMATICS2D_H
#define SYSTEMATICS2D_H

// Structure for holding 2D systematic error information
class Systematics2D
{
  public:

  TString sname;
  std::vector<TH2D*> universes;
  TH2D* mean_syst;
  TH2D* covariance;
  TH2D* frac_covariance;
  TH2D* correlation;

  //Systematics2D(){}

  // Constructor
  Systematics2D(TH2D* hist, TString name){
    mean_syst = (TH2D*)hist->Clone(TString(hist->GetName())+name);
    sname = TString(mean_syst->GetName());
    mean_syst->Reset();
    size_t nbins = mean_syst->GetNbinsX()*mean_syst->GetNbinsY();
    covariance = new TH2D(sname+"_covariance", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    frac_covariance = new TH2D(sname+"_frac_covariance", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    correlation = new TH2D(sname+"_correlation", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
  }
/*
  Systematics2D(TH2D* m, bool empty=false){
    mean_syst = m;
    if(empty) mean_syst->Reset();
  }
*/
/*
  Systematics2D(std::vector<TH2D*> u){
    universes = u;
    mean_syst = (TH2D*)universes[0]->Clone();
    Calculate();
  }
  */

  // Create 1D slice in Y bin
  Systematics* Slice(size_t i){
    int bin = i;
    std::vector<TH1D*> uni_s;
    TString slice_name = sname + Form("_%.1f_%.1f", mean_syst->GetYaxis()->GetBinLowEdge(i), mean_syst->GetYaxis()->GetBinLowEdge(i+1));
    for(size_t u = 0; u < universes.size(); u++){
      uni_s.push_back(universes[u]->ProjectionX(Form(slice_name+"_uni%i", (int)u), bin, bin));
    }
    TH1D* mean_s = mean_syst->ProjectionX(slice_name, bin, bin);
    Systematics* syst_s = new Systematics(sname, uni_s, mean_s);
    return syst_s;
  }

  // Create empty universes for variations
  void CreateUniverses(size_t nuni){
    universes.clear();
    for(size_t i = 0; i < nuni; i++){
      TH2D* uni = (TH2D*) mean_syst->Clone(Form(sname+"_uni%i",(int)i));
      uni->Reset();
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

    size_t nxbins = mean_syst->GetNbinsX();
    size_t nybins = mean_syst->GetNbinsY();
    size_t nbins = nxbins*nybins;

    // For constant errors covariance is just diagonal variance matrix
    if(universes.size()==0){
      for(size_t i = 1; i <= nbins; i++){
        size_t i_x = ceil((double)i/nybins); // number of multiples of y
        size_t i_y = i - nybins*(i_x-1);
        for(size_t j = 1; j <= nbins; j++){
          if(i==j){ 
            covariance->SetBinContent(i, j, pow(mean_syst->GetBinError(i_x, i_y),2));
            frac_covariance->SetBinContent(i, j, pow(mean_syst->GetBinError(i_x, i_y),2)/pow(mean_syst->GetBinContent(i_x, i_y),2));
            correlation->SetBinContent(i, j, 1.);
          }
        }
      }
      return;
    }

    mean_syst->Reset();

    // Calculate the mean and standard deviation for each bin over all universes
    std::vector<double> means;
    for(size_t x = 1; x <= mean_syst->GetNbinsX(); x++){
      for(size_t y = 1; y <= mean_syst->GetNbinsY(); y++){
        double mean = 0;
        for(size_t ns = 0; ns < universes.size(); ns++){
          mean += universes[ns]->GetBinContent(x, y);
        }
        mean /= universes.size();
        means.push_back(mean);
        double std_dev = 0;
        for(size_t ns = 0; ns < universes.size(); ns++){
          std_dev += std::pow(universes[ns]->GetBinContent(x, y) - mean, 2.);
        }
        std_dev = std::sqrt(std_dev/(universes.size()-1));
        mean_syst->SetBinContent(x, y, mean);
        mean_syst->SetBinError(x, y, std_dev);
      }
    }


    // Calculate the covariance and correlation over all universes
    for(size_t i = 1; i <= nbins; i++){
      double cv_i = means[i-1];
      size_t i_x = ceil((double)i/nybins); // number of multiples of y
      size_t i_y = i - nybins*(i_x-1);
      for(size_t j = 1; j <= nbins; j++){
        double cv_j = means[j-1];
        size_t j_x = ceil((double)j/nybins);
        size_t j_y = j - nybins*(j_x-1);
        double E_ij = 0;
        for(size_t ns = 0; ns < universes.size(); ns++){
          E_ij += (universes[ns]->GetBinContent(i_x, i_y)-cv_i)*(universes[ns]->GetBinContent(j_x, j_y)-cv_j);
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
