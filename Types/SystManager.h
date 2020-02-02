#ifndef SYSTMANAGER_H
#define SYSTMANAGER_H

#include "Configuration.h"
#include "PlotManager.h"

// Structure for holding interaction information
class SystManager
{
  public:

  Configuration *config;
  PlotManager *plotman;

  SystManager(Configuration *c, PlotManager *p)
  {
    config = c;
    plotman = p;
  }

  double BkgSubtractionError(double mid, double width, TH1D* bkg){

    int bin = bkg->GetXaxis()->FindBin(mid);
    double bin_width = bkg->GetBinWidth(bin);
    double scale = (width/bin_width) * (config->pot_scale/6.6e20);
    double percent_error = bkg->GetBinError(bin)/bkg->GetBinContent(bin);
    double subtracted = bkg->GetBinContent(bin) * scale;
    double subtraction_error = subtracted * percent_error;
    if(std::isnan(subtraction_error)) subtraction_error = 0.;

    return subtraction_error;
  }

  TH1D* GetBkgSystHist(TString name, std::vector<std::vector<double>> bin_edges, int tune_i, int i){

    double edges_array[bin_edges[i].size()];
    std::copy(bin_edges[i].begin(), bin_edges[i].end(), edges_array);
    TH1D *syst_hist = new TH1D("bkgsyst"+name+"_"+config->tune_name[tune_i], "", bin_edges[i].size()-1, edges_array);

    // Background templates should all be scales to 6.6e20
    TFile *bkg_file = new TFile("BackgroundTemplates.root");
    TH1D* hMomCos = (TH1D*)bkg_file->Get("hMomCosErr");
    TH1D* hCosThetaCos = (TH1D*)bkg_file->Get("hCosThetaCosErr");
    TH1D* hMomDirt = (TH1D*)bkg_file->Get("hMomDirtErr");
    TH1D* hCosThetaDirt = (TH1D*)bkg_file->Get("hCosThetaDirtErr");

    // Loop over the actual binning
    for(size_t n = 1; n <= syst_hist->GetNbinsX(); n++){
      double mid = syst_hist->GetBinCenter(n);
      double width = syst_hist->GetBinWidth(n);
      // Determine the plotting variable
      if(config->plot_variables[i]=="lep_mom"){
        // Determine the bin of the background template
        double cos_sub_err = BkgSubtractionError(mid, width, hMomCos);
        double dirt_sub_err = BkgSubtractionError(mid, width, hMomDirt);
        syst_hist->SetBinError(n, std::sqrt(std::pow(cos_sub_err, 2.)+std::pow(dirt_sub_err, 2.)));
      }
      else if(config->plot_variables[i]=="cos_lep_theta"){
        double cos_sub_err = BkgSubtractionError(mid, width, hCosThetaCos);
        double dirt_sub_err = BkgSubtractionError(mid, width, hCosThetaDirt);
        syst_hist->SetBinError(n, std::sqrt(std::pow(cos_sub_err, 2.)+std::pow(dirt_sub_err, 2.)));
      }
    }

    
    double error;
    double integral = syst_hist->IntegralAndError(0, syst_hist->GetNbinsX()+1, error);
    std::cout<<"External background = "<<integral<<" +/- "<<error<<"\n";

    return syst_hist;

  }

  TH1D* GetDetSystHist(TString name, std::vector<std::vector<double>> bin_edges, int tune_i, int i, int j = -1, int bin_j = -1, int k = -1, int bin_k = -1){

    int nsims = 50;

    double edges_array[bin_edges[i].size()];
    std::copy(bin_edges[i].begin(), bin_edges[i].end(), edges_array);

    TH1D *syst_hist = new TH1D("detsyst"+name+"_"+config->tune_name[tune_i], "", bin_edges[i].size()-1, edges_array);

    std::vector<TH1D*> syst_hists;
    for(size_t ns = 0; ns < nsims; ns++){
      TH1D *syst_hist = new TH1D(Form("detsyst"+name+"_"+config->tune_name[tune_i]+"%i",ns), "", bin_edges[i].size()-1, edges_array);
      syst_hists.push_back(syst_hist);
    }
    std::vector<double> total(nsims, 0);

    Selection sel(config);

    TFile data_file(config->input_file[tune_i]);

    TTreeReader tree_reader("XSecTree/detsyst", &data_file);
    TTreeReaderValue<double> vtx_x(tree_reader, "ds_vtx_x");
    TTreeReaderValue<double> vtx_y(tree_reader, "ds_vtx_y");
    TTreeReaderValue<double> vtx_z(tree_reader, "ds_vtx_z");
    TTreeReaderArray<bool> particles_contained(tree_reader, "ds_particles_contained");
    TTreeReaderArray<bool> lep_contained(tree_reader, "ds_lep_contained");
    TTreeReaderArray<int> cc(tree_reader, "ds_cc");
    TTreeReaderArray<int> nu_pdg(tree_reader, "ds_nu_pdg");
    TTreeReaderArray<double> lep_mom(tree_reader, "ds_lep_mom");
    TTreeReaderArray<double> lep_theta(tree_reader, "ds_lep_theta");

    while (tree_reader.Next()) {

      if(!sel.InFiducial(*vtx_x, *vtx_y, *vtx_z)) continue;

      for(size_t ns = 0; ns < nsims; ns++){
        if(!sel.IsSelected(nu_pdg[ns], cc[ns], lep_contained[ns], particles_contained[ns], 0, 0, 0, 0)) continue;
        if(config->plot_variables[i] == "lep_mom") syst_hists[ns]->Fill(lep_mom[ns]);
        if(config->plot_variables[i] == "lep_theta") syst_hists[ns]->Fill(lep_theta[ns]);
        if(config->plot_variables[i] == "cos_lep_theta") syst_hists[ns]->Fill(cos(lep_theta[ns]));
        total[ns]++;
      }

    }

    // Include scale factor
    for(size_t ns = 0; ns < nsims; ns++){
      syst_hists[ns]->Scale(config->pot_scale_fac[tune_i]);
    }

    // If plotting cross section convert from rate
    if(config->plot_xsec){
      double width = 1;
      if(j != -1) width = width * (bin_edges[j][bin_j+1] - bin_edges[j][bin_j]);
      if(k != -1) width = width * (bin_edges[k][bin_k+1] - bin_edges[k][bin_k]);
      double xsec_scale = 1e38/(width * config->flux[tune_i] * config->targets);
      for(size_t ns = 0; ns < nsims; ns++){
        syst_hists[ns]->Scale(xsec_scale, "width");
      }
    }
    // Else if max error used divide each bin by width
    else if (config->max_error > 0 || config->bin_edges[i].size()>1){
      for(size_t ns = 0; ns < nsims; ns++){
        syst_hists[ns]->Scale(1, "width");
      }
    }

    GetTotalRateSyst(total);

    // Average over the number of universes and get the mean and standard deviation for each bin
    std::vector<double> means = GetTotalSyst(syst_hists, syst_hist);

    if(config->plot_correlation){
      PlotCorrelation(syst_hists, means, "det_"+name, tune_i);
    }

    
    double error;
    double integral = syst_hist->IntegralAndError(0, syst_hist->GetNbinsX()+1, error);
    std::cout<<"Detector = "<<integral<<" +/- "<<error<<"\n";

    return syst_hist;

  }


  // Get the total (unstacked) histogram
  TH1D* GetSystHist(const std::vector<std::vector<double>> &data, const std::vector<bool> &used, TString name, std::vector<std::vector<double>> bin_edges, int tune_i, int i, int j = -1, int bin_j = -1, int k = -1, int bin_k = -1){
      
    int nsims = 100;

    double edges_array[bin_edges[i].size()];
    std::copy(bin_edges[i].begin(), bin_edges[i].end(), edges_array);

    TH1D *syst_hist = new TH1D("syst"+name+"_"+config->tune_name[tune_i], "", bin_edges[i].size()-1, edges_array);
    std::vector<TH1D*> syst_hists;
    for(size_t ns = 0; ns < nsims; ns++){
      TH1D *syst_hist = new TH1D(Form("syst"+name+"_"+config->tune_name[tune_i]+"%i",ns), "", bin_edges[i].size()-1, edges_array);
      syst_hists.push_back(syst_hist);
    }

    std::vector<double> total(nsims, 0);

    TFile data_file(config->input_file[tune_i]);

    //Read in TTree
    TTreeReader tree_reader(config->weight_path, &data_file);
    TTreeReaderArray<double> genie_weight(tree_reader, "genie_weights");
    TTreeReaderArray<double> flux_weight(tree_reader, "flux_weights");
    
    int index = 0;
    int data_i = 0;
    while (tree_reader.Next()) {
      if(!used[index]){ index++; continue; }
      if(j == -1 && k == -1){
        for(size_t ns = 0; ns < nsims; ns++){
          //double weight = genie_weight[ns]*flux_weight[ns];
          double weight = genie_weight[ns];
          if(weight < 0 || weight > 100) continue;
          total[ns] += weight;
          syst_hists[ns]->Fill(data[data_i][i], weight);
        }
      }
      /*else if(data[data_i][j] >= bin_edges[j][bin_j] && data[data_i][j] < bin_edges[j][bin_j+1]){
        if(k == -1){
          for(size_t ns = 0; ns < nsims; ns++){
            if(weight[ns] < 0 || weight[ns] > 100) continue;
            syst_hists[ns]->Fill(data[data_i][i], weight[ns]);
          }
        }
        else if(data[data_i][k] >= bin_edges[k][bin_k] && data[data_i][k] < bin_edges[k][bin_k+1]){
          for(size_t ns = 0; ns < nsims; ns++){
            if(weight[ns] < 0 || weight[ns] > 100) continue;
            syst_hists[ns]->Fill(data[data_i][i], weight[ns]);
          }
        }
      }*/
      index++;
      data_i++;
    }

    // Include scale factor
    for(size_t ns = 0; ns < nsims; ns++){
      syst_hists[ns]->Scale(config->pot_scale_fac[tune_i]);
    }

    // If plotting cross section convert from rate
    if(config->plot_xsec){
      double width = 1;
      if(j != -1) width = width * (bin_edges[j][bin_j+1] - bin_edges[j][bin_j]);
      if(k != -1) width = width * (bin_edges[k][bin_k+1] - bin_edges[k][bin_k]);
      double xsec_scale = 1e38/(width * config->flux[tune_i] * config->targets);
      for(size_t ns = 0; ns < nsims; ns++){
        syst_hists[ns]->Scale(xsec_scale, "width");
      }
    }
    // Else if max error used divide each bin by width
    else if (config->max_error > 0 || config->bin_edges[i].size()>1){
      for(size_t ns = 0; ns < nsims; ns++){
        syst_hists[ns]->Scale(1, "width");
      }
    }

    // Average over the number of universes and get the mean and standard deviation for each bin
    std::vector<double> means = GetTotalSyst(syst_hists, syst_hist);

    if(config->plot_correlation){
      PlotCorrelation(syst_hists, means, "rw_"+name, tune_i);
    }

    GetTotalRateSyst(total);

    double error;
    double integral = syst_hist->IntegralAndError(0, syst_hist->GetNbinsX()+1, error);
    std::cout<<"Genie + flux = "<<integral<<" +/- "<<error<<" "<<100*error/integral<<"\n";

    return syst_hist;
  }

  // Add histogram errors in quadrature, ignoring bin contents
  void AddErrors(TH1D* syst_hist, TH1D* hist){
    for(size_t i = 1; i <= syst_hist->GetNbinsX(); i++){
      double new_err = std::sqrt(std::pow(syst_hist->GetBinError(i),2)+std::pow(hist->GetBinError(i),2));
      syst_hist->SetBinError(i, new_err);
    }
  }

  // Plot correlation and covariance matrices
  void PlotCorrelation(std::vector<TH1D*> syst_hists, std::vector<double> means, TString name, int tune_i){
    if(syst_hists.size()==0) return;

    size_t nbins = syst_hists[0]->GetNbinsX();
    TH2D* covariance = new TH2D("covariance_"+name+"_"+config->tune_name[tune_i], "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    TH2D* fraccovar = new TH2D("frac_covariance_"+name+"_"+config->tune_name[tune_i], "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    TH2D* correlation = new TH2D("correlation_"+name+"_"+config->tune_name[tune_i], "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    for(size_t i = 1; i <= syst_hists[0]->GetNbinsX(); i++){
      double cv_i = means[i-1];
      for(size_t j = 1; j <= syst_hists[0]->GetNbinsX(); j++){
        double cv_j = means[j-1];
        double E_ij = 0;
        for(size_t ns = 0; ns < syst_hists.size(); ns++){
          E_ij += (syst_hists[ns]->GetBinContent(i)-cv_i)*(syst_hists[ns]->GetBinContent(j)-cv_j);
        }
        E_ij /= syst_hists.size();
        covariance->SetBinContent(i, j, E_ij);
        fraccovar->SetBinContent(i, j, E_ij/(cv_i*cv_j));
      }
    }

    for(size_t i = 1; i <= syst_hists[0]->GetNbinsX(); i++){
      for(size_t j = 1; j <= syst_hists[0]->GetNbinsX(); j++){
        double corr_bin = covariance->GetBinContent(i, j)/(std::sqrt(covariance->GetBinContent(i,i))*std::sqrt(covariance->GetBinContent(j,j)));
        correlation->SetBinContent(i, j, corr_bin);
      }
    }
    plotman->Plot2D(covariance, covariance->GetName(), "Bin i", "Bin j");
    plotman->Plot2D(fraccovar, fraccovar->GetName(), "Bin i", "Bin j");
    plotman->Plot2D(correlation, correlation->GetName(), "Bin i", "Bin j");
  }

  // Get total systematic errors from variations
  std::vector<double>  GetTotalSyst(std::vector<TH1D*> syst_hists, TH1D* syst_hist){
    std::vector<double> means;
    if(syst_hists.size() == 0) return means;
    //std::cout<<"number of hists = "<<syst_hists.size()<<"\n";

    for(size_t n = 1; n <= syst_hists[0]->GetNbinsX(); n++){
      double mean = 0;
      for(size_t ns = 0; ns < syst_hists.size(); ns++){
        mean += syst_hists[ns]->GetBinContent(n);
      }
      mean /= syst_hists.size();
      means.push_back(mean);
      double std_dev = 0;
      for(size_t ns = 0; ns < syst_hists.size(); ns++){
        std_dev += std::pow(syst_hists[ns]->GetBinContent(n) - mean, 2.);
      }
      std_dev = std::sqrt(std_dev/(syst_hists.size()-1));
      //std::cout<<"Bin "<<n<<" mean = "<<mean<<" stddev = "<<std_dev<<"\n";
      syst_hist->SetBinContent(n, mean);
      syst_hist->SetBinError(n, std_dev);
    }
    return means;
  }

  // Get total systematic errors from variations
  void  GetTotalRateSyst(std::vector<double> total){
    double mean = 0;
    for(size_t ns = 0; ns < total.size(); ns++){
      mean += total[ns];
    }
    mean /= total.size();
    double std_dev = 0;
    for(size_t ns = 0; ns < total.size(); ns++){
      std_dev += std::pow(total[ns] - mean, 2.);
    }
    std_dev = std::sqrt(std_dev/(total.size()-1));
    std::cout<<"Total rate = "<<mean<<" +/- "<<std_dev<<" "<<100*std_dev/mean<<"%\n";
  }

};

#endif
