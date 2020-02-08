#ifndef BINMANAGER_H
#define BINMANAGER_H

#include "Configuration.h"
#include "DataManager.h"

// Structure for holding binning information
class BinManager
{
  public:

  Configuration *config;
  DataManager *dataman;

  // Multiple files possible but binning must be same, use first
  // 1D binning - X bins
  std::map<TString, std::vector<double>> bin_edges_1D;
  // 2D binning - X and Y bins
  std::map<std::pair<TString, TString>, std::vector<std::vector<double>>> bin_edges_2D;

  //BinManager(){}

  // Constructor
  BinManager(Configuration *c, DataManager *d)
  {
    config = c;
    dataman = d;

    // Get the 1D binning for all plotting variables
    std::vector<std::vector<double>> bin_edges = GetBinning(dataman->total_data);
    for(size_t i = 0; i < config->plot_variables.size(); i++){
      bin_edges_1D[config->plot_variables[i]] = bin_edges[i];
      for(size_t j = 0; j < config->plot_variables.size(); j++){
        std::pair<TString, TString> key2D = std::make_pair(config->plot_variables[i], config->plot_variables[j]);
        std::vector<std::vector<double>> bins_2D = GetBinning2D(dataman->total_data, bin_edges, i, j);
        bin_edges_2D[key2D] = bins_2D;
      }
    }

  }

  // Accessor to 1D binning information
  std::vector<double> Get1DBinning(TString plot_var){
    std::vector<double> null;
    if(bin_edges_1D.find(plot_var) == bin_edges_1D.end()){
      std::cout<<"couldn't find binning for variables!\n";
      return null;
    }
    return bin_edges_1D[plot_var];
  }

  // Accessor to 2D binning information
  std::vector<std::vector<double>> Get2DBinning(TString plot_var1, TString plot_var2){
    std::vector<std::vector<double>> null;
    std::pair<TString, TString> key = std::make_pair(plot_var1, plot_var2);
    if(bin_edges_2D.find(key) == bin_edges_2D.end()){
      std::cout<<"couldn't find binning for variables!\n";
      return null;
    }
    return bin_edges_2D[key];
  }

  // Set the minimum bin value
  double GetMinBin(const std::vector<double>& data, int i){
    
    double min_bins = data[0];

    if(config->plot_variables[i]=="particles_contained"||
       config->plot_variables[i]=="lep_contained"||
       config->plot_variables[i]=="cc"){
      min_bins = 0;
    }
    if(config->plot_variables[i]=="nu_pdg") min_bins = -14;
    if(config->plot_variables[i]=="int_type") min_bins = 0;
    if(config->plot_variables[i]=="n_pr") min_bins = 0;
    if(config->plot_variables[i]=="n_pipm") min_bins = 0;
    if(config->plot_variables[i]=="n_pi0") min_bins = 0;

    return (min_bins);
   
  }

  // Set the maximum bin value
  double GetMaxBin(const std::vector<double> &data, int i){
    
    double max_bins = data[data.size()-1];

    if(config->plot_variables[i]=="particles_contained"||
       config->plot_variables[i]=="lep_contained"||
       config->plot_variables[i]=="cc"){
      max_bins = 2;
    }
    if(config->plot_variables[i]=="nu_pdg") max_bins = 15;
    if(config->plot_variables[i]=="int_type") max_bins = 11;
    if(config->plot_variables[i]=="n_pr") max_bins = 11;
    if(config->plot_variables[i]=="n_pipm") max_bins = 11;
    if(config->plot_variables[i]=="n_pi0") max_bins = 11;

    return (max_bins);
   
  }

  // Set the default number of bins
  int DefaultBins(int data_size, int dims, int i){
    
    int dflt_bins = (int)pow(2*pow(data_size,.33), 1./dims);

    if(config->plot_variables[i]=="particles_contained"||
       config->plot_variables[i]=="lep_contained"||
       config->plot_variables[i]=="cc"){
      dflt_bins = 2;
    }
    if(config->plot_variables[i]=="nu_pdg") dflt_bins = 29;
    if(config->plot_variables[i]=="int_type") dflt_bins = 11;
    if(config->plot_variables[i]=="n_pr") dflt_bins = 11;
    if(config->plot_variables[i]=="n_pipm") dflt_bins = 11;
    if(config->plot_variables[i]=="n_pi0") dflt_bins = 11;

    return (dflt_bins);
   
  }

  // Recursive function to merge bins to get up to required statistical error
  std::vector<double> ChangeBinning(TH1D* hist, double max){

    std::vector<double> bin_edges;
    for(size_t i = 0; i < hist->GetNbinsX(); i++){
      bin_edges.push_back(hist->GetBinLowEdge(i+1));
    }
    bin_edges.push_back(max);

    for(size_t i = 0; i < hist->GetNbinsX(); i++){
      if(hist->GetBinError(i+1)/hist->GetBinContent(i+1) > config->max_error || hist->GetBinContent(i+1) == 0){

        if(i != hist->GetNbinsX()-1){
          bin_edges.erase(bin_edges.begin()+i+1);
          double edges_array[bin_edges.size()];
          std::copy(bin_edges.begin(), bin_edges.end(), edges_array);
          TH1D* new_hist = (TH1D*)hist->Rebin(bin_edges.size()-1, "new", edges_array);
          return ChangeBinning(new_hist, max);
        }

        else{
          bin_edges.erase(bin_edges.begin()+(bin_edges.size()-2));
          return bin_edges;
        }

      }
    }
    return bin_edges;
  }

  // Function to get or calculate the number of bins, and min and max bins
  std::vector<std::vector<double>> GetBinning(const std::vector<std::vector<double>> &data_v){

    std::map<int, std::vector<double>> data_map;
    for(auto const& data : data_v){
      for(size_t i = 0; i < data.size(); i++){
        data_map[i].push_back(data[i]);
      }
    }
    for(size_t i = 0; i < data_map.size(); i++){
      std::sort(data_map[i].begin(), data_map[i].end());
    }
    std::vector<double> hist_min = config->min_value;
    std::vector<double> hist_max = config->max_value;
    std::vector<int> hist_bins = config->num_bins;
    for(size_t i = 0; i < config->plot_variables.size(); i++){
      if(config->min_value[i] < 0){
        hist_min[i] = GetMinBin(data_map[i], i);
      }
      if(config->max_value[i] <= 0){
        hist_max[i] = GetMaxBin(data_map[i], i);
      }
      if(config->num_bins[i] <= 0){
        int data_size = 0;
        for(auto const& data : data_v){
          if(data[i] >= hist_min[i] && data[i] <= hist_max[i]) data_size++;
        }
        hist_bins[i] = DefaultBins(data_size*config->pot_scale_fac[0], data_map.size(), i);
      }
    }
    
    std::vector<std::vector<double>> all_bin_edges;
    for(size_t i = 0; i < config->plot_variables.size(); i++){
      // If bin edges are set by the user
      if(config->bin_edges[i].size() > 1){
        all_bin_edges.push_back(config->bin_edges[i]);
      }
      else{
        TH1D *temp_hist = new TH1D("temp_hist", "", hist_bins[i], hist_min[i], hist_max[i]);
        std::vector<double> bin_edges;
        for(size_t j = 0; j < temp_hist->GetNbinsX(); j++){
          bin_edges.push_back(temp_hist->GetBinLowEdge(j+1));
        }
        delete temp_hist;
        bin_edges.push_back(hist_max[i]);
        all_bin_edges.push_back(bin_edges);
      }
    }

    // If a maximum bin error is set
    if(config->max_error <= 0) return all_bin_edges;


    for(size_t i = 0; i < config->plot_variables.size(); i++){
      // Get the optimal binning for the first variable
      TH1D *temp_hist = new TH1D("temp_hist", "", hist_bins[i], hist_min[i], hist_max[i]);
      // Loop over data
      for(auto const& data : data_v){
        // Fill temporary histogram
        temp_hist->Fill(data[i]);
      }
      // Include scale factor bin by bin as Scale() won't change errors
      for(size_t n = 1; n <= temp_hist->GetNbinsX(); n++){
        temp_hist->SetBinContent(n, temp_hist->GetBinContent(n)*config->pot_scale_fac[0]);
      }
      // Change the binning so that all bin errors below maximum
      std::vector<double> bin_edges = ChangeBinning(temp_hist, hist_max[i]);
      delete temp_hist;
      all_bin_edges[i] = bin_edges;
    }

    return all_bin_edges;
  }

  std::vector<std::vector<double>> GetBinning2D(const std::vector<std::vector<double>> &data_v, std::vector<std::vector<double>> bin_edges, size_t i, size_t j){

    // Perform two dimensional rebinning
    double xedges_array[bin_edges[i].size()];
    std::copy(bin_edges[i].begin(), bin_edges[i].end(), xedges_array);
    double yedges_array[bin_edges[j].size()];
    std::copy(bin_edges[j].begin(), bin_edges[j].end(), yedges_array);
    // Create 2D histogram with the calculated binning
    TH2D *temp_hist = new TH2D("temp_hist", "", bin_edges[i].size()-1, xedges_array, bin_edges[j].size()-1, yedges_array);
    // Fill with the data
    for(auto const& data : data_v){
      temp_hist->Fill(data[i], data[j]);
    }
    // Include scale factor bin by bin as Scale() won't change errors
    for(size_t i = 1; i <= temp_hist->GetNbinsX(); i++){
      for(size_t j = 1; j <= temp_hist->GetNbinsX(); j++){
        temp_hist->SetBinContent(i, j, temp_hist->GetBinContent(i, j)*config->pot_scale_fac[0]);
      }
    }
    std::vector<std::vector<double>> all_bin_edges = ChangeBinning2D(temp_hist, bin_edges[i]);
    delete temp_hist;
    
    return all_bin_edges;
  }

/*
  // Rebin to the maximum bin error for 2D histograms
  std::vector<double> ChangeBinning2D(const std::vector<std::vector<double>> &data, std::vector<std::vector<double>> bin_edges, int i, int j, int bin_j){

    double edges_array[bin_edges[i].size()];
    std::copy(bin_edges[i].begin(), bin_edges[i].end(), edges_array);
    TH1D *temp_hist = new TH1D("temp_hist", "", bin_edges[i].size()-1, edges_array);
    // Loop over data
    for(auto const& dat : data){
      // Fill temporary histogram
      if(dat[j] >= bin_edges[j][bin_j] && dat[j] < bin_edges[j][bin_j+1]){
        temp_hist->Fill(dat[i]);
      }
    }
    // Include scale factor bin by bin as Scale() won't change errors
    for(size_t n = 1; n <= temp_hist->GetNbinsX(); n++){
      temp_hist->SetBinContent(n, temp_hist->GetBinContent(n)*config->pot_scale_fac[0]);
    }
    // Change the binning so that all bin errors below maximum
    std::vector<double> bin_edges_new = ChangeBinning(temp_hist, bin_edges[i][bin_edges[i].size()-1]);
    delete temp_hist;
    return bin_edges_new;

  }
*/

  // Rebin to the maximum bin error for 2D histograms
  std::vector<std::vector<double>> ChangeBinning2D(TH2D* hist, std::vector<double> xbin_edges){
    
    std::vector<std::vector<double>> all_edges;
    // Loop over slices in one dimension
    for(size_t i = 1; i <= hist->GetNbinsY(); i++){
      TH1D* temp_hist = (TH1D*) hist->ProjectionX("temp_slice", i, i);
      // Rebin slice so errors below maximum
      std::vector<double> bin_edges_new = ChangeBinning(temp_hist, xbin_edges.back());
      // Push new bin edges to vector
      all_edges.push_back(bin_edges_new);
      delete temp_hist;
    }
    return all_edges;
  }

  /*
  // Rebin to the maximum bin error for 2D histograms
  std::pair<std::vector<double>, std::vector<double>> ChangeBinning2D(TH2D* hist, double max1, double max2){

    std::vector<double> xbin_edges;
    std::vector<double> ybin_edges;
    for(size_t i = 1; i <= hist->GetNbinsX(); i++){
      xbin_edges.push_back(hist->ProjectionX()->GetBinLowEdge(i));
    }
    xbin_edges.push_back(max1);
    for(size_t j = 1; j <= hist->GetNbinsY(); j++){
      ybin_edges.push_back(hist->ProjectionY()->GetBinLowEdge(j));
    }
    ybin_edges.push_back(max2);

    // Loop over all X and Y bins
    std::cout<<"n xbins = "<<hist->GetNbinsX()<<" n ybins = "<<hist->GetNbinsY()<<"\n";
    for(size_t i = 1; i <= hist->GetNbinsX(); i++){
      for(size_t j = 1; j <= hist->GetNbinsY(); j++){
        std::cout<<"i "<<i<<" j "<<j<<" err "<<hist->GetBinError(i,j)/hist->GetBinContent(i,j)<<" = "<<hist->GetBinError(i,j)<<"/"<<hist->GetBinContent(i,j)<<"\n";
        // See if the bin error is above the limit
        if(hist->GetBinError(i,j)/hist->GetBinContent(i,j) > config->max_error || hist->GetBinContent(i,j) == 0){
          std::cout<<"Above error \n";

          // if not in the last x and y bins
          if(i != hist->GetNbinsX() && j != hist->GetNbinsY()){
            std::cout<<"Both\n";
            // Take a look in the next bin over in x and y
            double xcontent = hist->GetBinContent(i+1, j);
            double ycontent = hist->GetBinContent(i, j+1);
            // Remove the bin edge in whichever direction has more entries
            if(xcontent > ycontent){
              xbin_edges.erase(xbin_edges.begin()+i);
            }
            else{
              ybin_edges.erase(ybin_edges.begin()+i);
            }
            Rebin2D(hist, xbin_edges, ybin_edges);
            return ChangeBinning2D(hist, max1, max2);
          }
          // If not in the last x bin
          else if(i != hist->GetNbinsX()){
            std::cout<<"i\n";
            xbin_edges.erase(xbin_edges.begin()+i);
            Rebin2D(hist, xbin_edges, ybin_edges);
            return ChangeBinning2D(hist, max1, max2);
          }
          // If not in the last y bin
          else if(j != hist->GetNbinsY()){
            std::cout<<"j\n";
            ybin_edges.erase(ybin_edges.begin()+i);
            Rebin2D(hist, xbin_edges, ybin_edges);
            return ChangeBinning2D(hist, max1, max2);
          }
          else{
            std::cout<<"neither\n";
            // Look at the previous bin in x and y
            double xcontent = hist->GetBinContent(i-1, j);
            double ycontent = hist->GetBinContent(i, j-1);
            // Remove the bin edge in whichever direction has least entries
            if(xcontent < ycontent){
              xbin_edges.erase(xbin_edges.begin()+(xbin_edges.size()-2));
            }
            else{
              ybin_edges.erase(ybin_edges.begin()+(ybin_edges.size()-2));
            }
            return std::make_pair(xbin_edges, ybin_edges);
          }
        }
      }
    }
    return std::make_pair(xbin_edges, ybin_edges);

  }

  // Rebin in 2D assuming new edges have corresponding edges in the old histogram
  void Rebin2D(TH2D* hist, std::vector<double> xbin_edges, std::vector<double> ybin_edges){
    double xedges_array[xbin_edges.size()];
    std::copy(xbin_edges.begin(), xbin_edges.end(), xedges_array);
    double yedges_array[ybin_edges.size()];
    std::copy(ybin_edges.begin(), ybin_edges.end(), yedges_array);
    TH2D* temp_hist = new TH2D("temp", "", xbin_edges.size()-1, xedges_array, ybin_edges.size()-1, yedges_array);
    for(size_t i = 1; i < hist->GetNbinsX(); i++){
      double xcenter = hist->GetXaxis()->GetBinCenter(i);
      for(size_t j = 1; j < hist->GetNbinsY(); j++){
        double ycenter = hist->GetYaxis()->GetBinCenter(j);
        double weight = hist->GetBinContent(i, j);
        temp_hist->Fill(xcenter, ycenter, weight);
      }
    }
    hist->Reset();
    hist->SetBins(xbin_edges.size()-1, xedges_array, ybin_edges.size()-1, yedges_array);
    for(size_t i = 1; i < hist->GetNbinsX(); i++){
      for(size_t j = 1; j < hist->GetNbinsY(); j++){
        hist->SetBinContent(i, j, temp_hist->GetBinContent(i, j));
      }
    }
    delete temp_hist;
    
  }
*/
};

#endif
