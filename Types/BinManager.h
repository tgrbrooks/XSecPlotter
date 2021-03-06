#ifndef BINMANAGER_H
#define BINMANAGER_H

#include "Configuration.h"
#include "DataManager.h"

// Structure for holding binning information
class BinManager
{
  public:

  // --- Service pointers ---
  Configuration *config; // Global configuration
  DataManager *dataman; // Stored data

  // Multiple files possible but binning must be same, use first
  // 1D binning - X bins
  std::map<TString, std::vector<double>> bin_edges_1D;
  // 2D binning - X and Y bins
  std::map<std::pair<TString, TString>, std::vector<std::vector<double>>> bin_edges_2D;

  // Constructor
  BinManager(Configuration *c, DataManager *d)
  {
    config = c;
    dataman = d;

    // Get the 1D binning for all plotting variables
    // Use a smaller error when getting the 1D binning on 2D slices
    if(config->plot_variables.size() == 2 && config->max_error > 0) config->max_error = config->max_error/3.;
    std::vector<std::vector<double>> bin_edges = GetBinning(dataman->total_data);
    // Switch the error back
    if(config->plot_variables.size() == 2 && config->max_error > 0) config->max_error = config->max_error*3.;

    // Get the 2D binning for all plotting variables
    for(size_t i = 0; i < config->plot_variables.size(); i++){
      bin_edges_1D[config->plot_variables[i]] = bin_edges[i];
      for(size_t j = 0; j < config->plot_variables.size(); j++){
        std::pair<TString, TString> key2D = std::make_pair(config->plot_variables[i], config->plot_variables[j]);
        std::vector<std::vector<double>> bins_2D = GetBinning2D(dataman->total_data, bin_edges, i, j);
        bin_edges_2D[key2D] = bins_2D;
      }
    }

  }

  // ------------------------------------------------------------------------------------------------------------------
  //                                                  ACCESSORS
  // ------------------------------------------------------------------------------------------------------------------

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

  // ------------------------------------------------------------------------------------------------------------------
  //                                          DEFAULT BINNING SETTINGS
  // ------------------------------------------------------------------------------------------------------------------
  
  // Set the minimum bin value
  double GetMinBin(const std::vector<double>& data, int i){
    
    double min_bins = data[0];

    // Some variables have special values
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

    // Some variables have special values
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

    // Some variables have special values
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

  // ------------------------------------------------------------------------------------------------------------------
  //                                          REBINNING FUNCTIONS
  // ------------------------------------------------------------------------------------------------------------------

  // Recursive function to merge bins to get up to required statistical error
  std::vector<double> ChangeBinning(TH1D* hist, double max){

    // Get the edges of the histogram
    std::vector<double> bin_edges;
    for(int i = 0; i < hist->GetNbinsX(); i++){
      bin_edges.push_back(hist->GetBinLowEdge(i+1));
    }
    bin_edges.push_back(max);

    // Loop over the bins
    for(int i = 0; i < hist->GetNbinsX(); i++){

      // If bin error is greater than max value
      if(hist->GetBinError(i+1)/hist->GetBinContent(i+1) > config->max_error || hist->GetBinContent(i+1) == 0){

        // If we're not in the last bin erase the next bin edges
        if(i != hist->GetNbinsX()-1){
          bin_edges.erase(bin_edges.begin()+i+1);
          double edges_array[bin_edges.size()];
          std::copy(bin_edges.begin(), bin_edges.end(), edges_array);
          TH1D* new_hist = (TH1D*)hist->Rebin(bin_edges.size()-1, "new", edges_array);
          return ChangeBinning(new_hist, max);
        }

        // Otherwise erase the bin edge before
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

    // Sort the data into separate vectors
    std::map<int, std::vector<double>> data_map;
    for(auto const& data : data_v){
      for(size_t i = 0; i < data.size(); i++){
        data_map[i].push_back(data[i]);
      }
    }
    // Sort by size
    for(size_t i = 0; i < data_map.size(); i++){
      std::sort(data_map[i].begin(), data_map[i].end());
    }

    // Get default values from configuration
    std::vector<double> hist_min = config->min_value;
    std::vector<double> hist_max = config->max_value;
    std::vector<int> hist_bins = config->num_bins;

    // Loop over the plotting variables
    for(size_t i = 0; i < config->plot_variables.size(); i++){
      // Set minimum and maximum values if not specified
      if(config->min_value[i] < 0){
        hist_min[i] = GetMinBin(data_map[i], i);
      }
      if(config->max_value[i] <= 0){
        hist_max[i] = GetMaxBin(data_map[i], i);
      }
      // Set the default number of bins if not specified
      if(config->num_bins[i] <= 0){
        int data_size = 0;
        for(auto const& data : data_v){
          if(data[i] >= hist_min[i] && data[i] <= hist_max[i]) data_size++;
        }
        hist_bins[i] = DefaultBins(data_size*config->pot_scale_fac[0], data_map.size(), i);
      }
    }
    
    // Store the bin edges
    std::vector<std::vector<double>> all_bin_edges;
    for(size_t i = 0; i < config->plot_variables.size(); i++){
      // If bin edges are set by the user
      if(config->bin_edges[i].size() > 1){
        all_bin_edges.push_back(config->bin_edges[i]);
      }
      // Otherwise use the calculated binning
      else{
        TH1D *temp_hist = new TH1D("temp_hist", "", hist_bins[i], hist_min[i], hist_max[i]);
        std::vector<double> bin_edges;
        for(int j = 0; j < temp_hist->GetNbinsX(); j++){
          bin_edges.push_back(temp_hist->GetBinLowEdge(j+1));
        }
        delete temp_hist;
        bin_edges.push_back(hist_max[i]);
        all_bin_edges.push_back(bin_edges);
      }
    }

    // If a maximum bin error is not set return the bin edges
    if(config->max_error <= 0) return all_bin_edges;

    // Loop over the plotting variables
    for(size_t i = 0; i < config->plot_variables.size(); i++){
      // Get the optimal binning for the first variable
      TH1D *temp_hist = new TH1D("temp_hist", "", hist_bins[i], hist_min[i], hist_max[i]);
      // Loop over data
      for(auto const& data : data_v){
        // Fill temporary histogram
        temp_hist->Fill(data[i]);
      }
      // Include scale factor bin by bin as Scale() won't change errors
      for(int n = 1; n <= temp_hist->GetNbinsX(); n++){
        temp_hist->SetBinContent(n, temp_hist->GetBinContent(n)*config->pot_scale_fac[0]);
      }
      // Change the binning so that all bin errors below maximum
      std::vector<double> bin_edges = ChangeBinning(temp_hist, hist_max[i]);
      delete temp_hist;
      all_bin_edges[i] = bin_edges;
    }

    return all_bin_edges;
  }

  // Get the two dimensional binning from slicing in one variable
  std::vector<std::vector<double>> GetBinning2D(const std::vector<std::vector<double>> &data_v, std::vector<std::vector<double>> bin_edges, size_t i, size_t j){

    // Get the bin edges as arrays
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
    for(int x = 1; x <= temp_hist->GetNbinsX(); x++){
      for(int y = 1; y <= temp_hist->GetNbinsY(); y++){
        temp_hist->SetBinContent(x, y, temp_hist->GetBinContent(x, y)*config->pot_scale_fac[0]);
      }
    }

    // If maximum error is not set return the bin edges
    if(config->max_error <= 0){
      std::vector<std::vector<double>> all_edges;
      for(int y = 1; y <= temp_hist->GetNbinsY(); y++){
        all_edges.push_back(bin_edges[i]);
      }
      delete temp_hist;
      return all_edges;
    }

    // Otherwise rebin
    std::vector<std::vector<double>> all_bin_edges = ChangeBinning2D(temp_hist, bin_edges[i]);
    delete temp_hist;
    
    return all_bin_edges;
  }

  // Rebin to the maximum bin error for 2D histograms
  std::vector<std::vector<double>> ChangeBinning2D(TH2D* hist, std::vector<double> xbin_edges){
    
    std::vector<std::vector<double>> all_edges;
    // Loop over slices in one dimension
    for(int i = 1; i <= hist->GetNbinsY(); i++){
      TH1D* temp_hist = (TH1D*) hist->ProjectionX("temp_slice", i, i);
      // Rebin slice so errors below maximum
      std::vector<double> bin_edges_new = ChangeBinning(temp_hist, xbin_edges.back());
      // Push new bin edges to vector
      all_edges.push_back(bin_edges_new);
      delete temp_hist;
    }
    return all_edges;
  }

  
};

#endif
