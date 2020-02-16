#ifndef CONFIGURATION_H
#define CONFIGURATION_H

// Structure for holding global configuration settings
class Configuration
{
  public:

  // File configurations
  std::vector<TString> input_file;
  std::vector<TString> tune_name;
  TString output_file;
  // Neutrino interaction configurations
  std::vector<int> nu_pdg;
  std::vector<int> is_cc;
  std::vector<bool> contained_lepton;
  std::vector<bool> contained_particles;
  std::vector<double> fiducial;
  bool plot_by_fsi;
  std::vector<int> num_protons;
  std::vector<int> num_pipm;
  std::vector<int> num_pi0;
  std::vector<int> interaction_type;
  // Plotting variable configurations
  TString stage;
  std::vector<TString> plot_variables;
  bool plot_slices;
  double pot_scale;
  // Plotting option configurations
  bool plot_stacked;
  TString stack_by;
  std::vector<double> min_value;
  std::vector<double> max_value;
  std::vector<int> num_bins;
  std::vector<std::vector<double>> bin_edges;
  double max_error;
  bool plot_xsec;
  bool plot_filled;
  // Optional extras
  bool show_info;
  bool show_error_band;
  bool show_stat_error;
  bool show_syst_error;
  std::vector<TString> systematics;
  double constant_syst;
  bool show_error_bars;
  bool plot_correlation;
  bool plot_universes;
  bool plot_eff_pur;
  bool plot_response; //TODO
  bool unfold; //TODO

  // Set by functions
  std::vector<double> pot;
  std::vector<double> pot_scale_fac;
  double fiducial_mass;
  std::vector<double> flux;
  double targets;
  // Constants
  // Integrated flux for each neutrino species
  std::map<int, double> nu_flux = {{14, 1.305e13}, {-14, 1.011e12}, {12, 7.924e10}, {-12, 8.4133e9}}; //[/6.6e20POT/cm^2]
  std::vector<int> cols = {46, 33, 38, 42, 40, 30, 49};

  // Constructor
  //Configuration(){}
  
  // Constructor
  Configuration(const std::string config_filename){
    Configure(config_filename);
  }

  // Retrieve the configuration from a file
  void Configure(const std::string config_filename){
    
    // Read in the configuration file
    std::ifstream config_file;
    config_file.open(config_filename);
    if(!config_file){
      std::cerr << "Unable to open configuration file!\n";
      exit(1);
    }
 
    std::string delim = ":";
    std::string line;
    while(std::getline(config_file, line)){
      std::string key = line.substr(0, line.find(delim));
      std::string value = line.erase(0, line.find(delim)+delim.length());
      value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
      // File configurations
      if(key.find("InputFile") != std::string::npos)    input_file = ToTStrings(value);
      if(key.find("TuneName") != std::string::npos)     tune_name = ToTStrings(value);
      if(key.find("OutputFile") != std::string::npos)   output_file = TString(value);
      // Neutrino configurations
      if(key.find("NuPdg") != std::string::npos)              nu_pdg = ToInts(value);
      if(key.find("IsCC") != std::string::npos)               is_cc = ToInts(value);
      if(key.find("ContainedLepton") != std::string::npos)    contained_lepton = ToBools(value);
      if(key.find("ContainedParticles") != std::string::npos) contained_particles = ToBools(value);
      if(key.find("Fiducial") != std::string::npos)           fiducial = ToDoubles(value);
      if(key.find("PlotByFsi") != std::string::npos)          plot_by_fsi = (value=="true");
      if(key.find("NumProtons") != std::string::npos)         num_protons = ToInts(value);
      if(key.find("NumPiPM") != std::string::npos)            num_pipm = ToInts(value);
      if(key.find("NumPi0") != std::string::npos)             num_pi0 = ToInts(value);
      if(key.find("InteractionType") != std::string::npos)    interaction_type = ToInts(value);
      // Plotting variable configurations
      if(key.find("Stage") != std::string::npos)        stage = TString(value);
      if(key.find("PlotVariable") != std::string::npos) plot_variables = ToTStrings(value);
      if(key.find("PlotSlices") != std::string::npos)   plot_slices = (value=="true");
      if(key.find("PotScale") != std::string::npos)     pot_scale = stod(value);
      // Plotting option configurations
      if(key.find("PlotStacked") != std::string::npos) plot_stacked = (value=="true");
      if(key.find("StackBy") != std::string::npos)     stack_by = TString(value);
      if(key.find("MinValue") != std::string::npos)    min_value = ToDoubles(value);
      if(key.find("MaxValue") != std::string::npos)    max_value = ToDoubles(value);
      if(key.find("NumBins") != std::string::npos)     num_bins = ToInts(value);
      if(key.find("BinEdges") != std::string::npos){
        for(auto const& val : ToVector(value, "],[")) bin_edges.push_back(ToDoubles(val));
      }
      if(key.find("MaxError") != std::string::npos)    max_error = stod(value);
      if(key.find("PlotXSec") != std::string::npos)    plot_xsec = (value=="true");
      if(key.find("PlotFilled") != std::string::npos)  plot_filled = (value=="true");
      // Optional extras
      if(key.find("ShowInfo") != std::string::npos)        show_info = (value=="true");
      if(key.find("ShowErrorBand") != std::string::npos)   show_error_band = (value=="true");
      if(key.find("ShowStatError") != std::string::npos)   show_stat_error = (value=="true");
      if(key.find("ShowSystError") != std::string::npos)   show_syst_error = (value=="true");
      if(key.find("Systematics") != std::string::npos)     systematics = ToTStrings(value);
      if(key.find("ConstantSyst") != std::string::npos)    constant_syst = stod(value);
      if(key.find("ShowErrorBars") != std::string::npos)   show_error_bars = (value=="true");
      if(key.find("PlotCorrelation") != std::string::npos) plot_correlation = (value=="true");
      if(key.find("PlotUniverses") != std::string::npos)   plot_universes = (value=="true");
      if(key.find("PlotEffPur") != std::string::npos)      plot_eff_pur = (value=="true");
      if(key.find("PlotResponse") != std::string::npos)    plot_response = (value=="true");
      if(key.find("Unfold") != std::string::npos)          unfold = (value=="true");
    }
 
    // Do some basic error checking
    if(plot_variables.size() != min_value.size()
       || plot_variables.size() != max_value.size()
       || plot_variables.size() != num_bins.size()
       || plot_variables.size() != bin_edges.size()){
      std::cout<<"Must have same number of binning parameters as plotting variables!\n";
      exit(1);
    }
    if(plot_variables.size() > 3){
      std::cout<<"Sorry, ROOT doesn't do >3D histograms...\n";
      exit(1);
    }
    if(plot_variables.size() < 1){
      std::cout<<"Need something to plot in.\n";
      exit(1);
    }
    if(!(stack_by=="fsi" || stack_by=="int" || stack_by=="nu")){
      std::cout<<"Unknown stack by parameter, not stacking.\n";
      plot_stacked = false;
    }
    if(input_file.size()>1 && input_file.size() != tune_name.size()){
      std::cout<<"Need to name inputs if using more than one.\n";
      exit(1);
    }
    if(input_file.size()>1 && pot_scale == -1){
      std::cout<<"Two inputs but no POT scale! scaling to POT of first file.\n";
    }
  }

  // Get the POT, flux and target number from the metadata
  void GetMetaData(){

    // Open the root tree file
    for(size_t i = 0; i < input_file.size(); i++){
      TFile data_file(input_file[i]);

      //Read in TTree
      TTreeReader tree_reader("XSecTree/metadata", &data_file);
      TTreeReaderValue<double> pot_val(tree_reader, "pot");
      double pot_count = 0;
    
      while (tree_reader.Next()) {
        pot_count += *pot_val;
      }
    
      pot.push_back(pot_count);

      if(pot_scale > 0){
        pot_scale_fac.push_back(pot_scale/pot[i]);
      }
      else if(input_file.size() == 1) pot_scale_fac.push_back(1);
      else pot_scale_fac.push_back(pot[0]/pot[i]);

      double flux_factor = 0;
      for(auto const& pdg : nu_pdg){
        if(nu_flux.find(pdg) == nu_flux.end()){
          std::cout<<"Unknown neutrino PDG code!\n";
          exit(1);
        }
        flux_factor += nu_flux[pdg];
      }
      flux.push_back(flux_factor * pot[i] * pot_scale_fac[i] / 6.6e20); // [cm^-2]
    }
    std::cout<<"Integrated flux = "<<flux[0]<<"\n";

    double volume = 400*400*500; // [cm^3]
    if(std::find(fiducial.begin(), fiducial.end(), -1) == fiducial.end() && fiducial.size() == 6){
      volume = (400-fiducial[0]-fiducial[3])*(400-fiducial[1]-fiducial[4])*(500-fiducial[2]-fiducial[5]); // [cm^3]
    }
    else if(std::find(fiducial.begin(), fiducial.end(), -1) == fiducial.end() && fiducial.size() == 8){
      volume = (400-fiducial[0]-fiducial[3]-2*fiducial[6])*(400-fiducial[1]-fiducial[4])*(500-fiducial[2]-fiducial[5]-2*fiducial[7]); // [cm^3]
    }
    std::cout<<"volume = "<<volume<<"\n";
    fiducial_mass = 1.3973*volume/1e6; //[tons]
    targets = 6.022e23 * fiducial_mass * 1e3 * 40/ (0.03995); // [/nucleon]
    std::cout<<"Number of targets = "<<targets<<"\n";

  }

  // Helper functions

  // Convert comma separated string to vector of strings
  std::vector<std::string> ToVector(std::string values, std::string delim = ","){

    // Expects comma separated string
    size_t pos = 0;
    std::vector<string> return_values;

    while ((pos = values.find(delim)) != std::string::npos) {
      std::string value = values.substr(0, pos);

      // Get rid of any []
      value.erase(std::remove(value.begin(), value.end(), '['), value.end());
      value.erase(std::remove(value.begin(), value.end(), ']'), value.end());

      return_values.push_back(value);
      values.erase(0, pos + delim.length());
    }
    values.erase(std::remove(values.begin(), values.end(), '['), values.end());
    values.erase(std::remove(values.begin(), values.end(), ']'), values.end());
    return_values.push_back(values);


    return return_values;
  }

  // Convert string to bool
  std::vector<bool> ToBools(std::string values){
    std::vector<bool> return_bools;
    for(auto const& value : ToVector(values)){
      return_bools.push_back(value=="true");
    }
    return return_bools;
  }

  // Convert string to int
  std::vector<int> ToInts(std::string values){
    std::vector<int> return_ints;
    for(auto const& value : ToVector(values)){
      return_ints.push_back(stoi(value));
    }
    return return_ints;
  }

  // Convert string to double
  std::vector<double> ToDoubles(std::string values){
    std::vector<double> return_doubles;
    for(auto const& value : ToVector(values)){
      return_doubles.push_back(stod(value));
    }
    return return_doubles;
  }

  // Convert string to TString
  std::vector<TString> ToTStrings(std::string values){
    std::vector<TString> return_doubles;
    for(auto const& value : ToVector(values)){
      return_doubles.push_back(TString(value));
    }
    return return_doubles;
  }

};

#endif
