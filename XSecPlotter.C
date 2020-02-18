#include <string>
#include <vector>
#include <TChain.h>
#include <TTree.h>

#include <map>
#include <cmath>
#include <iostream>
#include <ctime>
#include <algorithm>

#include <TFile.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLine.h>
#include <TVector3.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TH3D.h>
#include <THStack.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCanvas.h>
#include <sstream>
#include <fstream>
#include <math.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "Types/Configuration.h"
#include "Types/Titles.h"
#include "Types/BinManager.h"
#include "Types/Plotter.h"
#include "Types/HistManager.h"
#include "Types/DataManager.h"
#include "Types/SystCalculator.h"
#include "Types/Interaction.h"
#include "Types/Selection.h"
#include "Types/ChiSquare.h"
#include "Functions/StyleSetter.h"


// Main
void XSecPlotter(){

  // Setting the style
  SetStyle();

  //-------------------------------------------------------------------------
  //                            CONFIGURATION
  //-------------------------------------------------------------------------

  // Get the configuration file
  std::cout<<"Reading the config file...\n";
  std::string input_file = "config.txt";
  Configuration *config = new Configuration(input_file);
  config->GetMetaData();
  std::cout<<"...Finished.\n";

  // Get from configuration
  std::cout<<"Getting labels...\n";
  Titles *titles = new Titles(config);
  std::cout<<"...Finished.\n";

  //-------------------------------------------------------------------------
  //                          READING IN DATA
  //-------------------------------------------------------------------------

  // Get the plotting variables from input files
  std::cout<<"Reading from the tree...\n";
  // Loop over input files and create a data manager for each one
  std::vector<DataManager*> datamans;
  for(size_t file_i = 0; file_i < config->input_file.size(); file_i++){
    DataManager *dataman = new DataManager(config, file_i);
    datamans.push_back(dataman);
  }
  std::cout<<"...Finished.\n";

  if(datamans.size() < 1){
    std::cout<<"Not enough input files\n";
    exit(1);
  }

  //-------------------------------------------------------------------------
  //                          CREATING HISTOGRAMS
  //-------------------------------------------------------------------------

  // Get the binning using the data from the first input file
  // Only one bin manager instance as binning can't be different for the two files
  std::cout<<"Getting the correct binning...\n";
  BinManager *binman = new BinManager(config, datamans[0]);
  std::cout<<"...Finished.\n";

  // Make all of the histograms for both files
  std::cout<<"Creating all of the histograms...\n";
  std::vector<HistManager*> histmans;
  for(size_t file_i = 0; file_i < config->input_file.size(); file_i++){
    HistManager *histman = new HistManager(config, titles, datamans[file_i], binman, file_i);
    histmans.push_back(histman);
  }
  std::cout<<"...Finished.\n";

  //-------------------------------------------------------------------------
  //                  CALCULATING SYSTEMATIC UNCERTAINTIES
  //-------------------------------------------------------------------------
  
  // Calculate systematics and associate to histograms if requested
  if(config->show_syst_error){
    std::cout<<"Calculating the systematics...\n";
    for(size_t file_i = 0; file_i < config->input_file.size(); file_i++){
      SystCalculator(config, histmans[file_i], datamans[file_i], file_i);
    }
    std::cout<<"...Finished.\n";
  }
  
  // Calculate how many histograms the choice of variables and binning would produce
  int n_hists = histmans[0]->GetNHists();
  // Ask user if they want to make that many histograms if it's big
  std::string response = "y";
  if(n_hists>20 && config->plot_slices){
    std::cout<<"This will produce "<<n_hists<<" histograms(+ extras), continue (y/n)? ";
    std::cin>>response;
  }
  if(response=="n") exit(1);

  //-------------------------------------------------------------------------
  //                              PLOTTING
  //-------------------------------------------------------------------------

  // Initialise the plotter
  std::cout<<"Making the plots...\n";
  Plotter *plotter = new Plotter(config, titles);

  // Create a chi2 calculater for comparing models
  ChiSquare chsq;

  //-------------------------------------------------------------------------
  //                              1D PLOTS
  //-------------------------------------------------------------------------
  
  // Loop over the plotting variables
  for(size_t i = 0; i < config->plot_variables.size(); i++){
    std::vector<Histo1D*> histos_1D;

    // Fill a vector of 1D histograms with each input file
    for(size_t file_i = 0; file_i < config->input_file.size(); file_i++){
      histos_1D.push_back(histmans[file_i]->GetHisto1D(i));
    }

    // Plot the 1D histograms
    plotter->All1DPlots(histos_1D, i);

    // If the systematics are calculated make individual plots for each
    if(config->show_syst_error){
      for(auto const& systname : config->systematics){
        if(config->show_error_band) plotter->Plot1DWithErrors(histos_1D, i, systname);
        else plotter->Plot1D(histos_1D, i, systname);
      }
    }

    // Plot correlation, covariance, and universes for first file if selected
    if(config->plot_correlation && config->show_syst_error){
      plotter->PlotAllSysts(histmans[0]->GetHisto1D(i));
    }

    // Plot response matrix for first file if selected
    if(config->plot_response){
      plotter->PlotResponse(histmans[0]->GetHisto1D(i)->response);
    }

    // If two input files used calculate chi2 between models
    if(config->input_file.size() == 2){
      std::pair<double, int> chindof = chsq.Calculate(histmans[0]->GetHisto1D(i)->total_hist, histmans[1]->GetHisto1D(i));
      std::cout<<"1D chi^2 = "<<chindof.first<<", ndof = "<<chindof.second<<" chi^2/ndof = "<<chindof.first/chindof.second<<"\n";
    }
  }

  //-------------------------------------------------------------------------
  //                              2D PLOTS
  //-------------------------------------------------------------------------

  // Only plot slices in the second variable
  if(config->plot_variables.size() == 2){

    // Fill a vector of 2D histograms for each file
    std::vector<Histo2D*> histos_2D;
    for(size_t file_i = 0; file_i < config->input_file.size(); file_i++){
      histos_2D.push_back(histmans[file_i]->GetHisto2D(0, 1));
    }

    // Plot the 2D histogram
    plotter->Plot2DHisto(histmans[0]->GetHisto2D(0, 1), 0, 1);

    // If selected, plot all of the slices
    if(config->plot_slices){
      plotter->Plot2DSlices(histos_2D, 0);
    }

    // Plot correlation, covariance and universes for first file if selected
    if(config->plot_correlation && config->show_syst_error){
      plotter->PlotAllSysts(histmans[0]->GetHisto2D(0, 1));
    }

    // Plot response matrix for first file if selected
    if(config->plot_response){
      plotter->PlotResponse(histmans[0]->GetHisto2D(0, 1)->response);
    }

    // If more than two input files, calculate chi2
    if(config->input_file.size() == 2){
      std::pair<double, int> chindof = chsq.Calculate(histmans[0]->GetHisto2D(0, 1)->total_hist, histmans[1]->GetHisto2D(0, 1));
      std::cout<<"2D chi^2 = "<<chindof.first<<", ndof = "<<chindof.second<<" chi^2/ndof = "<<chindof.first/chindof.second<<"\n";
    }
  }
  std::cout<<"...Finished.\n";

}