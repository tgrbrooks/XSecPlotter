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
#include "Types/PlotManager.h"
#include "Types/Plotter.h"
#include "Types/HistManager.h"
#include "Types/DataManager.h"
#include "Types/SystManager.h"
#include "Types/SystCalculator.h"
#include "Types/Interaction.h"
#include "Types/Selection.h"
#include "Functions/StyleSetter.h"


// Main
void MultiPlots(){

  // Setting the style
  SetStyle();

  // Get the configuration file
  std::cout<<"Reading the config file...\n";
  std::string input_file = "multiconfig.txt";
  Configuration *config = new Configuration(input_file);
  config->GetMetaData();
  std::cout<<"...Finished.\n";

  // Get from configuration
  std::cout<<"Getting labels...\n";
  Titles *titles = new Titles(config);
  std::cout<<"...Finished.\n";

  // Get the plotting variable
  std::cout<<"Reading from the tree...\n";
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

  // Get the binning using the first input file
  std::cout<<"Getting the correct binning...\n";
  BinManager *binman = new BinManager(config, datamans[0]);
  std::cout<<"...Finished.\n";

  // Make every possible histogram
  std::cout<<"Creating all of the histograms...\n";
  HistManager *histman = new HistManager(config, titles, datamans[0], binman, 0);
  std::cout<<"...Finished.\n";

  // Calculate systematics and associate to histograms
  if(config->show_syst_error){
    std::cout<<"Calculating the systematics...\n";
    SystCalculator(config, histman, datamans[0], 0);
    std::cout<<"...Finished.\n";
  }
  
  // Calculate how many 1D histograms the choice of variables and binning would produce
  int n_hists = histman->GetNHists();
  // Ask user if they want to make that many histograms
  std::string response = "y";
  if(n_hists>20){
    std::cout<<"This will produce "<<n_hists<<" histograms, continue (y/n)? ";
    std::cin>>response;
  }
  if(response=="n") exit(1);

  // Initialise the plotter
  std::cout<<"Making the plots...\n";
  Plotter *plotter = new Plotter(config, titles);
  // Print summary
  // Plot all the 1D plots
  for(size_t i = 0; i < config->plot_variables.size(); i++){
    if(!config->show_plots[i]) continue;
    plotter->All1DPlots(histman->GetHisto1D(config->plot_variables[i]), i);
    if(config->plot_correlation){
      plotter->PlotAllSysts(histman->GetHisto1D(config->plot_variables[i]));
    }
  }
  // Only plot slices in the second variable
  if(config->plot_variables.size() == 2){
    std::pair<TString, TString> key = std::make_pair(config->plot_variables[0], config->plot_variables[1]);
    plotter->Plot2DHisto(histman->GetHisto2D(config->plot_variables[0], config->plot_variables[1]), 0, 1);
    plotter->Plot2DSlices(histman->GetHisto2D(config->plot_variables[0], config->plot_variables[1]), 0, 1);
    if(config->plot_correlation){
      plotter->PlotAllSysts(histman->GetHisto2D(config->plot_variables[0], config->plot_variables[1]));
    }
  }
  std::cout<<"...Finished.\n";

}
