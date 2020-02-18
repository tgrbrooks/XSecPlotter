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
void MultiPlots(){

  // Setting the style
  SetStyle();

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
  std::vector<HistManager*> histmans;
  for(size_t file_i = 0; file_i < config->input_file.size(); file_i++){
    HistManager *histman = new HistManager(config, titles, datamans[file_i], binman, file_i);
    histmans.push_back(histman);
  }
  std::cout<<"...Finished.\n";

  // Calculate systematics and associate to histograms
  if(config->show_syst_error){
    std::cout<<"Calculating the systematics...\n";
    for(size_t file_i = 0; file_i < config->input_file.size(); file_i++){
      SystCalculator(config, histmans[file_i], datamans[file_i], file_i);
    }
    std::cout<<"...Finished.\n";
  }
  
  // Calculate how many 1D histograms the choice of variables and binning would produce
  int n_hists = histmans[0]->GetNHists();
  // Ask user if they want to make that many histograms
  std::string response = "y";
  if(n_hists>20 && config->plot_slices){
    std::cout<<"This will produce "<<n_hists<<" histograms, continue (y/n)? ";
    std::cin>>response;
  }
  if(response=="n") exit(1);

  // Initialise the plotter
  std::cout<<"Making the plots...\n";
  Plotter *plotter = new Plotter(config, titles);
  ChiSquare chsq;
  // Print summary
  // Plot all the 1D plots
  for(size_t i = 0; i < config->plot_variables.size(); i++){
    std::vector<Histo1D*> histos_1D;
    for(size_t file_i = 0; file_i < config->input_file.size(); file_i++){
      histos_1D.push_back(histmans[file_i]->GetHisto1D(i));
    }
    plotter->All1DPlots(histos_1D, i);
    if(config->show_syst_error){
      for(auto const& systname : config->systematics){
        if(config->show_error_band) plotter->Plot1DWithErrors(histos_1D, i, systname);
        else plotter->Plot1D(histos_1D, i, systname);
      }
    }
    // Only plot correlation and response for first file
    if(config->plot_correlation && config->show_syst_error){
      plotter->PlotAllSysts(histmans[0]->GetHisto1D(i));
    }
    if(config->plot_response){
      plotter->PlotResponse(histmans[0]->GetHisto1D(i)->response);
    }
    if(config->input_file.size() == 2){
      std::pair<double, int> chindof = chsq.Calculate(histmans[0]->GetHisto1D(i)->total_hist, histmans[1]->GetHisto1D(i));
      std::cout<<"1D chi^2 = "<<chindof.first<<", ndof = "<<chindof.second<<" chi^2/ndof = "<<chindof.first/chindof.second<<"\n";
    }
  }
  // Only plot slices in the second variable
  if(config->plot_variables.size() == 2){
    std::vector<Histo2D*> histos_2D;
    for(size_t file_i = 0; file_i < config->input_file.size(); file_i++){
      histos_2D.push_back(histmans[file_i]->GetHisto2D(0, 1));
    }
    plotter->Plot2DHisto(histmans[0]->GetHisto2D(0, 1), 0, 1);
    if(config->plot_slices){
      plotter->Plot2DSlices(histos_2D, 0);
    }
    // Only plot correlation and response for first file
    if(config->plot_correlation && config->show_syst_error){
      plotter->PlotAllSysts(histmans[0]->GetHisto2D(0, 1));
    }
    if(config->plot_response){
      plotter->PlotResponse(histmans[0]->GetHisto2D(0, 1)->response);
    }
    if(config->input_file.size() == 2){
      std::pair<double, int> chindof = chsq.Calculate(histmans[0]->GetHisto2D(0, 1)->total_hist, histmans[1]->GetHisto2D(0, 1));
      std::cout<<"2D chi^2 = "<<chindof.first<<", ndof = "<<chindof.second<<" chi^2/ndof = "<<chindof.first/chindof.second<<"\n";
    }
  }
  std::cout<<"...Finished.\n";

}
