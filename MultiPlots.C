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
  HistManager *histman = new HistManager(config, datamans[0], binman, 0);
  std::cout<<"...Finished.\n";

  // Calculate systematics and associate to histograms
  std::cout<<"Calculating the systematics...\n";
  SystCalculator(config, histman, datamans[0], 0);
  std::cout<<"...Finished.\n";

  
  // Calculate how many 1D histograms the choice of variables and binning would produce
  int n_hists = histman->GetNHists();
  // Ask user if they want to make that many histograms
  std::string response = "y";
  if(n_hists>10){
    std::cout<<"This will produce "<<n_hists<<" histograms, continue (y/n)? ";
    std::cin>>response;
  }
  if(response=="n") exit(1);

  // Initialise the plotter
  std::cout<<"Making the plots...\n";
  Plotter *plotter = new Plotter(config, titles);
  for(size_t i = 0; i < config->plot_variables.size(); i++){
    if(!config->show_plots[i]) continue;
    if(config->show_stat_error) plotter->Plot1DWithErrors(histman->GetHisto1D(config->plot_variables[i]), i);
    else plotter->Plot1D(histman->GetHisto1D(config->plot_variables[i]), i);
    if(config->plot_eff_pur && config->stage == "reco"){
      plotter->PlotEffPur(histman->GetHisto1D(config->plot_variables[i]), i);
    }
    if(config->plot_correlation){
      plotter->PlotAllSysts(histman->GetHisto1D(config->plot_variables[i]));
    }
  }
  std::cout<<"...Finished.\n";


/*
  PlotManager *plotman = new PlotManager(config, titles);
  std::vector<std::vector<double>> bin_edges = binman->GetBinning(datamans[0]->total_data);
  //std::vector<std::vector<double>> bin_edges = binman->bin_edges;

  SystManager *systman = new SystManager(config, plotman);

  // Loop over the number of variables - this is how many sets of 1D hists we will have
  std::cout<<"Making the plots...\n";
  for(size_t d_i = 0; d_i < config->plot_variables.size(); d_i++){

    // Don't do anything if we don't want to see the plots for this variable
    if(!config->show_plots[d_i]) continue;

    // Get the file name and title of the histogram
    TString name_1D = config->plot_variables[d_i];
    TString title_1D = titles->hist_titles[d_i];

    // Get the statistical errors per bin
    std::vector<TH1D*> total_hist;
    std::vector<TH1D*> syst_hist;
    TH1D* error_band;
    std::pair<THStack*, TLegend*> stack;
    for(size_t i = 0; i < config->input_file.size(); i++){
      total_hist.push_back(histman->GetTotalHist(datamans[i]->total_data, name_1D, bin_edges, i, d_i));
      syst_hist.push_back(total_hist[i]);
      double stat_err;
      int integral = total_hist[i]->IntegralAndError(0, total_hist[i]->GetNbinsX()+1, stat_err);
      std::cout<<"Stat only = "<<integral<<" +/- "<<stat_err<<"\n";
      if(config->show_syst_error){
        systman->AddErrors(syst_hist[i], systman->GetSystHist(datamans[i]->total_data, datamans[i]->data_used, name_1D, bin_edges, i, d_i));
        //systman->AddErrors(syst_hist[i], systman->GetDetSystHist(name_1D, bin_edges, i, d_i));
        //systman->AddErrors(syst_hist[i], systman->GetBkgSystHist(name_1D, bin_edges, i, d_i));
      }
      if(config->input_file.size() == 1){
        error_band = histman->GetErrorBand(total_hist[0]);
        // Create a total 1D stacked histogram for each of the variables
        stack = histman->StackHist1D(datamans[0]->stack_data, name_1D, title_1D, bin_edges, d_i);
      }
    }

    // Draw the plots
    if(config->input_file.size() == 1){
      if(config->show_stat_error) plotman->Plot1DWithErrors(stack.first, stack.second, error_band, total_hist[0], syst_hist[0], d_i);
      else plotman->Plot1D(stack.first, stack.second, total_hist[0], syst_hist[0], d_i);
    }
    else plotman->PlotMulti1D(total_hist, syst_hist, d_i);

    // Plot efficiency and purity if option selected and reconstruction selected 
    if(config->input_file.size() == 1 && config->plot_eff_pur && config->stage == "reco"){
      plotman->PlotEffPur(datamans[0]->interactions, name_1D, bin_edges, d_i);
    }

    // Loop over the other variables - this is how many sets of 2D hists we will have
    for(size_t d_j = 0; d_j < config->plot_variables.size(); d_j++){
      if(d_j == d_i) continue;

      // Create a total 2D histogram for each combination of variables
      if(config->input_file.size() == 1){
        TH2D* hist_2D = histman->Get2DHist(datamans[0]->total_data, bin_edges, name_1D, d_i, d_j);
        plotman->Plot2D(hist_2D, config->plot_variables[d_i]+"_"+config->plot_variables[d_j], titles->names[d_i]+" ["+titles->units[d_i]+"]", titles->names[d_j]+" ["+titles->units[d_j]+"]");
      }

      // Loop over the bins for variable 2
      for(size_t bin_j = 0; bin_j < bin_edges[d_j].size()-1; bin_j++){
        // Only make these plots for 2 variables
        if(config->plot_variables.size() != 2) continue;

        // Rebin so every bin below maximum error if set
        std::vector<std::vector<double>> bin_edges_copy = bin_edges;
        if(config->max_error>0){
          std::vector<double> bin_edges_new = binman->ChangeBinning2D(datamans[0]->total_data, bin_edges, d_i, d_j, bin_j);
          bin_edges_copy[d_i] = bin_edges_new;
        }

        // Get the file name and title of the histogram
        TString name_2D = config->plot_variables[d_i] +"_"
                          + config->plot_variables[d_j] +"_"+ Form("%.1f", bin_edges_copy[d_j][bin_j]) +"_"+ Form("%.1f", bin_edges_copy[d_j][bin_j+1]);
        TString title_2D = titles->hist_titles[d_j] 
                           +": ["+ Form("%.2f", bin_edges_copy[d_j][bin_j]) +", "+ Form("%.2f", bin_edges_copy[d_j][bin_j+1]) +"]";

        // Get the statistical errors per bin
        std::vector<TH1D*> total_hist_2D;
        std::vector<TH1D*> syst_hist_2D;
        TH1D* error_band_2D;
        std::pair<THStack*, TLegend*> stack_2D;
        for(size_t i = 0; i < config->input_file.size(); i++){
          total_hist_2D.push_back(histman->GetTotalHist(datamans[i]->total_data, name_2D, bin_edges_copy, i, d_i, d_j, bin_j));
          syst_hist_2D.push_back(total_hist_2D[i]);
          if(config->show_syst_error){
            systman->AddErrors(syst_hist_2D[i], systman->GetSystHist(datamans[i]->total_data, datamans[i]->data_used, name_2D, bin_edges_copy, i, d_i, d_j, bin_j));
            systman->AddErrors(syst_hist_2D[i], systman->GetDetSystHist(name_2D, bin_edges_copy, i, d_i, d_j, bin_j));
            systman->AddErrors(syst_hist_2D[i], systman->GetBkgSystHist(name_2D, bin_edges_copy, i, d_i));
          }

          if(config->input_file.size() == 1){
            error_band_2D = histman->GetErrorBand(total_hist_2D[0]);
            //if(error_band_2D->Integral(0, error_band_2D->GetNbinsX()) == 0) continue;

            // Create a 1D stacked histogram for each of the bins
            stack_2D = histman->StackHist1D(datamans[0]->stack_data, name_2D, title_2D, bin_edges_copy, d_i, d_j, bin_j);
          }
        }

        // Draw the plots
        if(config->input_file.size() == 1){
          if(config->show_stat_error) plotman->Plot1DWithErrors(stack_2D.first, stack_2D.second, error_band_2D, total_hist_2D[0], syst_hist_2D[0], d_i, d_j);
          else plotman->Plot1D(stack_2D.first, stack_2D.second, total_hist_2D[0], syst_hist_2D[0], d_i, d_j);
        }
        else plotman->PlotMulti1D(total_hist_2D, syst_hist_2D, d_i, d_j);

        // Plot efficiency and purity if option selected and reconstruction selected
        if(config->input_file.size() == 1 && config->plot_eff_pur && config->stage == "reco"){
          plotman->PlotEffPur(datamans[0]->interactions, name_2D, bin_edges_copy, d_i, d_j, bin_j);
        }
        
      }

    }
  }
  std::cout<<"...Finished.\n";
  */

}
