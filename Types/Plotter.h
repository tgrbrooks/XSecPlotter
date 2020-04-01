#ifndef PLOTTER_H
#define PLOTTER_H

#include "Configuration.h"
#include "Titles.h"
#include "Interaction.h"
#include "Systematics.h"
#include "Histo1D.h"
#include "Histo2D.h"
#include "HistManager.h"
#include "XSecCalculator.h"

// Class for handling all plotting
class Plotter
{
  public:

  Configuration *config; // Global configuration
  Titles *titles;        // Histogram title manager

  // Constructor
  Plotter(Configuration *c, Titles *t)
  {
    config = c;
    titles = t;
  }

  // ------------------------------------------------------------------------------------------------------------------
  //                                        1D PLOTTING FUNCTIONS
  // ------------------------------------------------------------------------------------------------------------------
  
  // Create all 1D plots defined by configuration
  void All1DPlots(std::vector<Histo1D*> histos, size_t var_i, size_t var_j=-1){
    // Different plotters for cross section and rate
    if(config->plot_xsec){
      if(config->show_error_band) Plot1DWithErrorsXSec(histos, var_i, var_j);
      else Plot1DXSec(histos, var_i, var_j);
    }
    else{
      if(config->show_error_band) Plot1DWithErrors(histos, var_i);
      else Plot1D(histos, var_i);
    }
    if(config->plot_eff_pur && config->stage == "reco"){
      PlotEffPur(histos[0], var_i);
    }
    if(config->plot_universes && config->show_syst_error){
      PlotUniverses(histos[0], var_i);
    }
  }

  // Plot a 1D stacked hist with statistical errors on the bottom
  void Plot1DWithErrors(std::vector<Histo1D*> histos, size_t var_i, TString systname="total"){

    // Create the canvas
    TString name = histos[0]->stacked_hist->GetName();
    if(systname != "total") name = name + "_" + systname;
    TCanvas *canvas = new TCanvas("canvas"+name,"canvas");

    // Split the pad for histogram and error plot
    double pad_split = .3;
    TPad *upper_pad = new TPad("upper_pad", "" , 0., pad_split, 1.0, 1.0);
    upper_pad->SetMargin(0.12, 0.05, 0.075, 0.12); //LRBT

    TPad *lower_pad = new TPad("lower_pad", "", 0., 0., 1., pad_split);
    lower_pad->SetMargin(0.12, 0.05, 0.34, 0.04);

    upper_pad->Draw();
    lower_pad->Draw();

    // Fill the upper pad with histogram, info and legend
    upper_pad->cd();

    // Draw the stacked histogram and legend
    std::vector<TH1D*> error_hists = Draw1D(histos, systname);
    DrawLegends(histos);
    canvas->Update();
    canvas->Modified();
    
    double ymax = 0;
    for(size_t file_i = 0; file_i < config->input_file.size(); file_i++){
      int binmax = histos[file_i]->total_hist->GetMaximumBin();
      double max = (histos[file_i]->total_hist->GetBinContent(binmax) + error_hists[file_i]->GetBinError(binmax));
      if(max > ymax) ymax = max;
    }
    AxisConfig(histos[0], ymax, var_i, -1);
    canvas->Modified();

    // Info text
    // Text position and content
    double width = 0.7*(histos[0]->total_hist->GetXaxis()->GetXmax()
                        - histos[0]->total_hist->GetXaxis()->GetXmin())
                        + histos[0]->total_hist->GetXaxis()->GetXmin();
    double upper_text_size = 0.7*histos[0]->total_hist->GetYaxis()->GetTitleSize();
    if(config->show_info) DrawInfo(width, ymax, upper_text_size);
   
    // Fill the lower pad with percentage error per bin
    lower_pad->cd();
    lower_pad->SetTickx();
    lower_pad->SetTicky();
   
    if(config->show_syst_error && config->show_stat_error) histos[0]->SystErrorBand(systname, true);
    else if(config->show_syst_error) histos[0]->SystErrorBand(systname, false);
    double size_ratio = upper_pad->GetAbsHNDC()/lower_pad->GetAbsHNDC();
    DrawErrorBand(histos[0]->error_band, var_i, size_ratio, systname);
    
    TString output_file = config->output_file;
    output_file.ReplaceAll(".","_"+name+".");
    canvas->SaveAs(output_file);
  }

  // Plot a 1D stacked hist with statistical errors on the bottom
  void Plot1DWithErrorsXSec(std::vector<Histo1D*> histos, size_t var_i, size_t var_j=-1, TString systname="total"){

    // Create the canvas
    TString name = histos[0]->stacked_hist->GetName();
    if(systname != "total") name = name + "_" + systname;
    TCanvas *canvas = new TCanvas("canvas"+name,"canvas");

    histos[0]->xsec_hist->SetTitle(histos[0]->stacked_hist->GetTitle());

    // Split the pad for histogram and error plot
    double pad_split = .3;
    TPad *upper_pad = new TPad("upper_pad", "" , 0., pad_split, 1.0, 1.0);
    upper_pad->SetMargin(0.12, 0.05, 0.075, 0.12); //LRBT

    TPad *lower_pad = new TPad("lower_pad", "", 0., 0., 1., pad_split);
    lower_pad->SetMargin(0.12, 0.05, 0.34, 0.04);

    upper_pad->Draw();
    lower_pad->Draw();

    // Fill the upper pad with histogram, info and legend
    upper_pad->cd();

    // Draw the stacked histogram and legend
    std::vector<TH1D*> error_hists = Draw1D(histos, systname);
    DrawLegends(histos);
    canvas->Update();
    canvas->Modified();
    
    double ymax = 0;
    for(size_t file_i = 0; file_i < config->input_file.size(); file_i++){
      int binmax = histos[file_i]->xsec_hist->GetMaximumBin();
      double max = (histos[file_i]->xsec_hist->GetBinContent(binmax) + error_hists[file_i]->GetBinError(binmax));
      if(max > ymax) ymax = max;
    }
    AxisConfig(histos[0], ymax, var_i, var_j);
    canvas->Modified();

    // Info text
    // Text position and content
    double width = 0.7*(histos[0]->xsec_hist->GetXaxis()->GetXmax() 
                        - histos[0]->xsec_hist->GetXaxis()->GetXmin()) 
                        + histos[0]->xsec_hist->GetXaxis()->GetXmin();
    double upper_text_size = 0.7*histos[0]->xsec_hist->GetYaxis()->GetTitleSize();
    if(config->show_info) DrawInfo(width, ymax, upper_text_size);
   
    // Fill the lower pad with percentage error per bin
    lower_pad->cd();
    lower_pad->SetTickx();
    lower_pad->SetTicky();
   
    if(config->show_syst_error && config->show_stat_error) histos[0]->XSecErrorBand(systname, true);
    else if(config->show_syst_error) histos[0]->XSecErrorBand(systname, false);
    else if(config->show_stat_error) histos[0]->XSecErrorBand("", true);
    double size_ratio = upper_pad->GetAbsHNDC()/lower_pad->GetAbsHNDC();
    DrawErrorBand(histos[0]->error_band, var_i, size_ratio, systname);

    TString output_file = config->output_file;
    output_file.ReplaceAll(".","_"+name+".");
    canvas->SaveAs(output_file);
  }

  // Plot a 1D stacked hist
  void Plot1D(std::vector<Histo1D*> histos, size_t var_i, TString systname="total"){

    // Create the canvas
    TString name = histos[0]->stacked_hist->GetName();
    if(systname != "total") name = name + "_" + systname;
    TCanvas *canvas = new TCanvas("canvas"+name,"canvas");
    canvas->SetMargin(0.12, 0.04, 0.15, 0.15); //LRBT

    // Draw the stacked histogram and legend
    std::vector<TH1D*> error_hists = Draw1D(histos, systname);
    DrawLegends(histos);

    double ymax = 0;
    for(size_t file_i = 0; file_i < histos.size(); file_i++){
      int binmax = histos[file_i]->total_hist->GetMaximumBin();
      double max = (histos[file_i]->total_hist->GetBinContent(binmax) + error_hists[file_i]->GetBinError(binmax));
      if(max > ymax) ymax = max;
    }
    AxisConfig(histos[0], ymax, var_i, -1);
    canvas->Update();
    canvas->Modified();

    // Text position and content
    double width = 0.65*(histos[0]->total_hist->GetXaxis()->GetXmax() 
                         - histos[0]->total_hist->GetXaxis()->GetXmin())
                         + histos[0]->total_hist->GetXaxis()->GetXmin();
    double upper_text_size = 0.6*histos[0]->total_hist->GetYaxis()->GetTitleSize();
    if(config->show_info) DrawInfo(width, ymax, upper_text_size);
    
    TString output_file = config->output_file;
    output_file.ReplaceAll(".","_"+name+".");
    canvas->SaveAs(output_file);
  }

  // Plot a 1D stacked hist
  void Plot1DXSec(std::vector<Histo1D*> histos, size_t var_i, size_t var_j=-1, TString systname="total"){

    // Create the canvas
    TString name = histos[0]->stacked_hist->GetName();
    if(systname != "total") name = name + "_" + systname;
    TCanvas *canvas = new TCanvas("canvas"+name,"canvas");
    canvas->SetMargin(0.16, 0.04, 0.15, 0.15); //LRBT

    histos[0]->xsec_hist->SetTitle(histos[0]->stacked_hist->GetTitle());

    // Draw the stacked histogram and legend
    std::vector<TH1D*> error_hists = Draw1D(histos, systname);
    DrawLegends(histos);
    canvas->Update();
    canvas->Modified();

    double ymax = 0;
    for(size_t file_i = 0; file_i < config->input_file.size(); file_i++){
      int binmax = histos[file_i]->xsec_hist->GetMaximumBin();
      double max = (histos[file_i]->xsec_hist->GetBinContent(binmax) + error_hists[file_i]->GetBinError(binmax));
      if(max > ymax) ymax = max;
    }
    AxisConfig(histos[0], ymax, var_i, var_j);
    canvas->Modified();

    // Text position and content
    double width = 0.65*(histos[0]->xsec_hist->GetXaxis()->GetXmax() 
                         - histos[0]->xsec_hist->GetXaxis()->GetXmin()) 
                         + histos[0]->xsec_hist->GetXaxis()->GetXmin();
    double upper_text_size = 0.6 * 0.06;
    if(config->show_info) DrawInfo(width, ymax, upper_text_size);
    
    TString output_file = config->output_file;
    output_file.ReplaceAll(".","_"+name+".");
    canvas->SaveAs(output_file);
  }

  // ------------------------------------------------------------------------------------------------------------------
  //                                        2D PLOTTING FUNCTIONS
  // ------------------------------------------------------------------------------------------------------------------

  // Plot all slices in Y bins
  void Plot2DSlices(std::vector<Histo2D*> histos, size_t var_i, size_t var_j){
    // Loop over all bins in the second variable
    for(size_t j = 0; j < histos[0]->ybins.size()-1; j++){
      std::vector<Histo1D*> histos_1D;
      for(size_t i = 0; i < histos.size(); i++){
        histos_1D.push_back(histos[i]->Slice(j));
        histos_1D[i]->xsec_hist->SetFillStyle(config->fsty[i]);
        histos_1D[i]->xsec_hist->SetLineStyle(config->lsty[i]);
      }
      // Make all relevant plots
      All1DPlots(histos_1D, var_i, var_j);
    }
  }

  // Plot 2D histogram with correct axis labels
  void Plot2DHisto(Histo2D* histo, size_t var_i, size_t var_j, bool bins=false){

    TString name = config->plot_variables[var_i]+"_"+config->plot_variables[var_j]+"2D";
    if(bins) name = name + "_bins";

    TString xaxis = titles->names[var_i]+" ["+titles->units[var_i]+"]";
    if(titles->units[var_i]=="") xaxis = titles->names[var_i];

    TString yaxis = titles->names[var_j]+" ["+titles->units[var_j]+"]";
    if(titles->units[var_j]=="") yaxis = titles->names[var_j];

    TCanvas *canvas = new TCanvas(name, "", 900, 900);
    canvas->SetFrameLineWidth(4.);
    canvas->SetLineWidth(4.);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetMargin(0.16, 0.16, 0.16, 0.16);

    // Scale by width because this doesn't work with TH2Poly
    TH2Poly* hist = (TH2Poly*)histo->total_hist->Clone();
    for(auto const& obj : *hist->GetBins()){
      TH2PolyBin *bin = (TH2PolyBin*)obj;
      double wy = abs(bin->GetYMax() - bin->GetYMin());
      double wx = abs(bin->GetXMax() - bin->GetXMin());
      double width = wy*wx;
      int j = bin->GetBinNumber();
      if(config->plot_xsec){
        hist->SetBinContent(j, histo->xsec_hist->GetBinContent(j)/width);
      }
      else{
        hist->SetBinContent(j, histo->total_hist->GetBinContent(j)/width);
      }
      if(bins) hist->SetBinContent(j, j);
    }

    hist->SetMarkerSize(1.1);
    gStyle->SetPaintTextFormat("2.0f");
    hist->GetXaxis()->SetTitle(xaxis);
    hist->GetYaxis()->SetTitle(yaxis);
    // X axis config
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetXaxis()->SetTickLength(0.04);
    hist->GetXaxis()->SetTitleSize(1.1 * 0.06);
    hist->GetXaxis()->SetNdivisions(108);
    // Y axis config
    hist->GetYaxis()->SetTitleOffset(1.15);
    hist->GetYaxis()->SetTickLength(0.015);
    hist->GetYaxis()->SetTitleSize(1.1 * 0.06);
    hist->GetYaxis()->SetNdivisions(108);
    // Text size
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetZaxis()->SetLabelSize(0.035);

    if(!bins){
      hist->Draw("colz");
      hist->GetZaxis()->SetMaxDigits(3);
      /*TPaletteAxis *palette = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
      palette->SetMaxDigits(3);*/
      canvas->Modified();
      canvas->Update();
    }
    else hist->Draw("text");

    TString output_file = config->output_file;
    output_file.ReplaceAll(".","_"+name+".");
    canvas->SaveAs(output_file);
  }

  // Plot a 2D histogram
  void Plot2D(TH2D* hist, TString name, TString xaxis, TString yaxis){

    TCanvas *canvas = new TCanvas(name, "", 900, 900);
    canvas->SetFrameLineWidth(4.);
    canvas->SetLineWidth(4.);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetMargin(0.16, 0.16, 0.16, 0.16);

    hist->GetXaxis()->SetTitle(xaxis);
    hist->GetYaxis()->SetTitle(yaxis);
    // X axis config
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetXaxis()->SetTickLength(0.04);
    hist->GetXaxis()->SetTitleSize(1.1 * 0.06);
    hist->GetXaxis()->SetNdivisions(108);
    // Y axis config
    hist->GetYaxis()->SetTitleOffset(1.15);
    hist->GetYaxis()->SetTickLength(0.015);
    hist->GetYaxis()->SetTitleSize(1.1 * 0.06);
    hist->GetYaxis()->SetNdivisions(108);
    // Text size
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetZaxis()->SetLabelSize(0.035);

    hist->Draw("colz");

    /*TPaletteAxis *palette = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
    palette->SetMaxDigits(3);*/
    hist->GetZaxis()->SetMaxDigits(3);
    canvas->Modified();
    canvas->Update();

    TString output_file = config->output_file;
    output_file.ReplaceAll(".","_"+name+".");
    canvas->SaveAs(output_file);
  }

  // ------------------------------------------------------------------------------------------------------------------
  //                                      GENERAL HELPER FUNCTIONS
  // ------------------------------------------------------------------------------------------------------------------

  // Draw additional information on to hist
  void DrawInfo(double width, double height, double size){

    height = height - 0.05*height;
    TLatex *POT        = new TLatex(width, .97*height, titles->pot);
    TLatex *mass       = new TLatex(width, .91*height, titles->mass);
    TLatex *data_type  = new TLatex(width, .85*height, titles->data_type);
    TLatex *is_cc      = new TLatex(width, .79*height, titles->is_cc);
    TLatex *part_cont  = new TLatex(width, .72*height, titles->part_cont);
    TLatex *lep_cont   = new TLatex(width, .66*height, titles->lep_cont);
    TLatex *n_pr       = new TLatex(width, .60*height, titles->n_pr);
    TLatex *n_pipm     = new TLatex(width, .54*height, titles->n_pipm);
    TLatex *n_pi0      = new TLatex(width, .48*height, titles->n_pi0);
    TLatex *int_type   = new TLatex(width, .60*height, titles->int_type);

    // Set the text size
    POT->SetTextSize(size);
    mass->SetTextSize(size);
    data_type->SetTextSize(size);
    part_cont->SetTextSize(size);
    lep_cont->SetTextSize(size);
    is_cc->SetTextSize(size);
    n_pr->SetTextSize(size);
    n_pipm->SetTextSize(size);
    n_pi0->SetTextSize(size);
    int_type->SetTextSize(size);

    // Draw the info text
    POT->Draw("same");
    mass->Draw("same");
    data_type->Draw("same");
    part_cont->Draw("same");
    lep_cont->Draw("same");
    is_cc->Draw("same");
    if(config->plot_by_fsi){
      n_pr->Draw("same");
      n_pipm->Draw("same");
      n_pi0->Draw("same");
    }
    else{
      int_type->Draw("same");
    }

  }

  // Function for main 1D plot and errors
  std::vector<TH1D*> Draw1D(std::vector<Histo1D*> histos, TString systname="total"){

    bool first = true;
    std::vector<TH1D*> error_hists;
    // Loop over the input files
    for(size_t file_i = 0; file_i < histos.size(); file_i++){
      // Get a histogram for plotting the errors
      TH1D* error_hist;
      if(config->plot_xsec){
        error_hist = (TH1D*)histos[file_i]->xsec_hist->Clone();
      }
      else{
        error_hist = (TH1D*)histos[file_i]->total_hist->Clone();
      }
      error_hists.push_back(error_hist);

      // Plot the main histogram(s)
      if(first){
        if(config->plot_xsec){
          histos[file_i]->xsec_hist->Draw("HIST");
        }
        else{
          histos[file_i]->stacked_hist->Draw("HIST");
        }
        first = false;
      } else{ 
        if(config->plot_xsec){
          histos[file_i]->xsec_hist->Draw("HIST SAME");
        }
        else{
          histos[file_i]->stacked_hist->Draw("HIST SAME");
        }
      }

      // Set the errors if showing
      if(histos.size() > 1 && file_i == 0) continue;
      if(config->show_error_bars){
        // Add statistical and systematic errors bin by bin
        for(int n = 0; n <= error_hist->GetNbinsX(); n++){
          double stat_error = 0;
          if(config->show_stat_error) stat_error = error_hist->GetBinError(n);
          double syst_error = 0;
          if(config->show_syst_error) syst_error = histos[file_i]->systematics->GetSyst(systname)->mean_syst->GetBinError(n);
          error_hist->SetBinError(n, std::sqrt(std::pow(syst_error, 2) + std::pow(stat_error, 2)));
        }
        // Show the error bars
        error_hist->SetLineWidth(0);
        error_hist->SetMarkerStyle(0);
        error_hist->SetFillColor(15);
        error_hist->SetFillStyle(3001);
        if(config->show_syst_error) error_hist->Draw("E2 SAME");
        // Show stat error separately for rate
        if(!config->plot_xsec && config->show_stat_error){
          histos[file_i]->total_hist->SetLineWidth(0);
          histos[file_i]->total_hist->SetMarkerStyle(0);
          histos[file_i]->total_hist->SetFillColor(12);
          histos[file_i]->total_hist->SetFillStyle(3001);
          histos[file_i]->total_hist->Draw("E2 SAME");
        }
      }  
    }

    return error_hists;

  }

  // Draw the 1D legends
  void DrawLegends(std::vector<Histo1D*> histos){
    // Stacked histogram legend
    if(config->plot_stacked && !config->plot_xsec){
      histos[0]->legend->SetNColumns(histos[0]->legend->GetNRows()*histos[0]->legend->GetNColumns());
      
      histos[0]->legend->SetFillStyle(0);
      histos[0]->legend->Draw();
    }
    // Multiple input file legend
    if(config->input_file.size() > 1){
      //TLegend* legend = new TLegend(0.24, 0.85, 0.96, 0.91);
      TLegend* legend = new TLegend(0.24, 0.89, 0.96, 0.95);
      if(!config->show_error_band){
        legend->SetX1NDC(0.24);
        legend->SetY1NDC(0.91);
        legend->SetX2NDC(0.96);
        legend->SetY2NDC(0.97);
      }
      legend->SetFillStyle(0);
      legend->SetNColumns(config->input_file.size());
      for(size_t file_i = 0; file_i < config->input_file.size(); file_i++){
        TString lname = config->tune_name[file_i];
        lname.ReplaceAll("_", "");
        if(config->plot_xsec){
          histos[file_i]->xsec_hist->SetLineWidth(3);
          legend->AddEntry(histos[file_i]->xsec_hist, lname, "l");
        }
        else{
          histos[file_i]->total_hist->SetLineWidth(3);
          legend->AddEntry(histos[file_i]->total_hist, lname, "l");
        }
      }
      legend->Draw();
    }
  }

  void AxisConfig(Histo1D* hist, double max, int var_i, int var_j){

    double xtitle_offset = 1.8;
    if(!config->show_error_band) xtitle_offset = 1.;
    double xtick_length = 0.04;
    if(!config->show_error_band) xtick_length = 0.02;

    if(config->plot_xsec){
      double ytitle_offset = 1.0;
      if(!config->show_error_band) ytitle_offset = 1.2;
      // Set the titles
      hist->xsec_hist->GetYaxis()->SetTitle(titles->GetXSecTitle(var_i, var_j));
      if(!config->show_error_band){
        hist->xsec_hist->GetXaxis()->SetTitle(titles->names[var_i]+" ["+titles->units[var_i]+"]");
        if(titles->units[var_i]=="") hist->xsec_hist->GetXaxis()->SetTitle(titles->names[var_i]);
      }
      // X axis config
      if(config->show_error_band) hist->xsec_hist->GetXaxis()->SetLabelOffset(0.1);
      hist->xsec_hist->GetXaxis()->SetTitleOffset(xtitle_offset);
      hist->xsec_hist->GetXaxis()->SetTickLength(xtick_length);
      if(!config->show_error_band) hist->xsec_hist->GetXaxis()->SetTitleSize(1.1 * 0.06);
      // Y axis config
      hist->xsec_hist->GetYaxis()->SetTitleOffset(ytitle_offset);
      double title_size = 0.06;
      if(!config->show_error_band) title_size = 0.8*0.06;
      if(config->plot_variables.size()==2){ 
        title_size = 0.8*0.06;
        if(!config->show_error_band) title_size = 0.8*0.06;
        hist->xsec_hist->GetYaxis()->SetTitleOffset(1.1);
        if(!config->show_error_band) hist->xsec_hist->GetYaxis()->SetTitleOffset(1.6);
      }
      hist->xsec_hist->GetYaxis()->SetTitleSize(title_size);
      hist->xsec_hist->GetYaxis()->SetNdivisions(110);
      hist->xsec_hist->GetYaxis()->SetTickLength(0.015);
      hist->xsec_hist->GetYaxis()->SetRangeUser(0, max*1.1);
      hist->xsec_hist->GetYaxis()->SetMaxDigits(3);
    }

    else{
      double ytitle_offset = 0.8;
      if(!config->show_error_band) ytitle_offset = 0.9;
      // Set the titles
      hist->stacked_hist->GetYaxis()->SetTitle("Events (/bin width)");
      if(!config->show_error_band){
        hist->stacked_hist->GetXaxis()->SetTitle(titles->names[var_i]+" ["+titles->units[var_i]+"]");
        if(titles->units[var_i]=="") hist->stacked_hist->GetXaxis()->SetTitle(titles->names[var_i]);
      }
      // X axis config
      if(config->show_error_band) hist->stacked_hist->GetXaxis()->SetLabelOffset(0.1);
      hist->stacked_hist->GetXaxis()->SetTitleOffset(xtitle_offset);
      hist->stacked_hist->GetXaxis()->SetTickLength(xtick_length);
      if(!config->show_error_band) hist->stacked_hist->GetXaxis()->SetTitleSize(1.1 * 0.06);
      // Y axis config
      hist->stacked_hist->GetYaxis()->SetTitleOffset(ytitle_offset);
      double title_size = 1.1*0.06;
      hist->stacked_hist->GetYaxis()->SetTitleSize(title_size);
      hist->stacked_hist->GetYaxis()->SetNdivisions(110);
      hist->stacked_hist->GetYaxis()->SetTickLength(0.015);
      hist->stacked_hist->SetMaximum(max*1.1);
      hist->stacked_hist->GetYaxis()->SetMaxDigits(3);
    }
  }

  void DrawErrorBand(TH1D* error_band, int var_i, double size_ratio, TString systname){

    error_band->SetFillColor(38);
    error_band->SetLineColor(38);

    // Set axis titles
    TString err_title = "#sigma_{syst} (%)";
    if(config->show_stat_error && config->show_syst_error && systname=="total") err_title = "#sigma_{all} (%)";
    else if(config->show_stat_error && !config->show_syst_error) err_title = "#sigma_{stat} (%)";
    error_band->GetYaxis()->SetTitle(err_title);
    error_band->GetXaxis()->SetTitle(titles->names[var_i]+" ["+titles->units[var_i]+"]");
    if(titles->units[var_i]=="") error_band->GetXaxis()->SetTitle(titles->names[var_i]);

    // x axis config
    error_band->GetXaxis()->SetTitleSize(1.1 * size_ratio * 0.06);
    error_band->GetXaxis()->SetLabelSize(size_ratio * 0.06);
    error_band->GetXaxis()->SetLabelOffset(0.04);
    error_band->GetXaxis()->SetTickLength(size_ratio * 0.04);
    error_band->SetTitleOffset(1.0, "x");
    // y axis config
    error_band->GetYaxis()->SetTitleSize(1.1 * size_ratio * 0.06);
    error_band->GetYaxis()->SetLabelSize(size_ratio * 0.06);
    error_band->GetYaxis()->CenterTitle();
    error_band->GetYaxis()->SetTickLength(0.015);
    error_band->GetYaxis()->SetMaxDigits(3);
    error_band->SetNdivisions(105, "y");
    error_band->SetTitleOffset(0.32, "y");
    double emax = error_band->GetBinContent(error_band->GetMaximumBin());
    error_band->GetYaxis()->SetRangeUser(0, 1.1 * emax);

    // Draw the error bars
    if(error_band->GetNbinsX() < 50) error_band->Draw("B");
    else error_band->Draw("C");
  }

  // ------------------------------------------------------------------------------------------------------------------
  //                                EFFICIENCY AND RESPONSE PLOTTERS
  // ------------------------------------------------------------------------------------------------------------------

  // Make the efficiency and purity plots
  void PlotEffPur(Histo1D* histo, size_t var_i){

    // Efficiency: selected/total in true
    TString effname = TString(histo->total_hist->GetName())+"_efficiency";
    TString effxaxis = titles->names[var_i]+"^{true} ["+titles->units[var_i]+"]";
    PlotEfficiency(histo->efficiency.first, histo->efficiency.second, effname, effxaxis, "Efficiency [%/100]");

    // Purity: correct selected/total selected
    TString purname = TString(histo->total_hist->GetName())+"_purity";
    TString purxaxis = titles->names[var_i]+"^{reco} ["+titles->units[var_i]+"]";
    PlotEfficiency(histo->purity.first, histo->purity.second, purname, purxaxis, "Purity [%/100]");

  }

  // Draw efficiency/purity as function of some variable
  void PlotEfficiency(TH1D* select, TH1D* total, TString name, TString xaxis, TString yaxis){

    TCanvas *canvas = new TCanvas(name+yaxis, "");
    canvas->SetMargin(0.14, 0.04, 0.16, 0.1);

    TGraphAsymmErrors *graph = new TGraphAsymmErrors();

    graph->SetMarkerStyle(8);
    graph->SetMarkerSize(1);
    graph->SetLineWidth(3);
    graph->SetMarkerColor(46);
    graph->SetLineColor(46);
    graph->GetXaxis()->SetTitle(xaxis);
    graph->GetYaxis()->SetTitle(yaxis);
    // X axis config
    graph->GetXaxis()->SetTitleOffset(1.1);
    graph->GetXaxis()->SetTickLength(0.04);
    graph->GetXaxis()->SetTitleSize(1.1 * 0.06);
    // Y axis config
    graph->GetYaxis()->SetTitleOffset(0.95);
    graph->GetYaxis()->SetTickLength(0.015);
    graph->GetYaxis()->SetTitleSize(1.1 * 0.06);
    graph->GetYaxis()->SetNdivisions(108);
    graph->GetYaxis()->SetMaxDigits(3);

    graph->BayesDivide(select, total);
    graph->Draw("ap");
    graph->GetYaxis()->SetRangeUser(0, 1); 
    canvas->Modified();

    TString output_file = config->output_file;
    output_file.ReplaceAll(".","_"+name+".");
    canvas->SaveAs(output_file);

  }

  // Plot systematic matrices where appropriate
  void PlotResponse(TH2D* response){

    TString name = response->GetName();
    TCanvas *canvas = new TCanvas(name, "", 900, 900);
    bool show_text = false;
    if(response->GetNbinsX() <= 10 && response->GetNbinsY() <= 10) show_text = true;
    canvas->SetFrameLineWidth(4.);
    canvas->SetLineWidth(4.);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetMargin(0.16, 0.16, 0.16, 0.16);
    if(show_text) canvas->SetRightMargin(0.1);

    response->GetXaxis()->SetTitle("True bin i");
    response->GetYaxis()->SetTitle("Reco bin j");
    gStyle->SetPaintTextFormat("4.2f");
    response->SetMarkerSize(1.5);
    response->SetMarkerColor(10);
    response->GetXaxis()->SetNdivisions(110);
    response->GetXaxis()->SetTitleOffset(1.1);
    response->GetXaxis()->SetTickLength(0.04);
    response->GetXaxis()->SetTitleSize(1.1 * 0.06);
    response->GetXaxis()->CenterLabels();
    response->GetYaxis()->SetNdivisions(110);
    response->GetYaxis()->SetTitleOffset(1.15);
    response->GetYaxis()->SetTickLength(0.015);
    response->GetYaxis()->SetTitleSize(1.1 * 0.06);
    response->GetYaxis()->CenterLabels();
    // Text size
    response->GetXaxis()->SetTitleSize(0.05);
    response->GetYaxis()->SetTitleSize(0.05);
    response->GetXaxis()->SetLabelSize(0.04);
    response->GetYaxis()->SetLabelSize(0.04);
    response->GetZaxis()->SetLabelSize(0.035);

    if(show_text) response->Draw("COL TEXT");
    else response->Draw("COLZ");

    TString output_file = config->output_file;
    output_file.ReplaceAll(".","_"+name+".");
    canvas->SaveAs(output_file);
  }

  // ------------------------------------------------------------------------------------------------------------------
  //                                       SYSTEMATICS PLOTTERS
  // ------------------------------------------------------------------------------------------------------------------

  // Plot systematic matrices where appropriate
  void PlotAllSysts(Histo1D* histo){
    for(auto const& syst : config->systematics){
      PlotSyst(histo->systematics->GetSyst(syst));
    }
    PlotSyst(histo->systematics->total);
  }

  // Plot covariance, fractional covariance and correlation matrices
  void PlotSyst(Systematics* syst){
    Plot2D(syst->covariance, syst->covariance->GetName(), "Bin i", "Bin j");
    Plot2D(syst->frac_covariance, syst->frac_covariance->GetName(), "Bin i", "Bin j");
    Plot2D(syst->correlation, syst->correlation->GetName(), "Bin i", "Bin j");
  }

  // Plot systematic matrices where appropriate for 2D histograms
  void PlotAllSysts(Histo2D* histo){
    for(auto const& syst : config->systematics){
      PlotSyst(histo->systematics->GetSyst(syst));
    }
    PlotSyst(histo->systematics->total);
  }

  // Plot covariance, fractional covariance and correlation matrices for 2D histograms
  void PlotSyst(Systematics2D* syst){
    Plot2D(syst->covariance, syst->covariance->GetName(), "Bin i", "Bin j");
    Plot2D(syst->frac_covariance, syst->frac_covariance->GetName(), "Bin i", "Bin j");
    Plot2D(syst->correlation, syst->correlation->GetName(), "Bin i", "Bin j");
  }

  // Plot 1D universe variations where appropriate
  void PlotUniverses(Histo1D* histo, size_t var_i){
    for(auto const& syst : config->systematics){
      if(config->plot_xsec){
        if(syst=="genie") PlotUni(histo->xsec_hist, histo->systematics->genie->universes, var_i);
        if(syst=="flux") PlotUni(histo->xsec_hist, histo->systematics->flux->universes, var_i);
        if(syst=="detector") PlotUni(histo->xsec_hist, histo->systematics->detector->universes, var_i);
      }
      else{
        if(syst=="genie") PlotUni(histo->total_hist, histo->systematics->genie->universes, var_i);
        if(syst=="flux") PlotUni(histo->total_hist, histo->systematics->flux->universes, var_i);
        if(syst=="detector") PlotUni(histo->total_hist, histo->systematics->detector->universes, var_i);
      }
    }
  }

  // Plot a 1D hist with multiple universes drawn in background
  void PlotUni(TH1D* hist, std::vector<TH1D*> unis, size_t var_i){

    // Create the canvas
    TString name = unis[0]->GetName();
    name.ReplaceAll("_uni0","");
    TCanvas *canvas = new TCanvas("canvas"+name+"_unis", "canvas");
    canvas->SetMargin(0.15, 0.04, 0.15, 0.15);

    // Draw the stacked histogram and legend
    hist->SetFillStyle(0);
    hist->SetLineWidth(3);
    hist->SetLineColor(46);
    hist->Draw("HIST");

    double ymax = 0;
    for(size_t i = 0; i < unis.size(); i++){
      unis[i]->SetFillStyle(0);
      unis[i]->SetLineWidth(1);
      unis[i]->SetLineColor(44);
      unis[i]->Draw("HIST SAME");
      double max = unis[i]->GetBinContent(unis[i]->GetMaximumBin());
      if(max > ymax) ymax = max;
    }
    hist->Draw("HIST SAME");

    // Set the titles
    if(config->plot_xsec){
      hist->GetYaxis()->SetTitle(titles->GetXSecTitle(var_i));
    }
    else{
      hist->GetYaxis()->SetTitle("Events (/Bin width)");
    }
    hist->GetXaxis()->SetTitle(titles->names[var_i]+" ["+titles->units[var_i]+"]");
    if(titles->units[var_i]=="") hist->GetXaxis()->SetTitle(titles->names[var_i]);
    // X axis config
    hist->GetXaxis()->SetTitleOffset(1.);
    hist->GetXaxis()->SetLabelOffset(0.02);
    hist->GetXaxis()->SetTickLength(0.02);
    // Y axis config
    hist->GetYaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTickLength(0.015);
    double title_size = 1.1 * 0.06;
    if(config->plot_xsec && config->plot_variables.size()==1){ 
      title_size = 1.0 * 0.06;
      hist->GetYaxis()->SetTitleOffset(1.25);
    }
    if(config->plot_xsec && config->plot_variables.size()==2){ 
      title_size = 0.8 * 0.06;
      hist->GetYaxis()->SetTitleOffset(1.4);
    }
    hist->GetYaxis()->SetTitleSize(title_size);   
    hist->GetYaxis()->SetNdivisions(110);
    hist->GetYaxis()->SetRangeUser(0, ymax*1.1);
    hist->GetYaxis()->SetMaxDigits(3);
    canvas->Modified();

    // Text position and content
    double width = 0.65*(hist->GetXaxis()->GetXmax() 
                         - hist->GetXaxis()->GetXmin()) 
                        + hist->GetXaxis()->GetXmin();
    double height = hist->GetMaximum();
    double upper_text_size = 0.6*hist->GetYaxis()->GetTitleSize();
    if(config->show_info) DrawInfo(width, height, upper_text_size);
    
    TString output_file = config->output_file;
    output_file.ReplaceAll(".","_"+name+"_universes.");
    canvas->SaveAs(output_file);
  }

  
};

#endif
