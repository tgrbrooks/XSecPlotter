#ifndef PLOTTER_H
#define PLOTTER_H

#include "Configuration.h"
#include "Titles.h"
#include "Interaction.h"
#include "Systematics.h"
#include "Histo1D.h"
#include "Histo2D.h"

// Structure for holding interaction information
class Plotter
{
  public:

  Configuration *config;
  Titles *titles;

  Plotter(){}

  Plotter(Configuration *c, Titles *t)
  {
    config = c;
    titles = t;
  }

  // Draw additional information on to hist
  void DrawInfo(double width, double height, double size){

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

  // Plot a 1D stacked hist with statistical errors on the bottom
  void Plot1DWithErrors(Histo1D* histo, size_t var_i){

    // Create the canvas
    TString name = histo->stacked_hist->GetName();
    TCanvas *canvas = new TCanvas("canvas"+name,"canvas");

    // Split the pad for histogram and error plot
    double pad_split = .3;
    TPad *upper_pad = new TPad("upper_pad", "" , 0., pad_split, 1.0, 1.0);
    upper_pad->SetTopMargin(0.12);
    upper_pad->SetBottomMargin(0.075);
    upper_pad->SetLeftMargin(0.12);
    upper_pad->SetRightMargin(0.05);

    TPad *lower_pad = new TPad("lower_pad", "", 0., 0., 1., pad_split);
    lower_pad->SetTopMargin(0.01);
    lower_pad->SetBottomMargin(0.34);
    lower_pad->SetLeftMargin(0.12);
    lower_pad->SetRightMargin(0.05);

    upper_pad->Draw();
    lower_pad->Draw();

    // Fill the upper pad with histogram, info and legend
    upper_pad->cd();

    // Draw the stacked histogram and legend
    histo->stacked_hist->Draw("HIST");
    if(config->show_error_bars){
      TH1D* error_hist = (TH1D*)histo->total_hist->Clone();
      if(config->show_syst_error){
        for(size_t n = 0; n <= histo->systematics->total->GetNbinsX(); n++){
          error_hist->SetBinError(n, std::sqrt(std::pow(histo->systematics->total->GetBinError(n), 2)
                                            +std::pow(error_hist->GetBinError(n), 2)));
        }
      }
      error_hist->SetLineWidth(0);
      error_hist->SetMarkerStyle(0);
      error_hist->SetFillColor(15);
      error_hist->SetFillStyle(3001);
      if(config->show_syst_error) error_hist->Draw("E2 SAME");
      histo->total_hist->SetLineWidth(0);
      histo->total_hist->SetMarkerStyle(0);
      histo->total_hist->SetFillColor(12);
      histo->total_hist->SetFillStyle(3001);
      histo->total_hist->Draw("E2 SAME");
    }

    if(config->plot_stacked){
      histo->legend->SetNColumns(histo->legend->GetNRows());
      histo->legend->SetFillStyle(0);
      histo->legend->Draw();
    }
    // Set the titles
    if(config->plot_xsec){
      histo->stacked_hist->GetYaxis()->SetTitle(titles->GetXSecTitle(var_i));
    }
    else if(config->max_error > 0){
      histo->stacked_hist->GetYaxis()->SetTitle("Events (/bin width)");
    }
    else{
      histo->stacked_hist->GetYaxis()->SetTitle("Events");
    }
    // X axis config
    histo->stacked_hist->GetXaxis()->SetLabelOffset(0.1);
    histo->stacked_hist->GetXaxis()->SetTitleOffset(1.8);
    histo->stacked_hist->GetXaxis()->SetTickLength(0.04);
    // Y axis config
    histo->stacked_hist->GetYaxis()->SetTitleOffset(0.8);
    double title_size = 1.1*histo->stacked_hist->GetYaxis()->GetTitleSize();
    if(config->plot_xsec && config->plot_variables.size()==1){ 
      title_size = 1.0*histo->stacked_hist->GetYaxis()->GetTitleSize();
      histo->stacked_hist->GetYaxis()->SetTitleOffset(0.9);
    }
    if(config->plot_xsec && config->plot_variables.size()==2){ 
      title_size = 0.8*histo->stacked_hist->GetYaxis()->GetTitleSize();
      histo->stacked_hist->GetYaxis()->SetTitleOffset(1.0);
    }
    if(config->plot_xsec && config->plot_variables.size()==3){ 
      title_size = 0.6*histo->stacked_hist->GetYaxis()->GetTitleSize();
      histo->stacked_hist->GetYaxis()->SetTitleOffset(1.1);
    }
    histo->stacked_hist->GetYaxis()->SetTitleSize(title_size);
    histo->stacked_hist->GetYaxis()->SetNdivisions(110);
    histo->stacked_hist->GetYaxis()->SetTickLength(0.015);
    canvas->Modified();

    // Info text
    // Text position and content
    double width = 0.7*(histo->stacked_hist->GetXaxis()->GetXmax()-histo->stacked_hist->GetXaxis()->GetXmin())+histo->stacked_hist->GetXaxis()->GetXmin();
    double height = histo->stacked_hist->GetMaximum();
    double upper_text_size = 0.7*histo->stacked_hist->GetYaxis()->GetTitleSize();
    if(config->show_info) DrawInfo(width, height, upper_text_size);
   
    // Fill the lower pad with percentage error per bin
    lower_pad->cd();
    lower_pad->SetTickx();
    lower_pad->SetTicky();
   
    // Set axis titles
    histo->error_band->SetFillColor(38);
    histo->error_band->SetLineColor(38);
    histo->error_band->GetYaxis()->SetTitle("#sigma_{stat} (%)");
    histo->error_band->GetXaxis()->SetTitle(titles->names[var_i]+" ["+titles->units[var_i]+"]");

    double size_ratio = upper_pad->GetAbsHNDC()/lower_pad->GetAbsHNDC();
    // x axis config
    histo->error_band->GetXaxis()->SetTitleSize(1.1*size_ratio*histo->error_band->GetXaxis()->GetTitleSize());
    histo->error_band->GetXaxis()->SetLabelSize(size_ratio*histo->error_band->GetXaxis()->GetLabelSize());
    histo->error_band->GetXaxis()->SetLabelOffset(0.04);
    histo->error_band->GetXaxis()->SetTickLength(size_ratio*0.04);
    histo->error_band->SetTitleOffset(1.0, "x");
    // y axis config
    histo->error_band->GetYaxis()->SetTitleSize(1.1*size_ratio*histo->error_band->GetYaxis()->GetTitleSize());
    histo->error_band->GetYaxis()->SetLabelSize(size_ratio*histo->error_band->GetYaxis()->GetLabelSize());
    histo->error_band->GetYaxis()->CenterTitle();
    histo->error_band->GetYaxis()->SetTickLength(0.015);
    histo->error_band->SetNdivisions(105, "y");
    histo->error_band->SetTitleOffset(0.3, "y");

    // Draw the error bars
    if(histo->error_band->GetNbinsX() < 40) histo->error_band->Draw("B");
    else histo->error_band->Draw("C");
    
    TString output_file = config->output_file;
    output_file.ReplaceAll(".","_"+name+".");
    canvas->SaveAs(output_file);
  }

  // Plot a 1D stacked hist
  void Plot1D(Histo1D* histo, size_t var_i){

    // Create the canvas
    TString name = histo->stacked_hist->GetName();
    TCanvas *canvas = new TCanvas("canvas"+name,"canvas");

    // Split the pad for histogram and error plot
    canvas->SetTopMargin(0.15);
    canvas->SetBottomMargin(0.15);
    canvas->SetLeftMargin(0.12);
    canvas->SetRightMargin(0.04);

    // Draw the stacked histogram and legend
    histo->stacked_hist->Draw("HIST");

    TH1D* error_hist = (TH1D*)histo->total_hist->Clone();
    if(config->show_error_bars){
      if(config->show_syst_error){
        for(size_t n = 0; n <= error_hist->GetNbinsX(); n++){
          error_hist->SetBinError(n, std::sqrt(std::pow(histo->systematics->total->GetBinError(n), 2)
                                            +std::pow(error_hist->GetBinError(n), 2)));
        }
      }
      error_hist->SetLineWidth(0);
      error_hist->SetMarkerStyle(0);
      error_hist->SetFillColor(15);
      error_hist->SetFillStyle(3001);
      if(config->show_syst_error) error_hist->Draw("E2 SAME");
      histo->total_hist->SetLineWidth(0);
      histo->total_hist->SetMarkerStyle(0);
      histo->total_hist->SetFillColor(12);
      histo->total_hist->SetFillStyle(3001);
      histo->total_hist->Draw("E2 SAME");
    }

    if(config->plot_stacked){
      histo->legend->SetNColumns(histo->legend->GetNRows());
      histo->legend->SetFillStyle(0);
      histo->legend->Draw();
      canvas->Update();
      histo->legend->SetX1NDC(0.24);
      histo->legend->SetY1NDC(0.85);
      histo->legend->SetX2NDC(0.96);
      histo->legend->SetY2NDC(0.91);
      canvas->Modified();
    }
    // Set the titles
    if(config->plot_xsec){
      histo->stacked_hist->GetYaxis()->SetTitle(titles->GetXSecTitle(var_i));
    }
    else if(config->max_error > 0){
      histo->stacked_hist->GetYaxis()->SetTitle("Events (/Bin width)");
    }
    else{
      histo->stacked_hist->GetYaxis()->SetTitle("Events");
    }
    histo->stacked_hist->GetXaxis()->SetTitle(titles->names[var_i]+" ["+titles->units[var_i]+"]");
    // X axis config
    histo->stacked_hist->GetXaxis()->SetTitleOffset(1.);
    histo->stacked_hist->GetXaxis()->SetTickLength(0.02);
    histo->stacked_hist->GetXaxis()->SetTitleSize(1.1*histo->stacked_hist->GetXaxis()->GetTitleSize());
    // Y axis config
    histo->stacked_hist->GetYaxis()->SetTitleOffset(0.95);
    histo->stacked_hist->GetYaxis()->SetTickLength(0.015);
    double title_size = 1.1*histo->stacked_hist->GetYaxis()->GetTitleSize();
    if(config->plot_xsec && config->plot_variables.size()==1){ 
      title_size = 1.0*histo->stacked_hist->GetYaxis()->GetTitleSize();
      histo->stacked_hist->GetYaxis()->SetTitleOffset(1.05);
    }
    if(config->plot_xsec && config->plot_variables.size()==2){ 
      title_size = 0.8*histo->stacked_hist->GetYaxis()->GetTitleSize();
      histo->stacked_hist->GetYaxis()->SetTitleOffset(1.15);
    }
    if(config->plot_xsec && config->plot_variables.size()==3){ 
      title_size = 0.6*histo->stacked_hist->GetYaxis()->GetTitleSize();
      histo->stacked_hist->GetYaxis()->SetTitleOffset(1.25);
    }
      
    histo->stacked_hist->GetYaxis()->SetTitleSize(title_size);
    histo->stacked_hist->GetYaxis()->SetNdivisions(110);
    int binmax = histo->total_hist->GetMaximumBin();
    double ymax = (histo->total_hist->GetBinContent(binmax) + error_hist->GetBinError(binmax))*1.01;
    histo->stacked_hist->SetMaximum(ymax);
    canvas->Modified();

    // Text position and content
    double width = 0.65*(histo->stacked_hist->GetXaxis()->GetXmax()-histo->stacked_hist->GetXaxis()->GetXmin())+histo->stacked_hist->GetXaxis()->GetXmin();
    double height = histo->stacked_hist->GetMaximum();
    double upper_text_size = 0.6*histo->stacked_hist->GetYaxis()->GetTitleSize();
    if(config->show_info) DrawInfo(width, height, upper_text_size);
    
    TString output_file = config->output_file;
    output_file.ReplaceAll(".","_"+name+".");
    canvas->SaveAs(output_file);
  }

/*
  // Plot a 1D stacked hist
  void PlotMulti1D(std::vector<TH1D*> total_hist, std::vector<TH1D*> syst_hist, size_t i, size_t j = -1){

    // Create the canvas
    TString name = total_hist[0]->GetName();
    TCanvas *canvas = new TCanvas("canvas"+name,"canvas");
    TLegend *legend = new TLegend(0.14, 0., 0.94, 0.06);

    // Split the pad for histogram and error plot
    canvas->SetTopMargin(0.15);
    canvas->SetBottomMargin(0.15);
    canvas->SetLeftMargin(0.12);
    canvas->SetRightMargin(0.04);

    bool first = false;

    for(size_t tune_i = 0; tune_i < config->input_file.size(); tune_i++){

      total_hist[tune_i]->SetLineWidth(2);
      total_hist[tune_i]->SetMarkerStyle(1);
      total_hist[tune_i]->SetLineColor(config->cols[tune_i]);
      if(config->show_error_bars){
        if(first) total_hist[tune_i]->Draw("E2 HIST");
        else total_hist[tune_i]->Draw("E2 HIST SAME");
        for(size_t n = 1; n <= syst_hist[tune_i]->GetNbinsX(); n++){
          syst_hist[tune_i]->SetBinError(n, std::sqrt(std::pow(syst_hist[tune_i]->GetBinError(n), 2)
                                                      +std::pow(config->constant_syst*total_hist[tune_i]->GetBinContent(n)/100., 2)));
        }
        syst_hist[tune_i]->SetLineWidth(2);
        syst_hist[tune_i]->SetMarkerStyle(1);
        if(config->show_syst_error) syst_hist[tune_i]->Draw("E2 SAME");
      }
      else{
        if(first) total_hist[tune_i]->Draw("HIST");
        else total_hist[tune_i]->Draw("HIST SAME");
      }

      legend->AddEntry(total_hist[tune_i], config->tune_name[tune_i], "l");
      first = false;
    }

    legend->SetNColumns(legend->GetNRows());
    legend->SetFillStyle(0);
    legend->Draw();
    canvas->Update();
    legend->SetX1NDC(0.24);
    legend->SetY1NDC(0.85);
    legend->SetX2NDC(0.96);
    legend->SetY2NDC(0.91);
    canvas->Modified();

    // Set the titles
    if(config->plot_xsec){
      total_hist[0]->GetYaxis()->SetTitle(titles->GetXSecTitle(i, j));
    }
    else if(config->max_error > 0){
      total_hist[0]->GetYaxis()->SetTitle("Events (/Bin width)");
    }
    else{
      total_hist[0]->GetYaxis()->SetTitle("Events");
    }
    total_hist[0]->GetXaxis()->SetTitle(titles->names[i]+" ["+titles->units[i]+"]");
    // X axis config
    total_hist[0]->GetXaxis()->SetTitleOffset(1.);
    total_hist[0]->GetXaxis()->SetTickLength(0.02);
    total_hist[0]->GetXaxis()->SetTitleSize(1.1*total_hist[0]->GetXaxis()->GetTitleSize());
    // Y axis config
    total_hist[0]->GetYaxis()->SetTitleOffset(0.95);
    total_hist[0]->GetYaxis()->SetTickLength(0.015);
    double title_size = 1.1*total_hist[0]->GetYaxis()->GetTitleSize();
    if(config->plot_xsec && config->plot_variables.size()==1){ 
      title_size = 1.0*total_hist[0]->GetYaxis()->GetTitleSize();
      total_hist[0]->GetYaxis()->SetTitleOffset(1.05);
    }
    if(config->plot_xsec && config->plot_variables.size()==2){ 
      title_size = 0.8*total_hist[0]->GetYaxis()->GetTitleSize();
      total_hist[0]->GetYaxis()->SetTitleOffset(1.15);
    }
    if(config->plot_xsec && config->plot_variables.size()==3){ 
      title_size = 0.6*total_hist[0]->GetYaxis()->GetTitleSize();
      total_hist[0]->GetYaxis()->SetTitleOffset(1.25);
    }
      
    total_hist[0]->GetYaxis()->SetTitleSize(title_size);
    total_hist[0]->GetYaxis()->SetNdivisions(110);
    if(config->plot_xsec && config->plot_variables.size()==1)
    canvas->Modified();

    // Text position and content
    double width = 0.65*(total_hist[0]->GetXaxis()->GetXmax()-total_hist[0]->GetXaxis()->GetXmin())+total_hist[0]->GetXaxis()->GetXmin();
    double height = total_hist[0]->GetMaximum();
    double upper_text_size = 0.6*total_hist[0]->GetYaxis()->GetTitleSize();
    if(config->show_info) DrawInfo(width, height, upper_text_size);
    
    TString output_file = config->output_file;
    output_file.ReplaceAll(".","_"+name+".");
    canvas->SaveAs(output_file);
  }
*/

  
  // Make the efficiency and purity plots for up to 3 variables
  void PlotEffPur(Histo1D* histo, size_t var_i){

    // Efficiency: selected/total in true
    PlotEfficiency(histo->efficiency.first, histo->efficiency.second, config->plot_variables[var_i]+"_efficiency", titles->names[var_i]+"^{true} ["+titles->units[var_i]+"]", "Efficiency");
    // Purity: correct selected/total selected
    PlotEfficiency(histo->purity.first, histo->purity.second, config->plot_variables[var_i]+"_purity", titles->names[var_i]+"^{reco} ["+titles->units[var_i]+"]", "Purity");

  }

  // Draw efficiency/purity as function of some variable
  void PlotEfficiency(TH1D* select, TH1D* total, TString name, TString xaxis, TString yaxis){

    TCanvas *canvas = new TCanvas(name+yaxis, "");
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.16);
    canvas->SetLeftMargin(0.14);
    canvas->SetRightMargin(0.04);

    TGraphAsymmErrors *graph = new TGraphAsymmErrors();

    graph->SetMarkerColor(46);
    graph->SetLineColor(46);
    graph->GetXaxis()->SetTitle(xaxis);
    graph->GetYaxis()->SetTitle(yaxis);
    // X axis config
    graph->GetXaxis()->SetTitleOffset(1.1);
    graph->GetXaxis()->SetTickLength(0.04);
    graph->GetXaxis()->SetTitleSize(1.1*graph->GetXaxis()->GetTitleSize());
    // Y axis config
    graph->GetYaxis()->SetTitleOffset(0.95);
    graph->GetYaxis()->SetTickLength(0.015);
    graph->GetYaxis()->SetTitleSize(1.1*graph->GetYaxis()->GetTitleSize());
    graph->GetYaxis()->SetNdivisions(108);

    graph->BayesDivide(select, total);
    graph->Draw("ap");
    graph->GetYaxis()->SetRangeUser(0, 1); 
    canvas->Modified();

    TString output_file = config->output_file;
    output_file.ReplaceAll(".","_"+name+".");
    canvas->SaveAs(output_file);

  }

  void PlotAllSysts(Histo1D* histo){
    for(auto const& syst : config->systematics){
      if(syst=="genie") PlotSyst(histo->systematics->genie);
      if(syst=="flux") PlotSyst(histo->systematics->flux);
      if(syst=="detector") PlotSyst(histo->systematics->detector);
    }
  }

  void PlotSyst(Systematics* syst){
    Plot2D(syst->covariance, syst->covariance->GetName(), "Bin i", "Bin j");
    Plot2D(syst->frac_covariance, syst->frac_covariance->GetName(), "Bin i", "Bin j");
    Plot2D(syst->correlation, syst->correlation->GetName(), "Bin i", "Bin j");
  }

  // Plot a 2D histogram
  void Plot2D(TH2D* hist, TString name, TString xaxis, TString yaxis){

    TCanvas *canvas = new TCanvas(name, "", 900, 600);
    canvas->SetFrameLineWidth(4.);
    canvas->SetLineWidth(4.);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetBottomMargin(0.16);
    canvas->SetLeftMargin(0.16);
    canvas->SetRightMargin(0.16);

    hist->GetXaxis()->SetTitle(xaxis);
    hist->GetYaxis()->SetTitle(yaxis);
    // X axis config
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetXaxis()->SetTickLength(0.04);
    hist->GetXaxis()->SetTitleSize(1.1*hist->GetXaxis()->GetTitleSize());
    hist->GetXaxis()->SetNdivisions(108);
    // Y axis config
    hist->GetYaxis()->SetTitleOffset(0.95);
    hist->GetYaxis()->SetTickLength(0.015);
    hist->GetYaxis()->SetTitleSize(1.1*hist->GetYaxis()->GetTitleSize());
    hist->GetYaxis()->SetNdivisions(108);

    hist->Draw("colz");

    TString output_file = config->output_file;
    output_file.ReplaceAll(".","_"+name+".");
    canvas->SaveAs(output_file);
  }


  
};

#endif
