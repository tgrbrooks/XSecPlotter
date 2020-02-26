#ifndef FLUXMANAGER_H
#define FLUXMANAGER_H

#include "Configuration.h"

// Structure for holding interaction information
class FluxManager
{
  public:

  Configuration *config; // Global configurations

  double pot;            // Total simulated POT in flux file

  double numu_flux;      // Muon neutrino flux
  double antinumu_flux;  // Muon anti-neutrino flux
  double nue_flux;       // Electron neutrino flux
  double antinue_flux;   // Electron anti-neutrino flux

  std::vector<double> numu_uni;     // Reweighted muon neutrino flux
  std::vector<double> antinumu_uni; // Reweighted muon anti-neutrino flux
  std::vector<double> nue_uni;      // Reweighted electron neutrino flux
  std::vector<double> antinue_uni;  // Reweighted electron anti-neutrino flux

  // Constructor
  FluxManager(Configuration *c)
  {
    config = c;

    TString fname = "Flux/fluxreweight.root";

    // Get the simulated POT
    pot = GetPot(fname, "FluxReweight/metadata");

    double scale = 1e6/(pot*4*4);
 
    int nsims = config->reweight_nuni;
 
    // PLOTS
    TFile *flux_file = new TFile(fname, "READ");

    // Simulated histograms
    TH1D* hNuMu = (TH1D*) flux_file->Get("FluxReweight/hNuMu");
    TH1D* hAntiNuMu = (TH1D*) flux_file->Get("FluxReweight/hAntiNuMu");
    TH1D* hNuE = (TH1D*) flux_file->Get("FluxReweight/hNuE");
    TH1D* hAntiNuE = (TH1D*) flux_file->Get("FluxReweight/hAntiNuE");

    // Get the integrated fluxes
    numu_flux = hNuMu->Integral(1, hNuMu->GetNbinsX()+1);
    antinumu_flux = hAntiNuMu->Integral(1, hAntiNuMu->GetNbinsX()+1);
    nue_flux = hNuE->Integral(1, hNuE->GetNbinsX()+1);
    antinue_flux = hAntiNuE->Integral(1, hAntiNuE->GetNbinsX()+1);

    // Get the universe variations
    for(int i = 0; i < nsims; i++){
      TH1D* hNuMu_tmp = (TH1D*) flux_file->Get(Form("FluxReweight/hNuMu_rw%i", i));
      TH1D* hAntiNuMu_tmp = (TH1D*) flux_file->Get(Form("FluxReweight/hAntiNuMu_rw%i", i));
      TH1D* hNuE_tmp = (TH1D*) flux_file->Get(Form("FluxReweight/hNuE_rw%i", i));
      TH1D* hAntiNuE_tmp = (TH1D*) flux_file->Get(Form("FluxReweight/hAntiNuE_rw%i", i));

      double numu_integral = hNuMu_tmp->Integral(1, hNuMu_tmp->GetNbinsX()+1);
      double antinumu_integral = hAntiNuMu_tmp->Integral(1, hAntiNuMu_tmp->GetNbinsX()+1);
      double nue_integral = hAntiNuE_tmp->Integral(1, hNuE_tmp->GetNbinsX()+1);
      double antinue_integral = hAntiNuE_tmp->Integral(1, hAntiNuE_tmp->GetNbinsX()+1);

      numu_uni.push_back(numu_integral);
      antinumu_uni.push_back(antinumu_integral);
      nue_uni.push_back(nue_integral);
      antinue_uni.push_back(antinue_integral);
    }

  }

  // Get the integrated flux
  double IntegratedFlux(){
    double int_flux = 0;
    for(auto const& pdg : config->nu_pdg){
      if(pdg==14) int_flux += numu_flux;
      if(pdg==-14) int_flux += antinumu_flux;
      if(pdg==12) int_flux += nue_flux;
      if(pdg==-12) int_flux += antinue_flux;
    }
    // Scale to correct POT
    double pot_scale = config->pot[0];
    if(config->pot_scale > 0) pot_scale = config->pot_scale;
    // Return in units of cm^-2
    return int_flux * pot_scale / (pot * 400 * 400);
  }

  // Get the integrated flux from universe variation
  double IntegratedFlux(int uni){
    double int_flux = 0;
    if(uni >= (int)numu_uni.size()) return int_flux;
    for(auto const& pdg : config->nu_pdg){
      if(pdg==14) int_flux += numu_uni[uni];
      if(pdg==-14) int_flux += antinumu_uni[uni];
      if(pdg==12) int_flux += nue_uni[uni];
      if(pdg==-12) int_flux += antinue_uni[uni];
    }
    // Scale to correct POT
    double pot_scale = config->pot[0];
    if(config->pot_scale > 0) pot_scale = config->pot_scale;
    // Return in units of cm^-2
    return int_flux * pot_scale / (pot * 400 * 400);
  }

  // Get the POT from the flux file
  double GetPot(TString name, TString path){

    TFile f(name);

    // Get the simulated POT
    double total_pot = 0;
    TTreeReader potReader(path, &f);
    TTreeReaderValue<double> fpot (potReader, "pot");
    while(potReader.Next()) total_pot += 1e7; //*fpot;

    return total_pot;
  }

};


#endif
