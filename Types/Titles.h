#ifndef TITLES_H
#define TITLES_H

#include "Configuration.h"

// Structure for holding plot titles
class Titles
{
  public:

  std::vector<TString> hist_titles;
  std::vector<TString> names;
  std::vector<TString> units;
  TString data_type;
  TString part_cont;
  TString lep_cont;
  TString is_cc;
  TString n_pr;
  TString n_pipm;
  TString n_pi0;
  TString int_type;
  TString pot;
  TString mass;

  Titles(std::vector<TString> ht, std::vector<TString> n, std::vector<TString> u, TString dt, TString pc, TString lc, 
        TString ic, TString npr, TString npi, TString npi0, TString it, TString p, TString m)
  {
    hist_titles = ht;
    names = n;
    units = u;
    data_type = dt;
    part_cont = pc;
    lep_cont = lc;
    is_cc = ic;
    n_pr = npr;
    n_pipm = npi;
    n_pi0 = npi0;
    int_type = it;
    pot = p;
    mass = m;
  }

  Titles(){}

  Titles(Configuration *config)
  {
    // Set the units and histogram titles based on plotting variable
    int index = 0;
    for(auto const& var : config->plot_variables){
      bool apply_cos = false;
      TString plot_var = var;
      if(var(0, 4) == "cos_"){ 
        apply_cos = true;
        plot_var = var(4, plot_var.Length());
      }

      if (plot_var == "nu_energy"){
        names.push_back("E_{#nu}");
        units.push_back("GeV");
        hist_titles.push_back("neutrino energy");
      } 
      else if (plot_var == "lep_mom"){
        names.push_back("P_{lep}");
        units.push_back("GeV");
        hist_titles.push_back("lepton momentum");
      }
      else if (plot_var == "pr1_mom"){
        names.push_back("P_{p}");
        units.push_back("GeV");
        hist_titles.push_back("leading proton momentum");
      }
      else if (plot_var == "pipm1_mom"){
        names.push_back("P_{#pi}");
        units.push_back("GeV");
        hist_titles.push_back("leading #pi^{#pm} momentum");  
      }
      else if (plot_var == "lep_theta"){
        names.push_back("#theta_{lep}");
        units.push_back("rad");
        hist_titles.push_back("lepton #theta");
      }
      else if (plot_var == "pr1_theta"){
        names.push_back("#theta_{p}");
        units.push_back("rad");
        hist_titles.push_back("leading proton #theta");
      }
      else if (plot_var == "pipm1_theta"){
        names.push_back("#theta_{#pi}");
        units.push_back("rad");
        hist_titles.push_back("leading #pi^{#pm} #theta");
      }
      else if (plot_var == "lep_pr1_angle"){
        names.push_back("#theta_{lep,p}");
        units.push_back("rad");
        hist_titles.push_back("angle between lepton and proton");
      }
      else if (plot_var == "lep_pipm1_angle"){
        names.push_back("#theta_{lep,#pi}");
        units.push_back("rad");
        hist_titles.push_back("angle between lepton and #pi^{#pm}");
      }
      else if (plot_var == "nu_pdg"){
        names.push_back("PDG code");
        units.push_back("");
        hist_titles.push_back("neutrino PDG");
      }
      else if (plot_var == "lep_contained"){
        names.push_back("lepton contained?");
        units.push_back("");
        hist_titles.push_back("lepton containment");
      }
      else if (plot_var == "particles_contained"){
        names.push_back("particles contained?");
        units.push_back("");
        hist_titles.push_back("secondary particle containment");
      }
      else if (plot_var == "cc"){
        names.push_back("Is CC?");
        units.push_back("");
        hist_titles.push_back("charged or neutral current");
      }
      else if (plot_var == "int_type"){
        names.push_back("interaction code");
        units.push_back("");
        hist_titles.push_back("QE(0) RES(1) DIS(2) COH(3) MEC(10)");
      }
      else if (plot_var == "n_pr"){
        names.push_back("N_{p}");
        units.push_back("");
        hist_titles.push_back("number of protons");
      }
      else if (plot_var == "n_pipm"){
        names.push_back("N_{#pi}");
        units.push_back("");
        hist_titles.push_back("number of #pi^{#pm}");
      }
      else if (plot_var == "n_pi0"){
        names.push_back("N_{#pi^{0}}");
        units.push_back("");
        hist_titles.push_back("number of #pi^{0}");
      }
      else if (plot_var == "delta_pt"){
        names.push_back("#delta p_{T}");
        units.push_back("GeV");
        hist_titles.push_back("#delta p_{T}");
      }
      else if (plot_var == "delta_alphat"){
        names.push_back("#delta #alpha_{T}");
        units.push_back("deg");
        hist_titles.push_back("#delta #alpha_{T}");
      }
      else if (plot_var == "delta_phit"){
        names.push_back("#delta #phi_{T}");
        units.push_back("deg");
        hist_titles.push_back("#delta #phi_{T}");
      }
      else{
        std::cout<<"Invalid plotting variable!\n";
        exit(1);
      }

      // TODO remove units
      if(apply_cos){
        names[index] = "cos"+ names[index];
        units[index] = "";
        hist_titles[index] = "cos "+hist_titles[index];
      }
      index++;
    }

    // True or reco
    if (config->stage == "reco")
      data_type = "Reconstructed";
    else if (config->stage == "true")
      data_type = "Truth";
    else if (config->stage == "smeareff")
      data_type = "Smearing + Efficiency";
    else{
      std::cout<<"Unrecognised stage!\n";
      exit(1);
    }

    // Particle containment
    if (config->contained_particles.size() == 1){
      if (config->contained_particles[0])
        part_cont = "Contained Particles";
      else
        part_cont = "Exiting Particles";
    }
    else
      part_cont = "Cont+Exit Particles";

    // Lepton containment
    if (config->contained_lepton.size() == 1){
      if (config->contained_lepton[0])
        lep_cont = "Contained Lepton";
      else
        lep_cont = "Exiting Lepton";
    }
    else
      lep_cont = "Cont+Exit Lepton";
    
    // Neutrino pdg
    TString nu_type;
    for(auto const& pdg : config->nu_pdg){
      if(pdg == 12) nu_type += "#nu_{e} ";
      if(pdg == 14) nu_type += "#nu_{#mu} ";
      if(pdg == -12) nu_type += "#bar{#nu}_{e} ";
      if(pdg == -14) nu_type += "#bar{#nu}_{#mu} ";
    }

    // Charged or neutral current
    if (config->is_cc.size() == 1){
      if (config->is_cc[0])
        is_cc = "CC";
      else
        is_cc = "NC";
    }
    else
      is_cc = "CC+NC";

    is_cc = nu_type + is_cc;

    // Number of protons
    if(std::find(config->num_protons.begin(), config->num_protons.end(), -1) != config->num_protons.end())
      n_pr = "All ";
    else{
      for(auto const& num : config->num_protons)
        n_pr += std::to_string(num) + " ";
    }
    n_pr += "Proton(s)";

    // Number of charged pions
    if(std::find(config->num_pipm.begin(), config->num_pipm.end(), -1) != config->num_pipm.end())
      n_pipm = "All ";
    else{
      for(auto const& num : config->num_pipm)
        n_pipm += std::to_string(num) + " ";
    }
    n_pipm += "Charged Pion(s)";

    // Number of neutral pions
    if(std::find(config->num_pi0.begin(), config->num_pi0.end(), -1) != config->num_pi0.end())
      n_pi0 = "All ";
    else{
      for(auto const& num : config->num_pi0)
        n_pi0 += std::to_string(num) + " ";
    }
    n_pi0 += "Neutral Pion(s)";

    // Interaction type
    if(std::find(config->interaction_type.begin(), config->interaction_type.end(), -1) != config->interaction_type.end())
      int_type = "All Interactions";
    else{
      for(auto const& i_type : config->interaction_type){
        if(i_type == 0) int_type += "QE ";
        if(i_type == 1) int_type += "RES ";
        if(i_type == 2) int_type += "DIS ";
        if(i_type == 3) int_type += "COH ";
        if(i_type == 10) int_type += "MEC ";
      }
    }

    // POT
    std::stringstream pot_stream;
    double pot_d = config->pot_scale;
    if(config->pot_scale <= 0) pot_d = config->pot[0];
    pot_stream << std::setprecision(3) << "POT = " << pot_d << "}";
    std::string pot_string = pot_stream.str();
    pot_string.replace(pot_string.find("e"), 1, "#times10^{");
    pot_string.replace(pot_string.find("+"), 1, "");
    pot = TString(pot_string);

    // Fiducial mass
    std::stringstream mass_stream;
    mass_stream << std::setprecision(3) << "Fid Mass = " << config->fiducial_mass << " t";
    std::string mass_string = mass_stream.str();
    mass = TString(mass_string);

  }

  // Get the Y axis title when plotting cross sections
  TString GetXSecTitle(int i, int j = -1){

    TString xsec_title = "d#sigma/d"+names[i]+" [10^{-38}#frac{cm^{2}}{"+units[i]+" n}]";

    if(j != -1){
      xsec_title = "d^{2}#sigma/d"+names[i]+"d"+names[j]+" [10^{-38}#frac{cm^{2}}{"+units[i]+" "+units[j]+" n}]";
    }
    return TString(xsec_title);
  }

};

#endif
