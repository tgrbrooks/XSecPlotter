#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include "Configuration.h"
#include "Interaction.h"
#include "Selection.h"

// Structure for holding interaction information
class DataManager
{
  public:

  Configuration *config;

  size_t file_i;

  std::vector<Interaction> interactions;
  std::map<std::string, std::vector<std::vector<double>>> stack_data;
  std::vector<std::vector<double>> total_data;
  std::vector<bool> data_used;

  DataManager(){}

  DataManager(Configuration *c, size_t f)
  {
    config = c;
    file_i = f;

    // Read in all of the data
    interactions = ReadData(file_i);

    // Loop over all of the interactions and sort by type
    for(auto const& in : interactions){
      if(!in.selected){data_used.push_back(false); continue;}

      total_data.push_back(in.variables);

      if(!config->plot_stacked){
        stack_data["all"].push_back(in.variables);
      }
      else if(config->stack_by == "fsi"){
        stack_data[in.fsi].push_back(in.variables);
      }
      else if(config->stack_by == "int"){
        stack_data[in.int_type].push_back(in.variables);
      }
      else if(config->stack_by == "nu"){
        stack_data[in.nu_type].push_back(in.variables);
      }

      data_used.push_back(true);
    }
  }

  // Read in true variables
  std::vector<Interaction> ReadData(int i){

    // Open the root tree
    TFile data_file(config->input_file[i]);
    if(!data_file.IsOpen()){
      std::cout<<"Could not read input file!\n";
      exit(1);
    }

    TString prefix = config->stage+"_";

    Selection sel(config);

    //Read in TTree
    TTreeReader tree_reader(config->tree_path, &data_file);

    // True vertex
    TTreeReaderValue<double>       vtx_x(tree_reader, "vtx_x");
    TTreeReaderValue<double>       vtx_y(tree_reader, "vtx_y");
    TTreeReaderValue<double>       vtx_z(tree_reader, "vtx_z");

    // True quantities for stacked hists
    TTreeReaderValue<bool>         true_lep_contained(tree_reader, "true_lep_contained");
    TTreeReaderValue<bool>         true_particles_contained(tree_reader, "true_particles_contained");
    TTreeReaderValue<bool>         true_cc(tree_reader, "true_cc");
    TTreeReaderValue<int>          true_nu_pdg(tree_reader, "true_nu_pdg");
    TTreeReaderValue<int>          true_int_type(tree_reader, "true_int_type");
    TTreeReaderValue<unsigned int> true_n_pipm(tree_reader, "true_n_pipm");
    TTreeReaderValue<unsigned int> true_n_pi0(tree_reader, "true_n_pi0");
    TTreeReaderValue<unsigned int> true_n_pr(tree_reader, "true_n_pr");

    // Need true values for efficiency/purity/response
    TTreeReaderValue<double>       true_nu_energy(tree_reader, "true_nu_energy");
    TTreeReaderValue<double>       true_lep_mom(tree_reader, "true_lep_mom");
    TTreeReaderValue<double>       true_lep_theta(tree_reader, "true_lep_theta");
    TTreeReaderValue<double>       true_pr1_mom(tree_reader, "true_pr1_mom");
    TTreeReaderValue<double>       true_pr1_theta(tree_reader, "true_pr1_theta");
    TTreeReaderValue<double>       true_lep_pr1_angle(tree_reader, "true_lep_pr1_angle");
    TTreeReaderValue<double>       true_pipm1_mom(tree_reader, "true_pipm1_mom");
    TTreeReaderValue<double>       true_pipm1_theta(tree_reader, "true_pipm1_theta");
    TTreeReaderValue<double>       true_lep_pipm1_angle(tree_reader, "true_lep_pipm1_angle");
    TTreeReaderValue<double>       true_delta_pt(tree_reader, "true_delta_pt");
    TTreeReaderValue<double>       true_delta_alphat(tree_reader, "true_delta_alphat");
    TTreeReaderValue<double>       true_delta_phit(tree_reader, "true_delta_phit");

    //Associate TTree values with variables
    TTreeReaderValue<bool>         lep_contained(tree_reader, prefix+"lep_contained");
    TTreeReaderValue<bool>         particles_contained(tree_reader, prefix+"particles_contained");
    TTreeReaderValue<bool>         cc(tree_reader, prefix+"cc");
    TTreeReaderValue<int>          nu_pdg(tree_reader, prefix+"nu_pdg");
    TTreeReaderValue<int>          int_type(tree_reader, prefix+"int_type");
    TTreeReaderValue<unsigned int> n_pipm(tree_reader, prefix+"n_pipm");
    TTreeReaderValue<unsigned int> n_pi0(tree_reader, prefix+"n_pi0");
    TTreeReaderValue<unsigned int> n_pr(tree_reader, prefix+"n_pr");

    TTreeReaderValue<double>       nu_energy(tree_reader, prefix+"nu_energy");
    TTreeReaderValue<double>       lep_mom(tree_reader, prefix+"lep_mom");
    TTreeReaderValue<double>       lep_theta(tree_reader, prefix+"lep_theta");
    TTreeReaderValue<double>       pr1_mom(tree_reader, prefix+"pr1_mom");
    TTreeReaderValue<double>       pr1_theta(tree_reader, prefix+"pr1_theta");
    TTreeReaderValue<double>       lep_pr1_angle(tree_reader, prefix+"lep_pr1_angle");
    TTreeReaderValue<double>       pipm1_mom(tree_reader, prefix+"pipm1_mom");
    TTreeReaderValue<double>       pipm1_theta(tree_reader, prefix+"pipm1_theta");
    TTreeReaderValue<double>       lep_pipm1_angle(tree_reader, prefix+"lep_pipm1_angle");
    TTreeReaderValue<double>       delta_pt(tree_reader, prefix+"delta_pt");
    TTreeReaderValue<double>       delta_alphat(tree_reader, prefix+"delta_alphat");
    TTreeReaderValue<double>       delta_phit(tree_reader, prefix+"delta_phit");

    // Vector to fill
    std::vector<Interaction> interactions;
   
    // Loop over all the interactions
    while (tree_reader.Next()) {
     
      bool selected = sel.IsSelected(*nu_pdg, *cc, *lep_contained, *particles_contained, 
                                 *n_pr, *n_pipm, *n_pi0, *int_type);
      bool true_selected = sel.IsSelected(*true_nu_pdg, *true_cc, *true_lep_contained, 
                                      *true_particles_contained, *true_n_pr, *true_n_pipm, *true_n_pi0, *true_int_type);

      // Check true vertex inside fiducial volume
      if(!sel.InFiducial(*vtx_x, *vtx_y, *vtx_z)){
        selected = false;
        true_selected = false;
      }

      std::vector<double> variables;
      std::vector<double> true_variables;
      int index = 0;
      for(auto const& var : config->plot_variables){
        bool apply_cos = false;
        TString plot_var = var;
        if(var(0, 4) == "cos_"){ 
          apply_cos = true;
          plot_var = var(4, plot_var.Length());
        }

        if (plot_var == "lep_contained"){ 
          variables.push_back((double)(*lep_contained));
          true_variables.push_back((double)(*true_lep_contained));
        }
        if (plot_var == "particles_contained"){ 
          variables.push_back((double)(*particles_contained));
          true_variables.push_back((double)(*true_particles_contained));
        }
        if (plot_var == "cc"){ 
          variables.push_back((double)(*cc));
          true_variables.push_back((double)(*true_cc));
        }
        if (plot_var == "nu_pdg"){ 
          variables.push_back((double)(*nu_pdg));
          true_variables.push_back((double)(*true_nu_pdg));
        }
        if (plot_var == "int_type"){ 
          variables.push_back((double)(*int_type));
          true_variables.push_back((double)(*true_int_type));
        }
        if (plot_var == "n_pr"){ 
          variables.push_back((double)(*n_pr));
          true_variables.push_back((double)(*true_n_pr));
        }
        if (plot_var == "n_pipm"){ 
          variables.push_back((double)(*n_pipm));
          true_variables.push_back((double)(*true_n_pipm));
        }
        if (plot_var == "n_pi0"){ 
          variables.push_back((double)(*n_pi0));
          true_variables.push_back((double)(*true_n_pi0));
        }
        if (plot_var == "nu_energy"){ 
          variables.push_back(*nu_energy);
          true_variables.push_back(*true_nu_energy);
        }
        if (plot_var == "lep_mom"){ 
          variables.push_back(*lep_mom);
          true_variables.push_back(*true_lep_mom);
        }
        if (plot_var == "lep_theta"){ 
          variables.push_back(*lep_theta);
          true_variables.push_back(*true_lep_theta);
        }
        if (plot_var == "pr1_mom"){ 
          variables.push_back(*pr1_mom);
          true_variables.push_back(*true_pr1_mom);
        }
        if (plot_var == "pr1_theta"){ 
          variables.push_back(*pr1_theta);
          true_variables.push_back(*true_pr1_theta);
        }
        if (plot_var == "lep_pr1_angle"){ 
          variables.push_back(*lep_pr1_angle);
          true_variables.push_back(*true_lep_pr1_angle);
        }
        if (plot_var == "pipm1_mom"){ 
          variables.push_back(*pipm1_mom);
          true_variables.push_back(*true_pipm1_mom);
        }
        if (plot_var == "pipm1_theta"){ 
          variables.push_back(*pipm1_theta);
          true_variables.push_back(*true_pipm1_theta);
        }
        if (plot_var == "lep_pipm1_angle"){ 
          variables.push_back(*lep_pipm1_angle);
          true_variables.push_back(*true_lep_pipm1_angle);
        }
        if (plot_var == "delta_pt"){ 
          variables.push_back(*delta_pt);
          true_variables.push_back(*true_delta_pt);
        }
        if (plot_var == "delta_alphat"){ 
          variables.push_back(*delta_alphat);
          true_variables.push_back(*true_delta_alphat);
        }
        if (plot_var == "delta_phit"){ 
          variables.push_back(*delta_phit);
          true_variables.push_back(*true_delta_phit);
        }

        if(apply_cos){ 
          variables[index] = cos(variables[index]);
          true_variables[index] = cos(true_variables[index]);
        }
        index++;
      }

      // FSI: 0pi0p, 0pi1p, 0pi2+p, 1pi, 2+pi, 1+pi0
      std::string fsi_string = "other";
      if(*true_n_pi0 >= 1) fsi_string = "#geq1#pi^{0}";
      else if(*true_n_pipm == 1) fsi_string = "1#pi^{#pm}";
      else if(*true_n_pipm >= 2) fsi_string = "#geq2#pi^{#pm}";
      else if(*true_n_pr == 0) fsi_string = "0#pi0p";
      else if(*true_n_pr == 1) fsi_string = "0#pi1p";
      else if(*true_n_pr >= 2) fsi_string = "0#pi#geq2p";

      // Int: QE, RES, COH, DIS, MEC
      std::string int_string = "other";
      if(*true_int_type == 0) int_string = "QE";
      else if(*true_int_type == 1) int_string = "RES";
      else if(*true_int_type == 2) int_string = "DIS";
      else if(*true_int_type == 3) int_string = "COH";
      else if(*true_int_type == 10) int_string = "MEC";

      std::string nu_string = "other";
      if(*true_nu_pdg == -12 && *true_cc) nu_string = "#bar{#nu}_{e} CC";
      if(*true_nu_pdg == -12 && !*true_cc) nu_string = "#bar{#nu}_{e} NC";
      if(*true_nu_pdg == -14 && *true_cc) nu_string = "#bar{#nu}_{#mu} CC";
      if(*true_nu_pdg == -14 && !*true_cc) nu_string = "#bar{#nu}_{#mu} NC";
      if(*true_nu_pdg == 12 && *true_cc) nu_string = "#nu_{e} CC";
      if(*true_nu_pdg == 12 && !*true_cc) nu_string = "#nu_{e} NC";
      if(*true_nu_pdg == 14 && *true_cc) nu_string = "#nu_{#mu} CC";
      if(*true_nu_pdg == 14 && !*true_cc) nu_string = "#nu_{#mu} NC";

      if(std::find(true_variables.begin(), true_variables.end(), -99999) != true_variables.end()){ 
        selected = false;
        true_selected = false;
      }
      
      Interaction interaction(selected, true_selected, fsi_string, int_string, nu_string, variables, true_variables);
      interactions.push_back(interaction);
      
    }

    // Throw exception if nothing to plot
    if (interactions.size() == 0){
      std::cout << "No events match your input parameters." << std::endl;
      throw std::exception();
    }

    return interactions;
   
  }

};

#endif
