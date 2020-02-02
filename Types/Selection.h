#ifndef SELECTION_H
#define SELECTION_H

#include "Configuration.h"

// Structure for holding interaction information
class Selection
{
  public:

  Configuration *config;

  Selection(Configuration *c)
  {
    config = c;
  }

  // Check if an interaction matches the selection criteria
  bool IsSelected(int nu_pdg, bool cc, bool lep_contained, bool particles_contained, int n_pr, int n_pipm, int n_pi0, int int_type){
    bool selected = true;

    if (std::find(config->nu_pdg.begin(), config->nu_pdg.end(), nu_pdg) == config->nu_pdg.end())
      selected = false;
    if (std::find(config->is_cc.begin(), config->is_cc.end(), cc) == config->is_cc.end())
      selected = false;
    if (std::find(config->contained_lepton.begin(), config->contained_lepton.end(), lep_contained) == config->contained_lepton.end())
      selected = false;
    if (std::find(config->contained_particles.begin(), config->contained_particles.end(), particles_contained) == config->contained_particles.end())
      selected = false;
    if (config->plot_by_fsi){
      if (std::find(config->num_protons.begin(), config->num_protons.end(), -1) == config->num_protons.end()){
        if (std::find(config->num_protons.begin(), config->num_protons.end(), n_pr) == config->num_protons.end())
          selected = false;
      }
      if (std::find(config->num_pipm.begin(), config->num_pipm.end(), -1) == config->num_pipm.end()){
        if (std::find(config->num_pipm.begin(), config->num_pipm.end(), n_pipm) == config->num_pipm.end())
          selected = false;
      }
      if (std::find(config->num_pi0.begin(), config->num_pi0.end(), -1) == config->num_pi0.end()){
        if (std::find(config->num_pi0.begin(), config->num_pi0.end(), n_pi0) == config->num_pi0.end())
          selected = false;
      }
    }
    else{
      if (std::find(config->interaction_type.begin(), config->interaction_type.end(), -1) == config->interaction_type.end()){
        if (std::find(config->interaction_type.begin(), config->interaction_type.end(), int_type) == config->interaction_type.end())
          selected = false;
      }
    }

    return selected;
  }

  // Check true vertex inside fiducial volume
  bool InFiducial(double vtx_x, double vtx_y, double vtx_z){
    bool in_fiducial = true;
    if(std::find(config->fiducial.begin(), config->fiducial.end(), -1) == config->fiducial.end() && config->fiducial.size() == 6){
      if(vtx_x < -200+config->fiducial[0] || vtx_x > 200-config->fiducial[3] ||
         vtx_y < -200+config->fiducial[1] || vtx_x > 200-config->fiducial[4] ||
         vtx_z < 0+config->fiducial[2] || vtx_x > 500-config->fiducial[5]){ 
        in_fiducial = false;
      }
    }
    // If fiducial volume has APA and CPA cuts
    else if(std::find(config->fiducial.begin(), config->fiducial.end(), -1) == config->fiducial.end() && config->fiducial.size() == 8){
      if(vtx_x < -200+config->fiducial[0] || vtx_x > 200-config->fiducial[3] ||
         vtx_y < -200+config->fiducial[1] || vtx_x > 200-config->fiducial[4] ||
         vtx_z < 0+config->fiducial[2] || vtx_x > 500-config->fiducial[5] ||
         (vtx_x > -config->fiducial[6] && vtx_x < config->fiducial[6]) ||
         (vtx_z > 250-config->fiducial[7] && vtx_z < 250+config->fiducial[7])){ 
        in_fiducial = false;
      }
    }
    return in_fiducial;
  }
  
};

#endif
