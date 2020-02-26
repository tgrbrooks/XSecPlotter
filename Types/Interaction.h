#ifndef INTERACTION_H
#define INTERACTION_H

// Structure for holding interaction information
class Interaction
{
  public:

  bool selected;        // Is event selected
  bool true_selected;   // Does true interaction match selection
  bool in_fv;           // Is event in fiducial volume
  std::string fsi;      // Final state interaction of event
  std::string int_type; // Interaction type of event
  std::string nu_type;  // Neutrino type of event
  std::vector<double> variables;      // Plotting variables
  std::vector<double> true_variables; // Truth plotting variables

  // Constructor
  Interaction(bool s, bool ts, bool fv, std::string f, std::string i, std::string n, std::vector<double> v, std::vector<double> tv)
  {
    selected = s;
    true_selected = ts;
    in_fv = fv;
    fsi = f;
    int_type = i;
    nu_type = n;
    variables = v;
    true_variables = tv;
  }
};

#endif
