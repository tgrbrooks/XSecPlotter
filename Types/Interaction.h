#ifndef INTERACTION_H
#define INTERACTION_H

// Structure for holding interaction information
class Interaction
{
  public:

  bool selected;
  bool true_selected;
  std::string fsi;
  std::string int_type;
  std::string nu_type;
  std::vector<double> variables;
  std::vector<double> true_variables;

  Interaction(bool s, bool ts, std::string f, std::string i, std::string n, std::vector<double> v, std::vector<double> tv)
  {
    selected = s;
    true_selected = ts;
    fsi = f;
    int_type = i;
    nu_type = n;
    variables = v;
    true_variables = tv;
  }
};

#endif
