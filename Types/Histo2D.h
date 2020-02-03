#ifndef HISTO2D_H
#define HISTO2D_H

// Structure for holding 2D histogram information
class Histo2D
{
  public:

  TH2D* total_hist;
  std::vector<THStack*> stacked_hist;
  TLegend* legend;

  Histo2D(TH2D* hist, std::pair<std::vector<THStack*>, TLegend*> stack2D)
  {
    total_hist = hist;
    stacked_hist = stack2D.first;
    legend = stack2D.second;
  }

};

#endif
