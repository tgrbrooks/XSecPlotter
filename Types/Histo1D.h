#ifndef HISTO1D_H
#define HISTO1D_H

// Structure for holding 1D histogram information
class Histo1D
{
  public:

  TH1D* total_hist;
  THStack* stacked_hist;
  TLegend* legend;

  Histo1D()
  {
  }

};

#endif
