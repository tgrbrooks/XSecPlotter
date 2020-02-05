// Set some global style configurations here
void SetStyle(){

  Int_t font = 62;
  Double_t font_size = 0.06;
  Double_t line_width = 3.;
  Double_t label_offset = 0.01;

  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);
  gStyle->SetMarkerStyle(8);
  // Widths
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameLineWidth(line_width);
  gStyle->SetGridWidth(line_width);
  gStyle->SetLineWidth(line_width);
  // Fonts
  gStyle->SetTitleFont(font, "title");
  gStyle->SetTitleFont(font, "x");
  gStyle->SetTitleFont(font, "y");
  gStyle->SetTitleFont(font, "z");
  gStyle->SetLabelFont(font, "x");
  gStyle->SetLabelFont(font, "y");
  gStyle->SetLabelFont(font, "z");
  gStyle->SetTextFont(font);
  gStyle->SetLegendFont(font);
  // Sizes
  gStyle->SetTitleSize(1.3*font_size, "title");
  gStyle->SetTitleSize(font_size, "x");
  gStyle->SetTitleSize(font_size, "y");
  gStyle->SetTitleSize(font_size, "z");
  gStyle->SetLabelSize(font_size, "x");
  gStyle->SetLabelSize(font_size, "y");
  gStyle->SetLabelSize(font_size, "z");
  gStyle->SetMarkerSize(0.6);
  // Offsets
  gStyle->SetLabelOffset(label_offset, "x");
  gStyle->SetLabelOffset(label_offset, "y");
  gStyle->SetLabelOffset(label_offset, "z");
  // Legend
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

}

