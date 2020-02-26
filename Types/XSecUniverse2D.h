#ifndef XSECUNIVERSE2D_H
#define XSECUNIVERSE2D_H

// Structure for holding 2D systematic error information
class XSecUniverse2D
{
  public:

  TString hname;  // Base name

  TH2Poly* generated;   // True generated histogram
  TH2Poly* selected;    // True selected histogram
  TH2Poly* background;  // Selected background histogram
  TH2D* migration;      // Migration matrix

  // Constructor (default)
  XSecUniverse2D(TH2Poly* hist, TString name){
    hname = name;

    // Create empty histograms
    generated = (TH2Poly*)hist->Clone(name+"_gen");
    generated->ClearBinContents();
    selected = (TH2Poly*)hist->Clone(name+"_sel");
    selected->ClearBinContents();
    background = (TH2Poly*)hist->Clone(name+"_bkg");
    background->ClearBinContents();

    // Create empty migration matrix
    size_t nbins = generated->GetNumberOfBins();
    migration = new TH2D(name+"_resp", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
  }

  // Scale appropriate histogram
  void Scale(double scale){
    // Only background actually needs to be scaled
    background->Scale(scale);
    // Also scale by width if not total xsec
    size_t nbins = background->GetNbinsX();
    if(nbins > 1) background->Scale(1., "width");
  }

  // Calculate response matrix from migration matrix
  TH2D* Response(){
    size_t nbins = generated->GetNumberOfBins();
    TH2D* response = new TH2D(hname+"_resp", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    for(size_t i = 0; i <= nbins+1; i++){
      double total = 0;
      for(size_t j = 0; j <= nbins+1; j++){
        total += (double)migration->GetBinContent(i, j);
      }
      for(size_t j = 0; j <= nbins+1; j++){
        response->SetBinContent(i, j, (double)migration->GetBinContent(i, j)/total);
        if(std::isnan(response->GetBinContent(i, j))) response->SetBinContent(i, j, 0);
      }
    }

    return response;
  }

};

#endif
