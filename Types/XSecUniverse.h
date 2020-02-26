#ifndef XSECUNIVERSE_H
#define XSECUNIVERSE_H

// Structure for holding 1D systematic error information
class XSecUniverse
{
  public:

  TString hname;  // Base name

  TH1D* generated;  // True generated histogram
  TH1D* selected;   // True selected histogram
  TH1D* background; // Selected background histogram
  TH2D* migration;  // Migration matrix

  // Constructor (default)
  XSecUniverse(TH1D* hist, TString name){
    hname = name;

    // Create empty histograms
    generated = (TH1D*)hist->Clone(name+"_gen");
    generated->Reset();
    selected = (TH1D*)hist->Clone(name+"_sel");
    selected->Reset();
    background = (TH1D*)hist->Clone(name+"_bkg");
    background->Reset();

    // Create empty migration matrix
    size_t nbins = generated->GetNbinsX();
    double edges[nbins+1];
    generated->GetXaxis()->GetLowEdge(edges);
    edges[nbins] = generated->GetXaxis()->GetBinUpEdge(nbins);
    migration = new TH2D(name+"_mig", "", nbins, edges, nbins, edges);
  }

  // Scale the appropriate histograms
  void Scale(double scale){
    // Only background actually needs to be scaled
    background->Scale(scale);
    // Also scale by width if not total xsec
    size_t nbins = background->GetNbinsX();
    if(nbins > 1) background->Scale(1., "width");
  }

  // Calculate the response matrix from the migration matrix
  TH2D* Response(){
    size_t nbins = generated->GetNbinsX();
    TH2D* response = new TH2D(hname+"_resp", "", nbins, 1, nbins+1, nbins, 1, nbins+1);
    for(size_t i = 0; i <= nbins+1; i++){
      double total = 0;
      for(size_t j = 0; j <= nbins+1; j++){
        total+=(double)migration->GetBinContent(i, j);
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
