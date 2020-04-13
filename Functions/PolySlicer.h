#ifndef POLYSLICER_H
#define POLYSLICER_H

// Slice 2D poly hist in Y assuming y binning is constant
TH1D* SlicePoly(TH2Poly* poly, size_t i, TString name, std::vector<double> ybins, std::vector<std::vector<double>> xbins, bool xsec=false){

  // Get the upper and lower y edge that corresponds to slice
  double ylow = ybins[i];
  double yup = ybins[i+1];
  double width = (yup - ylow);
  // Get the x binning that corresponds to slice
  std::vector<double> xbin = xbins[i];
  double xbin_array[xbin.size()];
  std::copy(xbin.begin(), xbin.end(), xbin_array);
  // Create a 1D histogram with correct binning
  TH1D* hist = new TH1D(name, "", xbin.size()-1, xbin_array);
  // Loop over the poly bins
  for(auto const& obj : *poly->GetBins()){
    // If bin inside 
    TH2PolyBin *bin = (TH2PolyBin*)obj;
    double meany = (bin->GetYMax() + bin->GetYMin())/2.;
    double meanx = (bin->GetXMax() + bin->GetXMin())/2.;
    if(meany >= ylow && meany <= yup){
      int bin_i = hist->GetXaxis()->FindBin(meanx);
      // This doesn't preserve scaling
      hist->SetBinContent(bin_i, bin->GetContent());
    }
  }
  // FIXME work out how to scale by xsec
  if(!xsec) hist->Scale(1/width, "width");
  return hist;


}

#endif
