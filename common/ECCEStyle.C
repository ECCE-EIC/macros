//
// ECCE Style, based on a style file from BaBar, v0.1
//

#include <iostream>

#include "ECCEStyle.h"

#include "TROOT.h"

void SetECCEStyle ()
{
  static TStyle* ecceStyle = 0;
  std::cout << "ECCEStyle: Applying nominal settings." << std::endl ;
  if ( ecceStyle==0 ) ecceStyle = ECCEStyle();
  gROOT->SetStyle("ECCE");
  gROOT->ForceStyle();
}

TStyle* ECCEStyle() 
{
  TStyle *ECCEStyle = new TStyle("ECCE","ECCE style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  ECCEStyle->SetFrameBorderMode(icol);
  ECCEStyle->SetFrameFillColor(icol);
  ECCEStyle->SetCanvasBorderMode(icol);
  ECCEStyle->SetCanvasColor(icol);
  ECCEStyle->SetPadBorderMode(icol);
  ECCEStyle->SetPadColor(icol);
  ECCEStyle->SetStatColor(icol);
  //ECCEStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  ECCEStyle->SetPaperSize(20,26);

  // set margin sizes
  ECCEStyle->SetPadTopMargin(0.05);
  ECCEStyle->SetPadRightMargin(0.05);
  ECCEStyle->SetPadBottomMargin(0.16);
  ECCEStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  ECCEStyle->SetTitleXOffset(1.4);
  ECCEStyle->SetTitleYOffset(1.4);

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.05;
  ECCEStyle->SetTextFont(font);

  ECCEStyle->SetTextSize(tsize);
  ECCEStyle->SetLabelFont(font,"x");
  ECCEStyle->SetTitleFont(font,"x");
  ECCEStyle->SetLabelFont(font,"y");
  ECCEStyle->SetTitleFont(font,"y");
  ECCEStyle->SetLabelFont(font,"z");
  ECCEStyle->SetTitleFont(font,"z");
  
  ECCEStyle->SetLabelSize(tsize,"x");
  ECCEStyle->SetTitleSize(tsize,"x");
  ECCEStyle->SetLabelSize(tsize,"y");
  ECCEStyle->SetTitleSize(tsize,"y");
  ECCEStyle->SetLabelSize(tsize,"z");
  ECCEStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  ECCEStyle->SetMarkerStyle(20);
  ECCEStyle->SetMarkerSize(1.2);
  ECCEStyle->SetHistLineWidth(2.);
  ECCEStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  //ECCEStyle->SetErrorX(0.001);
  // get rid of error bar caps
  ECCEStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  ECCEStyle->SetOptTitle(0);
  //ECCEStyle->SetOptStat(1111);
  ECCEStyle->SetOptStat(0);
  //ECCEStyle->SetOptFit(1111);
  ECCEStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  ECCEStyle->SetPadTickX(1);
  ECCEStyle->SetPadTickY(1);

  // legend modificatin
  ECCEStyle->SetLegendBorderSize(0);
  ECCEStyle->SetLegendFillColor(0);
  ECCEStyle->SetLegendFont(font);


#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
  std::cout << "ECCEStyle: ROOT6 mode" << std::endl;
  ECCEStyle->SetLegendTextSize(tsize);
  ECCEStyle->SetPalette(kBird);
#else
  std::cout << "ECCEStyle: ROOT5 mode" << std::endl;
  // color palette - manually define 'kBird' palette only available in ROOT 6
  Int_t alpha = 0;
  Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
  Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
  TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
#endif

  ECCEStyle->SetNumberContours(80);

  return ECCEStyle;

}
