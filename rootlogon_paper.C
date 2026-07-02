{
  // ===== Global style for all ROOT canvases =====
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);

  // Default canvas size
  gStyle->SetCanvasDefW(800);
  gStyle->SetCanvasDefH(600);
  gStyle->SetCanvasDefX(50);
  gStyle->SetCanvasDefY(50);
  
  // Times-like font
  gStyle->SetTextFont(132);
  gStyle->SetLabelFont(132, "XYZ");
  gStyle->SetTitleFont(132, "XYZ");
  gStyle->SetLegendFont(132);
  gStyle->SetStatFont(132);

  // Canvas / pad
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);

  // Margins
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.17);

  // Axis fonts and sizes
  gStyle->SetTitleSize(0.055, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetTitleOffset(1.0, "X");
  gStyle->SetTitleOffset(1.15, "Y");

  // Axis ticks
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTickLength(0.03, "X");
  gStyle->SetTickLength(0.03, "Y");
  gStyle->SetNdivisions(510, "XYZ");

  // Lines / markers
  gStyle->SetFrameLineWidth(2);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineWidth(2);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetEndErrorSize(4);

  // Legend
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(132);

  // Palette
  gStyle->SetPalette(kBird);

  // Apply
  gROOT->ForceStyle();
}
