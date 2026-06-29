const int runnumber = 2447;

void inv_mass(){
  string dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/physics-735";

  TFile *file = new TFile(Form("%s/run%05d_DstTPCHelixTracking.root",
                               dir.c_str(), runnumber));
  auto hist_inv_L = (TH1D*)file->Get("Lambda_Mass");

  string outpdf = Form("result/inv-mass-run%05d.pdf", runnumber);

  gStyle->SetOptFit(1111);

  TCanvas *c1 = new TCanvas("c1","c1",900,700);
  hist_inv_L->SetLineColor(kBlue);
  hist_inv_L->Draw("hist");

  // fit range: Lambda peak 주변만
  double fit_min = 1.085;
  double fit_max = 1.145;

  // Lambda mass 근처 Gaussian + 2nd order polynomial background
  TF1 *fit = new TF1("fit",
    "[0]*TMath::Gaus(x,[1],[2],0) + [3] + [4]*x + [5]*x*x",
    fit_min, fit_max);

  fit->SetParNames("Amp", "Mean", "Sigma", "p0", "p1", "p2");

  fit->SetParameters(
    hist_inv_L->GetMaximum(),  // Amp
    1.115,                     // Mean
    0.004,                     // Sigma
    100, -100, 10              // background initial guess
  );

  fit->SetParLimits(1, 1.108, 1.122);
  fit->SetParLimits(2, 0.001, 0.015);

  hist_inv_L->Fit(fit, "R");

  // background only
  TF1 *bkg = new TF1("bkg", "[0] + [1]*x + [2]*x*x", fit_min, fit_max);
  bkg->SetParameters(fit->GetParameter(3),
                     fit->GetParameter(4),
                     fit->GetParameter(5));
  bkg->SetLineColor(kGreen+2);
  bkg->SetLineStyle(2);
  bkg->Draw("same");

  // peak only
  TF1 *sig = new TF1("sig", "[0]*TMath::Gaus(x,[1],[2],0)", fit_min, fit_max);
  sig->SetParameters(fit->GetParameter(0),
                     fit->GetParameter(1),
                     fit->GetParameter(2));
  sig->SetLineColor(kRed);
  sig->Draw("same");

  fit->SetLineColor(kMagenta);
  fit->Draw("same");

  // Yield 계산: mean ± 3 sigma 안의 signal integral
  double mean  = fit->GetParameter(1);
  double sigma = fit->GetParameter(2);

  double yield_min = mean - 3.0*sigma;
  double yield_max = mean + 3.0*sigma;

  double binw = hist_inv_L->GetBinWidth(1);

  double sig_yield = sig->Integral(yield_min, yield_max) / binw;
  double bkg_yield = bkg->Integral(yield_min, yield_max) / binw;
  double total_yield = fit->Integral(yield_min, yield_max) / binw;

  cout << "Lambda mean  = " << mean << endl;
  cout << "Lambda sigma = " << sigma << endl;
  cout << "Yield range  = " << yield_min << " - " << yield_max << endl;
  cout << "Signal yield = " << sig_yield << endl;
  cout << "Bkg yield    = " << bkg_yield << endl;
  cout << "Total yield  = " << total_yield << endl;

  TLegend *leg = new TLegend(0.55,0.60,0.88,0.82);
  leg->AddEntry(hist_inv_L, "Data", "l");
  leg->AddEntry(fit, "Gaussian + pol2", "l");
  leg->AddEntry(sig, "Lambda peak", "l");
  leg->AddEntry(bkg, "Background", "l");
  leg->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.55,0.54,Form("Yield = %.0f", sig_yield));
  lat.DrawLatex(0.55,0.49,Form("#mu = %.5f", mean));
  lat.DrawLatex(0.55,0.44,Form("#sigma = %.5f", sigma));

  c1->SaveAs(outpdf.c_str());
}
