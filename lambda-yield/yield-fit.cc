//const vector<int> runnumbers = {2447};
const vector<int> runnumbers = {2447, 2449, 2450, 2451, 2452, 2453, 2454, 2456, 2457, 2458, 2459, 2460, 2462, 2463, 2465, 2466, 2468};



double asymExpBkgPlusGauss(double *x, double *par)
{
  double xx = x[0];

  // background parameters
  double A_bkg = par[0];
  double x0    = par[1];
  double tauL  = par[2];
  double tauR  = par[3];
  double c0    = par[4];

  // Gaussian signal parameters
  double A_sig = par[5];
  double mean  = par[6];
  double sigma = par[7];
  
  double t = xx - x0;
  double bkg = c0;
  if (t >= 0) {
    bkg += A_bkg * (TMath::Exp(-t / tauR)- TMath::Exp(-t / tauL));
  }
    /*
  double bkg;
  if (xx < x0) {
    bkg = c0 + A_bkg * TMath::Exp((xx - x0) / tauL);
  } else {
    bkg = c0 + A_bkg * TMath::Exp(-(xx - x0) / tauR);
  }
    */

  double sig = A_sig * TMath::Exp(-0.5 * TMath::Power((xx - mean) / sigma, 2));

  return bkg + sig;
}

double asymExpBkgOnly(double *x, double *par)
{
  double xx = x[0];

  double A_bkg = par[0];
  double x0    = par[1];
  double tauL  = par[2];
  double tauR  = par[3];
  double c0    = par[4];

  double t = xx - x0;
  double bkg = c0;
  if (t >= 0) {
    bkg += A_bkg * (TMath::Exp(-t / tauR)- TMath::Exp(-t / tauL));
  }

  return bkg;
}

double gaussOnly(double *x, double *par)
{
  double xx = x[0];

  double A_sig = par[0];
  double mean  = par[1];
  double sigma = par[2];

  return A_sig * TMath::Exp(-0.5 * TMath::Power((xx - mean) / sigma, 2));
}


void yield_fit(){

  const double xmin = 1.07;
  const double xmax = 1.20;
  string inroot;
  string outpdf;
  if (runnumbers.size() == 1){
    outpdf = Form("result/yield-fit-run%05d.pdf", runnumbers.front());
    inroot = Form("result/cut-condition-run%05d.root", runnumbers.front());
  }
  else{
    outpdf = Form("result/yield-fit-run%05d-%05d-n%zu.pdf",
                  runnumbers.front(), runnumbers.back(), runnumbers.size());
    inroot = Form("result/cut-condition-run%05d-%05d-n%zu.root",
		   runnumbers.front(), runnumbers.back(), runnumbers.size());
  }


  TCanvas *c1 = new TCanvas("c1", "c1", 900, 700);
  TPaveText *p = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
  p->AddText("yield-fit.cc");
  if (runnumbers.size() == 1)
    p->AddText(Form("run%05d", runnumbers.front()));
  else
    p->AddText(Form("runs %05d-%05d, n=%zu", runnumbers.front(), runnumbers.back(), runnumbers.size()));
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",
                  now.GetYear(), now.GetMonth(), now.GetDay(),
                  now.GetHour(), now.GetMinute(), now.GetSecond()));
  p->Draw();
  c1->Print((outpdf + "(").c_str());
  c1->Clear();
  
  TFile *file = TFile::Open(inroot.c_str(),"READ");
  auto hist = (TH1D*)file->Get(Form("hist_mass_tight"));
  hist->GetListOfFunctions()->Clear();
  
  TF1 *fTot = new TF1("fTot", asymExpBkgPlusGauss,xmin,xmax,8);
  fTot->SetParameters(913*17,1.0775,0.0225,0.045,6,714*17,1.109,0.006);


  hist->Fit(fTot,"R");

  TF1 *fBkg = new TF1("fBkg", asymExpBkgOnly, xmin, xmax, 5);

  fBkg->SetParameters(
		      fTot->GetParameter(0),
		      fTot->GetParameter(1),
		      fTot->GetParameter(2),
		      fTot->GetParameter(3),
		      fTot->GetParameter(4)
		      );


  TF1 *fSig = new TF1("fSig", gaussOnly, xmin, xmax, 3);

  fSig->SetParameters(
		      fTot->GetParameter(5),
		      fTot->GetParameter(6),
		      fTot->GetParameter(7)
		      );
  

  double xmin_cut = 1.09;
  double xmax_cut = 1.13;


  double binw = hist->GetBinWidth(1);
  double NSig = fSig->Integral(xmin_cut,xmax_cut) / binw;
  double NBkg = fBkg->Integral(xmin_cut,xmax_cut) / binw;
  
  double Mean = fTot->GetParameter(6);
  double Sig = fTot->GetParameter(7);
  
  

  hist->SetFillColorAlpha(kBlue, 0.1);
  hist->GetYaxis()->SetTitle("Counts/(0.001 GeV/#it{c}^{2})");
  hist->GetXaxis()->SetTitle("#it{m}_{p#pi^{-}} [GeV/#it{c}^{2}]");
  hist->Draw();
  fTot->SetLineColor(kGreen+2);
  fBkg->SetLineColor(kBlack);
  fSig->SetLineColor(kRed);
  fTot->SetLineWidth(3);
  fBkg->SetLineWidth(3);
  fSig->SetLineWidth(3);
  fTot->Draw("same");
  fBkg->SetLineStyle(2);
  
  fBkg->Draw("same");
  //fSig->SetLineStyle(3);
  fSig->Draw("same");

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.55, 0.54, Form("Yield = %.0f", NSig));
  lat.DrawLatex(0.55, 0.49, Form("#mu = %.5f", Mean));
  lat.DrawLatex(0.55, 0.44, Form("#sigma = %.5f", Sig));
  lat.DrawLatex(0.55, 0.39, Form("Bkg = %.0f", NBkg));
  
  c1->Print((outpdf + ")").c_str());




}
